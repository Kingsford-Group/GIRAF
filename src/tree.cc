#include "tree.h"
#include <algorithm>

//=========================================================================
// Tree Output
//=========================================================================

//
// Convert spaces to _
//
string
FixSpaces(const string & s)
{
  string r;
  for(unsigned i=0;i<s.length();i++) r += (s[i]==' ')?'_':s[i];
  return r;
}

//
// Indent the given number of spaces
//
void
Indent(
  ostream & s,
  int i)
{
  for(;i>0;i--) s << ' ';
}

//
// Do most of the work to output the tree
//
void
WriteTree_Recurse(
  ostream & out, 
  const TreeNode * T,
  int indent,
  int options
  )
{
  assert(T!=0);
  bool len = !(options & NoLengths);
  bool ids = options & InternalIDs;
  bool sameline = options & SameLine;

  // for leaves, write name
  if(T->children.size() == 0) 
  {
    if(!sameline) Indent(out, indent);
    out << FixSpaces(T->id);
    if(len) out << ":" << T->length;
  }
  else
  {
    if(!sameline) Indent(out, indent);
    out << "(";
    if(!sameline) out << endl;
    for(TreeNode::const_child_iterator child = T->children.begin();
        child != T->children.end();
        ++child)
    {
      WriteTree_Recurse(out, *child, indent+2, options);
      if(child != T->children.end()-1) out << ",";
      if(!sameline) out << endl;
    }
    if(!sameline) Indent(out, indent);
    out << ")";
    if(ids && T->id != "") 
    {
      out << "\"" << T->id<< "\"";
    }
    if(len) out << ":" << T->length;
  }
}

//
// Write a tree in NH format; The output options are
// encoded as an | of the above enum
//
void
WriteTree(
  ostream & out,
  const TreeNode * T,
  int options)
{
  WriteTree_Recurse(out, T, 0, options);
  out << ";";
}

//
// Output the tree as a dot file
//
void
WriteTreeAsDot_Recurse(
  ostream & out,
  TreeNode * T
  )
{
  // write out my node
  out << T->id << ";" << endl;

  // write out edges to children
  for(TreeNode::child_iterator C = T->children.begin();
      C != T->children.end();
      ++C)
  {
    out << T->id << ";" << endl;
    WriteTreeAsDot_Recurse(out, (*C));
  }
}

//
// Output the tree as a dot file
//
void
WriteTreeAsDot(
  ostream & out,
  TreeNode * T,
  int index)
{
  out << "digraph G {" << endl;
  out << "graph [rotate=90,size=\"30,8\",page=\"8.5,11\"]" << endl;
  out << "node [shape=record,width=0,height=0]" << endl;
  WriteTreeAsDot_Recurse(out, T);
  out << "}" << endl;
}

//
// Caclulate the depth of each node & assign it to depth
//
void
CalcDepth(
  TreeNode * T,
  int depth)
{
  assert(T!=0);
  T->depth = depth;
  for_each(T->children.begin(), T->children.end(), 
           bind2nd(ptr_fun(CalcDepth), depth+1));
}


//=========================================================================
// Tree Input
//=========================================================================

void
ReadNexTranslate(
    istream & inp, 
    NodeNameMapping * mapping
    )
{
    vector<string> fields;
    string line;
    while(getline(inp, line)) 
    {
        // if the line starts with a translate command
        line = Trim(line);
        if(Upcase(line.substr(0, 9)) == "TRANSLATE") 
        {
            while(getline(inp, line))
            {
                SplitString(Trim(line), ' ', fields);
                if (fields.size() == 0) continue;
                if (fields[0] == ";") break; // semi on line by itself
                DIE_IF(fields.size() < 2, "Bad NEXUS translate command");

                unsigned namelen = fields[1].length();
                if (fields[1].at(namelen-1) == ',' || fields[1].at(namelen-1)==';')
                {
                    (*mapping)[fields[0]] = fields[1].substr(0, namelen-1);
                }
                else
                {
                    (*mapping)[fields[0]] = fields[1].substr(0, namelen);
                }

                if(fields[1].at(fields[1].length()-1) == ';') break;
            }

            break;
        }

        // stop once you see a tree command
        if(Upcase(line.substr(0,5)) == "TREE ") break;
    }
}


// apply the giving mapping to the leaves of the tree
void
TranslateLeaves(
    TreeNode * T, 
    NodeNameMapping & mapping
    )
{
    if (T->children.empty())
    {
        if(mapping.count(T->id) == 1) 
        {
            T->id = mapping[T->id];
        }
        else
        {
            WARN("missing NEXUS translate value for a leaf");
            cerr << "ID=" << T->id << std::endl;
        }
    }
    else
    {
        for(TreeNode::child_iterator C = T->children.begin();
            C != T->children.end();
            ++C)
        {
            TranslateLeaves(*C, mapping);
        }
    }
}


//
// Assign valid pointers to the leaves members of the tree
//
void
AssignLeavesLists(
  TreeNode * T)
{
  assert(T);
  T->leaves.clear();  

  // leaves just contain themselves
  if(T->children.empty())
  {
    T->leaves.push_back(T);
    return;
  }

  // leaves of T is the union of the leaves of its children
  for(TreeNode::child_iterator C = T->children.begin();
      C != T->children.end();
      ++C)
  {
    AssignLeavesLists(*C);
    copy((*C)->leaves.begin(), (*C)->leaves.end(), back_inserter(T->leaves));
  }
}


//
// Give each non-labeled node an id.
//
void
AssignIDs(
  TreeNode * T,
  int * i)
{
  if(T->id.empty())
  {
    ostringstream ids;
    ids << "n" << (*i)++;
    T->id = ids.str();
  }
  for_each(T->children.begin(), T->children.end(), 
           bind2nd(ptr_fun(AssignIDs), i));
}


//
// replace all 0 length edges with 'fix' & warn if we find any such edges
// 
void
CheckFixLengths(TreeNode * T, double fix)
{
    assert(T);
    if(T->length == 0) 
    {
        cerr << "warning: edge " << T->id 
             << " is 0; setting to " << fix << endl;
        T->length = fix;
    }
    DIE_IF(T->length < 0, "Tree has negative edge lenghts; can't use.");
    for_each(T->children.begin(), T->children.end(), 
        bind2nd(ptr_fun(CheckFixLengths), fix));
}


//
// Read a double from a stream; Used by ReadTree_Recurse
//
double
ReadLength(istream & in)
{
  string len_str = "";
  char ch;
  while(in >> ch && (isdigit(ch) || ch == '-' || ch == '.' || ch == 'e' || ch == 'E' || ch == '+'))
  {
    len_str += ch;
  }
  assert(ch == ')' || ch == ',' || ch == ';' || ch == '[' || isspace(ch));
  in.unget();
  return (double)atof(len_str.c_str());
}

//
// Actually read most of the tree from a stream. Called from ReadTree.
//
void
ReadTree_Recurse(
  istream & in, 
  TreeNode * P)
{
  TreeNode * C = 0;
  bool done = false;

  int line_number = 0;
  char ch;
  while(!done && in >> ch)
  {
    switch(ch)
    {
      // skip all whitespace
      case '\n': line_number++; continue;
      case ' ': case '\t': continue;

      // start of a new internal child
      case '(': 
        assert(C==0);
        C = new TreeNode;
        ReadTree_Recurse(in, C); 
        break;

      // finished reading a child, so add it
      case ')': done = true;  // fall through
      case ',': 
        ERROR_IF(C==0,line_number,"empty leaf name or internal node has no children.");
        P->children.push_back(C);
        C->parent = P;
        C=0;
        break;

      // read edge label e.g. :5.20
      case ':':
        ERROR_IF(C==0,line_number,"unexpected colon (:) in tree file.");
        C->length = ReadLength(in);
        break;

      // comments
      case '[':
        while(in >> ch && ch != ']') {}
        ERROR_IF(!in, line_number, "missing end comment (])");
        break;

      // any normal character starts a label
      default: 
        if(C==0) C = new TreeNode;
        C->id += ch;
        break;
    }
  }
}


//
// Read a tree that is saved in (a simplified variant)
// of the new hampshire format
//
TreeNode *
ReadTree(istream & in)
{
  char ch;
  in >> ch;
  DIE_IF(ch != '(', "Tree file must start with '('");
  TreeNode * P = new TreeNode();
  ReadTree_Recurse(in, P);

  bool done = false;
  while(in >> ch && !done)
  {
    switch(ch)
    {
      case ' ': case '\t': case '\n': continue;
      case ';': done = true; break;
      case ':': P->length = ReadLength(in); break;
      default: DIE("bad NH format");
    }
  }
  DIE_IF(!done, "missing ; in tree file");
  CalcDepth(P, 0);
  return P;
}

// remove things between [] from the line and return a new line
// with them removed.
std::string
RemoveNEXComments(const string & line)
{
    std::ostringstream oss;

    // for every char in the input string
    for (std::size_t i = 0; i < line.length(); i++) {
        // if this is a start of a comment
        if (line[i] == '[') {
            while (line[i] != ']') i++;  // skip chars until we get past the comment
        } else {
            oss << line[i]; 
        }
    }
    return oss.str();
}


// Read a .nex file that contains a collection of trees
void
ReadNexTrees(
    istream & inp,
    NodeNameMapping * mapping, 
    vector<TreeNode *> & list_of_trees,
    int burnin
    )
{
    int internal_node_count = 0;
    int tree_count = 0;
    string line;
    while(getline(inp, line)) 
    {
        // if the line starts with a tree command
        line = Trim(line);
        if(Upcase(line.substr(0, 5)) == "TREE ") 
        {
            // skip the first burnin trees
            tree_count++;
            if(tree_count <= burnin) continue;

            line = RemoveNEXComments(line);

            // find the start of the tree, if it exists
            size_t paren = line.find('(');
            if (paren != string::npos) 
            {
                // actually read the tree
                istringstream ios(line.substr(paren));
                TreeNode * T = ReadTree(ios);

                // assign internal nodes to have ids
                internal_node_count = 0;
                AssignIDs(T, &internal_node_count);

                // assign leaves their real names
                if (mapping) TranslateLeaves(T, *mapping);

                // store leaf sets at each internal node
                AssignLeavesLists(T);

                // save the tree in the tree list
                list_of_trees.push_back(T);
            }
            else
            {
                WARN("Skipping tree-like line.");
            }
        }
    }
}


void
AppendAllTreeNodes(
    TreeNode * Root,
    vector<TreeNode *> & nodes
    )
{
    nodes.push_back(Root);
    for (TreeNode::child_iterator C = Root->children.begin();
        C != Root->children.end();
        ++C)
    {
        AppendAllTreeNodes(*C, nodes);
    }
}


// free the memory associated with a tree
void
DeleteTree(TreeNode * T)
{
    vector<TreeNode *> nodes;
    AppendAllTreeNodes(T, nodes);
    for (vector<TreeNode *>::iterator N = nodes.begin();
        N != nodes.end();
        ++N)
    {
        delete (*N);
    }
}

//=========================================================================
// Tree Properties
//=========================================================================

// compute the sum of all the branch lengths (including
// the root branch leading to teh root node).
double
TotalTreeLength(
    TreeNode * T
    )
{
    double sum = T->length;

    for (TreeNode::child_iterator C = T->children.begin();
         C != T->children.end();
         ++C)
    {
        sum += TotalTreeLength(*C);
    }

    return sum;
}


// divides all the distances in the tree by the total tree length
void
ScaleDistanceMatrix(
    TreeNode * T,
    DistanceMatrix & M
    )
{
    double total_length = TotalTreeLength(T);

    for (DistanceMatrix::iterator I = M.begin();
         I != M.end();
         ++I)
    {
        for(map<string, double>::iterator J = I->second.begin();
            J != I->second.end();
            ++J)
        {
            J->second /= total_length;
        }
    }
}

void
ComputeDistancesFrom(
    TreeNode * From,
    DistanceMatrix & M,
    TreeNode *T,
    TreeNode *Last,
    float dist
    )
{
    // for leaves, we store the current distance
    if (T->children.empty())
    {
        if(T != From) {
            string a = T->id;
            string b = From->id;
            if (a > b) swap(a,b);
            M[a][b] = dist;
        }
    } 
    else
    {
        // recurse to the neighbors we haven't already visited
        for (TreeNode::child_iterator C = T->children.begin();
             C != T->children.end();
             ++C)
        {
            if (*C != Last) 
            {
                ComputeDistancesFrom(From, M, *C, T, dist + (*C)->length);
            }
        }
    }

    if (T->parent && T->parent != Last) 
    {
        ComputeDistancesFrom(From, M, T->parent, T, dist + T->length);
    }
}


void
ComputeDistanceMatrix(
    TreeNode *T, 
    DistanceMatrix & M
    )
{
    if (T->children.empty())
    {
        //if(T->id == "F0137") cout << "DEBUG: " << T->id << endl; 
        ComputeDistancesFrom(T, M, T, T, 0.0);
    }
    else
    {
        for (TreeNode::child_iterator C = T->children.begin();
             C != T->children.end();
             ++C)
        {
            ComputeDistanceMatrix(*C, M);
        }
    }
}


vector<DistanceMatrix> *
AllDistanceMatrices(
    vector<TreeNode *> & trees
    )
{
    vector<DistanceMatrix> * all_dist = new vector<DistanceMatrix>(trees.size());
    DIE_IF(all_dist == 0, "Not enough memory...");

    int i = 0;
    for(vector<TreeNode*>::iterator T = trees.begin();
        T != trees.end();
        ++T)
    {
        WriteStatusNumber(cout, i);
        ComputeDistanceMatrix(*T, (*all_dist)[i]);
        ScaleDistanceMatrix(*T, (*all_dist)[i]);
        i++;
    }
    return all_dist;
}


void
PrintDistances(
    ofstream & out,
    vector<DistanceMatrix> & matrices)
{
    if (matrices.size() == 0) return;

    streamsize pp = out.precision();
    out.precision(21);

    for (DistanceMatrix::iterator I = matrices[0].begin();
        I != matrices[0].end();
        ++I)
    {
        for (map<string,double>::iterator J = I->second.begin(); 
            J != I->second.end();
            ++J)
        {
            if (I->first < J->first) 
            {
                out << I->first << " " << J->first;
                for (unsigned i = 0; i < matrices.size(); i++)
                {
                    out << " " << matrices[i][I->first][J->first];
                }

                out << endl;
            }
        }
    }
    out.precision(pp);
}


//
// Return the set of ids for the leaves of this node
//
set<string>
LeavesOf(
  const TreeNode * T)
{
  assert(T);
  set<string> leaves;
  for(vector<TreeNode*>::const_iterator C = T->leaves.begin();
      C != T->leaves.end();
      ++C)
  {
    leaves.insert((*C)->id);
  }
  return leaves;
}

