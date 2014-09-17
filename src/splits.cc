#include "splits.h"
#include <algorithm>


string
SplitAsString(const Split & s)
{
    // construct the full set of taxa in sorted order
    set<string> taxa;
    set_union(s.first().begin(), s.first().end(), 
              s.second().begin(), s.second().end(), 
              inserter(taxa, taxa.begin())
              );

    // construct a signature string; note that the lexicagraphically first
    // taxa is always in the "*" set
    string sig = "";
    string Asymb = "*";
    string Bsymb = ".";
    for (set<string>::iterator I = taxa.begin();
        I != taxa.end();
        ++I)
    {
        if (I == taxa.begin() && s.first().find(*I) == s.first().end()) 
        {
            swap(Asymb, Bsymb);
        }
        sig += (s.first().find(*I) != s.first().end()) ?  Asymb : Bsymb;
    }
    return sig;
}


Split::Split(const set<string> & f, const set<string> & s)
    : _first(f), _second(s) 
{ 
    _rep = SplitAsString(*this);
}


Split::Split(const set<string> & f, const set<string> & s, int index)
    : _first(f), _second(s) 
{ 
    _rep = SplitAsString(*this);
    id = index;
}


// compare two splits
bool
operator==(const Split & a, const Split & b)
{
    return (a.first() == b.first() && a.second() == b.second()) || 
           (a.first() == b.second() && a.second() == b.first());
}


// Return the size of the intersection of two sets
int
SetIntersectionSize(
    const set<string> & a,
    const set<string> & b
    )
{
    set<string> I;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(I, I.begin()));
    return I.size();
}

// Return the intersection of two split sets
void
SetIntersection(
    const set<string> & a,
    const set<string> & b,
    set<string> & I
    )
{
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(I, I.begin()));
}

void
SetDifference(
    const set<string> & a,
    const set<string> & b,
    set<string> & D
    )
{
    set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(D, D.begin()));
}


// Return true iff the two splits are incompatible.
bool
SplitsAreIncompatible(
    const Split & a,
    const Split & b)
{
    return (SetIntersectionSize(a.first(), b.first()) > 0 &&
            SetIntersectionSize(a.first(), b.second()) > 0 &&
            SetIntersectionSize(a.second(), b.first()) > 0 &&
            SetIntersectionSize(a.second(), b.second()) > 0);
}


void
AddSplit(
    SplitDatabase & splits, 
    const set<string> & A,
    const set<string> & B,
    int tree_index
    )
{
    splits[Split(A,B)].insert(tree_index);
}


// walk down the tree recursively, adding the splits
void
AddSplits_Recurse(
    TreeNode * T,
    int tree_index,
    const set<string> & taxa,
    SplitDatabase & splits
    )
{
    // ignore trivial splits
    if (!T->children.empty())
    {
        // construct the split
        set<string> A = LeavesOf(T);
        set<string> B;
        set_difference(taxa.begin(), taxa.end(), A.begin(), A.end(), inserter(B, B.begin()));
        if (A.size() > B.size()) swap(A,B);

        // add it to the database
        AddSplit(splits, A, B, tree_index);

        // recurse on the children
        for(TreeNode::child_iterator C = T->children.begin();
            C != T->children.end();
            ++C)
        {
            AddSplits_Recurse(*C, tree_index, taxa, splits);
        }
    }
}


void
AllSplits(vector<TreeNode *> & trees, SplitDatabase & splits)
{
    // for every tree
    int tree_index = 0;
    for (vector<TreeNode*>::iterator Tree = trees.begin();
        Tree != trees.end();
        ++Tree)
    {
        // add its splits
        AddSplits_Recurse(*Tree, tree_index, LeavesOf(*Tree), splits);
        tree_index++;
    }
}


void
CullSplits(
    SplitDatabase & splits, 
    unsigned required_trees)
{
    SplitDatabase tmp;

    for (SplitDatabase::iterator S = splits.begin();
        S != splits.end();
        ++S)
    {
        if (S->second.size() >= required_trees)
        {
            tmp[S->first] = S->second;
        }
    } 
    splits = tmp;
}

//===================================================================================
// Split Printing and Reading
//===================================================================================

void
PrintSplitsReadable(
    ostream & out, 
    SplitDatabase & splits,
    unsigned min_occur
    )
{
    // for every split
    for (SplitDatabase::iterator S = splits.begin();
        S != splits.end();
        ++S)
    {
        if (S->second.size() > min_occur)
        {
            out << S->first.as_string() << " " << S->second.size() << endl;
        }
    }
}


void
PrintSplitsMapping(
    ostream & out,
    SplitDatabase & splits
    )
{
    // for every split
    int split = 0;
    for (SplitDatabase::iterator S = splits.begin();
        S != splits.end();
        ++S)
    {
        set<string> a = S->first.first();
        set<string> b = S->first.second();
        if (a.size() < b.size()) swap(a,b);

        out << split << " {" << SetAsString(a) << "} {" << SetAsString(b) << "}" << endl;
        split++;
    } 
}


// Will read a _splits file produced by WriteSplitsMapping
// Produces a SplitsDatabase with empty tree lists.
void
ReadSplitsMapping(
    istream & in,
    SplitDatabase & splits
    )
{
    vector<string> fields;

    string line;
    while(getline(in, line))
    {
        SplitString(line, ' ', fields);

        if(fields.size() >= 3)
        {
            int index = atoi(fields[0].c_str());

            set<string> A;
            set<string> B;

            bool in_first_set = true;
            for (unsigned i = 1; i < fields.size(); i++)
            {
                string taxon = fields[i];
                if (taxon[0] == '{') taxon = taxon.substr(1);
                if (taxon[taxon.length() - 1] == '}') 
                {
                    taxon = taxon.substr(0, taxon.length()-1);
                    (in_first_set?A:B).insert(taxon);
                    in_first_set = false;
                }
                else
                {
                    (in_first_set?A:B).insert(taxon);
                }
            }

            // add the split to the map...
            splits[Split(A,B,index)];
        }
    }
}


void
PrintTreesForSplits(
    ostream & out,
    int num_trees, 
    SplitDatabase & splits
    )
{
    out << ">> " << num_trees << " " << splits.size() << endl;

    int split = 0;
    for (SplitDatabase::iterator S = splits.begin();
        S != splits.end();
        ++S)
    {
        out << split << " " << S->second.size();
        for(set<int>::iterator I = S->second.begin();
            I != S->second.end();
            ++I)
        {
            out << " " << *I;
        }
        out << endl;
        split++;
    } 
}

