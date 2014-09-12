#include <cmath>
#include <list>
#include <map>
#include "tree.h"
#include "splits.h"
#include "dist.h"
#include "options.h"


#define PROG_NAME "build_incompat_graph"

typedef vector<pair<int, int> > IntEdgeList;

bool use_dist_opt = true;
bool out_pairs_opt = false;
bool out_unlabeled_opt = false;
bool all4tests_opt = false;
bool max_perl_compat_opt = false; //try to be like the perl version
double evalue_threshold = 0.01;

//===========================================================================
// Incompatability Graph
//===========================================================================

bool
ByLex(const pair<int,int> & a, const pair<int,int> & b)
{
    return (a.first < b.first || (a.first == b.first && a.second < b.second));
}

/*
 * Will create an edge list of incompatible splits between the two sets of
 * splits.
 */
void
CreateIncompatEdgeList(
    SplitDatabase & left_splits, 
    SplitDatabase & right_splits,
    IntEdgeList & E
    )
{
    // for every pair of splits
    for (SplitDatabase::iterator L = left_splits.begin();
         L != left_splits.end();
         ++L)
    {
        for (SplitDatabase::iterator R = right_splits.begin();
             R != right_splits.end();
             ++R)
        {
            // if the splits are incompatible, add the edge
            if(SplitsAreIncompatible(L->first, R->first)) 
            {
                E.push_back(make_pair(L->first.id, R->first.id));
            } 
        }
    }
    sort(E.begin(), E.end(), ByLex);
}


void
PrintIncompatGraph(
    ostream & out,
    IntEdgeList & E
    )
{
    for(IntEdgeList::iterator e = E.begin();
        e != E.end();
        ++e)
    {
        out << e->first << " " << e->second << endl;
    }
}


//===============================================================
// Graph Labeling Tests
//===============================================================

void
ComputeMovedMatrix(
    DistanceMatrix & pair_distances, 
    DistanceMatrix & is_greater,  // inout
    double & ge_freq, // out
    double & le_freq  // out
    )
{

    // compute the number of entries in the matrix
    long num_of_mw_tests = 0;
    for(DistanceMatrix::iterator M = is_greater.begin();
        M != is_greater.end();
        ++M)
    {
        num_of_mw_tests += M->second.size();
    }

    // for every entry, convert is_greater to -1, 0, or 1
    // 1 if is greater and passes the test
    // -1 if is not greater and passes the test
    // 0 if fails the test
    long ge_count, le_count;
    ge_count = le_count = 0;

    for (DistanceMatrix::iterator M = is_greater.begin();
         M != is_greater.end();
         ++M)
    {
        for(map<string, double>::iterator I = M->second.begin();
            I != M->second.end();
            ++I)
        {
            assert(pair_distances.find(M->first)!= pair_distances.end());
            assert(pair_distances[M->first].find(I->first) != 
                        pair_distances[M->first].end());
            double pval = expm1(pair_distances[M->first][I->first])+1;
            if(I->second && num_of_mw_tests * pval < evalue_threshold)
            {
                I->second = 1;
                ge_count++;
            }
            else if(!I->second && num_of_mw_tests * pval < evalue_threshold)
            {
                I->second = -1;
                le_count++;
            }
            else
            {
                I->second = 0;
            }
        }
    }

    ge_freq = ((float)ge_count) / num_of_mw_tests;
    le_freq = ((float)le_count) / num_of_mw_tests;
}


// create a mapping of split ids to split objects
void
IndexSplitsById(
    SplitDatabase & splits,     // in
    map<int, Split *> & index   // out
    )
{
    for (SplitDatabase::iterator S = splits.begin();
         S != splits.end();
         ++S)
    {
        // make sure id doesn't already exist
        assert(index.find(S->first.id) == index.end());

        index[S->first.id] = const_cast<Split*>(&(S->first));
    }
}


void
AddLabel(
    map<string, int> & labels,
    set<string> & I, 
    int & index
    )
{
    string str = SetAsString(I);
    if(labels.find(str) == labels.end()) 
    {
        labels[str] = index;
        index++;
    }
}

bool 
BySize(const set<string> & a, const set<string> & b)
{
    return a.size() < b.size();
}


struct CandidateSets
{
    CandidateSets(
        const set<string> & A, 
        const set<string> & B, 
        const set<string> & C, 
        const set<string> & D,
        map<string, int> & labels
        )
    {
        // sort the candidate sets by size
        vector<set<string> > Tmp;
        Tmp.push_back(A);
        Tmp.push_back(B);
        Tmp.push_back(C);
        Tmp.push_back(D);
        sort(Tmp.begin(), Tmp.end(), BySize);

        a = Tmp[0]; b = Tmp[1]; c = Tmp[2]; d = Tmp[3];

        ai = labels[SetAsString(a)];
        bi = labels[SetAsString(b)];
        ci = labels[SetAsString(c)];
        di = labels[SetAsString(d)];
    }

    set<string> a,b,c,d;
    long ai, bi, ci, di;
};


// for every edge in the incompatibility graph we construct the 4 candidate
// sets (stored in a list of CandidateSets structs in the same order as the
// edges in the IntEdgeList list). We also construct a mapping from each set
// (as a string) to an integer that we will use as an id in future output.
//
void
ConstructCandidateSets(
    SplitDatabase & left_splits,
    SplitDatabase & right_splits,
    IntEdgeList & IG,
    map<string, int> & labels, // out
    vector<CandidateSets> & sets // out
    )
{
    map<int, Split*> left_index;
    IndexSplitsById(left_splits, left_index);

    map<int, Split*> right_index;
    IndexSplitsById(right_splits, right_index);

    int curr_index = 0;

    // for every edge in the incompatibility graph
    for (IntEdgeList::iterator E = IG.begin();
         E != IG.end();
         ++E)
    {
        // get the splits corresponding to this graph edge
        Split * L = left_index[E->first];
        Split * R = right_index[E->second];

        // compute the intersection & complement
        set<string> I1;
        set<string> D1;
        I1.clear(); D1.clear();
        SetIntersection(L->first(), R->first(), I1);
        SetDifference(L->first(), R->first(), D1);

        AddLabel(labels, I1, curr_index);
        AddLabel(labels, D1, curr_index);

        set<string> I2;
        set<string> D2;
        I2.clear(); D2.clear();
        SetIntersection(L->second(), R->first(), I2);
        SetDifference(L->second(), R->first(), D2);

        AddLabel(labels, I2, curr_index);
        AddLabel(labels, D2, curr_index);

        sets.push_back(CandidateSets(I1, D1, I2, D2, labels));
    }
}


// Print labels, sorted by label index
void
PrintLabelMapping(
    ostream & out,
    map<string, int> & labels
    )
{
    map<int, string> reverse;
    for (map<string, int>::iterator L = labels.begin();
         L != labels.end();
         ++L)
    {
        reverse[L->second] = L->first;
    }

    for (map<int, string>::iterator R = reverse.begin();
         R != reverse.end();
         ++R)
    {
        out << R->first << " " << R->second << endl << endl;
    }
}


void
PrintAllCandidates(
    ostream & out,
    IntEdgeList & IG,
    vector<CandidateSets> & candidates
    )
{
    // for every edge in the incompatibility graph
    int i = 0;
    for (IntEdgeList::iterator E = IG.begin();
         E != IG.end();
         ++E, ++i)
    {
        out << E->first << " " << E->second << " ";
        CandidateSets abcd = candidates[i];

        out << abcd.ai << " " << abcd.bi << " " 
            << abcd.ci << " " << abcd.di << " " << endl;
    }
}


// return an item from the distance matrix,
// handling the symmetry of (a,b) and (b,a)
double
MovedItem(
    DistanceMatrix & moved_matrix,
    const string & a,
    const string & b
    )
{
    if (moved_matrix.find(a) != moved_matrix.end()) 
    {
        if (moved_matrix[a].find(b) != moved_matrix[a].end()) 
        {
            return moved_matrix[a][b];
        }
    }

    if (moved_matrix.find(b) != moved_matrix.end())
    {
        if (moved_matrix[b].find(a) != moved_matrix[b].end())
        {
            return moved_matrix[b][a];
        }
    }
    assert(false);
}

extern "C" { float betai(float, float, float); }

double
GetBinPval(long count, long n, double freq)
{
    if( ((double)count) / n > freq) 
    {
        return betai(count+1, n - count, freq);
    }
    return 1.0;
}


void
CompareSets(
    set<string> & a,
    set<string> & b,
    DistanceMatrix & moved_matrix,
    double ge_freq,
    double le_freq,

    double & ge_pval,  // out
    double & le_pval   // out
    )
{
    long ge_count, le_count;
    ge_count = le_count = 0;

    for (set<string>::iterator A = a.begin();
         A != a.end();
         ++A)
    {
        for(set<string>::iterator B = b.begin();
            B != b.end();
            ++B)
        {
            double mi = MovedItem(moved_matrix, *A, *B);
            if (mi < 0)
            {
                le_count += 1;
            }
            if (mi > 0)
            {
                ge_count += 1;
            }
        }
    }

    long n = a.size() * b.size();

    ge_pval = GetBinPval(ge_count, n, ge_freq);
    le_pval = GetBinPval(le_count, n, le_freq);

    //char tmp[1024];
    //sprintf(tmp, "(%ld, %ld, %f, %f ,,, %ld, %ld, %f, %f)\n",
    //        ge_count, n, ge_freq, ge_pval, 
    //        le_count, n, le_freq, le_pval);
    //cout << tmp;
}


// tests whether set a has moved relative to one of the other sets
bool
TestCandidate(
    set<string> & a,
    set<string> & b,
    set<string> & c,
    set<string> & d,

    DistanceMatrix & moved_matrix,
    double ge_freq,
    double le_freq
    )
{
    vector<set<string>*> others;
    others.push_back(&b);
    others.push_back(&c);
    others.push_back(&d);

    double greater, lesser;
    greater = lesser = 0;

    for (vector<set<string>*>::iterator I = others.begin();
         I != others.end();
         ++I)
    {
        double ge_pval, le_pval;
        CompareSets(a, **I, moved_matrix, ge_freq, le_freq, ge_pval, le_pval);

        if (ge_pval < evalue_threshold) {
            greater += (le_pval < evalue_threshold) ? 0.5 : 1.0;
        }
        if (le_pval < evalue_threshold) {
            lesser += (ge_pval < evalue_threshold) ? 0.5 : 1.0;
        }
    }

    return greater*lesser >= 0.5;
}


void
PrintFilteredLabeledGraph(
    ostream & out,
    IntEdgeList & IG,
    vector<CandidateSets> & candidates,
    DistanceMatrix & moved_matrix,
    double ge_freq,
    double le_freq
    )
{
    // for every edge in the incompatibility graph
    int i = 0;
    for (IntEdgeList::iterator E = IG.begin();
         E != IG.end();
         ++E, ++i)
    {
        out << E->first << " " << E->second << " ";
        CandidateSets abcd = candidates[i];

        if(TestCandidate(abcd.a, abcd.b, abcd.c, abcd.d, moved_matrix, ge_freq, le_freq)) 
        {
            out << abcd.ai << " ";
        }
        if(TestCandidate(abcd.b, abcd.a, abcd.c, abcd.d, moved_matrix, ge_freq, le_freq)) 
        {
            out << abcd.bi << " ";
        }
        if(TestCandidate(abcd.c, abcd.a, abcd.b, abcd.d, moved_matrix, ge_freq, le_freq)) 
        {
            out << abcd.ci << " ";
        }

        // we do the test on the largest set if requested or if the "largest" set
        // is the same size as the 3rd largest.
        
        if(!max_perl_compat_opt && (all4tests_opt || abcd.d.size() == abcd.c.size())) 
        {
            if(TestCandidate(abcd.d, abcd.a, abcd.b, abcd.c, moved_matrix, ge_freq, le_freq)) 
            {
                out << abcd.di;
            }
        }
        out << endl;
    }
}


//===================================================================================
// Main Program
//===================================================================================

// Options for build_incompat_graph

const char *GRAPH_OPTIONS = "h";

enum {GRAPH_DIST_OPT=1, GRAPH_BAD_OPT, OUT_PAIRS_OPT, OUT_UNLABELED_OPT, ALL4TESTS_OPT, VER09_OPT};

static struct option MAYBE_UNUSED graph_long_options[] = {
    {"use-dist", 1, 0, GRAPH_DIST_OPT},
    {"ignore-bad-options", 0, 0, GRAPH_BAD_OPT},
    {"debug-out-pairs", 0, 0, OUT_PAIRS_OPT},
    {"debug-out-unlabeled", 0, 0, OUT_UNLABELED_OPT},
    {"test-all-candidates", 1, 0, ALL4TESTS_OPT},
    {"version-0.9-compat", 0, 0, VER09_OPT},
    {0,0,0,0}
};

void
GraphUsage(bool show_cmd=true)
{
    if(show_cmd)
    {
        cerr << "build_incompat_graph [options] base1 base2 outbase" 
             << endl 
             << endl
             << "OPTIONS" 
             << endl;
    }

    if(show_cmd) 
    {
        cerr << "   --use-dist=[0,1]   : if 0 ignore distances (default 1)" << endl;
    }

    cerr << "   --test-all-candidates=[0,1] : if 1, test even large candidate sets (default 0)" << endl 
         << "   --debug-out-pairs     : output result of statistical tests (debugging only)" << endl
         << "   --debug-out-unlabeled : output unlabeled incompat graph (debugging only)" << endl
         << endl;
    if(show_cmd) exit(3);
}

static
int
ProcessOptions(int argc, char *argv[])
{
    bool ignore_bad_opt = false;
    int a;
    
    opterr = 0;
    // GNU ONLY!!! By setting optind to 0, we force getopt to reinitialize
    // this is required if we want giraf to call this program
    optind = 0; 
    while ((a=getopt_long(argc, argv, GRAPH_OPTIONS, graph_long_options, 0)) != -1)
    {
        switch(a)
        {
            case 'h': GraphUsage(); break;
            case GRAPH_DIST_OPT: use_dist_opt = (bool)atoi(optarg); break;
            case GRAPH_BAD_OPT: ignore_bad_opt = true; opterr=0; break;
            case OUT_PAIRS_OPT: out_pairs_opt = true; break;
            case OUT_UNLABELED_OPT: out_unlabeled_opt = true; break;
            case ALL4TESTS_OPT: all4tests_opt = (bool)atoi(optarg); break;
            case VER09_OPT: max_perl_compat_opt = true; break;
            default:
                if(!ignore_bad_opt) {
                    cerr << "Unknown option." << endl;
                    GraphUsage();
                }
        }
    }
    if (optind+2 >= argc) GraphUsage();
    return optind;
}


void
CheckInFile(istream & in, const string & name)
{
    if (!in) {
        cerr << PROG_NAME ": Can't find file " << name << endl;
        exit(3);
    }
}


int
main_build_incompat_graph(int argc, char * argv[])
{
    cout << PROG_NAME ": built on " << __DATE__ << endl;

    // read the command line
    int first_base_index = ProcessOptions(argc, argv);
    string base1 = argv[first_base_index];
    string base2 = argv[first_base_index + 1];
    string outbase = argv[first_base_index + 2];

    cout << PROG_NAME ": Use Distance = " << use_dist_opt << endl;
    cout << PROG_NAME ": In Files = " << base1 << " " << base2 << endl;
    cout << PROG_NAME ": Out Base = " << outbase << endl;
    if(max_perl_compat_opt) 
    {
        cout << PROG_NAME 
             << ": Trying to be as simlar to GIRAF Version 0.9 as possible." 
             << endl;
    }

    // read the splits
    cout << PROG_NAME ": reading left splits." << endl;
    string tmp;
    tmp = base1 + "_splits";
    ifstream left_splits_file(tmp.c_str());
    CheckInFile(left_splits_file, tmp);
    SplitDatabase left_splits; 
    ReadSplitsMapping(left_splits_file, left_splits); 
    left_splits_file.close();
    cout << PROG_NAME ": found " << left_splits.size() << " left splits." 
         << endl;

    cout << PROG_NAME ": reading right splits." << endl;
    tmp = base2 + "_splits";
    ifstream right_splits_file(tmp.c_str());
    CheckInFile(right_splits_file, tmp);
    SplitDatabase right_splits;
    ReadSplitsMapping(right_splits_file, right_splits);
    right_splits_file.close();
    cout << PROG_NAME ": found " << right_splits.size() << " right splits." 
         << endl;

    // Construct the incompatible split list
    cout << PROG_NAME ": finding incompatible splits." << endl;
    IntEdgeList IG;
    CreateIncompatEdgeList(left_splits, right_splits, IG);

    // Write out the incompatibility graph
    if (out_unlabeled_opt) 
    {
        cout << PROG_NAME ": writing incompatibility graph." << endl;
        tmp = outbase + "_graph";
        ofstream outgraph(tmp.c_str());
        PrintIncompatGraph(outgraph, IG);
        outgraph.close();
    }
    
    // get the candidates implied by the incompatibile splits
    cout << PROG_NAME ": getting candidate taxa sets." << endl;
    map<string, int> labels;
    vector<CandidateSets> candidates;
    ConstructCandidateSets(left_splits, right_splits, IG, labels, candidates);
    tmp = outbase + "_graph.labels";
    ofstream graph_labels(tmp.c_str());
    PrintLabelMapping(graph_labels, labels);
    graph_labels.close();

    // compute the labels for every edge
    cout << PROG_NAME ": computing labelled graph." << endl;
    tmp = outbase + "_graph.labelled";
    ofstream new_graph(tmp.c_str());
    if(!use_dist_opt)
    {
        // if user asked to not filter
        PrintAllCandidates(new_graph, IG, candidates);
    }
    else
    {
        // read the distances & compute the pair-test results
        cout << PROG_NAME ": computing pair distances." << endl;
        tmp = base1 + "_dist";
        ifstream left_dist_file(tmp.c_str());
        CheckInFile(left_dist_file, tmp);
        tmp = base2 + "_dist";
        ifstream right_dist_file(tmp.c_str());
        CheckInFile(right_dist_file, tmp);

        // output pair_test_file if requested (only for backward compat)
        ofstream *pair_test_file = 0;
        if (out_pairs_opt)
        {
            tmp = outbase + "_pair_test_results"; 
            pair_test_file = new ofstream(tmp.c_str());
        }
        DistanceMatrix pair_distances;
        DistanceMatrix is_greater;
        cout << PROG_NAME ": ";
        ComputePairDistances(left_dist_file, right_dist_file, max_perl_compat_opt,
            pair_test_file, pair_distances, is_greater);
        cout << endl;

        // close up the files
        left_dist_file.close();
        right_dist_file.close();
        if(pair_test_file) 
        {
            pair_test_file->close();
            delete pair_test_file;
        }

        // compute the pairs that seemed to have moved
        cout << PROG_NAME ": calculating distance statistics." << endl;
        double ge_freq, le_freq;
        ComputeMovedMatrix(pair_distances, is_greater, ge_freq, le_freq); 

        cout << PROG_NAME ": writing graph." << endl;
        PrintFilteredLabeledGraph(new_graph, IG, candidates, 
            is_greater, ge_freq, le_freq);
    }
    new_graph.close();
    return 0;
}
