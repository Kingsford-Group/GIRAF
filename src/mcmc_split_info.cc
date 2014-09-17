#include "tree.h"
#include "splits.h"
#include "options.h"

#define PROG_NAME "mcmc_split_info"

int dist_opt = 1;
int burnin_opt = 500;
float cull_opt = 0.05;

// Options for mcmc_split_info
const char * SPLIT_OPTIONS = "h";

enum {DIST_OPT=1, BURNIN_OPT, CULL_OPT, SPLIT_BAD_OPT};

static struct option MAYBE_UNUSED split_long_options[] = {
    {"use-dist", 1, 0, DIST_OPT},
    {"burnin", 1, 0, BURNIN_OPT},
    {"cull", 1, 0, CULL_OPT},
    {"ignore-bad-options", 0, 0, SPLIT_BAD_OPT},
    {0,0,0,0}
};

void
SplitUsage(bool show_cmd = true)
{
    if(show_cmd) {
        cerr << PROG_NAME << " [options] nexfile.nex [file2.nex...]" << endl << endl;

        cerr << "OPTIONS" << endl;
    }
    cerr << "   --use-dist=[0,1] : if 1, compute the distances (default 1)" << endl
         << "   --burnin=N       : drop N trees" << endl
         << "   --cull=F         : drop splits that occur < F fraction time" << endl 
         << endl;
    if(show_cmd) exit(3);
}

static
int
ProcessOptions(int argc, char *argv[])
{
    bool ignore_bad_opt = false;
    opterr = 0;
    // GNU ONLY!!! By setting optind to 0, we force getopt to reinitialize
    // this is required if we want giraf to call this program
    optind = 0; 
    int a;
    while((a=getopt_long(argc, argv, SPLIT_OPTIONS, split_long_options, 0)) != -1)
    {
        switch(a)
        {
            case 'h': SplitUsage(); break;
            case DIST_OPT: dist_opt = atoi(optarg); break;
            case BURNIN_OPT: burnin_opt = atoi(optarg); break;
            case CULL_OPT: cull_opt = atof(optarg); break;
            case SPLIT_BAD_OPT: ignore_bad_opt = true; break;
            default:
                if(!ignore_bad_opt) {
                    cerr << "Unknown option." << endl;
                    SplitUsage();
                }
        }
    }
    if(optind+1 >= argc) SplitUsage();
    return optind;
}


int
main_mcmc_split_info(int argc, char *argv[])
{
    cout << PROG_NAME ": built on " << __DATE__ << endl;

    int first_file_index = ProcessOptions(argc, argv);
    string basename = argv[first_file_index];

    cout << PROG_NAME ": Burn-in = " << burnin_opt << endl;
    cout << PROG_NAME ": Distance = " << dist_opt << endl;
    cout << PROG_NAME ": Cull = " << cull_opt << endl;

    vector<TreeNode *> trees;

    for (int i = first_file_index+1; i < argc; i++)
    {
        cout << PROG_NAME ": Reading " << argv[i] << " ..." << endl;
        ifstream nexus(argv[i]);
        if(!nexus) {
            DIE("Couldn't read tree file.");
        }

        NodeNameMapping leafs, *ptr;
        ReadNexTranslate(nexus, &leafs);
        ptr = (leafs.empty()?0:&leafs);
        if (!ptr)
        {
            WARN("No translate table found!");
        }

        // write what we found
        cout << PROG_NAME << ": Found " << leafs.size() 
             << " mapping entries in " << argv[i] << endl;
       
        // reset the file
        nexus.seekg(0, ios::beg);

        // read the tree collection
        ReadNexTrees(nexus, ptr, trees, burnin_opt);
    } 

    cout << PROG_NAME ": Read " << trees.size() << " trees total." << endl;

    // find all the splits
    SplitDatabase splits;
    AllSplits(trees, splits);
    cout << PROG_NAME ": Extracted " << splits.size() << " splits." << endl;

    // Remove the splits that don't occur very often
    CullSplits(splits, (int)(trees.size()*cull_opt)); 

    cout << PROG_NAME ": Found " << splits.size() << " splits total." << endl;

    string tmp;
    tmp = basename + "_splits";
    ofstream outsplits(tmp.c_str());
    PrintSplitsMapping(outsplits, splits);
    outsplits.close();

    tmp = basename + "_trees";
    ofstream outtrees(tmp.c_str());
    PrintTreesForSplits(outtrees, trees.size(), splits);

    if (dist_opt > 0)
    {
        tmp = basename + "_dist";
        ofstream outdist(tmp.c_str());
        cout << PROG_NAME ": Computing matrix:";
        vector<DistanceMatrix> * matrices = AllDistanceMatrices(trees);
        cout << endl;
        cout << PROG_NAME ": Finished computing distance matrices..." << endl;
        PrintDistances(outdist, *matrices);
        delete matrices;
        outdist.close();
    }

    // delete all the trees
    for_each(trees.begin(), trees.end(), DeleteTree);

    return 0; 
}

