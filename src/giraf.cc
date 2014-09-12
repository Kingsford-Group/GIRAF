#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "util.h"
#include "catalog.h"
#include "timer.h"
#include "options.h"


using namespace std;

#define PROG_NAME "giraf"

// the main functions we call
int main_mcmc_split_info(int, char**);
int main_build_incompat_graph(int, char**);
vector<string> main_extract_reassortments(int, char**);

// the global options and data
vector<string> names;       // holds the tree names
vector<vector<string>*> tree_files; // holds the tree files
vector<string> options; // holds all other options

string catalog_filename = "catalog";
int arch_threshold = 0;

const int MAX_CMD_LINE = 4048;

const char *GIRAF_OPTIONS = "h";

enum {ARCH_THRESH_OPT=1, CATFILE_OPT};

static struct option MAYBE_UNUSED giraf_long_options[] = {
    {"arch-threshold", 1, 0, ARCH_THRESH_OPT},
    {"out-catalog", 1, 0, CATFILE_OPT},
    {0,0,0,0}
};


// create a mock commandline array
char **
MakeCmdArray(const vector<string> & options, char * buffer)
{
    char ** argv = new char*[options.size()];

    // create commandline buffer
    char * b = buffer;
    for(unsigned i = 0; i < options.size(); ++i) 
    {
        DIE_IF(MAX_CMD_LINE-(b-buffer) < (int)options[i].length()+1, 
            "Command line is too long");
        strncpy(b, options[i].c_str(), MAX_CMD_LINE-(b-buffer)-1);
        argv[i] = b;
        b += options[i].length()+1;
    }
    return argv;
}


void
ParseTreeFile(const string & treefile)
{
    ifstream in(treefile.c_str());
    string line;
    vector<string> words;
    while(getline(in, line)) 
    {
        line = Trim(line);
        if(line.length() == 0 || line[0]=='#') continue;
        
        // split line into words
        istringstream iss(line);
        words.clear();
        string w;
        while(iss >> w) 
        {
            words.push_back(w);
        }
        assert(iss.eof());
        DIE_IF(words.size() < 2, "Each line in " + treefile + 
           " must specify tree name and at least one tree file");

        names.push_back(words[0]);
        tree_files.push_back(new vector<string>());
        copy(words.begin()+1, words.end(), back_inserter(*tree_files.back()));

        // provide some output to let user know what is happening
        cout << PROG_NAME << ": tree name = " << words[0] << ";"
             << " tree files =";
        for(vector<string>::iterator I = tree_files.back()->begin();
            I != tree_files.back()->end();
            ++I)
        {
            cout << " " << *I;
        }
        cout << endl;
    }
}

//=======================================================================
// Command-line options
//=======================================================================

static
void
Usage()
{
    void SplitUsage(bool);
    void GraphUsage(bool);
    void ExtractUsage(bool);

    cerr << "Usage: " << PROG_NAME << " [options] in.giraf" << endl << endl;

    cerr << "OPTIONS" << endl
         << "   --out-catalog=filename : save catalog here (default \"catalog\")" << endl
         << "   --arch-threshold=N     : require N segment pairs to support a reassortment"  << endl
         << "                              (default max{3, #seg-2})" << endl 
         << endl; 
    SplitUsage(false);
    GraphUsage(false);
    ExtractUsage(false);
    
    exit(3);
}

static
int
ProcessGirafOptions(int argc, char *argv[])
{
    opterr = 0;
    // GNU ONLY!!! By setting optind to 0, we force getopt to reinitialize
    // this is required if we want giraf to call this program
    int a;
    while ((a=getopt_long(argc, argv, GIRAF_OPTIONS, giraf_long_options, 0)) != -1)
    {
        switch(a)
        {
            case 'h': Usage(); break;
            case ARCH_THRESH_OPT: 
                arch_threshold = atoi(optarg); 
                DIE_IF(arch_threshold < 0, "Argument to --arch-treshold must be >= 0"); 
                break;

            case CATFILE_OPT: catalog_filename = optarg; break;
            /*default:
                cerr << "Unknown option." << endl;
                Usage(); */
            }
    }
    if(optind >= argc) Usage();
    return optind;
}


static
void
ParseOptions(int argc, char * argv[])
{
    // do our own argument processing to break up the command into parts
    string treefile = "";
    string s;
    for(int i=1; i<argc; i++) 
    {
        s = argv[i];
        if(s.length() < 1) continue; // shouldn't happen

        // if there is a naked option except for the last, die
        if(s[0] != '-' && i != argc-1) 
        {
            cerr << "Unrecongnized option " << s << endl;
            DIE("You must use --opt=X syntax instead of --opt X");
        } 
        else if(s[0] != '-')
        {
            treefile = s;
        } 
        else 
        {
            options.push_back(s);
        }
    }
    if(treefile == "") {
        cerr << PROG_NAME 
             << ": Must specify file listing segments and their associated tree files" 
             << endl;
        Usage();
    }

    cerr << PROG_NAME << ": Reading " << treefile << " for tree files." << endl;

    ParseTreeFile(treefile);
    DIE_IF(names.size() < 2, "At least 2 segments must be present in " + treefile);

    // check if any options apply to us
    ProcessGirafOptions(argc, argv);
}


int
main(int argc, char *argv[])
{
    char cmdline[MAX_CMD_LINE];
    char **newargv;
    vector<string> vec;

    Timer T;
    cout << PROG_NAME << ": " << T.start();

    ParseOptions(argc, argv);

    // create the left mcmc_split_info command line
    for(unsigned i=0; i<names.size(); i++) 
    {
        vec.push_back("mcmc_split_info");
        vec.push_back("--ignore-bad-options");
        copy(options.begin(), options.end(), back_inserter(vec));    
        vec.push_back(names[i]);
        copy(tree_files[i]->begin(), tree_files[i]->end(), back_inserter(vec));

        // run the left mcmc_split_info command
        newargv = MakeCmdArray(vec, cmdline);
        main_mcmc_split_info(vec.size(), newargv);
        vec.clear();
    }

    ReassortDB allreassort;

    string pair_name;
    // for every pair of segments, run the comparison
    for(unsigned seg1 = 0; seg1 < names.size(); ++seg1)
    {
        for(unsigned seg2 = seg1+1; seg2 < names.size(); ++seg2)
        {
            vector<string> found;

            cout << PROG_NAME << ": Processing " 
                 << names[seg1] <<  " " << names[seg2] << endl;
            pair_name = names[seg1] + "-" + names[seg2];

            // create the build_incomp_graph command line
            vec.push_back("build_incompat_graph");
            vec.push_back("--ignore-bad-options");
            copy(options.begin(), options.end(), back_inserter(vec));    
            vec.push_back(names[seg1]);
            vec.push_back(names[seg2]);
            vec.push_back(pair_name);

            // run the build_incompat_graph command
            newargv = MakeCmdArray(vec, cmdline);
            //for(int j = 0; j < 100; j++) cout << cmdline[j] << endl;
            main_build_incompat_graph(vec.size(), newargv);
            vec.clear();
            
            // create the extract_reassortments command line
            vec.push_back("extract_reassortments");
            vec.push_back("--ignore-bad-options");
            copy(options.begin(), options.end(), back_inserter(vec));    
            vec.push_back(names[seg1]);
            vec.push_back(names[seg2]);
            vec.push_back(pair_name);
            vec.push_back(pair_name);

            // run the extract_reassortments command
            newargv = MakeCmdArray(vec, cmdline);
            found = main_extract_reassortments(vec.size(), newargv);
            for(vector<string>::iterator F = found.begin();
                F != found.end();
                ++F)
            {
                allreassort[*F].push_back(make_pair(names[seg1], names[seg2]));
            }
            vec.clear();
        }
    }
    // allreassort = map from reassort to vector of evidence

    // catalog reassortments if more than 2 segments given
    if(names.size() >= 3) {
        cout << PROG_NAME << ": Writing architectures to file: " 
             << catalog_filename << endl;
        ofstream catout(catalog_filename.c_str());
        DIE_IF(!catout, PROG_NAME ": Error opening catalog file.");
        int catthreshold = (arch_threshold<=0) ? (min(3, (int)names.size()-2)) : arch_threshold;
        cout << PROG_NAME << ": Pair threshold = " << catthreshold << endl;
        CatalogReassortments(allreassort, catout, catthreshold);
    }
    cout << PROG_NAME << ": " << T.stop();
    cout << PROG_NAME << ": " << T.report() << endl;
    cout << PROG_NAME << ": done." << endl;
}
