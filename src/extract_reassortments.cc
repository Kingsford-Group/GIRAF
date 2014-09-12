#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "label_types.h" 
#include "bigraph.h"
#include "scored_set.h"
int non_star = 1;
#include "mica.h"
#include "timer.h"
#include "options.h"

#define PROG_NAME "extract_reassortments"

using namespace std;

double threshold_opt = 0.70;
bool require_multiple_cliques = true;

// Options for extract_reassortments

const char *EXTRACT_OPTIONS = "ht:";

enum {TYPE_OPT=1, THRESHOLD_OPT, SINGLE_OPT, EXTRACT_BAD_OPT};

static struct option MAYBE_UNUSED extract_long_options[] = {
    {"non-star", 1, 0, TYPE_OPT},
    {"threshold", 1, 0, THRESHOLD_OPT},
    {"single", 0, 0, SINGLE_OPT},
    {"ignore-bad-options", 0, 0, EXTRACT_BAD_OPT},
    {0,0,0,0}
};

void
ExtractUsage(bool show_cmd=true)
{
    if(show_cmd)
    {
        cerr << PROG_NAME << " [options] base1 base2 jointbase outbase"
             << endl 
             << endl
             << "OPTIONS" << endl;;
    }
    cerr << "   --non-star=[0,1] : If 1, only consider non-star bicliques (default)" << endl
         //<< "          2 for just edges" << endl
         << "   --threshold=F    : confidence cutoff (default F=0.50)" << endl
         << "   --single         : allow single cliques" << endl;
    if(show_cmd) exit(3);
}


int
ProcessOptions(int argc, char *argv[])
{
    bool ignore_bad_opt = false;
    // these are required to resuse getopt
    opterr = 0;
    // GNU ONLY!!! By setting optind to 0, we force getopt to reinitialize
    // this is required if we want giraf to call this program
    optind = 0; 
    int a;
    while ((a=getopt_long(argc, argv, EXTRACT_OPTIONS, extract_long_options, 0)) != -1)
    {
        switch(a)
        {
            case 'h': ExtractUsage(); break;
            case TYPE_OPT: non_star = atoi(optarg); break;
            case 't': case THRESHOLD_OPT: threshold_opt = sqrt(atof(optarg)); break;
            case SINGLE_OPT: require_multiple_cliques = false; break;
            case EXTRACT_BAD_OPT: ignore_bad_opt = true; opterr = 0; break;
            default:
                if(!ignore_bad_opt) {
                    cerr << "Unknown option " << optopt << endl;
                    ExtractUsage();
                }
        }
    }
    if (optind+3 >= argc) ExtractUsage();
    return optind;
}


double
GetConfidenceScore(
    const string & results,
    int & num_cliques
    )
{
    istringstream iss(results);

    double pval = 1.0;
    num_cliques = 0;

    string line;
    while(iss)
    {
        getline(iss, line);
        getline(iss, line);
        getline(iss, line);

        if(iss) {
            istringstream liness(line);
            double pa, pb;
            liness >> pa >> pb;

            pval *= (1-pa*pb);
            num_cliques++;

            getline(iss, line);
        }
    }
    return 1.0 - pval;
}


void
ReadLabelMapping(
    istream & in,
    map<edge_label_t, string> & Labels
    )
{
    string line;
    while(getline(in, line))
    {
        if (line.length() > 0) {
            unsigned space = line.find(' ');
            edge_label_t id = atoi(line.substr(0, space).c_str());
            string taxa = line.substr(space+1);
            Labels[id] = taxa;
        }
    }
}


void
initialize_global_variables()
{
    num_of_left_trees = 0;
    num_of_right_trees = 0;
    left_trees.clear();
    right_trees.clear();
}


vector<string>
main_extract_reassortments(int argc, char** argv)
{
    // parse the command line
    cout << PROG_NAME ": built on " << __DATE__ << endl;
    initialize_global_variables();
    int base_index = ProcessOptions(argc, argv);

    string base1 = argv[base_index];
    string base2 = argv[base_index+1];
    string jointbase = argv[base_index+2];
    string outbase = argv[base_index+3];

    // list the values of the options
    cout << PROG_NAME ": Biclique Type = " << non_star << endl;
    cout << PROG_NAME ": SquareRoot(Threshold) = " << threshold_opt << endl;
    cout << PROG_NAME ": MultipleCliques = " << require_multiple_cliques << endl;

    // Open the output files
    ofstream results((outbase + "_bicliques").c_str());
    if(!results) 
    {
        cerr << PROG_NAME ": Error in opening results file." << endl;
        exit(3);
    }

    ofstream report((outbase + "_report").c_str());
    if(!report)
    {
        cerr << PROG_NAME ": Error in opening report file." << endl;
        exit(3);
    }
    report << "<<GiRaF Report>>" << endl << endl;
    report.precision(15);

    // read the trees & graph
    Labelled_Bigraphs<left_label_t, right_label_t, edge_label_t> graphs;

    read_trees((base1 + "_trees").c_str(), num_of_left_trees, left_trees);
    read_trees((base2 + "_trees").c_str(), num_of_right_trees, right_trees);
    graphs.read((jointbase + "_graph.labelled").c_str());

    // read a mapping between labels and sets
    map<edge_label_t, string> Labels;
    ifstream label_map_file((jointbase + "_graph.labels").c_str());
    ReadLabelMapping(label_map_file, Labels);
    label_map_file.close();

    // for every label, find the maximal bicliques

    // tracks the sets of taxa we output
    vector<string> found_sets;

    for(set<edge_label_t>::iterator it = graphs.getEdgeTypes().begin();
        it != graphs.getEdgeTypes().end(); it++) {

        ostringstream output;

        if(non_star == 2)
        {
            good_edges<left_label_t, right_label_t, tree_label_t>(
                graphs.getGraph(*it), threshold_opt, results
            );
        }
        else
        {
            mica<left_label_t, right_label_t, tree_label_t>(
                graphs.getGraph(*it), threshold_opt, output, cerr
            );
        }

        if(!output.str().empty()) 
        {
            results << "Label: " << *it << endl << endl << output.str() << endl;
            results << "------------" << endl;

            int num_cliques;
            double conf = GetConfidenceScore(output.str(), num_cliques);

            if ((!require_multiple_cliques) || num_cliques > 1)
            {
                report << "Candidate = ID: " << *it 
                       << ", Conf: " << conf
                       << ", Taxa: {" << Labels[*it] << "}"
                       << endl;
                found_sets.push_back(Labels[*it]);
            }
        }
    }

    results.close();
    report.close();

    return found_sets;
}

