#include <sstream>
#include <vector>
#include <cmath>
#include "dist.h"
#include "util.h"
#include "tree.h"


double
Average(const vector<double> & vec)
{
    double sum = 0.0;
    for(vector<double>::const_iterator i = vec.begin();
        i != vec.end();
        ++i)
    {
        sum += *i;
    }
    return sum / vec.size();
}


double
StdDev(const vector<double> & vec)
{
    double avg = Average(vec);

    double sum = 0.0;
    for(vector<double>::const_iterator i = vec.begin();
        i != vec.end();
        ++i)
    {
        sum += (*i - avg) * (*i - avg);
    }
    if (sum == 0.0) sum = 1; // XXX: check this
    return sqrt(sum / (vec.size()-1));
}


extern "C" { double ln_gamma_prob(double, double); }

double
normal_invcdf(double x)
{
    #define LN2 0.6931471805599453094172321214581765680755
    return -LN2 + ln_gamma_prob(0.5, x*x/2.0);
}

void
DebugPrintVector(const vector<double> & vec)
{
    cout.precision(13);
    for(vector<double>::const_iterator i = vec.begin();
        i != vec.end();
        ++i)
    {
        cout << *i << " ";
    }
}

double
ComputePValue(
    const vector<double> & dist1, 
    const vector<double> & dist2,
    bool & is_greater,
    bool asymetric
    )
{
    double z_score;

    // asymetric=true is for is for perl compatibility
    if(!asymetric) 
    {
        z_score = (Average(dist1)-Average(dist2)) / 
                   max(StdDev(dist1),StdDev(dist2));
    }
    else
    {
        z_score = (Average(dist1)-Average(dist2)) / StdDev(dist2);
    }
    double corr_z_score = 2.0 * sqrt(z_score*z_score); 

    is_greater = (z_score > 0);

    return normal_invcdf(corr_z_score); 
}


// parse a line formated as "T1 T2 dist1 dist2 ..."
void
ParseDistLine(
    const string & line, 
    string & taxon1, 
    string & taxon2,
    vector<double> & distvec
    )
{
    istringstream iss(line);
    iss >> taxon1 >> taxon2;
    distvec.clear();

    double tmp;
    while(iss >> tmp)
    {
        distvec.push_back(tmp);
    }
    assert(iss.eof()); 
}


// scan in1 and in2 _dist files and produce a pair_test_results file
void
ComputePairDistances(
    istream & in1,
    istream & in2,
    bool asymetric, // FALSE for normal; TRUE for perl compat
    ostream * out, // 0 if no output file
    DistanceMatrix & D,   // out
    DistanceMatrix & G    // out
    )
{
    string line1, line2;
    string taxon1, taxon2;
    string check1, check2;
    vector<double> distvec1, distvec2;

    long count = 0;
    unsigned long linenum = 0;
    unsigned long size1=0, size2=0;

    while(getline(in1, line1))
    {
        getline(in2, line2);
        linenum++;

        ParseDistLine(line1, taxon1, taxon2, distvec1);
        ParseDistLine(line2, check1, check2, distvec2);

        assert(taxon1 == check1 && taxon2 == check2);
        if(linenum == 1)
        {
            // save the size to check later
            size1 = distvec1.size();
            size2 = distvec2.size();

            // the two files can have different length vectors, but 
            // this is suspicious
            if (size1 != size2)
            {
                cout << "warning: Different # of trees were sampled for the two segments." << endl;
            }
        }
        else
        {
            // make sure vectors in each file are the same size
            DIE_IF(distvec1.size() != size1 || distvec2.size() != size2,
                "Vectors in _dist files are not all the same length!");
        }

        bool is_greater;
        double log_pvalue = ComputePValue(distvec1, distvec2, is_greater, asymetric);
        if (out)
        {
            (*out) << taxon1 << " " << taxon2 << " " << (is_greater?1:0) 
                   << " " << log_pvalue << endl;
                   //<< " " << StdDev(dist1) << " " << StdDev(dist2) << endl;
        }

        D[taxon1][taxon2] = log_pvalue;
        G[taxon1][taxon2] = is_greater;

        if(count++ % 1000 == 0) cout << "." << flush;
    }
}
