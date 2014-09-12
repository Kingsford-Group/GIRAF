#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>

#include "catalog.h"
#include "util.h"

bool
VectorContains(const vector<string> & vec, const string & item)
{
    return find(vec.begin(), vec.end(), item) != vec.end();
}


bool
GetArchitecture(
    const vector<pair<string,string> >::const_iterator & ebegin, 
    const vector<pair<string,string> >::const_iterator & eend, 
    vector<string> & leftseg, 
    vector<string> & rightseg)
{
    // for every pair in evidence

    vector<string> orig_left;
    vector<string> orig_right;

    bool error = false;
    for(vector<pair<string,string> >::const_iterator E = ebegin;
        E != eend;
        ++E)
    {
        string a,b;
        a = E->first;
        b = E->second;
        //cerr << "a,b=" << a  << " " << b <<endl;
        bool lefta = VectorContains(leftseg, a);
        bool leftb = VectorContains(leftseg, b);
        bool righta = VectorContains(rightseg, a);
        bool rightb = VectorContains(rightseg, b);

        // if both already assigned to left or both assigned to right, return error
        if ((lefta && leftb) || (righta && rightb)) 
        {
            return true;
        } 

        // if neither assigned to left, 

        if (!lefta && !leftb)
        {
            // but one is assigned to right, assign the other to left
            if (rightb)
            {
                leftseg.push_back(a);
            }
            else if(righta)
            {
                leftseg.push_back(b);
            }
            // and neither is assigned to the right 
            else
            {
                orig_left = leftseg;
                orig_right = rightseg;
                // try assigning a to left and b to right
                leftseg.push_back(a);
                rightseg.push_back(b);
                error = GetArchitecture(E+1, eend, leftseg, rightseg);

                if (error)
                {
                    leftseg = orig_left;
                    rightseg = orig_right;
                    leftseg.push_back(b);
                    rightseg.push_back(a);
                    error = GetArchitecture(E+1, eend, leftseg, rightseg);
                }
            }
        }
        else if(lefta && !rightb)
        {
            rightseg.push_back(b);
        }
        else if(leftb && !righta)
        {
            rightseg.push_back(a);
        }
    }

    return error;
}

string
EvidenceAsString(
    const vector<pair<string,string> > & evidence
    )
{
    ostringstream oss;
    for(vector<pair<string,string> >::const_iterator E = evidence.begin();
        E != evidence.end();
        ++E)
    {
        if(E != evidence.begin()) oss << " ";
        oss << E->first << "-v-" << E->second;
    }
    return oss.str();
}

struct ByEvidence {
    ReassortDB & _res;
    ByEvidence(ReassortDB & allres)
        : _res(allres) {} 
    bool operator()(const string & a, const string & b)
    {
        return _res[a].size() > _res[b].size();
    }
};


void
CatalogReassortments(
    ReassortDB & allreassort,
    ostream & out,
    unsigned min_evidence = 3
    )
{
    // create a sorted list of found sets
    set<string> allsets;
    vector<string> found;

    for(ReassortDB::iterator R = allreassort.begin();
        R != allreassort.end();
        ++R)
    {
        allsets.insert(R->first);
    }

    copy(allsets.begin(), allsets.end(), back_inserter(found));
    sort(found.begin(), found.end(), ByEvidence(allreassort));

    // for every set found, sorted by decreasing evidence
    for(vector<string>::iterator F = found.begin();
        F != found.end();
        ++F)
    {
        if(allreassort[*F].size() >= min_evidence)
        {
            bool error;
            vector<string> leftseg;
            vector<string> rightseg;

            // compute the arch
            out << *F << " :  " << EvidenceAsString(allreassort[*F]) << " :  ";
            leftseg.clear();
            rightseg.clear();
            error = GetArchitecture(allreassort[*F].begin(), allreassort[*F].end(), 
                        leftseg, rightseg);
            sort(leftseg.begin(), leftseg.end());
            sort(rightseg.begin(), rightseg.end());

            out << VectorAsString(leftseg) << " |  " 
                << VectorAsString(rightseg) 
                << (error ? "**" : "") << endl;
        }
    }
}

