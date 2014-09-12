#ifndef SPLITS_H
#define SPLITS_H
#include "tree.h"
#include <set>
#include <map>
#include <string>


struct Split
{
    int id;
    int x,y;
public:
    Split(const set<string> & f, const set<string> & s);
    Split(const set<string> & f, const set<string> & s, int index);

    const set<string> & first() const { return _first; }
    const set<string> & second() const { return _second; }

    const set<string> & smaller() const {
        return (first().size() < second().size()) ? first() : second();
    }
    
    string as_string() const { return _rep; }

    bool operator<(const Split & s) const { return _rep < s._rep; }


private:
    set<string> _first;
    set<string> _second;

    string _rep;
};

typedef map<Split, set<int> > SplitDatabase;

void AllSplits(vector<TreeNode *> &, SplitDatabase &);
void CullSplits(SplitDatabase &, unsigned); 
bool SplitsAreIncompatible(const Split &, const Split &);

void SetDifference(const set<string> &, const set<string> &, set<string> &);
void SetIntersection(const set<string> &, const set<string> &, set<string> &);


// split printing 
void PrintSplitsReadable(ostream & , SplitDatabase &, unsigned = 1);
void PrintSplitsMapping(ostream &, SplitDatabase &);
void ReadSplitsMapping(istream &, SplitDatabase &);
void PrintTreesForSplits(ostream & , int , SplitDatabase & );
#endif
