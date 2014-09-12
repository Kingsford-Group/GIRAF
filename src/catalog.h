#ifndef CATALOG_H
#define CATALOG_H

#include <ostream>

using namespace std;

typedef map<string, vector<pair<string, string> > > ReassortDB;

void
CatalogReassortments( ReassortDB &, ostream &, unsigned);

#endif
