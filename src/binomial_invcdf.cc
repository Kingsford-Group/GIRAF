#include <iostream>
#include "gamma-prob.c"

using namespace std;

int main(int argc, char** argv) {

  int count = atoi(argv[1]);
  int total = atoi(argv[2]);
  double p = atof(argv[3]);

  double pval = betai(count+1, total-count, p);
  cout << pval << endl;

  return 0;
}

