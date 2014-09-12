#include <iostream>
#include "gamma-prob.c"

using namespace std;

#define LN2 0.6931471805599453094172321214581765680755

int main(int argc, char** argv) {

  double x = atof(argv[1]);

  double pval = -LN2 + ln_gamma_prob(0.5, x*x/2.0);
  cout << pval << endl;

  return 0;
}

