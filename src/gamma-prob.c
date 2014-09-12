#ifndef GAMMA_PROB_C
#define GAMMA_PROB_C

/* Copyright 1991, 1992, 1993, 1994, 1998 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   Dept. of Molecular, Cellular, and Developmental Biology
 *   University of Colorado
 *   Campus Box 347
 *   Boulder, CO  80309-0347
 *
 *   hertz@colorado.edu
 */

/* These function determine from the log-likelihood ratio the natural
 * logarithm of the gamma/chi square probabilities.

 * The functions' parameters are:
 * the alpha parameter or the degrees of freedom,
 * negative log-likelihood ratio.

 * Also included are functions for calculating the complete and incomplete psi
 * functions, and for summing large numbers by working with their logarithms.
 */

#include <stdio.h>
#include <math.h>
#define S_ERR 1.0e-12       /* A smaller Error for numerical calculations. */
#define L_ERR 1.0e-6        /* A larger Error for numerical calculations. */

/* Natural logarithm of gamma function for values > 0. */
double ln_gamma(double alpha) {

  /* Magic coefficients for calculating the gamma function. */
  static double coef[6] = {76.18009173, -86.50532033, 24.01409822,
			     -1.231739516, 0.120858003e-2, -0.536382e-5};
  double ln_factor;
  double series;
  double a;
  int i;

  /* Use the magic coefficients to calculate ln(GAMMA). */
  if (alpha > 1.0)
    {
      ln_factor = alpha + 4.5;
      ln_factor = (alpha - 0.5) * log(ln_factor) - (ln_factor);

      for (i = 0, a = alpha, series = 1.0; i < 6; ++i, a = a + 1.0)
	series = series + coef[i] / a;

      return(ln_factor + log(2.50662827465 * series));
    }

  /* Use the relationship of GAMMA(z) = GAMMA(z+1) / z since
   * the magic coefficients work less well between 0 and 1. */
  else if (alpha > 0.0) return(ln_gamma(alpha + 1.0) - log(alpha));

  /* Magic coefficients do not work with values <= 0. */
  else
    {
      fprintf(stderr,
	      "The function \"ln_gamma()\" only works with values > 0\n");
      return -1;
    }
}

/* Determine the natural logarithm of the incomplete gamma function's
 * complement using the continued fraction approximation. */
double ln_gamma_cf(double alpha, double ln_likelihood) {

  double frac_old = 0.0;        /* Old continued fraction. */

  /* Variables A and B for calculating the continued fraction. */
  double a_0 = 0.0, a_1 = 1.0;
  double b_0 = 1.0, b_1 = ln_likelihood;

  /* Variables for calculating A and B. */
  double i;
  double i_alpha;     /* i - alpha. */

  /* Do 2 cycles of the reiteration for each index. */
  for (i = 1.0; ; i = i + 1.0)
    {
      /* Cycle one of reiteration. */
      i_alpha = i - alpha;
      a_0 = a_1 + i_alpha * a_0;
      b_0 = b_1 + i_alpha * b_0;

      /* Cycle two of the reiteration. */
      a_1 = ln_likelihood * a_0 + i * a_1;
      b_1 = ln_likelihood * b_0 + i * b_1;

      /* Determine the continued fraction if the denominator is not 0. */
      if (b_1)
	{
	  a_0 = a_0 / b_1;
	  a_1 = a_1 / b_1;
	  b_0 = b_0 / b_1;
	  b_1 = 1.0;
	  if (fabs((a_1 - frac_old)/a_1) < L_ERR)
	    return(-ln_likelihood + alpha*log(ln_likelihood) + log(a_1)
		   - ln_gamma(alpha));
	  frac_old = a_1;
	}
    }
}

/* Determine the incomplete gamma function using the series approximation. */
double gamma_ser(double alpha, double ln_likelihood) {
 
  double sum;                   /* The sum of the series. */
  double alpha_n;               /* "alpha" + "n". */
  double term;                  /* Current term in the series. */
  
  for (alpha_n = alpha + 1.0, sum = term = 1.0/alpha;
       term >= sum * L_ERR; alpha_n = alpha_n + 1.0)
    {
      term = term * ln_likelihood / alpha_n;
      sum = sum + term;
    }
  return(exp(-ln_likelihood + alpha * log(ln_likelihood)
	     + log(sum) - ln_gamma(alpha)));
}

/* Returns the natural logarithm of the gamma probability,
 * integrated from ln_like to infinity, assuming beta = 1.*/
double ln_gamma_prob(double alpha, double ln_like) {

  if ((ln_like < -L_ERR) || (alpha < 0.0))
    {
      fprintf(stdout, "PROGRAM BUG:\n");
      fprintf(stdout, "Cannot determine the gamma probability when ");
      fprintf(stdout, "the negative log-likelihood\nratio is less than 0 or ");
      fprintf(stdout, "the alpha parameter is less than 0.\n\n");

      fprintf(stderr, "log-likelihood=%g  alpha=%g\n", ln_like, alpha);
      return -1;
    }

  else if ((ln_like <= 0.0) || (alpha == 0.0)) return(0.0);

  /* Use the continued fraction approximation of the
   * incomplete gamma function's complement. */
  else if (ln_like >= alpha + 1.0) return(ln_gamma_cf(alpha, ln_like));

  /* Use the series approximation of the incomplete gamma function. */
  else if (ln_like > 0.0) return(log(1.0 - gamma_ser(alpha, ln_like)));
  fprintf(stdout, "PROGRAM BUG!");
  return -1;
}


double log_chi_inv_cdf(double alpha, double ln_like) {

  return ln_gamma_prob(alpha/2.0, ln_like/2.0)-ln_gamma(alpha/2.0);
}


/* Returns the natural logarithm of the one-dimensional
 * normal distribution, integrated from x to infinity. */
double integrate_normal(double x, double avg, double var) {

#define LN2 0.6931471805599453094172321214581765680755   /* ln(2) */

  double ln_prob;   /* The ln(integral) */
  double x_avg = x - avg;

  ln_prob = -LN2 + ln_gamma_prob(0.5, (x_avg * x_avg) / (2.0 * var));

  if (x_avg < 0.0) ln_prob = log(1.0 - exp(ln_prob));

  return(ln_prob);
}


#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

float betacf(float a, float b, float x) {

  int m,m2;
  float aa,c,d,del,h,qab,qam,qap;
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;

  if (fabs(d) < FPMIN) d=FPMIN;

  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {

    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;

    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; 

    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;

    d=1.0/d;
    del=d*c;
    h *= del;

    if (fabs(del-1.0) < EPS) break;
  }

  if (m > MAXIT) 
    fprintf(stdout, "a or b too big, or MAXIT too small in betacf");

  return h;
}

float gammln(float xx) {
											  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

float betai(float a, float b, float x) {

  float bt;
  if (x < 0.0 || x > 1.0) 
    fprintf(stdout, "Bad x in routine betai");

  if (x == 0.0 || x == 1.0) bt=0.0;
  else 
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

#endif
