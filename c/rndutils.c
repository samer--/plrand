/*
 * Pseudorandom generators for SWI Prolog
 * Samer Abdallah (2009)
 * Some samplers adapted from code in Tom Minka's lightspeed Matlab toolbox.
 * Other bits of code culled from SWI Prolog distribution random.c
 *
 * To do:
 *    work out lower bound for skew stable distributions
 * 	stream splitting with jumps
 * 	reduce blob copying
 * 	fast discrete distributions
 * 	ziggurat for normals
 */

#define _USE_MATH_DEFINES 1

#include <math.h>
#include <float.h>
#include "rndutils.h"

#if HAVE_GSL
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>

void initialise() { gsl_set_error_handler_off(); }
#else
double gs_ran_exponential_pdf(double x, double a) { return NAN; }
double gs_ran_gaussian_pdf(double x, double a) { return NAN; }
double gs_ran_beta_pdf(double x, double a, double b) { return NAN; }
double gs_ran_dirichlet_pdf(size_t n, double *a, double *x) { return NAN; }
double gs_ran_dirichlet_lnpdf(size_t n, double *a, double *x) { return NAN; }
double gs_ran_gamma_pdf(double x, double a, double b) { return NAN; }
double gs_ran_poisson_pdf(unsigned int x, double a) { return NAN; }
double gs_ran_binomial_pdf(unsigned int x, double a, unsigned int b) { return NAN; }
double gsl_sf_zeta(const double s) { return NAN; }
double gsl_sf_hzeta(const double s, const double k) { return NAN; }
double gsl_sf_psi(double x) { return NAN; }

void initialise() {}
#endif

double pdf_Normal(double x)                        { return gsl_ran_gaussian_pdf(x,1); }
double pdf_Beta(double a, double b, double x)      { return gsl_ran_beta_pdf(x,a,b); }
double pdf_Gamma(double a, double x)               { return gsl_ran_gamma_pdf(x,a,1); }
double pdf_Exponential(double x)                   { return gsl_ran_exponential_pdf(x,1); }
double pdf_Zeta(double s, long x)                  { return pow(x,-s)/gsl_sf_zeta(s); }
double pdf_Poisson(double m, long x)               { return gsl_ran_poisson_pdf(x,m); }
double pdf_Binomial(double p, long n, long x)      { return gsl_ran_binomial_pdf(x,p,n); }
double pdf_Dirichlet(long n, double *a, double *x) { return gsl_ran_dirichlet_pdf(n, a, x); }
double pdf_Discrete(long n, double *p, double tot, long x) { return p[x]/tot; }
double pdf_Uniform(double x)                       { return 1.0; }

double logpdf_Dirichlet(long n, double *a, double *x) { return gsl_ran_dirichlet_lnpdf(n, a, x); }
double logpdf_Dirichlet_log(long n, double *a, double *lx) {
	double s=0;
	long i;
	for (i=0; i<n; i++) { s += (a[i]-1)*lx[i]; }
	return s - vector_betaln(n,a);
}

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif


// -----------------------------------------------------

double psi(double x) { return gsl_sf_psi(x); }
double lngamma(double x) { return gsl_sf_lngamma(x); }

/* Computes psi (digamma) vector function */
int vector_psi(long n, double *a, double *x)
{
	long i;
	double tot=0, psi_tot;
	for (i=0; i<n; i++) { tot += a[i]; }
	psi_tot = gsl_sf_psi(tot);
	for (i=0; i<n; i++) { x[i] = gsl_sf_psi(a[i]) - psi_tot; }
	return 1;
}

/* Computes log of vector Beta function */
double vector_betaln(long n, double *a)
{
	long i;
	double sum_a=0, sum_ga=0;
	for (i=0; i<n; i++) { sum_a += a[i]; sum_ga += gsl_sf_lngamma(a[i]); }
	return sum_ga - gsl_sf_lngamma(sum_a);
}

// Returns a sample from Normal(0,1)
double Normal(RndState *S)
{
	double x=S->prev_normal;

	if(!isnan(x)) {
		S->prev_normal=NAN;
	} else {
		double y,radius;

		/* Generate a random point inside the unit circle */
		do {
			x = 2*Uniform(S)-1;
			y = 2*Uniform(S)-1;
			radius = (x*x)+(y*y);
		} while((radius >= 1.0) || (radius == 0.0));
		/* Box-Muller formula */
		radius = sqrt(-2*log(radius)/radius);
		x *= radius;
		y *= radius;
		S->prev_normal=y;
	}
	return x;
}


/* Returns a sample from Gamma(a, 1).
 * For Gamma(a,b), scale the result by b.
 */
double Gamma(RndState *S, double a)
{
  /* Algorithm:
   * G. Marsaglia and W.W. Tsang, A simple method for generating gamma
   * variables, ACM Transactions on Mathematical Software, Vol. 26, No. 3,
   * Pages 363-372, September, 2000.
   * http://portal.acm.org/citation.cfm?id=358414
   */
  double boost, d, c, v;
  if(a < 1) {
    /* boost using Marsaglia's (1961) method: gam(a) = gam(a+1)*U^(1/a) */
    boost = exp(log(Uniform(S))/a);
    a++;
  } 
  else boost = 1;
  d = a-1.0/3; c = 1.0/sqrt(9*d);
  while(1) {
    double x,u;
    do {
      x = Normal(S);
      v = 1+c*x;
    } while(v <= 0);
    v = v*v*v;
    x = x*x;
    u = Uniform(S);
    if((u < 1-.0331*x*x) || 
       (log(u) < 0.5*x + d*(1-v+log(v)))) break;
  }
  return( boost*d*v );
}

/* Returns a sample from Dirichlet(n,a) where n is dimensionality */
int Dirichlet(RndState *S, long n, double *a, double *x)
{
	long i;
	double tot=0, z;
	for (i=0; i<n; i++) { z=Gamma(S,a[i]); tot+=z; x[i]=z; }
	for (i=0; i<n; i++) { x[i]/=tot; }
	return 1;
}

/* Returns a sample from Beta(a,b) */
double Beta(RndState *S, double a, double b)
{
  double g = Gamma(S,a);
  return g/(g + Gamma(S,b));
}

double Hyperbolic(RndState *S, double a)
{
	return floor(pow(1-Uniform(S),-1/a));
}

/* Sample from zeta distribution
 *
 * This is an inverse CDF method. It works by using  an 
 * approximation of the Hurwitz Zeta function and then
 * walking left or right until we get to the correct point.
 * The approximation is good for large values so usually we
 * only have to walk left or right a few steps for small 
 * values.
 */
long Zeta(RndState *S, double s)
{
	double v=(1-Uniform(S))*gsl_sf_zeta(s);
	double k0, k=round(pow(v*(s-1),1/(1-s))); // approx inverse of hzeta
	double zk0, zk=gsl_sf_hzeta(s,k);

	if (zk>=v) {
		do {
			k0=k; zk0=zk;
			k=k0+1; zk=zk0-pow(k0,-s);
		} while (zk>=v && zk!=zk0);
		return (long)k0;
	} else {
		do {
			k0=k; zk0=zk;
			k=k0-1; zk=zk0+pow(k,-s);
		} while (zk<v && zk!=zk0);
		return (long)k;
	}
}

/* Very fast binomial sampler. 
 * Returns the number of successes out of n trials, with success probability p.
 */
int Binomial(RndState *S, double p, int n)
{
  int r = 0;
  if(isnan(p)) return 0;
  if(p < DBL_EPSILON) return 0;
  if(p >= 1-DBL_EPSILON) return n;
  if((p > 0.5) && (n < 15)) {
    /* Coin flip method. This takes O(n) time. */
    int i;
    for(i=0;i<n;i++) {
      if(Uniform(S) < p) r++;
    }
    return r;
  }
  if(n*p < 10) {
    /* Waiting time method.  This takes O(np) time. */
    double q = -log(1-p), e = -log(Uniform(S)), s;
    r = n;
    for(s = e/r; s <= q; s += e/r) {
      r--;
      if(r == 0) break;
      e = -log(Uniform(S));
    }
    r = n-r;
    return r;
  }
  if (1) {
    /* Recursive method.  This makes O(log(log(n))) recursive calls. */
    int i = (int)floor(p*(n+1));
    double b = Beta(S,i, n+1-i);
    if(b <= p) r = i + Binomial(S,(p-b)/(1-b), n-i);
    else r = i - 1 - Binomial(S,(b-p)/b, i-1);
    return r;
  }
}

double Exponential(RndState *S) { return -log(1-Uniform(S)); }

long Poisson(RndState *S, double lambda)
{
	long r;
	if (lambda>=15) {
		double m=floor(lambda*7/8);
		double x=Gamma(S,m);

		if (x>lambda) r=Binomial(S,lambda/x,m-1);
		else          r=m+Poisson(S,lambda-x);
	} else {
		double p, elambda = exp(-lambda);
		for (p=1, r=-1; p>elambda; p*=Uniform(S)) r++; 
	}
	return r;
}

/* Returns a sample from stable distribution */
double Stable(RndState *S, int param, double alpha, double beta)
{
	double X, theta=M_PI*(Uniform(S)-0.5);

	// These methods come from John Nolan's book about
	// stable distributions. They return standardised stable random 
	// numbers in the S(alpha,beta;0) parameterisation

	if (beta==0) {
		if (alpha==1) {
			X=tan(theta);
		} else {
			double W=Exponential(S);
			X=sin(alpha*theta)/cos(theta)*pow(cos((alpha-1)*theta)/(W*cos(theta)),1/alpha-1);
		}
	} else {
		double W=Exponential(S);
		double zeta=beta*tan(alpha*M_PI/2);
		if (alpha==1) {
			double b=2*beta/M_PI;
			double bt=1+b*theta;
			X=bt*tan(theta) - b*log(W*cos(theta)/bt);
		} else {
			double ath = alpha*theta; 
			double a1th = theta-ath; 
			X = ((sin(ath) + zeta*cos(ath))/cos(theta))
				 * pow((cos(a1th) + zeta*sin(a1th))/(W*cos(theta)),1/alpha-1);
			if (param==0) X=X-zeta;
		}
	}
	return X;
}

/* Returns a sample from discrete distribution */
long Discrete(RndState *S, long n, double *p, double tot)
{
	double u;
	long i;

	for (i=0, u=tot*Uniform(S)-p[0]; u>=0; u-=p[++i]); 
	return i;
}

double Uniform(RndState *S0)
{
	return RngStream_Double(&S0->g,&S0->g);
}

double Raw(RndState *S0)
{
	return RngStream_Raw(&S0->g,&S0->g);
}

void InitRndState(RndState *S)
{
	RngStream_InitState(&S->g);
	S->prev_normal=NAN;
}


int spawn_gen(RndState *S0, RndState *S1)
{
	unsigned seeds[6];
	int i;

	for (i=0; i<6; i++) seeds[i]=(unsigned)RngStream_Raw(&S0->g,&S0->g);
	RngStream_SetState(&S1->g,seeds);
	for (i=0; i<3; i++) RngStream_Raw(&S1->g,&S1->g);
	S1->prev_normal=NAN;
	return 1;
}

int randomise_from(FILE *p, RndState *S)
{
	unsigned seeds[6];
	size_t r=fread((void *)seeds,sizeof(unsigned),6,p);

	if (r==6) {
		int i;
		if (RngStream_SetState(&S->g,seeds)<0) return 0;
		for (i=0; i<3; i++) RngStream_Raw(&S->g,&S->g);
		S->prev_normal=NAN;
		return 1;
	} else return 0;
}

