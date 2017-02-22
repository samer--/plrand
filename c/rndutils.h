#include <stdio.h>
#include "RngStream.h"

// random generator state
typedef struct {
	struct rng_state g;
	double prev_normal;
} RndState;

typedef struct {
	struct rng_jump j;
	int e;
} RndJump;

int randomise_from(FILE *p, RndState *S);
int spawn_gen(RndState *S0, RndState *S1);
void InitRndState(RndState *S);

double Raw(RndState *);
double Uniform(RndState *);
double Normal(RndState *);
double Exponential(RndState *);
double Gamma(RndState *,double a);
long   Poisson(RndState *,double a);
int    Dirichlet(RndState *,long n, double *a, double *x);
double Beta(RndState *,double a, double b);
long   Zeta(RndState *,double a);
int 	 Binomial(RndState *,double p, int n);
double Stable(RndState *,int param,double a, double b);
long   Discrete(RndState *, long n, double *p, double tot);

double pdf_Normal(double x);
double pdf_Beta(double a, double b, double x);
double pdf_Gamma(double a, double x);
double pdf_Exponential(double x);
double pdf_Zeta(double s, long x);
double pdf_Poisson(double m, long x);
double pdf_Binomial(double p, long n, long int x);
double pdf_Dirichlet(long n, double *a, double *x);
double pdf_Discrete(long n, double *p, double tot, long x);
double pdf_Uniform(double x);

int  vector_psi(long n, double *a, double *x);
