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
double Zeta(RndState *,double a);
int 	 Binomial(RndState *,double p, int n);
double Stable(RndState *,int param,double a, double b);
long   Discrete(RndState *, long n, double *p, double tot);

