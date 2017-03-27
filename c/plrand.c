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

#include <SWI-Stream.h>
#include <SWI-Prolog.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "RngStream.h"
#include "rndutils.h"
#include "plutils.c"

static PL_blob_t rs_blob, rj_blob;
static functor_t functor_rs7;
static RndState CurrentState;

int rs_write(IOSTREAM *,atom_t,int);
int rj_write(IOSTREAM *,atom_t,int);

install_t install();

foreign_t init_rnd_state( term_t s); 
foreign_t get_rnd_state( term_t s); 
foreign_t set_rnd_state( term_t s); 
foreign_t is_rnd_state( term_t s); 
foreign_t randomise( term_t s); 
foreign_t spawn( term_t s0, term_t s1, term_t s2);
foreign_t jump( term_t j, term_t s1, term_t s2);
foreign_t init_jump( term_t n, term_t j);
foreign_t double_jump( term_t j1, term_t j2);

foreign_t rnd_state_to_term( term_t rs, term_t t);
foreign_t term_to_rnd_state( term_t t, term_t rs);

foreign_t sample_Raw( term_t x, term_t s0, term_t s1); 
foreign_t sample_Uniform01( term_t x, term_t s0, term_t s1); 
foreign_t sample_Normal( term_t x, term_t s0, term_t s1); 
foreign_t sample_Exponential( term_t x, term_t s0, term_t s1); 
foreign_t sample_Gamma( term_t a, term_t x, term_t s0, term_t s1);
foreign_t sample_Poisson( term_t a, term_t x, term_t s0, term_t s1);
foreign_t sample_Beta( term_t a, term_t b, term_t x, term_t s0, term_t s1);
foreign_t sample_Zeta( term_t a, term_t x, term_t s0, term_t s1);
foreign_t sample_Binomial( term_t p, term_t n, term_t x, term_t s0, term_t s1);
foreign_t sample_Stable( term_t a, term_t b, term_t x, term_t s0, term_t s1);
foreign_t sample_Dirichlet( term_t n, term_t a, term_t x, term_t s0, term_t s1);
foreign_t sample_DirichletF( term_t n, term_t a, term_t x, term_t s0, term_t s1);
foreign_t sample_Discrete( term_t n, term_t p, term_t x, term_t s0, term_t s1);
foreign_t sample_DiscreteF( term_t n, term_t p, term_t x, term_t s0, term_t s1);

foreign_t prob_Uniform01( term_t x, term_t p); 
foreign_t prob_Normal( term_t x, term_t p);
foreign_t prob_Exponential( term_t x, term_t p);
foreign_t prob_Gamma( term_t a, term_t x, term_t p);
foreign_t prob_Poisson( term_t a, term_t x, term_t p);
foreign_t prob_Beta( term_t a, term_t b, term_t x, term_t p);
foreign_t prob_Zeta( term_t a, term_t x, term_t p);
foreign_t prob_Binomial( term_t q, term_t n, term_t x, term_t p);
foreign_t prob_Dirichlet( term_t n, term_t a, term_t x, term_t p);
foreign_t prob_Discrete( term_t n, term_t q, term_t x, term_t p);

foreign_t log_prob_Dirichlet( term_t n, term_t a, term_t x, term_t p);
foreign_t mean_log_Dirichlet( term_t n, term_t a, term_t x);
foreign_t log_partition_Dirichlet( term_t n, term_t a, term_t z);
foreign_t kldiv_Dirichlet( term_t n, term_t a, term_t b, term_t d);

foreign_t sample_Single_( term_t x); 
foreign_t sample_Double_( term_t x); 

foreign_t crp_prob( term_t alpha, term_t classes, term_t x, term_t pprob, term_t p);
foreign_t crp_sample( term_t alpha, term_t classes, term_t action, term_t rnd1, term_t rnd2);
foreign_t crp_sample_obs( term_t alpha, term_t classes, term_t x, term_t probx, term_t act, term_t prob, term_t rnd1, term_t rnd2);
foreign_t crp_sample_rm( term_t classes, term_t x, term_t class, term_t rnd1, term_t rnd2);
foreign_t sample_dp_teh( term_t ApSumKX, term_t B, term_t NX, term_t p1, term_t p2, term_t rnd1, term_t rnd2);
foreign_t sample_py_teh( term_t ThPrior, term_t DPrior, term_t CountsX, term_t p1, term_t p2, term_t rnd1, term_t rnd2);

static atom_t atom_new;
static functor_t functor_old1, functor_old2;
static functor_t functor_dp1, functor_py2;

/*
static void test() 
{
	RndState S;

	clock_t t1,t2;
	int	i, N=1000000;
	double dt;

	InitRndState(&S);

	t1=clock();
	for (i=0; i<N; i++) Uniform(&S);
	t2=clock();
	dt=(t2-t1)/(double)CLOCKS_PER_SEC;
	printf("%d doubles in %lf s --> %lf k doubles/s\n",N,dt,N/(1000.0*dt));
}
*/

install_t install() { 
	PL_register_foreign("init_rnd_state", 1, (void *)init_rnd_state, 0);
	PL_register_foreign("get_rnd_state", 1, (void *)get_rnd_state, 0);
	PL_register_foreign("set_rnd_state", 1, (void *)set_rnd_state, 0);
	PL_register_foreign("is_rnd_state", 1, (void *)is_rnd_state, 0);
	PL_register_foreign("randomise", 1, (void *)randomise, 0);
	PL_register_foreign("spawn", 3, (void *)spawn, 0);
	PL_register_foreign("jump", 3, (void *)jump, 0);

	PL_register_foreign("init_jump", 2, (void *)init_jump, 0);
	PL_register_foreign("double_jump", 2, (void *)double_jump, 0);

	PL_register_foreign("rnd_state_to_term", 2, (void *)rnd_state_to_term, 0);
	PL_register_foreign("term_to_rnd_state", 2, (void *)term_to_rnd_state, 0);

	PL_register_foreign("sample_Single_", 1, (void *)sample_Single_, 0);
	PL_register_foreign("sample_Double_", 1, (void *)sample_Double_, 0);
	PL_register_foreign("sample_Raw", 3, (void *)sample_Raw, 0);

	PL_register_foreign("sample_Uniform01", 3, (void *)sample_Uniform01, 0);
	PL_register_foreign("sample_Normal",    3, (void *)sample_Normal, 0);
	PL_register_foreign("sample_Exponential",    3, (void *)sample_Exponential, 0);
	PL_register_foreign("sample_Gamma",     4, (void *)sample_Gamma, 0);
	PL_register_foreign("sample_Poisson",   4, (void *)sample_Poisson, 0);
	PL_register_foreign("sample_Beta",      5, (void *)sample_Beta, 0);
	PL_register_foreign("sample_Zeta",      4, (void *)sample_Zeta, 0);
	PL_register_foreign("sample_Dirichlet", 5, (void *)sample_Dirichlet, 0);
	PL_register_foreign("sample_DirichletF",5, (void *)sample_DirichletF, 0);
	PL_register_foreign("sample_Binomial",  5, (void *)sample_Binomial, 0);
	PL_register_foreign("sample_Stable",    5, (void *)sample_Stable, 0);
	PL_register_foreign("sample_Discrete",  5, (void *)sample_Discrete, 0);
	PL_register_foreign("sample_DiscreteF", 5, (void *)sample_DiscreteF, 0);

	PL_register_foreign("prob_Uniform01", 2, (void *)prob_Uniform01, 0);
	PL_register_foreign("prob_Normal",    2, (void *)prob_Normal, 0);
	PL_register_foreign("prob_Exponential",    2, (void *)prob_Exponential, 0);
	PL_register_foreign("prob_Gamma",     3, (void *)prob_Gamma, 0);
	PL_register_foreign("prob_Poisson",   3, (void *)prob_Poisson, 0);
	PL_register_foreign("prob_Beta",      4, (void *)prob_Beta, 0);
	PL_register_foreign("prob_Zeta",      3, (void *)prob_Zeta, 0);
	PL_register_foreign("prob_Dirichlet", 4, (void *)prob_Dirichlet, 0);
	PL_register_foreign("prob_Binomial",  4, (void *)prob_Binomial, 0);
	PL_register_foreign("prob_Discrete",  4, (void *)prob_Discrete, 0);
	PL_register_foreign("log_prob_Dirichlet", 4, (void *)log_prob_Dirichlet, 0);
	PL_register_foreign("mean_log_Dirichlet", 3, (void *)mean_log_Dirichlet, 0);
	PL_register_foreign("log_partition_Dirichlet", 3, (void *)log_partition_Dirichlet, 0);
	PL_register_foreign("kldiv_Dirichlet", 4, (void *)kldiv_Dirichlet, 0);

	PL_register_foreign("crp_prob", 5, (void *)crp_prob, 0);
	PL_register_foreign("crp_sample", 5, (void *)crp_sample, 0);
	PL_register_foreign("crp_sample_obs", 8, (void *)crp_sample_obs, 0);
	PL_register_foreign("crp_sample_rm", 5, (void *)crp_sample_rm, 0);
	PL_register_foreign("sample_dp_teh", 7, (void *)sample_dp_teh, 0);
	PL_register_foreign("sample_py_teh", 7, (void *)sample_py_teh, 0);

	rs_blob.magic = PL_BLOB_MAGIC;
	rs_blob.flags = PL_BLOB_UNIQUE; // !!! temporary fix to avoid assertion fail (2010/05)
	rs_blob.name = "rndstate";
	rs_blob.acquire = 0; // rs_acquire;
	rs_blob.release = 0; // rs_release;
	rs_blob.compare = 0; // rs_compare;
	rs_blob.write   = rs_write;

	rj_blob.magic = PL_BLOB_MAGIC;
	rj_blob.flags = PL_BLOB_UNIQUE;
	rj_blob.name = "rndjump";
	rj_blob.acquire = 0; 
	rj_blob.release = 0;
	rj_blob.compare = 0;
	rj_blob.write   = rj_write;

	functor_rs7 = PL_new_functor(PL_new_atom("rs"),7);
	functor_dp1  = PL_new_functor(PL_new_atom("dp"),1);
	functor_py2  = PL_new_functor(PL_new_atom("py"),2);
	functor_old1 = PL_new_functor(PL_new_atom("old"),1);
	functor_old2 = PL_new_functor(PL_new_atom("old"),2);
	atom_new     = PL_new_atom("new");

	InitRndState(&CurrentState);
	//test();
}


static int check(int cond, const char *expected, term_t got) {
	if (!cond) return PL_domain_error(expected, got);
	else return TRUE;
}

// unify Prolog BLOB with RndState structure
static int unify_state(term_t state,RndState *S) {
	return PL_unify_blob(state, S, sizeof(RndState), &rs_blob); 
}

static int unify_state_term(term_t t, RndState *S) {
	double	*Cg=&S->g.Cg[0];
	int64_t	*prev=(int64_t *)&S->prev_normal;
	return PL_unify_term(t, PL_FUNCTOR, functor_rs7, 
			PL_INT64, (int64_t)Cg[0],
			PL_INT64, (int64_t)Cg[1],
			PL_INT64, (int64_t)Cg[2],
			PL_INT64, (int64_t)Cg[3],
			PL_INT64, (int64_t)Cg[4],
			PL_INT64, (int64_t)Cg[5],
			PL_INT64, *prev);
}

static int get_state_term(term_t t, RndState *S) {
	term_t	a=PL_new_term_ref();
	double	*Cg=&S->g.Cg[0];
	int64_t	*prev=(int64_t *)&S->prev_normal;

	return PL_is_functor(t,functor_rs7)
		&&	(_PL_get_arg(1,t,a), PL_get_float(a, &Cg[0]))
		&&	(_PL_get_arg(2,t,a), PL_get_float(a, &Cg[1]))
		&&	(_PL_get_arg(3,t,a), PL_get_float(a, &Cg[2]))
		&&	(_PL_get_arg(4,t,a), PL_get_float(a, &Cg[3]))
		&&	(_PL_get_arg(5,t,a), PL_get_float(a, &Cg[4]))
		&&	(_PL_get_arg(6,t,a), PL_get_float(a, &Cg[5]))
		&& (_PL_get_arg(7,t,a), PL_get_int64(a, prev));
}

// extract RndState structure from Prolog BLOB
static int get_state(term_t state, RndState *S0)
{ 
	PL_blob_t *type;
	size_t    len;
	RndState *S;
  
	PL_get_blob(state, (void **)&S, &len, &type);
	if (type != &rs_blob) {
		return type_error(state, "rndstate");
	} else {
		*S0=*S;
		return TRUE;
	}
} 

int rs_write(IOSTREAM *s, atom_t a, int flags) 
{ 
	PL_blob_t *type;
	size_t    len;
	RndState *p=(RndState *)PL_blob_data(a,&len,&type);
	if (p) {
		int i;
		Sfprintf(s,"<rs");
		for (i=0; i<6; i++) Sfprintf(s," %.0lf",p->g.Cg[i]);
		Sfprintf(s,">");
	}
	return TRUE; 
}

// unify Prolog BLOB with RndJump structure
static int unify_jump(term_t jump,RndJump *J) {
	return PL_unify_blob(jump, J, sizeof(RndJump), &rj_blob); 
}

// extract RndState structure from Prolog BLOB
static int get_jump(term_t jump, RndJump *J0)
{ 
	PL_blob_t *type;
	size_t    len;
	RndJump *J;
  
	PL_get_blob(jump, (void **)&J, &len, &type);
	if (type != &rj_blob) {
		return type_error(jump, "rndjump");
	} else {
		*J0=*J;
		return TRUE;
	}
} 

int rj_write(IOSTREAM *s, atom_t a, int flags) 
{ 
	PL_blob_t *type;
	size_t    len;
	RndJump *p=(RndJump *)PL_blob_data(a,&len,&type);
	if (p) {
		Sfprintf(s,"<rj %d>",p->e);
	}
	return TRUE; 
}

// -----------------------------------------------------

// unify jump with a new blob representing a sequence jump
foreign_t init_jump(term_t pow, term_t jump) { 
	RndJump J;
	long e;
	
	if (get_long(pow,&e)) {
		J.e=e;
		RngStream_InitJump(&J.j,e);
		return unify_jump(jump,&J);
	} else return FALSE;
}

foreign_t double_jump(term_t j0, term_t j1) { 
	RndJump J;
	
	if (get_jump(j0,&J)) {
		J.e++;
		RngStream_DoubleJump(&J.j,&J.j);
		return unify_jump(j1,&J);
	} else return FALSE;
}

// unify Prolog term with the default random state 
foreign_t init_rnd_state(term_t state) { 
	RndState S0;
	InitRndState(&S0);
	return unify_state(state, &S0);
}

// unify Prolog term with current random state structure
foreign_t get_rnd_state(term_t state) { return unify_state(state, &CurrentState); }

// set current random state structure to values in Prolog term
foreign_t set_rnd_state(term_t state) { return get_state(state, &CurrentState); }

foreign_t rnd_state_to_term( term_t state, term_t t) 
{
	RndState S;
	return get_state(state, &S) && unify_state_term(t,&S);
}

foreign_t term_to_rnd_state( term_t t, term_t state)
{
	RndState S;
	return get_state_term(t, &S) && unify_state(state, &S);
}

// get truly random state from /dev/random
foreign_t randomise(term_t state) { 
	RndState S;
	FILE		*p;

	p=fopen("/dev/random","r");
	if (p!=NULL) {
		int rc=randomise_from(p,&S);
		fclose(p);
		return rc && unify_state(state, &S); 
	}
	{
		term_t ex = PL_new_term_ref();
		int rc = PL_unify_term(ex, PL_FUNCTOR_CHARS, "error", 2,
					PL_FUNCTOR_CHARS, "file_error", 1,
					  PL_CHARS, "fread /dev/random",
					PL_VARIABLE);

	  return rc && PL_raise_exception(ex);
	}
}
	
foreign_t spawn(term_t new, term_t orig, term_t next) { 
	RndState S0, S1;

	return get_state(orig, &S0)
			&& spawn_gen(&S0, &S1)
			&& unify_state(next, &S0) 
			&& unify_state(new, &S1); 
}

foreign_t jump(term_t j, term_t s0, term_t s1) { 
	RndState S;
	RndJump  J;

	if (get_jump(j, &J) && get_state(s0, &S)) {
		RngStream_Advance(&J.j, &S.g, &S.g);
		return unify_state(s1, &S); 
	} else return FALSE;
}

// set current random state structure to values in Prolog term
foreign_t is_rnd_state(term_t state) { 
	PL_blob_t *type;
	return PL_is_blob(state,&type) && type==&rs_blob;
}

// unify one term with numeric value and other with random state structure
static int unify_float_state(term_t x, double z, term_t s, RndState *S) {
	return  PL_unify_float(x, z) && unify_state(s, S);
}
	
// unify one term with numeric value and other with random state structure
//static int unify_integer(term_t x, long z, term_t s, RndState *S) {
//	return  PL_unify_integer(x, z) && unify_state(s,S);
//}
	
foreign_t sample_Single_(term_t x) {
	return PL_unify_float(x, RngStream_Float(&CurrentState.g,&CurrentState.g));
}

foreign_t sample_Double_(term_t x) {
	return PL_unify_float(x, RngStream_Double(&CurrentState.g,&CurrentState.g));
}


foreign_t sample_Raw(term_t x, term_t s0, term_t s1) {
	RndState S;
	return get_state(s0,&S) &&
	       unify_float_state(x, Raw(&S), s1, &S);
}


foreign_t sample_Uniform01(term_t x, term_t s0, term_t s1) {
	RndState S;
	return get_state(s0, &S) && PL_unify_float(x, Uniform(&S)) && unify_state(s1, &S);
}


foreign_t sample_Normal(term_t x, term_t s0, term_t s1) {
	RndState S;
	return get_state(s0,&S) &&
		    unify_float_state(x, Normal(&S), s1, &S);
}

foreign_t sample_Exponential(term_t x, term_t s0, term_t s1) {
	RndState S;
	return get_state(s0,&S) &&
		    unify_float_state(x, Exponential(&S), s1, &S);
}



foreign_t sample_Gamma(term_t a, term_t x, term_t s0, term_t s1) {
	RndState S;
	double A;

	return get_state(s0,&S) &&
		    get_double(a,&A) &&
			 unify_float_state(x, Gamma(&S,A), s1, &S);
}


foreign_t sample_Poisson(term_t a, term_t x, term_t s0, term_t s1) {
	RndState S;
	double A;

	return get_state(s0, &S) &&
		    get_double(a, &A) &&
	       PL_unify_integer(x, Poisson(&S,A)) &&
			 unify_state(s1, &S);
}

foreign_t sample_Dirichlet(term_t n, term_t a, term_t x, term_t s0, term_t s1)
{
	long N;

	if (get_long(n,&N)) {
		double   *A=NULL, *X=NULL;
		RndState S;

		int rc = alloc_array(N,sizeof(double),(void **)&A) 
				&& alloc_array(N,sizeof(double),(void **)&X) 
				&& get_list_doubles(a,A)
				&& get_state(s0,&S)
				&& Dirichlet(&S,N,A,X)
				&& unify_list_doubles(x,X,N)
				&& unify_state(s1, &S);

		if (X) free(X);
		if (A) free(A);
		return rc;
	} else return FALSE;
}

foreign_t sample_DirichletF(term_t n, term_t a, term_t x, term_t s0, term_t s1)
{
	long N;

	if (get_long(n,&N)) {
		double   *A=NULL, *X=NULL;
		RndState S;

		int rc = alloc_array(N,sizeof(double),(void **)&A) 
				&& alloc_array(N,sizeof(double),(void **)&X) 
				&& get_args_doubles(a,A,N)
				&& get_state(s0,&S)
				&& Dirichlet(&S,N,A,X)
				&& unify_args_doubles(x,X,N)
				&& unify_state(s1, &S);

		if (X) free(X);
		if (A) free(A);
		return rc;
	} else return FALSE;
}



foreign_t sample_Beta(term_t a, term_t b, term_t x, term_t s0, term_t s1) {
	RndState S;
	double A, B;

	return get_state(s0,&S) &&
		    get_double(a,&A) &&
	       get_double(b,&B) && 
			 unify_float_state(x, Beta(&S,A,B),s1,&S);
}

foreign_t sample_Zeta(term_t a, term_t x, term_t s0, term_t s1) {
	RndState S;
	double A;

	return get_state(s0,&S) &&
		    get_double(a,&A) &&
			 check(A>1.0, "Zeta parameter > 1.0", a) &&
			 PL_unify_integer(x, Zeta(&S,A)) &&
			 unify_state(s1, &S);
}

foreign_t sample_Binomial(term_t p, term_t n, term_t x, term_t s0, term_t s1) {
	RndState S;
	double P;
	long   N;

	return get_state(s0,&S) &&
		    get_double(p,&P) &&
	       get_long(n,&N) && 
	       PL_unify_integer(x, Binomial(&S,P,(int)N)) &&
			 unify_state(s1, &S);
}


foreign_t sample_Stable(term_t a, term_t b, term_t x, term_t s0, term_t s1) {
	RndState S;
	double A, B;

	return get_state(s0,&S) &&
		    get_double(a,&A) &&
	       get_double(b,&B) && 
			 unify_float_state(x, Stable(&S,1,A,B),s1,&S);
}

foreign_t sample_Discrete(term_t n, term_t p, term_t x, term_t s0, term_t s1)
{
	long N;

	if (get_long(n,&N)) {
		double   *P;
		RndState S;

		int rc = alloc_array(N,sizeof(double),(void **)&P) 
				&& get_list_doubles(p,P)
				&& get_state(s0,&S)
				&& PL_unify_integer(x, 1+Discrete(&S,N,P,sum_array(P,N)));
			   unify_state(s1, &S);

		if (P) free(P);
		return rc;
	} else return FALSE;
}

foreign_t sample_DiscreteF(term_t n, term_t p, term_t x, term_t s0, term_t s1)
{
	long N;

	if (get_long(n,&N)) {
		double   *P;
		RndState S;

		int rc = alloc_array(N,sizeof(double),(void **)&P) 
				&& get_args_doubles(p,P,N)
				&& get_state(s0,&S)
				&& PL_unify_integer(x, 1+Discrete(&S,N,P,sum_array(P,N)));
			   unify_state(s1, &S);

		if (P) free(P);
		return rc;
	} else return FALSE;
}

foreign_t prob_Uniform01(term_t x, term_t p)   { double X; return get_double(x,&X) && PL_unify_float(p, pdf_Uniform(X)); }
foreign_t prob_Normal(term_t x, term_t p)      { double X; return get_double(x,&X) && PL_unify_float(p, pdf_Normal(X)); }
foreign_t prob_Exponential(term_t x, term_t p) { double X; return get_double(x,&X) && PL_unify_float(p, pdf_Exponential(X)); }
foreign_t prob_Gamma(term_t a, term_t x, term_t p)   { double A, X; return get_double(a,&A) && get_double(x,&X) && PL_unify_float(p, pdf_Gamma(A,X)); }
foreign_t prob_Poisson(term_t a, term_t x, term_t p) { double A, X; return get_double(a,&A) && get_double(x,&X) && PL_unify_float(p, pdf_Poisson(A,X)); }
foreign_t prob_Zeta(term_t a, term_t x, term_t p)    { 
	double A, X; 
	return get_double(a,&A) && 
			 check(A>1.0, "Zeta parameter > 1.0", a) &&
	     	 get_double(x,&X) && PL_unify_float(p, pdf_Zeta(A,X)); }

foreign_t prob_Beta(term_t a, term_t b, term_t x, term_t p) { 
	double A, B, X; 
	return get_double(a,&A) && get_double(b,&B) && get_double(x,&X) 
		 && PL_unify_float(p, pdf_Gamma(A,X)); 
}

foreign_t prob_Binomial(term_t q, term_t n, term_t x, term_t p) { 
	double Q; 
	long N, X; 
	return get_double(q,&Q) && get_long(n,&N) && get_long(x,&X) 
		 && PL_unify_float(p, pdf_Binomial(Q,N,X));
}

static int gen_prob_Dirichlet(double (*prob_fn)(long, double*, double*), term_t n, term_t a, term_t x, term_t p) { 
	long N;
	double *A=NULL, *X=NULL; 
	int r = get_long(n,&N) 
	     && alloc_array(N,sizeof(double),(void **)&A) && get_list_doubles(a,A) 
	     && alloc_array(N,sizeof(double),(void **)&X) && get_list_doubles(x,X)
		  && PL_unify_float(p, prob_fn(N,A,X)); 
	if (A) free(A);
	if (X) free(X);
	return r;
}

foreign_t prob_Dirichlet(term_t n, term_t a, term_t x, term_t p) { return gen_prob_Dirichlet(pdf_Dirichlet,n,a,x,p); }
foreign_t log_prob_Dirichlet(term_t n, term_t a, term_t x, term_t p) { return gen_prob_Dirichlet(logpdf_Dirichlet,n,a,x,p); }

foreign_t prob_Discrete(term_t n, term_t q, term_t x, term_t p) { 
	long N, X;
	double *Q=NULL; 
   int r = get_long(n,&N) && get_long(x,&X)
	     && alloc_array(N,sizeof(double),(void **)&Q) && get_list_doubles(q,Q)
     	  && PL_unify_float(p, pdf_Discrete(N,Q,sum_array(Q,N),X)); 
	if (Q) free(Q);
	return r;
}

foreign_t mean_log_Dirichlet(term_t n, term_t a, term_t x) { 
	long N;
	double *A=NULL, *X=NULL; 
	int r = get_long(n,&N) 
	     && alloc_array(N,sizeof(double),(void **)&A) && get_list_doubles(a,A) 
	     && alloc_array(N,sizeof(double),(void **)&X) && vector_psi(N,A,X)
		  && unify_list_doubles(x,X,N);
	if (A) free(A);
	if (X) free(X);
	return r;
}

foreign_t log_partition_Dirichlet(term_t n, term_t a, term_t z) { 
	long N;
	double *A=NULL; 
	int r = get_long(n,&N) 
	     && alloc_array(N,sizeof(double),(void **)&A) && get_list_doubles(a,A) 
	     && PL_unify_float(z, vector_betaln(N,A));
	if (A) free(A);
	return r;
}

static double kldiv_dirichlet(int N, double *A, double *B) {
	double D = vector_betaln(N,B) - vector_betaln(N,A);
	double *PsiA;

	alloc_array(N,sizeof(double),(void **)&PsiA);
	vector_psi(N,A,PsiA);
	for (int i=0; i<N; i++) { D += (A[i]-B[i])*PsiA[i]; }
	free(PsiA);
	return D;
}

foreign_t kldiv_Dirichlet(term_t n, term_t a, term_t b, term_t d) { 
	long N;
	double *A=NULL, *B=NULL; 
	int r = get_long(n,&N) 
	     && alloc_array(N,sizeof(double),(void **)&A) && get_list_doubles(a,A) 
	     && alloc_array(N,sizeof(double),(void **)&B) && get_list_doubles(b,B) 
	     && PL_unify_float(d, kldiv_dirichlet(N,A,B));
	if (A) free(A);
	if (B) free(B);
	return r;
}

// ----------------------------------------------------------------------------
// Chinese Restaurant Processes

int counts_dist( term_t gem, term_t counts, size_t len, double *dist);
int get_classes(term_t Classes, term_t Counts, term_t Vals, long *len);
void stoch(double *x, size_t len);

/*
%% crp_prob( +GEM:gem_model, +Classes:classes(A), +X:A, +PProb:float, -Prob:float) is det.
%
%  Compute the probability Prob of observing X given a CRP
%  and a base probability of PProb.
crp_prob( Alpha, classes(Counts,Vals), X, PProb, P) :-
	counts_dist( Alpha, Counts, Counts1),
	stoch( Counts1, Probs, _),
	maplist( equal(X), Vals, Mask),
	maplist( mul, [PProb | Mask], Probs, PostProbs),
	sumlist( PostProbs, P).

*/

foreign_t crp_prob( term_t Alpha, term_t Classes, term_t X, term_t PProb, term_t Prob)
{
	term_t Counts=PL_new_term_ref();
	term_t Vals=PL_new_term_ref();
	double prob=0, pprob;
	double *dist=NULL;
	long len=0;

	int rc =	get_double(PProb, &pprob)
		&&	get_classes(Classes, Counts, Vals, &len)
		&& alloc_array(len+1, sizeof(double), (void **)&dist)
		&& counts_dist(Alpha, Counts, len, dist);

	if (rc) {
		term_t Val = PL_new_term_ref();
		int i;

		prob = pprob*dist[0];
		for (i=1; i<=len && PL_get_list(Vals,Val,Vals); i++) {
			if (PL_unify(Val,X)) prob += dist[i];
		}
		prob /= sum_array(dist,len+1);
	} else rc=0;
	if (dist) free(dist);
	return rc && PL_unify_float(Prob,prob);
}

/*


%% crp_sample( +GEM:gem_model, +Classes:classes(A), -A:action(A))// is det.
%
%  Sample a new value from CRP, Action A is either new, which means
%  that the user should sample a new value from the base distribtion,
%  or old(X,C), where X is an old value and C is the index of its class.
%  Operates in random state DCG.
crp_sample( Alpha, classes(Counts,Vals), Action, Prob, RS1, RS2) :-
	counts_dist(Alpha, Counts, Counts1),
	discrete(Counts1,Z,RS1,RS2),
	( Z>1 -> succ(C,Z), nth1(C,Vals,X), Action=old(X,C)
	; Action=new).

*/

foreign_t crp_sample( term_t Alpha, term_t Classes, term_t Action, term_t Rnd1, term_t Rnd2)
{
	term_t Counts=PL_new_term_ref();
	term_t Vals=PL_new_term_ref();
	double *dist=NULL;
	RndState rs;
	long len=0;

	int rc =	get_classes(Classes, Counts, Vals, &len)
		&& alloc_array(len+1, sizeof(double), (void **)&dist)
		&& counts_dist(Alpha, Counts, len, dist)
		&& get_state(Rnd1,&rs);

	if (rc) {
		int z=Discrete( &rs, len+1, dist, sum_array(dist,len+1));

		if (z==0) { rc = PL_unify_atom(Action,atom_new); }
		else { 
			term_t X=PL_new_term_ref();
			int 	i=0;
			while (i<z && PL_get_list(Vals,X,Vals)) i++;
			rc = (i==z) && PL_unify_term(Action, PL_FUNCTOR, functor_old2, PL_TERM, X, PL_INT, z);
		}
	}
	if (dist) free(dist); 
	return rc && unify_state(Rnd2, &rs);
}

/*

%% crp_sample_obs( +GEM:gem_model, +Classes:classes(A), +X:A, +PProb:float, -A:action, -Prob:float)// is det.
%
%  Sample class appropriate for observation of value X. PProb is the
%  base probability of X from the base distribution. Action A is new
%  or old(Class).
%  Operates in random state DCG.
crp_sample_obs( Alpha, classes(Counts,Vals), X, ProbX, A, Prob, RS1, RS2) :-
	counts_dist( Alpha, Counts, [CNew|Counts1]),	
	PNew is CNew*ProbX,
	maplist( post_count(X),Vals,Counts1,Counts2),
	discrete( [PNew|Counts2], Z, RS1, RS2),
	discrete( [PNew|Counts2], Z, p(1), p(Prob)),
	( Z=1 -> A=new; succ(C,Z), A=old(C)).

*/

foreign_t crp_sample_obs( term_t Alpha, term_t Classes, term_t X, term_t Probx, term_t Act, term_t Prob, term_t Rnd1, term_t Rnd2)
{
	term_t 		Counts=PL_new_term_ref();
	term_t 		Vals=PL_new_term_ref();
	double	probx=0;
	double 	*dist=NULL;
	long 		len=0;
	RndState rs;

	int rc =	get_double(Probx,&probx)
		&& get_classes(Classes, Counts, Vals, &len)
		&& alloc_array(len+1, sizeof(double), (void **)&dist)
		&& counts_dist(Alpha, Counts, len, dist)
		&& get_state(Rnd1,&rs);

	if (rc) {
		term_t 	Val=PL_new_term_ref();
		double   total = sum_array(dist,len+1), paccum;
		int		i, z;

		dist[0] *= probx;
		paccum = dist[0];
		for (i=1; i<=len && PL_get_list(Vals,Val,Vals); i++) {
			if (!PL_unify(Val,X)) dist[i]=0; else paccum+=dist[i];
		}

		z=Discrete( &rs, len+1, dist, sum_array(dist,len+1));
		if (z==0) { rc = PL_unify_atom(Act,atom_new); }
		else { 
			rc = PL_unify_term(Act, PL_FUNCTOR, functor_old1, PL_INT, z);
		}
		rc &= PL_unify_float(Prob, paccum/total);
	}
	if (dist) free(dist); 
	return rc && unify_state(Rnd2, &rs);
}

/*
%% crp_sample_rm( +Classes:classes(A), +X:A, -C:natural)// is det.
%
%  Sample appropriate class from which to remove value X.
%  Operates in random state DCG.
crp_sample_rm( classes(Counts,Vals), X, Class, RS1, RS2) :-
	maplist(post_count(X),Vals,Counts,Counts1),
	discrete( Counts1, Class, RS1, RS2).

*/

foreign_t crp_sample_rm( term_t Classes, term_t X, term_t Class, term_t Rnd1, term_t Rnd2)
{
	term_t 		Counts=PL_new_term_ref();
	term_t 		Vals=PL_new_term_ref();
	double 	*dist=NULL;
	long 		len=0;
	RndState rs;

	int rc =	get_classes(Classes, Counts, Vals, &len)
		&& alloc_array(len, sizeof(double), (void **)&dist)
		&& get_list_doubles(Counts, dist)
		&& get_state(Rnd1,&rs);

	if (rc) {
		term_t 	Val=PL_new_term_ref();
		int		i, z;

		for (i=0; i<len && PL_get_list(Vals,Val,Vals); i++) {
			if (!PL_unify(Val,X)) dist[i]=0;
		}

		z = Discrete( &rs, len, dist, sum_array(dist,len));
		rc = (z<len) && PL_unify_integer(Class, z+1);
	}
	if (dist) free(dist); 
	return rc && unify_state(Rnd2, &rs);
}

/*
post_count(X,Val,Count,PC) :- X=Val -> PC=Count; PC=0.

% -----------------------------------------------------------
% Dirichlet process and Pitman-Yor process
% pseudo-counts models.

counts_dist(_,[],0,[1]) :- !.
counts_dist(dp(Alpha),Counts,_,[Alpha|Counts]) :- !.
counts_dist(py(Alpha,Discount),Counts,K,[CNew|Counts1]) :- !,
	CNew is Alpha+Discount*K,
	maplist(sub(Discount),Counts,Counts1).

*/

int get_float_arg(int n,term_t Term, double *px)
{
	term_t X=PL_new_term_ref();
	return PL_get_arg(n,Term,X) && PL_get_float(X,px);
}

int counts_dist( term_t gem, term_t counts, size_t len, double *dist)
{
	if (len==0) { dist[0]=1; return TRUE; }
	else {
		if (PL_is_functor(gem, functor_dp1)) {
			double alpha;
			term_t head=PL_new_term_ref();
			int i, rc = get_float_arg(1,gem,&alpha);

			dist[0] = alpha;
			for(i=1; rc && i<=len && PL_get_list(counts,head,counts); i++) {
				rc = rc && PL_get_float(head,&dist[i]);
			}
			return rc;
		} else if (PL_is_functor(gem, functor_py2)) {
			double	theta, disc, c;
			term_t	head=PL_new_term_ref();

			int i, rc = get_float_arg(1,gem,&theta)
						&& get_float_arg(2,gem,&disc);

			dist[0] = theta + disc*len;
			for(i=1; rc && i<=len && PL_get_list(counts,head,counts); i++) {
				rc = rc && PL_get_float(head,&c);
				dist[i] = c-disc;
			}
			return rc;
		} else return FALSE;
	}
}

int get_classes(term_t Classes, term_t Counts, term_t Vals, long *len)
{
	term_t K=PL_new_term_ref();

	return 	PL_get_arg(1,Classes,K)
			&& PL_get_arg(2,Classes,Counts)
			&& PL_get_arg(3,Classes,Vals)
			&& PL_get_long(K,len);
}


void stoch(double *x, size_t len)
{
	int i;
	double total=0, *xp;
	for (i=0, xp=x; i<len; i++, xp++) total += *xp;
	for (i=0, xp=x; i<len; i++, xp++) *xp /= total; 
}

/*
sample_dp_teh( ApSumKX, B, NX, dp(Alpha1), dp(Alpha2)) -->
	{ Alpha1_1 is Alpha1+1  },
	seqmap(beta(Alpha1_1),NX,WX),
	seqmap(bernoulli(Alpha1),NX,SX),
	{	maplist(log,WX,LogWX),
		sumlist(SX,SumSX), 
		sumlist(LogWX,SumLogWX), 
		A1 is ApSumKX-SumSX, B1 is B-SumLogWX
	},
	gamma(A1,B1,Alpha2).

%	run_left( seqmap(accum_log_beta(Alpha1_1),NX), 0, SumLogWX), 
%	run_left( seqmap(accum_bernoulli(Alpha1),NX), 0, SumSX), 
%accum_log_beta(A,B) --> \> beta(A,B,X), { LogX is log(X) }, \< add(LogX).
%accum_bernoulli(A,B) --> \> bernoulli(A,B,X), \< add(X).  

 */

int Bernoulli(RndState *rs,double a,double b) {
	if ((a+b)*Uniform(rs)<b) return 1; else return 0;
}


foreign_t sample_dp_teh( term_t ApSumKX, term_t B, term_t NX, term_t p1, term_t p2, term_t rnd1, term_t rnd2)
{
	term_t	N=PL_new_term_ref();
	double	apsumkx, b, alphap1;
	double	alpha1=0, alpha2;
	double 	sum_log_wx, sum_sx;
	long		n=0;
	RndState rs;

	int rc = get_double(ApSumKX, &apsumkx)
			&& get_double(B, &b)
			&& get_float_arg(1,p1,&alpha1)
			&& get_state(rnd1,&rs);

	alphap1 = alpha1+1;
	sum_log_wx = sum_sx = 0;
	while (rc && PL_get_list(NX,N,NX)) {
		rc = get_long(N,&n);
		sum_log_wx += log(Beta(&rs,alphap1,n));
		sum_sx     += Bernoulli(&rs,alpha1,n);
	}
	alpha2 = Gamma(&rs, apsumkx-sum_sx)/(b-sum_log_wx);

	return rc && PL_unify_term(p2, PL_FUNCTOR, functor_dp1, PL_FLOAT, alpha2)
			&& unify_state(rnd2, &rs);
}

foreign_t sample_py_teh( term_t ThPrior, term_t DPrior, term_t CountsX, term_t p1, term_t p2, term_t rnd1, term_t rnd2)
{
	term_t	Counts = PL_new_term_ref();
	term_t	Count = PL_new_term_ref();
	double	theta_a,  disc_a;
	double	theta_b,  disc_b;
	double	theta1,   disc1;
	double	theta2,   disc2;
	double   theta1_1, disc1_1;
	double 	sum_log_wx, sum_sx, sum_nsx, sum_zx;
	RndState rs;

	int rc = get_float_arg(1,ThPrior,&theta_a)
			&& get_float_arg(2,ThPrior,&theta_b)
			&& get_float_arg(1,DPrior,&disc_a)
			&& get_float_arg(2,DPrior,&disc_b)
			&& get_float_arg(1,p1,&theta1)
			&& get_float_arg(2,p1,&disc1)
			&& get_state(rnd1,&rs);

	theta1_1 = theta1+1;
	disc1_1	= 1-disc1;
	sum_log_wx = sum_sx = sum_nsx = sum_zx = 0;
	while (rc && PL_get_list(CountsX,Counts,CountsX)) {
		int n, k, i;
		long c=0;

		for(k=0, n=0; rc && PL_get_list(Counts,Count,Counts); k++, n+=c) {
			rc = get_long(Count,&c);
			if (k>0) { if (Bernoulli(&rs, disc1*k, theta1)) sum_sx++; else sum_nsx++; }
			for (i=0; i<c-1; i++) sum_zx += Bernoulli(&rs, i, disc1_1);
		}
		if (n>1) sum_log_wx += log(Beta(&rs, theta1_1, n-1));
	}

	theta2 = Gamma(&rs, theta_a + sum_sx)/(theta_b-sum_log_wx);
	disc2  = Beta(&rs, disc_a + sum_nsx, disc_b + sum_zx);
	return rc && unify_state(rnd2, &rs)
			&& PL_unify_term(p2, PL_FUNCTOR, functor_py2, PL_FLOAT, theta2, PL_FLOAT, disc2);
}

