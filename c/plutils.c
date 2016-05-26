/*
 * Prolog-C utilities
 * Samer Abdallah (2009)
 */

//#define _USE_MATH_DEFINES 1
//#include "plutils.h"

//#include <math.h>
//#include <float.h>

// throws a Prolog exception to signal type error
int type_error(term_t actual, const char *expected)
{ 
	term_t ex = PL_new_term_ref();
	int rc = PL_unify_term(ex, PL_FUNCTOR_CHARS, "error", 2,
		      PL_FUNCTOR_CHARS, "type_error", 2,
		        PL_CHARS, expected,
		        PL_TERM, actual,
		      PL_VARIABLE);

  return rc && PL_raise_exception(ex);
}

double sum_array(double *p, int n) {
	double tot=0;
	int	i;
	for (i=0; i<n; i++) tot+=*p++;
	return tot;
}

int memory_error(size_t amount)
{ 
	term_t ex = PL_new_term_ref();
	int rc = PL_unify_term(ex, PL_FUNCTOR_CHARS, "error", 2,
		      PL_FUNCTOR_CHARS, "memory_error", 1,
		        PL_INTEGER, amount,
		      PL_VARIABLE);

	return rc && PL_raise_exception(ex);
}

// extract double from Prolog float
int get_double(term_t term, double *p)
{ 
	if (PL_get_float(term, p)) return TRUE; 
	else return type_error(term, "float");
}

// extract long from Prolog integer
int get_long(term_t term, long *p)
{ 
	if (PL_get_long(term, p)) return TRUE; 
	else return type_error(term, "integer");
}


// unify Prolog list of floats with array of doubles
int unify_list_doubles(term_t list, double *x, int n)
{
	int i;
	list=PL_copy_term_ref(list);

	for (i=0; i<n; i++) {
		term_t head=PL_new_term_ref();
		term_t tail=PL_new_term_ref();
		if (!PL_unify_list(list,head,tail)) PL_fail; 
		if (!PL_unify_float(head,x[i])) PL_fail;
		list=tail;
	}
	return PL_unify_nil(list);
}

// read list of floats from term and write to double array
int get_list_doubles(term_t list, double *vals)
{
	term_t  head=PL_new_term_ref();
	int 		n;

	// copy term ref so as not to modify original
	list=PL_copy_term_ref(list);
	for (n=0;PL_get_list(list,head,list);n++) {
			if (!PL_get_float(head,&vals[n])) return FALSE;
	}
	if (!PL_get_nil(list)) return FALSE; 
	return TRUE;
}

int unify_args_doubles(term_t state, double *p, int n)
{
	term_t arg = PL_new_term_ref();
	int i;

	for (i=0; i<n; i++) {
		if (!PL_put_float(arg,p[i])) PL_fail;
		if (!PL_unify_arg(i+1,state,arg)) PL_fail;
	}
	return TRUE;
}

int get_args_doubles(term_t state, double *p, int n)
{
	term_t arg = PL_new_term_ref();
	int i;

	for (i=0; i<n; i++) {
		_PL_get_arg(i+1, state, arg);
		if (!get_double(arg, p+i)) return FALSE;
	}
	return TRUE;
}

int alloc_array(size_t N, size_t SZ, void **PP) {
	*PP=calloc(N,SZ);
	if (*PP) return TRUE;
	else return memory_error(N*SZ);
}

