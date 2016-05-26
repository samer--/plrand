/*
 * Prolog-C utilities
 * Samer Abdallah (2009--2011)
 */
#include <SWI-Prolog.h>

int get_long(term_t term, long *p)
int get_double(term_t term, double *p)
int unify_list_doubles(term_t list, double *x, int n)
int get_list_doubles(term_t list, double *vals)
int alloc_array(size_t N, size_t SZ, void **PP) {
int get_args_doubles(term_t state, double *p, int n)
int unify_args_doubles(term_t state, double *p, int n)
double sum_array(double *p, int n) {

int type_error(term_t actual, const char *expected)
int memory_error(size_t amount)

