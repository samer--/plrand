/***********************************************************************\
 *
 * File:           RngStream.c for multiple streams of Random Numbers
 * Language:       ANSI C
 * Copyright:      Pierre L'Ecuyer, University of Montreal
 * Date:           14 August 2001
 * License: 	   GPL version 2 or later
 *
 * Notice:         Please contact P. L'Ecuyer at <lecuyer@iro.UMontreal.ca>
 *                 for commercial purposes.
 *
\***********************************************************************/


#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*---------------------------------------------------------------------*/
/* Private part.                                                       */
/*---------------------------------------------------------------------*/

typedef struct rng_state *RngStream;

#define norm  2.328306549295727688e-10
#define m1    4294967087.0
#define m2    4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

#define two17   131072.0
#define two53   9007199254740992.0
#define fact  5.9604644775390625e-8    /* 1 / 2^24 */



/* Default initial seed */
static double nextSeed[6] = { 12345, 12345, 12345, 12345, 12345, 12345 };


/* The following are the transition matrices of the two MRG components */
/* (in matrix form), raised to the powers 2^0, 2^76, and 2^127, resp.*/

static double A1p0[3][3] = {
          {       0.0,        1.0,       0.0 },
          {       0.0,        0.0,       1.0 },
          { -810728.0,  1403580.0,       0.0 }
          };

static double A2p0[3][3] = {
          {        0.0,        1.0,       0.0 },
          {        0.0,        0.0,       1.0 },
          { -1370589.0,        0.0,  527612.0 }
          };


static double A1p76[3][3] = {
          {      82758667.0, 1871391091.0, 4127413238.0 }, 
          {    3672831523.0,   69195019.0, 1871391091.0 }, 
          {    3672091415.0, 3528743235.0,   69195019.0 }
          };

static double A2p76[3][3] = {
          {    1511326704.0, 3759209742.0, 1610795712.0 }, 
          {    4292754251.0, 1511326704.0, 3889917532.0 }, 
          {    3859662829.0, 4292754251.0, 3708466080.0 }
          };

static double A1p127[3][3] = {
          {    2427906178.0, 3580155704.0,  949770784.0 }, 
          {     226153695.0, 1230515664.0, 3580155704.0 },
          {    1988835001.0,  986791581.0, 1230515664.0 }
          };

static double A2p127[3][3] = {
          {    1464411153.0,  277697599.0, 1610723613.0 },
          {      32183930.0, 1464411153.0, 1022607788.0 },
          {    2824425944.0,   32183930.0, 2093834863.0 }
          };





/*-------------------------------------------------------------------------*/

static double MultModM (double a, double s, double c, double m)
   /* Compute (a*s + c) % m. m must be < 2^35.  Works also for s, c < 0 */
{
   double v;
   long a1;
   v = a * s + c;
   if ((v >= two53) || (v <= -two53)) {
      a1 = (long) (a / two17);
      a -= a1 * two17;
      v = a1 * s;
      a1 = (long) (v / m);
      v -= a1 * m;
      v = v * two17 + a * s + c;
   }
   a1 = (long) (v / m);
   if ((v -= a1 * m) < 0.0)
      return v += m;
   else
      return v;
}


/*-------------------------------------------------------------------------*/


/* static void print_vec(double *p) */
/* { */
/* 	int i; */
/* 	for (i=0; i<6; i++) printf("%.0lf ",p[i]); */
/* 	printf("\n"); */
/* } */


static void MatVecModM (double A[3][3], double s[3], double v[3], double m)
   /* Returns v = A*s % m.  Assumes that -m < s[i] < m. */
   /* Works even if v = s. */
{
   int i;
   double x[3];
   for (i = 0; i < 3; ++i) {
      x[i] = MultModM (A[i][0], s[0], 0.0, m);
      x[i] = MultModM (A[i][1], s[1], x[i], m);
      x[i] = MultModM (A[i][2], s[2], x[i], m);
   }
	memcpy(v,x,3*sizeof(double));
}


static void MatMatModM (double A[3][3], double B[3][3], double C[3][3],
                        double m)
   /* Returns C = A*B % m. Work even if A = C or B = C or A = B = C. */
{
   int i, j;
   double V[3], W[3][3];
   for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) V[j] = B[j][i];
      MatVecModM (A, V, V, m);
      for (j = 0; j < 3; ++j) W[j][i] = V[j];
   }
	memcpy(C,W,9*sizeof(double));
}

static void MatTwoPowModM (double A[3][3], double B[3][3], double m, long e)
  /* Compute matrix B = (A^(2^e) % m);  works even if A = B */
{
   int i, j;

   /* initialize: B = A */
   if (A != B) {
      for (i = 0; i < 3; i++) {
         for (j = 0; j < 3; ++j)
            B[i][j] = A[i][j];
      }
   }
   /* Compute B = A^{2^e} */
   for (i = 0; i < e; i++)
      MatMatModM (B, B, B, m);
}


/* returns the next integer in the sequence 
 * new state goes in g1
 * g1 can be the same as g, resulting in 
 * an update in place.
 */
double next(RngStream g, RngStream g1)
{
   double p1, p2;

//	printf("before: "); print_vec(g->Cg);
   /* Component 1 */
   p1 = a12 * g->Cg[1] - a13n * g->Cg[0];
   p1 -= m1*floor(p1/m1);
   g1->Cg[0] = g->Cg[1];
   g1->Cg[1] = g->Cg[2];
   g1->Cg[2] = p1;

   /* Component 2 */
   p2 = a21 * g->Cg[5] - a23n * g->Cg[3];
   p2 -= m2*floor(p2/m2);
   g1->Cg[3] = g->Cg[4];
   g1->Cg[4] = g->Cg[5];
   g1->Cg[5] = p2;

	// printf("after: "); print_vec(g1->Cg);
   /* Combination */
   return (p1 > p2) ? (p1 - p2) : (p1 - p2 + m1);
}


/*-------------------------------------------------------------------------*/
static int CheckSeed (unsigned seed[6])
{
   /* Check that the seeds are legitimate values. Returns 0 if legal seeds,
	 * i:1..6 if seed[i] is bad
	 * 8 if seeds 1..3 are zero
	 * 9 if seeds 4..6 are zero
    */
   int i;

   for (i = 0; i < 3; ++i) {
      if (seed[i] >= m1) return i+1;
   }
   for (i = 3; i < 6; ++i) {
      if (seed[i] >= m2) return i+1;
   }
   if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0) return 8;
   if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0) return 9;

   return 0;
}


/*---------------------------------------------------------------------*/
/* Public part.                                                        */
/*---------------------------------------------------------------------*/


void RngStream_InitState(RngStream g) {
	memcpy(g->Cg,nextSeed,6*sizeof(double));
	//printf("init: "); print_vec(g->Cg);
}

void RngStream_InitJump(struct rng_jump *j, int e)
{
	void *a1, *a2;

	if (e>=127)     { a1=A1p127; a2=A2p127; e-=127; }
	else if (e==76) { a1=A1p76; a2=A2p76; e-=76; }
	else            { a1=A1p0; a2=A2p0; }

	memcpy(j->T1,a1,9*sizeof(double));
	memcpy(j->T2,a2,9*sizeof(double));
	MatTwoPowModM(j->T1,j->T1,m1,e);
	MatTwoPowModM(j->T2,j->T2,m2,e);
}
		
void RngStream_DoubleJump(struct rng_jump *j1, struct rng_jump *j2)
{
	MatMatModM(j1->T1,j1->T1,j2->T1,m1);
	MatMatModM(j1->T2,j1->T2,j2->T2,m2);
}

void RngStream_Advance(struct rng_jump *j, RngStream g1, RngStream g2)
{
   MatVecModM (j->T1, &g1->Cg[0], &g2->Cg[0], m1);
   MatVecModM (j->T2, &g1->Cg[3], &g2->Cg[3], m2);
}

double RngStream_Float(RngStream g, RngStream g1) {
	return norm*next(g,g1);
}

double RngStream_Double(RngStream g, RngStream g1) {
   double u = RngStream_Float(g,g1);
	u += fact*RngStream_Float(g1,g1);
	return (u < 1.0) ? u : (u - 1.0);
}

double RngStream_Raw(RngStream g, RngStream g1) { 
	return next(g,g1); 
}

/*-------------------------------------------------------------------------*/

int RngStream_SetState (RngStream g, unsigned seed[6])
{
   int i;
   if (CheckSeed (seed))
      return -1;                    /* FAILURE */
   for (i = 0; i < 6; ++i) 
      g->Cg[i] = (double)seed[i];
	
   return 0;                       /* SUCCESS */ 
}


void RngStream_GetState (RngStream g, unsigned seed[6])
{
   int i;
   for (i = 0; i < 6; ++i) seed[i] = (unsigned)g->Cg[i];
}

