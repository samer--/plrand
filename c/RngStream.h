/* RngStream.h for ANSI C */
#ifndef RNGSTREAM_H
#define RNGSTREAM_H

struct rng_state {
   double Cg[6]; // state of generator
};

struct rng_jump {
	double T1[3][3]; // transition matrix for a jump in first generator
	double T2[3][3]; // ditto for second
};


void RngStream_InitState(struct rng_state *g);
void RngStream_InitJump(struct rng_jump *j,int e);
void RngStream_DoubleJump(struct rng_jump *j, struct rng_jump *j2);
int  RngStream_SetState(struct rng_state *g, unsigned seed[6]);
void RngStream_GetState(struct rng_state *g, unsigned seed[6]);
void RngStream_Advance(struct rng_jump *j, struct rng_state *g1, struct rng_state *g2);
double RngStream_Float(struct rng_state *g1, struct rng_state *g2);
double RngStream_Double(struct rng_state *g1, struct rng_state *g2);
double RngStream_Raw(struct rng_state *g1, struct rng_state *g2);

#endif
 

