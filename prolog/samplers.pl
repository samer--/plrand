/*
 * Prolog part of random generator library
 * Samer Abdallah (2009)
*/
	  
:- module(plrand, [
		rv/2 					% -Spec, -Type
	,	sample/4				% +Dist, -Value, +StateIn, -StateOut
	,	sample/2				% +Dist, -Value
	,	dps_dist/3

	,  op(900,fx,\\)
	]).

/** <module> Random variable sampling interpreter

   This module provides a mini language for describing and sampling
   from random variables.
*/
	
:-	use_foreign_library(foreign(plrand)).
:- use_module(library(plrand))
:- multifile sample/4, rv/2.

%% rv(-Spec, -Type) is multi.
%
%  Unifies Spec and Type with specifications for all the distributions
%  known to the system.
%
%  @param Spec is a term representing a distribution that can be used
%         with sample/2 or sample/4. The head functor represents the
%         distribution to sampled from and any arguments represent the
%         types of the corresponding arguments to be used when sampling.
%  @param Type represents the type of the sampled value. 
rv( raw, natural).
rv( uniform01, nonneg).
rv( normal,    real).
rv( exponential, nonneg).
rv( gamma(nonneg), nonneg).
rv( students_t(nonneg), real).
rv( poisson(nonneg), natural).
rv( invgamma(nonneg), nonneg).
rv( beta(nonneg,nonneg), nonneg).
rv( zeta(nonneg), nonneg).
rv( binomial(nonneg,natural), natural).
rv( dirichlet(natural,list(nonneg)), list(nonneg)).
rv( dirichlet(list(nonneg)), list(nonneg)).
rv( dirichlet(tuple(nonneg)), tuple(nonneg)).
rv( stable(nonneg,real), real).
rv( bernoulli(nonneg), atom).
rv( bernoulli(nonneg), atom).
rv( discrete(natural,list(nonneg)),natural).
rv( discrete(list(nonneg)),natural).
rv( discrete(tuple(nonneg)),natural).

%% sample( +DistExpression, -Value) is det.
%% sample( +DistExpression, -Value, +StateIn, -StateOut) is det.
%
% 	sample/2 and sample/4 implement a small language for describing
%  and sampling from a distribution. 
%
%  sample/2 uses and modifies the global state.
%  sample/4 uses a given random generator state and returns the state
%  on completion, and is designed to compatible with DCG syntax to hide
%  the threading of the random state through several consecutive calls.
%
%  DistExpression is an expression describing a distribution.  The head
%  functor can be a distribution name (as listed by rv/2) or one of a number
%  of arithmetic operators or term constructors. The arguments (in almost all
%  cases) can then be further DistExpressions which are evaluated recursively. 
%  Valid non-distributional termsare:
%
%  	* X*Y
%		returns product of samples from X and Y
%  	* X/Y
%		returns ratio of samples from X and Y
%  	* X+Y
%		returns sum of samples from X and Y
%  	* X-Y
%		returns difference of samples from X and Y
%  	* -X
%		returns the negation of a sample from X
%		* sqrt(X)
%		returns square root of sample from X
%		* [X,Y,...]
%     return samples from X, Y etc in a list
%		* rep(N,X)
%     returns N independent sample from X in a list (N must be a constant)
%     * factorial(N,X)
%     returns N independent sample from X in an N argument term \(X1,...,XN)
%     * <any number>
%     returns itself as a constant
%
%     For example
%     ==
%     ?- sample(invgamma(1)*rep(8,discrete(dirichlet(rep(4,0.5)))),X).
%     X = [3, 2, 3, 3, 3, 1, 2, 2] .
%     ==

sample( raw, X)   --> sample_Raw(X1), {X is integer(X1)}.

sample( uniform01, X)   --> sample_Uniform01(X).
sample( normal, X)      --> sample_Normal(X).
sample( exponential, X) --> sample_Exponential(X).
sample( gamma(A), X)    --> sample(A,A1), sample_Gamma(A1,X).
sample( poisson(A), X)  --> sample(A,A1), sample_Poisson(A1,X).
sample( invgamma(A), X) --> sample(A,A1), sample_Gamma(A1,Y), {X is 1/Y}.
sample( beta(A,B), X)   --> sample(A,A1), sample(B,B1), sample_Beta(A1,B1,X).
sample( zeta(A), X)     --> sample(A,A1), sample_Zeta(A1,X1), {X is integer(X1)}.
sample( pareto(A), X)   --> sample(A,A1), sample_Uniform01(Y), {X is (1-Y)**(-1/A1) }.
sample( binomial(P,N), X)  --> sample(P,P1), sample(N,N1), sample_Binomial(P1,N1,X).

sample( stable(A,B), X)    --> sample(A,A1), sample(B,B1), sample_Stable(A1,B1,X).

sample( dirichlet(N,A), X) --> sample(A,A1), sample_Dirichlet(N,A1,X).
sample( dirichlet(A), X)   --> sample(A,A1), 
	(	{A1=[_|_]} 
	->	{length(A1,N)}, sample_Dirichlet(N,A1,X)
	;	{functor(A1,F,N), functor(X,F,N)}, sample_DirichletF(N,A1,X)
	).

sample( discrete(N,P), X)  --> sample(P,P1), sample_Discrete(N,P1,X).
sample( discrete(P), X)    --> 
	sample(P,P1), 
	(	{P1=[_|_]}
	->	{length(P1,N)}, sample_Discrete(N,P1,X)
	;	{functor(P1,_,N)}, sample_DiscreteF(N,P1,X)
	).

sample( bernoulli(P), X)   --> sample(P,P1), sample_Uniform01(U), {U<P1->X=1;X=0}.
sample( students_t(V), X)   --> sample(V/2,V1), sample(normal*sqrt(V1/gamma(V1)),X).

% dps(Vals) represents a countably infinite discrete distribution. It is an infinite
% stream of weight:value pairs.

sample(dps([B1:X1|Vals]),X) -->
   sample_Uniform01(U),
	(	{U<B1} -> {X=X1}
	;  sample(dps(Vals),X)
	).


sample( X*Y, Z) --> sample(X,X1), sample(Y,Y1), {Z is X1*Y1}.
sample( X/Y, Z) --> sample(X,X1), sample(Y,Y1), {Z is X1/Y1}.
sample( X+Y, Z) --> sample(X,X1), sample(Y,Y1), {Z is X1+Y1}.
sample( X-Y, Z) --> sample(X,X1), sample(Y,Y1), {Z is X1-Y1}.
sample( -X, Z)  --> sample(X,X1), {Z is -X1}.
sample( sqrt(X), Z) --> sample(X,X1), {Z is sqrt(X1)}.

sample( [], []) --> [].
sample( [X|XX], [Z|ZZ]) --> sample(X,Z), sample(XX,ZZ).
sample( rep(N,X), Z) --> {length(Z,N)}, seqmap(sample(X),Z).
sample( stream(X), Z) --> spawn(S), {freeze(Z,unfold_stream(S,X,Z))}.
sample( factorial(N,X), Z) --> {functor(Z,\,N)}, seqmapargs(sample(X),N,Z).
sample( \\(G), X) --> call(G,X).
sample( Tuple, Value) -->
	{functor(Tuple,F,N), functor(Value,F,N), tuple_functor(F)},
	seqmapargs(sample,N,Tuple,Value).

sample( N, N, S, S) :- number(N), !.

sample(M,X)     :- get_rnd_state(S1), sample(M,X,S1,S2), set_rnd_state(S2).

unfold_stream(S1,X,[Z1|ZX]) :- sample(X,Z1,S1,S2), freeze(ZX,unfold_stream(S2,X,ZX)).

% Return truncated infinite discrete distribution
dps_dist(dps(L),Probs,Vals) :- unfreeze_dps(1,L,Probs,Vals).

unfreeze_dps(_,_,[],[]).
unfreeze_dps(P0,[Q:V|T],[P|PX],[V|VX]) :- 
	P is P0*Q, P1 is P0*(1-Q),
	unfreeze_dps(P1,T,PX,VX).
	

tuple_functor(\).
tuple_functor(tuple).
tuple_functor(vec).

seqmapargs(P,N,X1) -->
	(	{N>0}
	->	{succ(M,N), arg(N,X1,X1N)},
		call(P,X1N),
		seqmapargs(P,M,X1)
	;	[]
	).
seqmapargs(P,N,X1,X2) -->
	(	{N>0}
	->	{succ(M,N), arg(N,X1,X1N), arg(N,X2,X2N)},
		call(P,X1N,X2N),
		seqmapargs(P,M,X1,X2)
	;	[]
	).

seqmap(_,[])             --> [].
seqmap(P,[A|AX])         --> call(P,A), seqmap(P,AX).
seqmap(_,[],[])          --> [].
seqmap(P,[A|AX],[B|BX])  --> call(P,A,B), seqmap(P,AX,BX).



