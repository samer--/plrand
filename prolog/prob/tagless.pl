:- module(prob_tagless,
          [ uniform01//1
          , uniform//2, uniformP//2
          , normal//1
          , gaussian//3
          , exponential//1
          , poisson//2
          , stable//3
          , dirichlet//2
          , discrete//2
          , discrete//3
          , binomial//3
          , beta//3
          , zeta//2
          , gamma//2
          , inv_gamma//2
          , bernoulli//2
          , students_t//2
          , mixture//3
          , pair//3
          ]).

/** <module> Random predicates, untagged RNG state

This module provides a set of predicates for sampling from various
distributions. The state of the random generator is threaded through
using the DCG idiom.

@author Samer Abdallah
 */

:- module_transparent stream/2.

:- use_module(library(dcg_core)).
:- use_module(library(dcg_pair)).
:- use_module(library(plrand),[]).

term_expansion(stub(Arity,Name,Pred), Head --> plrand:Body) :-
   length(Args, Arity),
   Head =.. [Name | Args],
   Body =.. [Pred | Args].

%% bernoulli( +A:prob, -X:oneof([0,1]))// is det.
%  Sample binary random variable.
bernoulli(P,X) --> plrand:sample_Uniform01(U), {U<P->X=1;X=0}.

%% binomial( +P:float, +N:natural, -X:natural)// is det.
%
%  Sample X from a binomial distribution, ie the number of
%  successful trials out of N trials where the probability
%  of success of each trial is P.
stub(3,binomial,sample_Binomial).


%% poisson( +A:nonneg, -X:float)// is det.
%  Sample from Poisson distribution of rate A.
stub(2,poisson,sample_Poisson).

%% discrete( +A:list(prob), -X:natural)// is det.
%  Sample from a discrete distribution over natural numbers.
discrete(Ps,I) --> {length(Ps,N)}, plrand:sample_Discrete(N,Ps,I).

%% discrete( +O:list(T), +A:list(prob), -X:T)// is det.
%  Sample from a discrete distribution over list of objects.
discrete(Xs,Ps,X) --> {length(Ps,N)}, plrand:sample_Discrete(N,Ps,I), {nth1(I,Xs,X)}.

%% uniform01( -X:float)// is det.
%
%  Sample X from uniform distribution on [0,1).
stub(1,uniform01,sample_Uniform01).

%% normal( -X:float)// is det.
%  Sample from zero-mean unit-variance Gaussian.
stub(1,normal,sample_Normal).

%% exponential( -X:float)// is det.
%  Sample from unit-mean exponential distribution.
stub(1,exponential,sample_Exponential).

%% stable( +A, +B, -X:float)// is det.
%  Sample from a Levy-stable distribution.
stub(3,stable,sample_Stable).

%% dirichlet( +A:list(nonneg), -X:list(prob))// is det.
%  Sample from a Dirichlet distribution.
dirichlet(A,X) --> {length(A,N)}, plrand:sample_Dirichlet(N,A,X).


%% uniform( +Items:list(A), -A)// is det.
%
%  Uniform distribution over a finite number of items.
%  uniform :: list(A) -> expr(A).

uniform(O,X) --> 
   {length(O,N)},
   plrand:sample_Uniform01(U),
   {I is 1+floor(N*U), nth1(I,O,X)}.

%% uniformP(+P:dcg(-A), -A)// is det.
%  Sample uniformly from all solutions to call(P,X).
:- meta_predicate uniformP(3,-,+,-).
uniformP(P,X) -->
   {findall(Y,call(P,Y),YY)},
   uniform(YY,X).

%% beta( +A:nonneg, +B:nonneg, -X:prob)// is det.
%  Sample from beta distribution.
stub(3,beta,sample_Beta).

%% zeta( +A:nonneg, -X:natural)// is det.
%  Sample from zeta (hyperbolic or power law) distribution over natural numbers.
%  NB: Must have A > 1.
stub(2,zeta,sample_Zeta).

%% gamma( +A:nonneg, -X:float)// is det.
%  Sample from gamma distribution with parameter A.
stub(2,gamma,sample_Gamma).


% ^ above use plrand samplers and need randstate
% ---------------------- DERIVED DISTRIBUTIONS ---------------------
% V below do not use state directly.

%% gaussian( +Mean:float, +Var:nonneg, -X:float)// is det.
%  gaussian :: \(float, nonneg) -> expr(float).
%  Sample from Gaussian with given mean and variance.
gaussian(Mean, Var, X) --> normal(U), {X is Mean + Var*U}.

%% inv_gamma( +A:nonneg, -X:float)// is det.
%  Sample from inverse gamma distribution with parameter A.
inv_gamma(A,X)  --> gamma(A,Y), {X is 1/Y}.

%% pareto( +A:nonneg, -X:float)// is det.
%  Sample from pareto (power-law) distribution over non-negative reals.
pareto(A,X)    --> uniform01(Y), { X is (1-Y)**(-1/A) }.

%% students_t( +V:nonneg, -X:float)// is det.
%  Sample from student's t distribution with V degrees of freedom.
students_t(V,X)--> {V1 is V/2}, normal(Z), gamma(V1,Y), {X is Z*sqrt(V1/Y)}.

%% pair(+F:dist(A), +G:dist(B), -X:pair(A,B))// is det.
%
%  Sample a pair from two independent distributions.
pair(F,G,X-Y) --> call(F,X), call(G,Y).


%% mixture( +Sources:list(expr(A)), +Probs:list(prob), -X:A)// is det.
%
%  Sample from discrete distribution over Sources with probabilities Probs
%  and then sample from the resulting distribution.
%
%  mixture :: \(list(expr(A)), list(prob)) -> expr(A).
mixture( Sources, Dist, X) -->
   discrete(Dist,I),
   {nth1(Sources,I,S)},
   call(S,X).
