:- module(prob_tagged,
          [ uniform01//1
          , uniform//2, uniformT//2, uniformP//2
          , normal//1
          , gaussian//3
          , exponential//1
          , poisson//2
          , stable//3
          , dirichlet//2
          , discrete//2
          , discrete//3, discreteT//3
          , binomial//3
          , beta//3
          , zeta//2
          , gamma//2
          , inv_gamma//2
          , bernoulli//2
          , students_t//2
          , mixture//3
          , product_pair//3
          , prob/3, prob/2
          , pdf/3, pdf/2
          ]).

/** <module> Random predicates

This module provides a set of predicates for sampling from various
distributions. The state of the random generator is threaded through
using the DCG idiom.

@author Samer Abdallah
 */

:- module_transparent stream/2.

:- use_module(library(dcg_core)).
:- use_module(library(dcg_pair)).
:- use_module(library(plrand),[]).

term_expansion(wrap_rs(Arity,Name,Pred), Head :- plrand:Body) :-
   length(A1, Arity),
   append(A1, [rs(S1),rs(S2)], HeadArgs), Head =.. [Name | HeadArgs],
   append(A1, [S1,S2],         PredArgs), Body =.. [Pred | PredArgs].

term_expansion(wrap_disc(Arity,Name,Pred), Head :- (plrand:Body, P2 is P1*P)) :-
   length(A1, Arity),
   append(A1, [p(P1),p(P2)], HeadArgs), Head =.. [Name | HeadArgs],
   append(A1, [P],           PredArgs), Body =.. [Pred | PredArgs].

term_expansion(wrap_cont(Arity,Name,Pred), [Head1, (Head :- (plrand:Body, P2 is P1*P))]) :-
   length(A1, Arity),
   append(A1, [p(_),p(0)],     Head1Args), Head1 =.. [Name | Head1Args],
   append(A1, [pd(P1),pd(P2)], HeadArgs), Head =.. [Name | HeadArgs],
   append(A1, [P],             PredArgs), Body =.. [Pred | PredArgs].

%% bernoulli( +A:prob, -X:oneof([0,1]))// is det.
%  Sample binary random variable.
bernoulli(P,X,rs(S1),rs(S2)) :- !, plrand:sample_Uniform01(U,S1,S2), (U<P->X=1;X=0).
bernoulli(P,0,P1,P2) :- !, P2 is (1-P)*P1.
bernoulli(P,1,P1,P2) :- P2 is P*P1.

%% binomial( +P:float, +N:natural, -X:natural)// is det.
%
%  Sample X from a binomial distribution, ie the number of
%  successful trials out of N trials where the probability
%  of success of each trial is P.
wrap_rs(3,binomial,sample_Binomial).
wrap_disc(3,binomial,prob_Binomial).


%% poisson( +A:nonneg, -X:float)// is det.
%  Sample from Poisson distribution of rate A.
wrap_rs(2,poisson,sample_Poisson).
wrap_disc(2,poisson,prob_Poisson).

%% discrete( +A:list(prob), -X:natural)// is det.
%  Sample from a discrete distribution over natural numbers.
discrete(Ps,I,p(P1),p(P2))   :- !, nth1(I,Ps,P), P2 is P1*P.
discrete(Ps,I,rs(S1),rs(S2)) :- length(Ps,N), plrand:sample_Discrete(N,Ps,I,S1,S2).


%% discrete( +O:list(T), +A:list(prob), -X:T)// is det.
%  Sample from a discrete distribution over list of objects.
discrete(Xs,Ps,X,rs(S1),rs(S2)) :- !, length(Ps,N), plrand:sample_Discrete(N,Ps,I,S1,S2), nth1(I,Xs,X).
discrete(Xs,Ps,X,p(P1),p(P2))   :-
   aggregate(sum(P), I^(nth1(I,Xs,X), nth1(I,Ps,P)), Prob),
   P2 is P1*Prob.

discreteT(Xs,Ps,X,rs(S1),rs(S2)) :- !, functor(Ps,_,N), plrand:sample_DiscreteF(N,Ps,I,S1,S2), arg(I,Xs,X).
discreteT(Xs,Ps,X,p(P1),p(P2))   :-
   aggregate(sum(P), I^(arg(I,Xs,X), arg(I,Ps,P)), Prob),
   P2 is P1*Prob.

%% uniform01( -X:float)// is det.
%
%  Sample X from uniform distribution on [0,1).
uniform01(_,p(_),p(0)) :- !.
uniform01(_,pd(P),pd(P)) :- !.
wrap_rs(1,uniform01,sample_Uniform01).

%% normal( -X:float)// is det.
%  Sample from zero-mean unit-variance Gaussian.
normal(_,p(_),p(0)) :- !.
wrap_rs(1,normal,sample_Normal).
wrap_cont(1,normal,prob_Normal).

%% exponential( -X:float)// is det.
%  Sample from unit-mean exponential distribution.
wrap_rs(1,exponential,sample_Exponential).
wrap_cont(1,exponential,prob_Exponential).


%% stable( +A, +B, -X:float)// is det.
%  Sample from a Levy-stable distribution.
wrap_rs(3,stable,sample_Stable).

%% dirichlet( +A:list(nonneg), -X:list(prob))// is det.
%  Sample from a Dirichlet distribution.
dirichlet(A,X,rs(S1),rs(S2)) :- !, length(A,N), plrand:sample_Dirichlet(N,A,X,S1,S2).
dirichlet(A,X,p(P1),p(P2))   :- length(A,N), plrand:prob_Dirichlet(N,A,X,P), P2 is P1*P.


%% uniform( +Items:list(A), -A)// is det.
%
%  Uniform distribution over a finite number of items.
%  uniform :: list(A) -> expr(A).

uniform(O,X,rs(S1),rs(S2)) :- !,
   length(O,N),
   plrand:sample_Uniform01(U,S1,S2),
   I is 1+floor(N*U), nth1(I,O,X).
uniform(O,X,p(P1),p(P2)) :-
   length(O,N),
   aggregate(count,member(X,O),K),
   P2 is P1*K/N.

uniformT(O,X,p(P1),p(P2)) :- !,
   functor(O,_,N),
   aggregate(count,I^arg(I,O,X),K),
   P2 is K*P1/N.

uniformT(O,X,rs(S1),rs(S2)) :-
   functor(O,_,N),
   plrand:sample_Uniform01(U,S1,S2),
   I is 1+floor(N*U), arg(I,O,X).

%% uniformP(+P:dcg(-A), -A)// is det.
%  Sample uniformly from all solutions to call(P,X).
:- meta_predicate uniformP(3,-,+,-).
uniformP(P,X) -->
   {findall(Y,call(P,Y),YY)},
   uniform(YY,X).

%% beta( +A:nonneg, +B:nonneg, -X:prob)// is det.
%  Sample from beta distribution.
wrap_rs(3,beta,sample_Beta).
wrap_cont(3,beta,prob_Beta).

%% zeta( +A:nonneg, -X:natural)// is det.
%  Sample from zeta (hyperbolic or power law) distribution over natural numbers.
%  NB: Must have A > 1.
wrap_rs(2,zeta,sample_Zeta).
wrap_disc(2,zeta,prob_Zeta).
% zeta(A,X,)      --> {A>1}, sample_Zeta(A,X1), { X is integer(X1) }.

%% gamma( +A:nonneg, -X:float)// is det.
%  Sample from gamma distribution with parameter A.
wrap_rs(2,gamma,sample_Gamma).
wrap_cont(2,gamma,prob_Gamma).


% ^ above use plrand samplers and need randstate
% ---------------------- DERIVED DISTRIBUTIONS ---------------------
% V below do not use state directly.

%% gaussian( +Mean:float, +Var:nonneg, -X:float)// is det.
%  gaussian :: \(float, nonneg) -> expr(float).
%  Sample from Gaussian with given mean and variance.
gaussian( Mean, Var, X) --> normal(U), {X is Mean + Var*U}.

%% inv_gamma( +A:nonneg, -X:float)// is det.
%  Sample from inverse gamma distribution with parameter A.
inv_gamma(A,X)  --> gamma(A,Y), {X is 1/Y}.

%% pareto( +A:nonneg, -X:float)// is det.
%  Sample from pareto (power-law) distribution over non-negative reals.
pareto(A,X)    --> uniform01(Y), { X is (1-Y)**(-1/A) }.

%% students_t( +V:nonneg, -X:float)// is det.
%  Sample from student's t distribution with V degrees of freedom.
students_t(V,X)--> {V1 is V/2}, normal(Z), gamma(V1,Y), {X is Z*sqrt(V1/Y)}.

%% product_pair(+F:dist(A), +G:dist(B), -X:pair(A,B))// is det.
%
%  Sample a pair from two independent distributions.
product_pair(F,G,X-Y) --> call(F,X), call(G,Y).


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

prob(Expr,Val,Prob) :- call(Expr,Val,p(1),p(Prob)).
pdf(Expr,Val,Prob) :- call(Expr,Val,pd(1),pd(Prob)).
prob(Goal,Prob) :- call_dcg(Goal,p(1),p(Prob)).
pdf(Goal,Prob)  :- call_dcg(Goal,pd(1),pd(Prob)).

