:- module(rmemo,
          [ gdpmem//3
          , gdp//3
          , memo//2
          , memo_lookup//3
          , gdp_info//3
          , gdpmem_info//4
          ]).

/** <module> Memoisation and stochastic memoisation of random predicates built on strand

   This module provides memoisation services for predicates that run in the =|strand|=
   store + random state DCG. Memoisation of sampling predicates can be used to do lazy
   sampling and implement random world semantics. Stochastic memoisation is used to
   implement Dirichlet Processes and Pitman Yor processes.
*/

:- use_module(library(dcg_pair)).
:- use_module(library(dcg_macros)).
:- use_module(library(data/store), [store_add//2, store_get//2, store_apply//2]).
:- use_module(library(prob/crp),   [empty_classes/1, crp_sample//3, add_class//2, inc_class//1]).
:- use_module(library(prob/strand), [(*)//4]).

%% memo(+F:dcg(A,B,strand), -G:dcg(A,B,strand))// is det.
%
%  Runs in strand DCG. Produces G, a callable which is a memoised version of F,
%  both of which are binary DCG goals.

:- meta_predicate memo(4,-,+,-).
memo(F,rmemo:memf(TabRef,F)) -->
   {empty_assoc(Tab)},
   \< store_add(Tab,TabRef).

memf(TabRef,F,X,Y) -->
   \< store_get(TabRef,T1),
   (  {get_assoc(X,T1,Y)} -> []
   ;  call(F,X,Y),
      \< store_apply(TabRef,assoc_put(X,Y))
   ).

assoc_put(K,V,T1,T2) :- put_assoc(K,T1,V,T2).

%% memo_lookup(+G:dcg(A,B,strand), -X:A, -Y:B)// is nondet.
%
%  Looks up previously memoised computations of G.
memo_lookup(rmemo:memf(TabRef,_), X, Y) -->
   \< store_get(TabRef,T),
   {gen_assoc(X,T,Y)}.

%% gdp(+GEM:gem_param, +Base:dcg(-A,strand), -CRP:dcg(-A,strand))// is det.
%
%  Generalised Dirichlet process: builds a CRP (so-called Chinese Restaurant
%  Process) representing the distribution obtained by sampling from a Dirichlet
%  or Pitman-Yor process with given parameter and base distribution.
:- meta_predicate gdp(+,3,-,+,-).
gdp(GEM,H,rmemo:crp(Ref,GEM,H)) -->
   {empty_classes(Classes)},
   \< store_add(Classes,Ref).

crp(Ref,GEM,H,X) -->
   \< store_get(Ref,Classes),
   \> crp_sample(GEM,Classes,Action),
   crp_action(Action,Ref,H,X).

crp_action(new,Ref,H,X) -->
   call(H,X),
   \< store_apply(Ref,add_class(X,_)).

crp_action(old(X,Idx),Ref,_,X) -->
   \< store_apply(Ref,inc_class(Idx)).

gdp_info(rmemo:crp(Ref,GEM,_),GEM, Classes) -->
   \< store_get(Ref,Classes).

%% gdpmem(+GEM:gem_param, +Base:dcg(+B,-A,strand), -CRP:dcg(+B,-A,strand))// is det.
%
%  Stochastic memosation based on generalised Dirichlet processes: Base is binary DCG goal
%  in strand such that call(Base,X,Y) samples Y from some distribtion depending on X.
%  gdpmem//3 produces a stochastically memoised version of Base that defines a GDP for
%  each distinct value of X, created on demand the first time that X is supplied.
%  See gdp//3 for more info.
:- meta_predicate gdpmem(+,4,-,+,-).
gdpmem(GEM,F,rmemo:(call*G)) --> memo(new_gdp(GEM,F),G).
new_gdp(GEM,F,X,GX) --> gdp(GEM,call(F,X),GX).

gdpmem_info(rmemo:(_*G), X, GEM, Classes) -->
   memo_lookup(G, X, CRP),
   gdp_info(CRP, GEM, Classes).
