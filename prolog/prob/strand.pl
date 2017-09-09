:- module(strand, [ strand/0
                  , strand/1
                  , strand//1
                  , clear//0
                  , hold_store//1
                  , pure//2
                  , pure1//2
                  , marginal_prob//2
                  , marginal_prob//3
                  ]).

/** <module> Stateful random generation DCG
   
   This module provides tools for working in a sort of random generator plus
   store of mutable references monad (using DCG lanugage to manage state threading).

   @tbd I'm not happy with the state of this module (no pun intended... no wait a
   minute, that's a pretty good pun actually: the state representation _is_ one
   of the problems). The use of an untagged union is not good, and the stuff trying 
   to compute marginals doesn't really belong in a low level module like this.
*/

:- meta_predicate strand(//), strand(//,+,-)
                , hold_store(//,?,?)
                , pure(3,-,+,-)
                .

:- module_transparent strand/0.

:- use_module(library(plrand)).
:- use_module(library(dcg_shell)).
:- use_module(library(dcg_core)).
:- use_module(library(dcg_pair)).
:- use_module(library(callutils)).
:- use_module(library(data/store)).



%% strand is det.
%  Start a DCG shell (see dcg_shell:dcgshell//0) with an empty store
%  and the current state of the plrand random generator.
strand :- 
   context_module(M),
   strand:strand(strand:shell_in(M)).

shell_in(M, S1, S2) :- @(dcgshell(strand, S1, S2), M).

%% strand(+Cmd:dcg(strand)) is det.
%  Run Cmd in a DCG where the state is of type =|strand = pair(store,plrand:state)|=.
strand(Cmd) :- with_rnd_state(strand(Cmd)).

%% strand(+Cmd:dcg(strand),+RS1:(plrand:state),-RS2:(plrand:state)) is det.
%  Run Cmd in a DCG where the state is of type =|strand = pair(store,plrand:state)|=,
%  taking and returning initial and final states of plrand generator.
strand(Cmd) --> {store_new(H0)}, run_left(Cmd,H0,_).

%% clear// is det.
%  Clear everything out of the store. Runs in strand DCG.
clear --> \< set_with(store_new).

%% hold_store(+Cmd:dcg(strand))// is det.
%  Runs Cmd leaving the store unchanged.
hold_store(Cmd) --> \< get(H), \> run_left(Cmd,H,_).

%% pure(+Dist:tagged_dist(A), -X:A, +S1:pair(store,number), -S2:pair(store,number)) is det.
%% pure(+Dist:tagged_dist(A), -X:A, +S1:pair(store,rndstate), -S2:pair(store,rndstate)) is det.
%  Samples or gets probability of value X with tagged distribution Dist. NB: 2nd half of DCG
%  state type is an UNTAGGED union of numbers (probabilities) and random states.
pure(Base,X,H-P1,H-P2) :- number(P1), !, call(Base,X,p(P1),p(P2)).
pure(Base,X,H-S1,H-S2) :- call(Base,X,rs(S1),rs(S2)).

:- meta_predicate pure1(3,?,?,?).
pure1(Dist, X) --> \> call(Dist,X).

:- meta_predicate marginal_prob(//,-,+,-).
:- meta_predicate marginal_prob(3,?,-,+,-).
marginal_prob(G,Prob,S1-P1,S1-P1) :-
   aggregate(sum(P), S2^call_dcg(G, S1-p(1), S2-p(P)), Prob).

marginal_prob(G,X,Prob,S1-P1,S1-P1) :-
   aggregate(sum(P), S2^call(G, X, S1-p(1), S2-p(P)), Prob).
