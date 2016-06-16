:- module(strand, [ strand/0
                  , strand/1
                  , clear//0
                  , hold_store//1
                  , (*)/4
                  , (*)//4
                  , constf//3
                  , pairf//3
                  , (*:)//3
                  , op(600,yfx,*:)
                  , pure//2
                  ]).

/** <module> Stateful random generation DCG
   
   This module provides tools for working in a sort of random generator plus
   store of mutable references monad (using DCG lanugage to manage state threading).
*/

:- meta_predicate strand(//)
                , hold_store(//,?,?)
                , *(2,2,?,?)
                , *(4,4,?,?,?,?)
                , constf(3,?,?,?,?)
                , pairf(3,3,?,?,?)
                , pure(3,-,+,-)
                .

:- module_transparent strand/0.

:- use_module(library(plrand)).
:- use_module(library(dcg_shell)).
:- use_module(library(dcg_core)).
:- use_module(library(dcg_pair)).
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
strand(Cmd) :- 
   store_new(H0),
   with_rnd_state(run_left(Cmd,H0,_)).

%% clear// is det.
%  Clear everything out of the store. Runs in strand DCG.
clear --> \< set_with(store_new).


%% hold_store(+Cmd:dcg(strand))// is det.
%  Runs Cmd leaving the store unchanged.
hold_store(Cmd) --> \< get(H), \> run_left(Cmd,H,_).

% pure and stateful function composition
*(P,Q,X,Z) --> call(Q,X,Y), call(P,Y,Z).
*(P,Q,X,Z) :- call(Q,X,Y), call(P,Y,Z).

% stateful piping of generator G into function P
*:(P,G,Y) --> call(G,X), call(P,X,Y).

pairf(F,G,X-Y) --> call(F,X), call(G,Y).
constf(F,_,X) --> call(F,X).

%% pure(+Dist:tagged_dist(A), +X:A, +P1:tagged_prob, -P2:tagged_prob) is det.
%% pure(+Dist:tagged_dist(A), -X:A, +S1:pair(store, rndstate), -S2:pair(store,rndstate)) is det.
%  Samples or gets probability of value X with tagged distribution Dist.
pure(Base,X,p(P1),p(P2)) :- call(Base,X,p(P1),p(P2)).
pure(Base,X,H-S1,H-S2) :- call(Base,X,rs(S1),rs(S2)).

