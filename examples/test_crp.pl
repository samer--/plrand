:- use_module(library(plrand)).
:- use_module(library(dcg_shell)).
:- use_module(library(dcg_pair)).
:- use_module(library(dcg_core)).
:- use_module(library(crp)).

run :- 
   empty_classes(Classes),
   with_rnd_state(run_left(dcgshell,Classes,_)).

sample(P,Base,X) -->
   \< get(Classes),
   \> crp_sample(P,Classes,X).


