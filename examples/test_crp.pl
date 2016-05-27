:- use_module(library(plrand)).
:- use_module(library(dcg_shell)).
:- use_module(library(dcg_pair)).
:- use_module(library(dcg_core)).
:- use_module(library(crp)).

run :- 
   empty_classes(Classes),
   with_rnd_state(run_left(dcgshell, Classes, _)).

sample(GEM, Base, X) -->
   \< get(Classes),
   \> crp_sample(GEM, Classes, Action),
   sample_action(Action, Base, X).

prob(GEM, Base, X, P) -->
   \< get(Classes),
   { call(Base, X, P0), crp_prob(GEM, Classes, X, P0, P) }.


sample_action(new, Base, X) --> 
   \> call(Base, X),
   \< add_class(X, _).

sample_action(old(X,N), _, X) -->
   \< inc_class(N).

observe(GEM, Base, X) -->
   \< get(Classes),
   { call(Base, X, P) },
   \> crp_sample_obs(GEM, Classes, X, P, Action),
   \< obs_action(Action, X).

obs_action(new, X) --> add_class(X, _).
obs_action(old(N), _) --> inc_class(N).

