:- use_module(library(dcg_shell)).
:- use_module(library(dcg_pair)).
:- use_module(library(dcg_core)).
:- use_module(library(data/store)).
:- use_module(library(plrand)).
:- use_module(library(crp)).

/* 
   This is an example of how to build Dirichlet processes using the
   tools in plrand and crp. It works but setting up a DCG to manage
   the states of the CRPs and the state of the random generator.
   Once in the DCG (using dcgshell for an interactive DCG prompt) 
   you can add observations and sample conditionally on those observations. 
   For example: 

   ==
   ?- run.

	user: call_dcg (dcg) >> sample(dp(3,binomial(0.5,100)), DP), seqmap_n(9, sample(DP), X).                                                                                                                            

	  DP = crp(dp(3), binomial(0.5, 100), 13) 
	  X = [48, 54, 49, 42, 48, 48, 45, 52, 42] 

	user: call_dcg (dcg) >> sample(counter,C), sample(py(2,0.5,C), DP), seqmap_n(9, sample(DP), X).                                                                                                                     

	  C = counter(21) 
	  DP = crp(py(2, 0.5), counter(21), 22) 
	  X = [0, 1, 2, 1, 3, 4, 5, 5, 1] 
	==

   NB. it requires the dcgutils and genutils packs to be installed.
*/

run :- 
   store_new(Store), 
   with_rnd_state(run_left(dcgshell, Store, _)).

sample(dp(Conc, Base), crp(dp(Conc),Base,Ref)) -->
   {empty_classes(Classes)},
   \< store_add(Classes, Ref).

sample(py(Conc, Disc, Base), crp(py(Conc,Disc),Base,Ref)) -->
   {empty_classes(Classes)},
   \< store_add(Classes, Ref).

sample(counter, counter(Ref)) --> \< store_add(0, Ref).
sample(normal, X)        --> \> plrand:sample_Normal(X).
sample(exponential, X)   --> \> plrand:sample_Exponential(X).
sample(gamma(A), X)      --> \> plrand:sample_Gamma(A,X).
sample(zeta(S), X)       --> \> plrand:sample_Zeta(S,X).
sample(binomial(P,N), X) --> \> plrand:sample_Binomial(P,N,X).
sample(counter(Ref), X)  --> \< store_get(Ref, X), \< store_apply(Ref, succ).

sample(crp(GEM,Base,Ref), X) -->
   \< store_get(Ref, Classes),
   \> crp_sample(GEM, Classes, Action),
   sample_action(Action, Base, Ref, X).

sample_action(new, Base, Ref, X) --> 
   sample(Base, X),
   \< store_apply(Ref, add_class(X,_)).

sample_action(old(X,N), _, Ref, X) -->
   \< store_apply(Ref, inc_class(N)).

prob(normal,  _, 0) --> nop.
prob(exponential, _, 0) --> nop.
prob(gamma(_), _, 0) --> nop.
prob(zeta(S), X, P) --> {plrand:prob_Zeta(S,X,P)}.
prob(binomial(Q,N), X, P) --> {plrand:prob_Binomial(Q,N,X,P)}.
prob(counter(Ref), X, P) --> \< store_get(Ref,Y), {delta(X,Y,P)}. 
prob(crp(GEM,Base,Ref), X, P) -->
   prob(Base, X, P0), 
   \< store_get(Ref, Classes),
   { crp_prob(GEM, Classes, X, P0, P) }.

delta(X,X,1).
delta(X,Y,0) :- dif(X,Y).

observe(crp(GEM,Base,Ref), X) -->
   prob(Base, X, P),
   \< store_get(Ref, Classes),
   \> crp_sample_obs(GEM, Classes, X, P, Action),
	{writeln(observe_action(Action))},
   \< store_apply(Ref, obs_action(Action, X)).

obs_action(new, X) --> add_class(X, _).
obs_action(old(N), _) --> inc_class(N).

resample(CRP) -->
	{CRP=crp(_,_,Ref)},
	\< store_get(Ref, classes(_,Counts,Values)),
	seqmap(resample_class(CRP),Counts,Values).

resample_class(CRP, Count, X) --> rep(Count, resample_value(CRP, X)).
resample_value(CRP, X) --> {writeln(reseating(X))}, unobserve(CRP, X), observe(CRP, X).
unobserve(crp(_,_,Ref), X) -->
	\< store_get(Ref, Classes), 
	\> crp_sample_rm(Classes, X, N),
	\< store_apply(Ref, remove_from_class(N)).

remove_from_class(N) -->
	dec_class(N,C,_),
	({C>0}; {C=0}, remove_class(N)).
	
gem_samples(gamma(A,B), crp(dp(Conc), _, Ref), N, Concs) -->
   \< store_get(Ref, classes(_,Counts,_)), 
   {dp_sampler_teh(gamma(A,B), [Counts], Sampler)},
   \> collect_dp_samples(N, Sampler, Conc, Concs).

collect_dp_samples(0, _, _, []) --> [].
collect_dp_samples(N, Sampler, Conc0, [Conc1|Concs]) -->
   {succ(M,N)},
   call(Sampler, dp(Conc0), dp(Conc1)),
   collect_dp_samples(M, Sampler, Conc1, Concs).

