:- module(crp,
		[	empty_classes/1
		,	dec_class//3
		,	inc_class//1
		,	remove_class//1
		,	add_class//2

		,	crp_prob/5
		,	crp_sample/5
		,	crp_sample/6
		,	crp_sample_obs/7
		,	crp_sample_rm/5

		,	dp_sampler_teh/3
		,	py_sampler_teh/4
		]).

/**	<module> Chinese Restaurant Process utilities
   
   This module provides some building blocks for implementing a family of random processes
   related to Dirichlet processes, including Pitman Yor processe, the Chinese Restaurant
   process, and the stick breaking model (GEM). The Dirichlet processes takes a single
   concentration parameter, representated as =|dp(Conc)|=, while the Pitman Yor process
   takes a concentration parameter and a discount parameter, representated as =|py(Conc,Disc)|=.

	==
	gem_param   ---> dp(nonneg) ; py(nonneg,0--1).
	gamma_prior ---> gamma(nonneg, nonneg).
	beta_prior  ---> beta(nonneg, nonneg).
   classes(A)  ---> classes(natural, list(nonneg), list(A)).
   action(A)   ---> new ; old(A, class_idx).
   action      ---> new ; old(class_idx).

   rndstate  == plrand:state
   class_idx == natural
   prob      == 0--1 

	param_sampler == pred(+gem_param, -gem_param, +rndstate, -rndstate).
	==

   This may seems like a very low-level library for building CRPs, leaving a lot
   for the implemeenter to do, but this is intential, to allow the implementer
   freedom to decide how to manage the states (terms of type classes(_)) of one
   or more CRPs, as well as the state of the random generator, in whatever way
   is most appropriate. See the the example implementation of test_crp.pl for
   one way to do this.
*/

% :- use_module(library(dcg_core)).
% :- use_module(library(dcg_macros)).
:- use_module(library(apply_macros)).
:- use_module(library(plrand),   [spawn/3, crp_prob/5, crp_sample/5, crp_sample/6, crp_sample_obs/7, crp_sample_rm/5]).


%% crp_prob( +GEM:gem_param, +Classes:classes(A), +X:A, +PBase:prob, -Prob:prob) is det.
%
%  Compute the probability Prob of observing X given a CRP with already observed
%  values in Classes if the probability of drawing X from the base distribution is PBase.


%% crp_sample( +GEM:gem_param, +Classes:classes(A), -A:action(A), -P:prob)// is det.
%% crp_sample( +GEM:gem_param, +Classes:classes(A), -A:action(A))// is det.
%
%  Sample a new value from CRP, Action A is either new, which means
%  that the user should sample a new value from the base distribtion,
%  or old(X,ID), where X is an old value and C is the class index.
%  Operates in random state DCG. crp_sample//4 additionally returns the probability
%  of the action choosen.
crp_sample(GEM, Classes, Action) --> crp_sample(GEM, Classes, Action, _).


%% crp_sample_obs( +GEM:gem_param, +Classes:classes(A), +X:A, +PBase:prob, -A:action, -P:prob)// is det.
%% crp_sample_obs( +GEM:gem_param, +Classes:classes(A), +X:A, +PBase:prob, -A:action)// is det.
%
%  Sample action appropriate for observation of value X. PBase is the
%  probability of X from the base distribution. Action A is new
%  or old(N) where N is the class index. crp_sample_obs//6 additionally returns the
%  probability of the observation, equivalent to calling crp_prob with X BEFORE
%  calling crp_sample_obs//5.
%  Operates in random state DCG.
crp_sample_obs(GEM, Classes, X, PBase, Action) --> 
   crp_sample(GEM, Classes, X, PBase, Action, _).


%% crp_sample_rm( +Classes:classes(A), +X:A, -N:class_idx)// is det.
%
%  Sample appropriate class index N from which to remove value X.
%  Operates in random state DCG.



% --------------------------------------------------------------------------------
% classes data structure (basic CRP stuff)


%% empty_classes( -Classes:classes(_)) is det.
%
%  Unify Classes with an empty classes structure.
empty_classes(classes(0,[],[])).


%% dec_class( +N:class_idx, -C:natural, -X:A, +C1:classes(A), -C2:classes(A)) is det.
%
%  Decrement count associated with class id N. C is the count after
%  decrementing and X is the value associated with the Nth class.
dec_class(N,CI,X,classes(K,C1,Vs),classes(K,C2,Vs)) :- dec_nth(N,CI,C1,C2), nth1(N,Vs,X).
dec_nth(1,Y,[X|T],[Y|T]) :- succ(Y,X).
dec_nth(N,Y,[X|T1],[X|T2]) :- succ(M,N), dec_nth(M,Y,T1,T2).

%% inc_class( +N:class_idx, +C1:classes(A), -C2:classes(A)) is det.
%
%  Increment count associated with class N.
inc_class(N,classes(K,C1,V),classes(K,C2,V)) :- inc_nth(N,C1,C2).
inc_nth(1,[X|T],[Y|T]) :- succ(X,Y).
inc_nth(N,[X|T1],[X|T2]) :- succ(M,N), inc_nth(M,T1,T2).


%% remove_class( +N:class_idx, +C1:classes(A), -C2:classes(A)) is det.
%
%  Removes Nth class.
remove_class(N,classes(K1,C1,V1),classes(K2,C2,V2)) :-
	remove_nth(N,C1,C2),
	remove_nth(N,V1,V2),
	succ(K2,K1).

%% add_class( +X:A, -ID:class_idx, +C1:classes(A), -C2:classes(A)) is det.
%
%  Add a class associated with value X. N is the index of the new class.
add_class(X,N2,classes(N1,C1,V1),classes(N2,C2,V2)) :-
	succ(N1,N2),
	append(C1,[1],C2),
	append(V1,[X],V2). 


remove_nth(1,[_|T],T).
remove_nth(N,[Y|T1],[Y|T2]) :- 
	(	var(N) 
	->	remove_nth(M,T1,T2), succ(M,N)
	;	succ(M,N), remove_nth(M,T1,T2)
	).


% ---------------------------------------------------------------
% PARAMETER SAMPLING
% Initialisers in Prolog, samplers written in C.

%% dp_sampler_teh( +Prior:gamma_prior, +Counts:list(natural), -S:param_sampler) is det.
%
%	Prepares a predicate for sampling the concentration parameter of a Dirichlet process.
%	The sampler's =|gem_prior|= arguments must be of the form =|dp(_)|=.
%	Prior specifies the Gamma distribution prior for the concentration parameter,
%	as =|gamma(a,b)|=, where a is the shape parameter and b is the rate parameter
%	(ie the inverse of the scale parameter).
dp_sampler_teh( gamma(A,B), CX, plrand:sample_dp_teh(ApSumKX,B,NX)) :-
	maplist(sumlist,CX,NX),
	maplist(length,CX,KX), 
	sumlist(KX,SumKX), 
	ApSumKX is A+SumKX.

%% py_sampler_teh( +ConcPrior:gamma_prior, +DiscPr:beta_prior, +Counts:list(natural), -S:param_sampler) is det.
%
%	Prepares a predicate for sampling the concentration and discount 
%	parameters of a Pitman-Yor process.
%	The sampler's =|gem_prior|= arguments must be of the form =|py(_,_)|=.
%	See dp_sampler_teh/3 for tha description of the gamma_prior type. DiscPr is a 
%	Beta distribution prior for the concentration parameter.
py_sampler_teh( ThPrior, DiscPrior, CountsX, Sampler) :-
	Sampler = plrand:sample_py_teh( ThPrior, DiscPrior, CountsX).


/*  ---------------------------------------------------------------------------------------i
    PURE PROLOG VERSIONS (not including MCMC parameter samplers)

:- use_module(library(math), [sub/3, equal/3, stoch/3, mul/3]).

sample_discrete(O,P,Y) --> {length(P,N)}, plrand:sample_Discrete(N,P,X), {nth1(X,O,Y)}.

crp_prob( Alpha, classes(_,Counts,Vals), X, PProb, P) :-
	counts_dist( Alpha, Counts, Counts1),
	stoch( Counts1, Probs, _),
	maplist( equal(X), Vals, Mask),
	maplist( mul, [PProb | Mask], Probs, PostProbs),
	sumlist( PostProbs, P).

%% crp_sample( +GEM:gem_model, +Classes:classes(A), -A:action(A))// is det.
%
%  Sample a new value from CRP, Action A is either new, which means
%  that the user should sample a new value from the base distribtion,
%  or old(X,C), where X is an old value and C is the index of its class.
%  Operates in random state DCG.
crp_sample( Alpha, classes(_,Counts,Vals), Action, RS1, RS2) :-
	counts_dist(Alpha, Counts, Counts1),
	sample_discrete(Counts1,Z,RS1,RS2),
	( Z>1 -> succ(C,Z), nth1(C,Vals,X), Action=old(X,C)
	; Action=new).


%% crp_sample_obs( +GEM:gem_model, +Classes:classes(A), +X:A, +PProb:float, -A:action)// is det.
%
%  Sample class appropriate for observation of value X. PProb is the
%  base probability of X from the base distribution. Action A is new
%  or old(Class).
%  Operates in random state DCG.
crp_sample_obs( Alpha, classes(_,Counts,Vals), X, ProbX, A, RS1, RS2) :-
	counts_dist( Alpha, Counts, [CNew|Counts1]),	
	PNew is CNew*ProbX,
	maplist( posterior_count(X),Vals,Counts1,Counts2),
	sample_discrete( [PNew|Counts2], Z, RS1, RS2),
	(Z=1 -> A=new; succ(C,Z), A=old(C)).


%% crp_sample_rm( +Classes:classes(A), +X:A, -C:natural)// is det.
%
%  Sample appropriate class from which to remove value X.
%  Operates in random state DCG.
crp_sample_rm( classes(_,Counts,Vals), X, Class, RS1, RS2) :-
	maplist(posterior_count(X),Vals,Counts,Counts1),
	sample_discrete( Counts1, Class, RS1, RS2).


posterior_count(X,Val,Count,PC) :- X=Val -> PC=Count; PC=0.

% -----------------------------------------------------------
% Dirichlet process and Pitman-Yor process
% pseudo-counts models.

counts_dist(dp(Alpha),Counts,[Alpha|Counts]) :- !.
counts_dist(py(_,_),[],[1]) :- !.
counts_dist(py(Alpha,Discount),Counts,[CNew|Counts1]) :- !,
	length(Counts,K),
	CNew is Alpha+Discount*K,
	maplist(sub(Discount),Counts,Counts1).

*/
