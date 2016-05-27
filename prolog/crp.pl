:- module(crp,
		[	empty_classes/1
		,	dec_class//3
		,	inc_class//1
		,	remove_class//1
		,	add_class//2

		,	crp_prob/5
		,	crp_sample/5
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
	gem_model   ---> dp(nonneg) ; py(nonneg,0--1).
	gamma_prior ---> gamma(nonneg, nonneg).
	beta_prior  ---> beta(nonneg, nonneg).
   classes(A)  ---> classes(natural, list(nonneg), list(A)).
   action(A)   ---> new ; old(A, class_idx).
   action      ---> new ; old(class_idx).

   rndstate  == plrand:state
   class_idx == natural
   prob      == 0--1 

	param_sampler == pred(+gem_model, -gem_model, +rndstate, -rndstate).
	==

*/

% :- use_module(library(dcg_core)).
% :- use_module(library(dcg_macros)).
:- use_module(library(apply_macros)).
:- use_module(library(plrand),   [spawn/3, crp_prob/5, crp_sample/5, crp_sample_obs/7, crp_sample_rm/5]).


%% crp_prob( +GEM:gem_model, +Classes:classes(A), +X:A, +PBase:prob, -Prob:prob) is det.
%
%  Compute the probability Prob of observing X given a CRP with already observed
%  values in Classes if the probability of drawing X from the base distribution is PBase.


%% crp_sample( +GEM:gem_model, +Classes:classes(A), -A:action(A))// is det.
%
%  Sample a new value from CRP, Action A is either new, which means
%  that the user should sample a new value from the base distribtion,
%  or old(X,ID), where X is an old value and C is the class index.
%  Operates in random state DCG.


%% crp_sample_obs( +GEM:gem_model, +Classes:classes(A), +X:A, +PBase:prob, -A:action)// is det.
%
%  Sample action appropriate for observation of value X. PBase is the
%  probability of X from the base distribution. Action A is new
%  or old(N) where N is the class index.
%  Operates in random state DCG.


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
dp_sampler_teh( gamma(A,B), CX, crp:sample_dp_teh(ApSumKX,B,NX)) :-
	maplist(sumlist,CX,NX),
	maplist(length,CX,KX), 
	sumlist(KX,SumKX), 
	ApSumKX is A+SumKX.

%% py_sampler_teh( +ThPrior:gamma_prior, +DiscPr:beta_prior, +Counts:list(natural), -S:param_sampler) is det.
%
%	Prepares a predicate for sampling the concentration and discount 
%	parameters of a Pitman-Yor process.
%	The sampler's =|gem_prior|= arguments must be of the form =|py(_,_)|=.
py_sampler_teh( ThPrior, DiscPrior, CountsX, crp:Sampler) :-
	Sampler = sample_py_teh( ThPrior, DiscPrior, CountsX).

