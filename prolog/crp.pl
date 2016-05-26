:- module(crp,
		[	empty_classes/1
		,	classes_value/2
		,	classes_counts/2
		,	classes_update/3
		,	seqmap_classes//2
		,	dec_class//3
		,	inc_class//1
		,	remove_class//1
		,	add_class//2

		,	crp_prob/5
		,	crp_sample/5
		,	crp_sample_obs/7
		,	crp_sample_rm/5
		,	crp_dist/6

		,	dp_sampler_teh/3
		,	py_sampler_teh/4
		]).

/**	<module> Chinese Restaurant Process utilities

	==
	gem_model ---> dp(Alpha:nonneg)
	             ; py(Alpha:nonneg,Discount:nonneg).

	gamma_prior ---> gamma(nonneg,nonneg).
	beta_prior  ---> beta(nonneg,nonneg).
	param_sampler == pred(+gem_model,-gem_model,+rndstate,-rndstate).
	==

*/

:- meta_predicate seqmap_classes(4,+,?,?).

:- use_module(library(dcg_core)).
:- use_module(library(dcg_pair)).
:- use_module(library(dcg_macros)).
:- use_module(library(apply_macros)).
:- use_module(library(plrand),   [spawn/3, crp_prob/5, crp_sample/5, crp_sample_obs/7, crp_sample_rm/5]).
:- use_module(library(lazy),     [lazy_unfold/4, lazy_unfold/5]).
:- use_module(library(randpred), [dirichlet//2, beta//3]).
:- add_import_module(randpred, plrand, start).


%% crp_prob( +GEM:gem_model, +Classes:classes(A), +X:A, +PProb:float, -Prob:float) is det.
%
%  Compute the probability Prob of observing X given a CRP
%  and a base probability of PProb.


%% crp_sample( +GEM:gem_model, +Classes:classes(A), -A:action(A))// is det.
%
%  Sample a new value from CRP, Action A is either new, which means
%  that the user should sample a new value from the base distribtion,
%  or old(X,ID), where X is an old value and C is the class ID.
%  Operates in random state DCG.


%% crp_sample_obs( +GEM:gem_model, +Classes:classes(A), +X:A, +PProb:float, -A:action)// is det.
%
%  Sample class appropriate for observation of value X. PProb is the
%  base probability of X from the base distribution. Action A is new
%  or old(ID) where ID is the class id.
%  Operates in random state DCG.


%% crp_sample_rm( +Classes:classes(A), +X:A, -C:class_id)// is det.
%
%  Sample appropriate class from which to remove value X. C is the
%  class id of the chosen class.
%  Operates in random state DCG.


%% crp_dist( +GEM:gem_model, +Classes:classes(A), +Base:dist(A), -Dist:dist(A))// is det.
%
%  Get posterior distribution associated with node using stick breaking method.
%  Operates in random state DCG.
crp_dist( dp(Alpha), classes(_,Counts,Values), Base, Dist, RS1, RS3) :-
	sumlist(Counts,Total), 
	Norm is Total+Alpha,

	(	Total>0
	-> dirichlet(Counts,Probs1, RS1, RS2),
		lazy_dp(Alpha,Base,Alpha,ValuesT,ProbsT, RS2, RS3),
		maplist(mul(Total),Probs1,Probs2),
		append(Probs2,ProbsT,ProbsA),
		append(Values,ValuesT,ValuesA),
		Dist=lazy_discrete(ValuesA,ProbsA,Norm)
	;	lazy_dp(Alpha, Base, 1, ValuesT, ProbsT, RS1, RS3),
		Dist=lazy_discrete(ValuesT,ProbsT,1)
	).
	

% --------------------------------------------------------------------------------
% classes data structure (basic CRP stuff)

user:portray(classes(_,Counts,Vals)) :- format('<crp|~p:~p>',[Counts,Vals]).


%% empty_classes( -Classes:classes(_)) is det.
%
%  Unify Classes with an empty classes structure.
empty_classes(classes(0,[],[])).


%% classes_value( +Classes:classes(A), +X:A) is semidet.
%% classes_value( +Classes:classes(A), -X:A) is multi.
%
%  Check that X is one of the values represented in Classes.
%  If X is unbound on entry, it is unified with all values on backtracking.
classes_value(classes(_,_,Vals),X) :- member(X,Vals).


%% classes_counts( +Classes:classes(A), -Counts:list(natural)) is det.
%
%  Gets the list of counts, one per class.
classes_counts( classes(_,Counts,_), Counts).

%% seqmap_classes( +P:pred(natural,A,T,T), +Classes:classes(A), +S1:T, -S2:T) is multi.
%
%  Sequentiall apply phrase P to all classes. Arguments to P are the number of items
%  in the class and the value (of type A) associated with it.
seqmap_classes(P, classes(_,Counts,Vals)) --> seqmap( P, Counts, Vals).

user:goal_expansion(seqmap_classes(P,CX,S1,S2), (CX=classes(_,Counts,Vals), seqmap(P, Counts,Vals,S1,S2))).

%% dec_class( +ID:class_id, -C:natural, -X:A, +C1:classes(A), -C2:classes(A)) is det.
%
%  Decrement count associated with class id N. C is the count after
%  decrementing and X is the value associated with the Nth class.
dec_class(N,CI,X,classes(K,C1,V),classes(K,C2,V)) :- dec_nth(N,_,CI,C1,C2), nth1(N,V,X).
dec_nth(1,X,Y,[X|T],[Y|T]) :- succ(Y,X).
dec_nth(N,A,B,[X|T1],[X|T2]) :- succ(M,N), dec_nth(M,A,B,T1,T2).

%% inc_class( +ID:class_id, +C1:classes(A), -C2:classes(A)) is det.
%
%  Increment count associated with class N.
inc_class(C,classes(K,C1,V),classes(K,C2,V)) :- inc_nth(C,C1,C2).
inc_nth(1,[X|T],[Y|T]) :- succ(X,Y).
inc_nth(N,[X|T1],[X|T2]) :- succ(M,N), inc_nth(M,T1,T2).


%% remove_class( +ID:class_id, +C1:classes(A), -C2:classes(A)) is det.
%
%  Removes class N.
remove_class(I,classes(K1,C1,V1),classes(K2,C2,V2)) :-
	remove_from_list(I,_,C1,C2),
	remove_from_list(I,_,V1,V2),
	succ(K2,K1).

%% add_class( +X:A, -ID:class_id, +C1:classes(A), -C2:classes(A)) is det.
%
%  Add a class associated with value X. N is the id of the new class.
add_class(X,K2,classes(K1,C1,V1),classes(K2,C2,V2)) :-
	succ(K1,K2),
	append(C1,[1],C2),
	append(V1,[X],V2). 


remove_from_list(1,X,[X|T],T).
remove_from_list(N,X,[Y|T1],[Y|T2]) :- 
	(	var(N) 
	->	remove_from_list(M,X,T1,T2), succ(M,N)
	;	succ(M,N), remove_from_list(M,X,T1,T2)
	).


%------------------------------------------------------------------
% Get posterior distribution at node using stick-breaking
% construction.

lazy_dp(A,H,P0,Vals,Probs) -->
	spawn(S0), { lazy_unfold(unfold_dp(A,H),Vals,Probs,P0-S0,_) }.

lazy_dp_paired(A,H,P0,ValsProbs) -->
	spawn(S0), { lazy_unfold(unfold_dp(A,H),ValsProbs,P0-S0,_) }.

unfold_dp(A,H,V,X) --> \> call(H,V), unfold_gem(A,X).
unfold_dp(A,H,V:X) --> \> call(H,V), unfold_gem(A,X).

% lazy_gem(A,Probs) --> spawn(S0), { lazy_unfold(unfold_gem(A),Probs,(1,S0),_) }.

unfold_gem(A,X) -->
	\> beta(1,A,P),
	\< trans(P0,P1),
	{ X is P*P0, P1 is P0-X }.

%% classes_update( +Action:action(A), +C1:classes(A), -C2:classes(A)) is det.
%
%	Update classes structure with a new observation.
classes_update(old(_,ID),C1,C2) :- inc_class(ID,C1,C2).
classes_update(new(X,ID),C1,C2) :- add_class(X,ID,C1,C2).




% PARAMETER SAMPLING



% ---------------------------------------------------------------
% Initialisers
% Samplers written in C.

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
%	The sampler's =|gem_prior|= arguments must be of the form =|dp(_)|=.
py_sampler_teh( ThPrior, DiscPrior, CountsX, crp:Sampler) :-
	Sampler = sample_py_teh( ThPrior, DiscPrior, CountsX).

/*
slow_sample_py_teh( gamma(A,B), beta(DA,DB), CountsX, py(Theta1,Disc1), py(Theta2,Disc2)) -->
	% do several lots of sampling auxillary variables, one per client node
	%	seqmap( py_sample_s_z_w(Theta1,Disc1), CountsX, SX, NSX, ZX, WX),
	seqmap( py_sample_s_z_log_w(Theta1,Disc1), CountsX, SX, NSX, ZX, LogWX),
	{	% maplist(log,WX,LogWX),
		sumlist(SX,SumSX),
		sumlist(NSX,SumNSX),
		sumlist(ZX,SumZX),
		sumlist(LogWX,SumLogWX),
		A1  is A+SumSX,   B1 is B-SumLogWX,
		DA1 is DA+SumNSX, DB1 is DB+SumZX },
	gamma(A1, B1, Theta2), 
	beta(DA1, DB1, Disc2).

py_sample_s_z_w(Theta,Disc,Counts,S,NS,Z,W) -->
	py_sample_bern_z(Disc,Counts,Z),
	py_sample_bern_s(Theta,Disc,Counts,S,NS),
	py_sample_beta_w(Theta,Counts,W).

py_sample_s_z_log_w(Theta,Disc,Counts,S,NS,Z,LogW) -->
	py_sample_bern_z(Disc,Counts,Z),
	py_sample_bern_s(Theta,Disc,Counts,S,NS),
	py_sample_beta_log_w(Theta,Counts,LogW).

py_sample_beta_w(_, [], 1) --> !.
py_sample_beta_w(Theta, Counts, W) -->
	{sumlist(Counts,N), Th1 is Theta+1, N1 is N-1},
	beta( Th1, N1, W).

py_sample_beta_log_w(_, [], 0) --> !.
py_sample_beta_log_w(Theta, Counts, LogW) -->
	{sumlist(Counts,N), Th1 is Theta+1, N1 is N-1},
	beta( Th1, N1, W), { LogW is log(W) }.

py_sample_bern_s(Theta,Disc,Counts,SumS,SumNS) -->
	(	{Counts=[_|Cm1], length(Cm1,Kminus1), numlist(1,Kminus1,KX)}
	->	{maplist(mul(Disc),KX,KDX)},
		sum_bernoulli(KDX, Theta, SumS),
		{SumNS is Kminus1 - SumS}
	;	{SumS=0,SumNS=0}
	).

py_sample_bern_z(Disc,Counts,Z) --> 
	{Disc1 is 1-Disc}, 
	seqmap( sample_bern_z(Disc1), Counts, ZX),
	{sumlist(ZX,Z)}.

sample_bern_z(Disc1,Count,SumZ) --> 
	{CountM2 is Count-2},
	(	{CountM2<0} -> {SumZ=0}
	;	{numlist(0,CountM2,I)}, 
		sum_bernoulli(I, Disc1, SumZ)
	).

sum_bernoulli(AX,B,T,S1,S2) :- sum_bernoulli(AX,B,0,T,S1,S2).
sum_bernoulli([],_,T,T,S,S) :- !.
sum_bernoulli([A|AX],B,T1,T3,S1,S3) :- 
	bernoulli(A,B,X,S1,S2), T2 is T1+X,
	sum_bernoulli(AX,B,T2,T3,S2,S3).

% Gamma distribution with rate parameter B.
:- procedure gamma(1,1).
gamma(A,B,X) --> gamma(A,U), {X is U/B}.

% Bernoulli with unnormalised weights for 0 and 1.
:- procedure bernoulli(1,1).
bernoulli(A,B,X) -->
	uniform01(U),
	({(A+B)*U<B} -> {X=1}; {X=0} ). 
*/

mul(X,Y,Z)   :- Z is X*Y.

