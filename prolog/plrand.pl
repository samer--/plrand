/*
 * Prolog part of random generator library
 * Samer Abdallah (2009)
*/
	  
:- module(plrand, [
	 	get_rnd_state/1	% -state
	,	set_rnd_state/1	% +state
	,	is_rnd_state/1		% +state
	,	init_rnd_state/1	% -state
	,	with_rnd_state/1	% +phrase(state)
	,	randomise/1			% -state
	,	init_jump/2       % +integer, -jump
	,	double_jump/2		% +jump, -jump
	,	spawn/3           % -state, +state, -state
	,	jump/3				% +jump, +state, -state

	,	rnd_state_term/2  
	,	run_rnd_state/3	

	,	stream_init/3     % +state, +integer, -stream
	,	stream_split/3    % +stream, -stream, -stream

	,	sample_Single_/1  % -float
	,	sample_Double_/1  % -float
	]).
/** <module> Skippable, splittable psuedorandom generator

   This module provides an interface to Pierre Ecuyer's pseudorandom generator
   library RngStream, which is a high-quality, long period pseudorandom generator
   that can do the neat trick of quickly skipping ahead an arbitrary number of samples
   in the stream without actually generating all the samples in between. This
   capability is used to implement a stream of random bits that can be split 
   recursively into 'pseudoindependent' substreams.

   ---+++ Types

      * state
         The state of a random generator, as a blob atom.
      * state_term
         The state represented as an ordinary Prolog term that can be stored as text.
      * jump
         A blob representing the operator to a skip ahead.
      * stream
         A spittable stream, containing a =|state|= and a =|jump|=.

*/

:-	use_foreign_library(foreign(plrand)).
:- meta_predicate with_rnd_state(:).
:- meta_predicate run_rnd_state(:,+,-).


%% get_rnd_state(-State:state) is det.
%
%  Unifies State with the current global RNG state.
%  @see set_rnd_state/1


%% set_rnd_state(+State:state) is det.
%
%  Sets the globab RNG state to State.
%  @see get_rnd_state/1

%% init_rnd_state(-State:state) is det.
%
%  Unifies State with the state that was set at load time.

%% randomise(-State:state) is det.
%
%  Returns a truly random state by reading from /dev/random.

%% is_rnd_state(+State:state) is semidet.
%
%  Succeeds if State is a BLOB atom representing a random generator state.


%% with_rnd_state(Cmd:phrase(state)) is nondet.
%
%  Runs DCG goal Cmd using the current global RNG state as the initial state
%  the global RNG state afterwards.
with_rnd_state(P) :- get_rnd_state(S1), call_dcg(P,S1,S2), set_rnd_state(S2).


%% run_rnd_state(Cmd:phrase(state), +State1:state, -State2:state) is nondet.
%% run_rnd_state(Cmd:phrase(state), +State1:state, -State2:state) is nondet.
%
%  Runs DCG phrase Cmd as with call_dcg/3, except that, if State1 is in the term representation 
%  (type =|state|=), it is converted to the blob representation (type =|state|=) and the final
%  state is converted back from term representation to get State2.
run_rnd_state(P,S1,S2) :- functor(S1,rs,7), !,
	term_to_rnd_state(S1,R1),
	call_dcg(P,R1,R2),
	rnd_state_to_term(R2,S2).

run_rnd_state(P,S1,S2) :- is_rnd_state(S1), !, call_dcg(P,S1,S2).


%% init_jump(+E:integer, -Jump:jump) is det.
%
% Unifies Jump with a BLOB atom representing an operator to jump forwards
% 2^E steps in the stream of random numbers. The generator has a period
% of about 2^191. Operators for jumping ahead by 2^76 and 2^127 are
% precomputed and so can be returned faster. The resulting BLOB can be
% used with jump/3 to advance any given random generator state.

%% double_jump( +Jump1:jump, -Jump2:jump) is det.
%
% Unifies Jump2 with an operator that jumps ahead twice as far Jump1,
% ie if Jump1 jumps by 2^E, then Jump2 jumps by 2^(E+1). Jump operators
% can be created by init_jump/2 and applied by jump/3.

%% jump( +Jump:jump, +State1:state, -State2:state) is det.
%
% Advances random generator state represented in State1 by the number
% of steps represented by Jump, unifying State2 with the result.
%
% @see double_jump/2
% @see init_jump/2

%% spawn( -New:state, +Orig:state, -Next:state) is det.
%
% Samples enough bits from Orig to specify a new generator state New,
% leaving the original generator in Next. If generator states represent
% streams of random numbers, then you can think of it as sampling a whole
% stream of values instead of just one value.
% Note: New is likely to point to a new point in the original stream
% far away from Orig and Next, simply because the period of the generator
% is so large (about 2^191) but there is no guarantee of this. Therefore, it's 
% possible (but unlikely) that New might produce a stream that overlaps
% significantly with samples drawn from Next. If you need to be sure
% that New is a long way from Orig and Next, then use jump/3 instead.
%
% @param New is the state of the newly created generator.
% @param Orig is the state of the source generator, the original stream.
% @param Next is the state of the source generator after extracting New.


%% stream_init(+State:state, +E:integer, -Stream:stream) is det.
%
% ==
% stream ---> stream(state, integer).
% ==
%
% [EXPERIMENTAL] Used to initialise a splittable stream. The idea is that
% a stream of random number can be recursively split into many substreams,
% each of which will not overlap unless very many samples are drawn.
%
% For example =| init_rnd_state(S), stream_init(S,76,Z) |=
% produces a stream such that the first split jumps 2^76 steps along.
% Thus, 2^76 samples can be drawn from either before overlap occurs.
% Since the period of the generator is more than 2^190, there can be up to
% 114 levels of recursive splitting, resulting in up to 2^114 substreams.
%
% @param State is the initial state of the generator.
% @param E is a nonnegative integer such that, when the stream is first split,
%        the new substream begins 2^E steps along the sequence. Subsequent splits
%        will jump twice as far each time.
% @param Stream is the new stream (a term)
stream_init(State,E,stream(State,Jump)) :- init_jump(E,Jump).

%% stream_split(+Stream1,-Stream2,-Stream3) is det.
%
%  [EXPERIMENTAL] Stream1 is split into independent streams Stream2 and Stream3.
%  Stream2 is actually the same as Stream1. Stream3 is a window on the same strem,
%  but starting a long way down the sequence.
stream_split(stream(S0,J0),stream(S0,J1),stream(S1,J1)) :- !,
	jump(J0,S0,S1), double_jump(J0,J1).


%% sample_Single_(-Float) is det.
%
% Samples a single precision floating point value in [0,1) using the 
% internal (global) random generator state. It consumes one 32 bit value
% from the generator. This is the fasted way to generate a random value 
% using this library.

%% sample_Double_(-Float) is det.
%
% Samples a double precision floating point value in [0,1) using the 
% internal (global) random generator state. It consumes two 32 bit values
% from the generator.


%% rnd_state_term( +State:state, -Term:state_term) is det.
%% rnd_state_term( -State:state, +Term:state_term) is det.
%
%  Convert between blob and term representation of random state.
rnd_state_term(RS,T) :-
	(	nonvar(RS)	-> rnd_state_to_term(RS,T)
	;	ground(T)	->	term_to_rnd_state(T,RS)
	).

