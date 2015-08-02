(** A simple implementation of the Gillespie algorithm, with csv output. *)

(** The time variable is represented by floats *)
type time = float

(** Specification of continuous time markov processes *)
module type SYSTEM =
  sig
    (** type of the values taken by the process *)
    type state
    (** type of the auxiliary statistics on the process *)
    type aux
    (** header for the csv output (column names) *)
    val csv_init : unit -> string list list
    (** representation for the state of the system to be output in the csv.
      * Should be compatible with the output of csv_init *)
    val csv_line : time -> aux -> state -> string list
    (** list of rate functions of the process *)
    (** computes the auxiliary statistics *)
    val aux_fun : time -> state -> aux
    val fl : (time -> state -> float) list
    (** list of modification functions of the process.
        Should be compatible with fl *)
    val ml : (state -> state) list
  end

(** Parameters controlling the algorithm *)
module type ALGPARAMS =
  sig
    (** maximum time resolution wanted *)
    val min_step : float
  end;;

(** Integration of continuous time markov processes *)
module type INTEGR =
  sig
    (** type of the values taken by the process *)
    type state
    (** type of the auxiliary statistics on the process *)
    type aux
    (** type of the process *)
    type process = {past : (time * aux * state) list ;
                    present : time * aux * state ; 
                    future : (time * aux * state) list}
    (** go back one step in (known) time. Fails if at the start of time. *)
    val backward : process -> process
    (** go forward one step in (known) time. Fails if at the end of time. *)
    val forward : process -> process
    (** csv representation of the process. Built from csv_init and csv_line *)
    val csv_of_proc : process -> string list list
    (** simulate the process.

        simulate fname tf x0 integrates the system from time 0 and state x0 to
        time tf, stores the result in the file designated by fname, 
        and returns the results. *)
    val simulate :
          out_channel ->
          time ->
          state ->
          time * aux * state
  end

(** functor to create integrators from system descriptions *)
module Integrator : functor (Sys : SYSTEM) -> 
                    functor (Algp : ALGPARAMS) -> 
                            INTEGR with type state = Sys.state
                                    and type aux = Sys.aux
