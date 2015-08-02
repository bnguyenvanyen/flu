(** Implementation of the Dormand-Prince integrator in Ocaml and BLAS, 
    with csv output. *)

(** ODE system description *)
module type SYSTEM =
  sig
    (** dimension of the system *)
    val n : int
    (** number of auxiliary variables *)
    val m : int
    (** vector field of the system. Updates z. *)
    val f : 
            ?z:Lacaml_float64.vec -> 
            float -> 
            Lacaml_float64.vec -> 
            Lacaml_float64.vec
    (** aux computes auxiliary statistics on the system. *)
    val aux : 
            ?z:Lacaml_float64.vec ->
            float ->
            Lacaml_float64.vec ->
            Lacaml_float64.vec  
    (** check_in_domain y returns true if y is in the domain, false otherwise *)
    val check_in_domain : 
            Lacaml_float64.vec ->
            bool
    (** shift_in_domain y updates and returns y.
        Called if y is not in the domain, it should shift it back in. *)
    val shift_in_domain : 
            ?z:Lacaml_float64.vec ->
            Lacaml_float64.vec ->
            Lacaml_float64.vec
    val csv_init :
            unit ->
            string list
  end

(** Parameter values for the algorithm *)
module type ALGPARAMS =
  sig
    (** initial step size *)
    val h0 : float
    (** maximum tolerated local error. Must be > 0. *)
    val delta : float
    (** maximum resolution in time. Must be > 0. *)
    val min_step : float
    (** minimum resolution in time. Must be > 0. *)
    val max_step : float
  end

(** Integrator of an ODE system *)
module type INTEGR =
  sig
   (*
    val update_k : 
          float -> 
          float -> 
          Lacaml_float64.mat -> 
          Lacaml_float64.vec -> 
          int -> 
          Lacaml_float64.mat
    val compute_k :
          float ->
          float ->
          Lacaml_float64.mat ->
          Lacaml_float64.vec ->
          Lacaml_float64.mat
    val runge_kutta :
          float ->
          Lacaml_float64.mat -> 
          Lacaml_float64.vec -> 
          Lacaml_float64.vec -> 
          Lacaml_float64.vec ->
          Lacaml_float64.vec
    val loop : 
          float -> 
          float -> 
          Lacaml_float64.mat ->
          Lacaml_float64.vec ->
          Lacaml_float64.vec ->
          Lacaml_float64.vec ->
          (float * float * Lacaml_float64.vec)
   *)
    (** simulate the process.

        simulate chan tf h0 x0 integrates from time 0 to time tf, 
        with initial stepsize of h0 and initial conditions x0 ;
        the resulting time series is output to the channel designated by chan and
        the final state is returned. *)
    val simulate : 
          out_channel ->
          float -> 
          Lacaml_float64.vec -> 
          (float * Lacaml_float64.vec)
  end

(** Functor building an integrator from an ODE system description and
    algorithmic parameters *)
module Integrator : functor (Sys : SYSTEM) -> 
                    functor (Algp : ALGPARAMS) -> INTEGR
