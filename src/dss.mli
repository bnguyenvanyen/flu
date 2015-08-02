(** Double-strain stochastic simulation module *)

(** Integration uses a very simple Gillespie algorithm.
    We assume that events happen much faster than the infectivity beta varies *)

(** Module for parameter values *)
module type PARS =
  sig
    val n : float
    val r0 : float
    val e : float
    val etaN1 : float
    val etaN2 : float
    val g1 : float
    val g2 : float
    val nu : float
    val q : float
  end
(** Functor initializing the module associated to parameter values *)
module Sys : functor (Pars : PARS) -> 
                     (Gill.SYSTEM with type state = int array
                                   and type aux = float * float)

(** Reasonable algorithm parameter values *)
module Default_Algp : Gill.ALGPARAMS
