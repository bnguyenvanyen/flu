(** Single strain stochastic simulation module *)

(** Integration uses a very simple Gillespie algorithm.
    We assume that events happen much faster than the infectivity beta varies *)

(** Module for parameter values *)
module type PARS =
  sig
    val n : float
    val r0 : float
    val e : float
    val etaN : float
    val g : float
    val nu : float
  end
(** Functor initializing the module associated to parameter values *)
module Sys : functor (Pars : PARS) -> 
                     (Gill.SYSTEM with type state = int * int * int 
                                   and type aux = float)

(** Reasonable algorithm parameter values *)
module Default_Algp : Gill.ALGPARAMS
