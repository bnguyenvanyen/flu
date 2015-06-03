(** Single strain deterministic simulation module *)

(** Module for parameter values *)
module type PARS =
  sig
    val size : float
    val r0 : float
    val e : float
    val etaN : float
    val g : float
    val nu : float
    val init_perturb : float
    val dilat_bound : float
    val contract_bound : float
  end;;

(** Functor initializing the module associated to parameter values *)
module Sys : functor (Pars : PARS) -> Dopri5.SYSTEM

(** Reasonable algorithm parameter values *)
module Default_Algp : Dopri5.ALGPARAMS


