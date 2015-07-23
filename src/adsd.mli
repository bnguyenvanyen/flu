(** Double strain age-structured deterministic simulation module *)

(** Module for parameter values *)
module type PARS =
  sig
    val size : float
    val r0 : float
    val e : float
    val etaN1 : float
    val etaN2 : float
    val g1 : float
    val g2 : float
    val nu : float
    val q : float
    val sensi_v : Lacaml_float64.vec
    val prop_v : Lacaml_float64.vec
    val cont_m : Lacaml_float64.mat
    val init_perturb : float
    val dilat_bound : float
  end;;

(** Functor initializing the module associated to parameter values *)
module Sys : functor (Pars : PARS) -> Dopri5.SYSTEM

(** Reasonable algorithm parameter values *)
module Default_Algp : Dopri5.ALGPARAMS


