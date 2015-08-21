(** Single strain deterministic simulation module *)

(** Module for parameter values *)
module type PARS =
  sig
    (** Number of hosts *)
    val size : float
    (** Basic reproductive ratio *)
    val r0 : float
    (** Strength of the seasonal forcing *)
    val e : float
    (** Strength of the seasonal immigration forcing *)
    val b : float
    (** Base immigration rate (per host) *)
    val etaN : float
    (** Phase shift between the two forcings *)
    val phi : float
    (** Immunity loss rate *)
    val g : float
    (** Recovery rate *)
    val nu : float
    val init_perturb : float
    val dilat_bound : float
  end;;

(** Functor initializing the module associated to parameter values *)
module Sys : functor (Pars : PARS) -> Dopri5.SYSTEM

(** Reasonable algorithm parameter values *)
module Default_Algp : Dopri5.ALGPARAMS


