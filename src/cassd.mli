(** Single strain age-structured city-structured deterministic simulation module *)

(** Module for parameter values *)
module type PARS =
  sig
    (** Number of age classes *)
    val a : int
    (** Number of cities *)
    val c : int
    (** Number of hosts *)
    val size : float
    (** Basic reproductive ratio *)
    val r0 : float
    (** Strength of the seasonal forcing *)
    val e : float
    (** Immunity loss rate *)
    val g : float
    (** Recovery rate *)
    val nu : float
    (** Full cross-immunity loss rate *)
    val sensi_base_v : Lacaml_float64.vec
    (** Proportions of the age classes *)
    val age_prop_v : Lacaml_float64.vec
    (** Contact matrix *)
    val cont_base_m : Lacaml_float64.mat
    (** Transport rate between cities *)
    val eta_base_m : Lacaml_float64.mat
    (** Proportions of the cities *)
    val city_prop_v : Lacaml_float64.vec
    val init_perturb : float
    val dilat_bound : float
  end;;

(** Functor initializing the module associated to parameter values *)
module Sys : functor (Pars : PARS) -> Dopri5.SYSTEM

(** Reasonable algorithm parameter values *)
module Default_Algp : Dopri5.ALGPARAMS


