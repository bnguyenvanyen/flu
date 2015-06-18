(* Stochastic single strain simulation *)
(* Hypothesis : events are fast compared to the forcing *)

(* internal type for simulations *)

let pi = 4. *. atan 1.;;

module type PARS =
  sig
    val n : float
    val r0 : float
    val e : float
    val etaN : float
    val g : float
    val nu : float
  end

module Sys (Pars : PARS) : (Gill.SYSTEM with type state = int * int * int) =
  struct
    open Pars
    type state = (int * int * int)
    let eta = etaN *. n;;
    let bet0 = r0 *. nu;;
    
    let infectivity t =
      bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.));;
    
    let immu_loss_rate t st =
      let _, _, r = st in
      let fr = float_of_int r in
      g *. fr;;
    
    let infection_rate t st =
      let s, i, _ = st in
      let beta = infectivity t in
      let fs = float_of_int s in
      let fi = float_of_int i in
      beta *. fs /. n *. (fi +. eta);;
       
    let recovery_rate t st =
      let _, i, _ = st in
      let fi = float_of_int i in
      nu *. fi;;
    
    let immu_loss_modif st =
      let s, i, r = st in
      (s + 1, i, r - 1);;
    
    let infection_modif st =
      let s, i, r = st in
      (s - 1, i + 1, r);;
    
    let recovery_modif st =
      let s, i, r = st in
      (s, i - 1, r + 1);;

    let min_step = 1.
    let fl = [immu_loss_rate ; infection_rate ; recovery_rate]
    let ml = [immu_loss_modif ; infection_modif ; recovery_modif]
    
    let csv_init () =
      ["t" ; "S" ; "I" ; "R"]
    
    let csv_line t st =
      let s, i, r = st in
      [string_of_float t ; 
       string_of_int s ; string_of_int i ; string_of_int r]
  end;;

module Default_Algp =
  struct
    let step_size = 0.25
  end
