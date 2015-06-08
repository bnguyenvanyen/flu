(* Simulate the single-strain ODE system using Dopri5 *)
open Lacaml.D;;

let pi = 4. *. atan 1.;;

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

module Sys (Pars : PARS) : Dopri5.SYSTEM =
  struct
    open Pars
    let eta = etaN *. size;;
    let bet0 = r0 *. nu;;
    let a = Mat.make0 3 3;; (* to simplify computations *)
    a.{1,3} <- g ;;
    a.{2,2} <- ~-. nu ;;
    a.{3,2} <- nu ;;
    a.{3,3} <- ~-. g ;;
    let j = Mat.make0 3 3;;(* jacobian *)
    j.{1,3} <- g ;;
    j.{3,2} <- nu ;;
    j.{3,3} <- ~-. g ;;

    let x = Vec.make0 3;; (* the s, i, r values *)
    let dx = Vec.make0 3;; (* the ds, di, dr values *)
    let tmp1 = Vec.make0 3;;
    let tmp2 = Vec.make0 3;;

    let n = 6

    let f t y ~z =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      let beta_i = beta *. (y.{2} +. eta) /. size in
      let beta_s = beta *. y.{1} /. size in
      a.{1,1} <- ~-. beta_i ;
      a.{2,1} <- beta_i ;

      j.{1,1} <- ~-. beta_i ;
      j.{1,2} <- ~-. beta_s ;
      j.{2,1} <- beta_i ;
      j.{2,2} <- beta_s -. nu ;

      let x = copy ~y:x ~n:3 ~ofsx:1 y in
      let dx = copy ~y:dx ~n:3 ~ofsx:4 y in
      let tmp1 = gemv ~alpha:1. ~beta:0. ~y:tmp1 a x in
      let tmp2 = gemv ~alpha:1. ~beta:0. ~y:tmp2 j dx in
      let z = copy ~y:z ~ofsy:1 tmp1 in
      let z = copy ~y:z ~ofsy:4 tmp2 in
      z

    let norm1_var y =
      let dx = copy ~y:dx ~n:3 ~ofsx:4 y in
      amax dx

    let norm2_var y =
      let dx = copy ~y:dx ~n:3 ~ofsx:4 y in
      sqrt (dot dx dx)
    
    (* FIXME breaks if two values are < 0 at once ... *)
    (* FIXME also, useless manipulations of y ... *)
    let check_in_domain y =
      if Vec.min y < 0. then false else 
      if norm2_var y > init_perturb *. dilat_bound then false else
      true

    let shift_in_domain y ~z =   
      let s, i, r =
        if y.{1} < 0. then 
          (0., y.{2} +. y.{1} /. 2., y.{3} +. y.{1} /. 2.)
        else 
        if y.{2} < 0. then 
          (y.{1} +. y.{2} /. 2., 0., y.{3} +. y.{2} /. 2.)
        else 
        if y.{3} < 0. then 
          (y.{1} +. y.{3} /. 2., y.{2} +. y.{3} /. 2., 0.)
        else (y.{1}, y.{2}, y.{3})
      in let ampl =
        let n = norm2_var y in
        if n > init_perturb *. dilat_bound then
          n /. init_perturb
        else if n < init_perturb /. contract_bound then
          n /. init_perturb
        else 
          1.
      in 
      z.{1} <- s ;
      z.{2} <- i ;
      z.{3} <- r ;
      z.{4} <- y.{4} /. ampl ;
      z.{5} <- y.{5} /. ampl ;
      z.{6} <- y.{6} /. ampl ;
      z

    let csv_init () =
      ["t" ; "h" ; "S" ; "I" ; "R"; "dS" ; "dI" ; "dR"]
  end;;

module Default_Algp =
  struct
    let h0 = 1. /. (24. *. 60.);;
    let delta = 0.1;; (* an error of 0.1 person seems ok *)
    let min_step = 1. /. (24. *. 3600.);; (* more than one step a second looks overkill *)
    let max_step = 1.;; (* we always want at least a resolution of a day *)
  end;;


