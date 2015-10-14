(* Simulate the single-strain ODE system using Dopri5 *)
open Lacaml.D;;

let pi = 4. *. atan 1.;;

module type PARS =
  sig
    val size : float
    val r0 : float
    val e : float
    val b : float
    val etaN : float
    val phi : float
    val g : float
    val nu : float
    val nr : int
    val init_perturb : float
    val dilat_bound : float
  end;;

module Sys (Pars : PARS) : Dopri5.SYSTEM =
  struct
    open Pars
    let et0 = etaN *. size
    let bet0 = r0 *. nu
    let gnr = g *. float_of_int nr
    let n = 2 * (2 + nr)
    let m = 1
    let a = Mat.make0 (n / 2) (n / 2);; (* to simplify computations *)
    a.{1, 2 + nr} <- gnr ;
    a.{2, 2} <- ~-. nu ;
    a.{3, 2} <- nu ;;
    for i = 3 to 2 + nr do
      a.{i, i} <- ~-. gnr
    done ;;
    for i = 4 to 2 + nr do
      a.{i, i - 1} <- gnr
    done 
    let j = lacpy a;;(* jacobian *)

    let x = Vec.make0 (n / 2) (* the s, i, r values *)
    let dx = Vec.make0 (n / 2) (* the ds, di, dr values *)
    let tmp1 = Vec.make0 (n / 2)
    let tmp2 = Vec.make0 (n / 2);;

    let f ?(z=Vec.make0 n) t y =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      let eta = et0 *. (1. +. b *. cos (2. *. pi *. (t +. phi ) /. 365.)) in
      let beta_i = beta *. (y.{2} +. eta) /. size in
      let beta_s = beta *. y.{1} /. size in
      a.{1,1} <- ~-. beta_i ;
      a.{2,1} <- beta_i ;

      j.{1,1} <- ~-. beta_i ;
      j.{1,2} <- ~-. beta_s ;
      j.{2,1} <- beta_i ;
      j.{2,2} <- beta_s -. nu ;

      let x = copy ~y:x ~n:(n / 2) ~ofsx:1 y in
      let dx = copy ~y:dx ~n:(n / 2) ~ofsx:(n / 2 + 1) y in
      let tmp1 = gemv ~alpha:1. ~beta:0. ~y:tmp1 a x in
      let tmp2 = gemv ~alpha:1. ~beta:0. ~y:tmp2 j dx in
      let z = copy ~y:z ~ofsy:1 tmp1 in
      let z = copy ~y:z ~ofsy:(n / 2 + 1) tmp2 in
      z

    let aux ?(z=Vec.make0 m) t y =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      let eta = et0 *. (1. +. b *. cos (2. *. pi *. (t +. phi ) /. 365.)) in
      let infct = beta *. y.{1} *. (y.{2} +. eta) /. size in
      let inc = infct *. 7. *. 100000. /. size in
      (* Weekly incidence for 100 000 *)
      z.{1} <- inc ;
      z

    let norm1_var y =
      let dx = copy ~y:dx ~n:(n / 2) ~ofsx:(n / 2 + 1) y in
      amax dx

    let norm2_var y =
      let dx = copy ~y:dx ~n:(n / 2) ~ofsx:(n / 2 + 1) y in
      sqrt (dot dx dx)
    
    (* FIXME breaks if two values are < 0 at once ... *)
    (* FIXME also, useless manipulations of y ... *)
    let check_in_domain y =
      if Vec.min ~n:(n / 2) y < 0. then false else 
      if norm2_var y > init_perturb *. dilat_bound then false else
      true

    let shift_in_domain ?(z=Vec.make0 n) y =   
      for i = 1 to (n / 2) do
        if y.{i} < 0. then
          let a = y.{i} /. float_of_int (n / 2 - 1) in
          for j = 1 to (i - 1) do
            z.{j} <- y.{j} +. a
          done;
          for j = (i + 1) to (n / 2) do
            z.{j} <- y.{j} +.a
          done;
          z.{i} <- 0.;
      done;  
      let ampl =
        let nrm = norm2_var y in
        if nrm > init_perturb *. dilat_bound then
          nrm /. init_perturb
        else 
          1.
      in 
      for i = (n / 2 + 1) to n do
        z.{i} <- y.{i} /. ampl
      done;
      z

    let csv_init () =
      let rec f s n s_l =
        match n > 0 with
        | true -> 
            f s (n - 1) ((s ^ string_of_int n) :: s_l)
        | false -> 
            s_l
      in
        ["t" ; "h" ; "inc" ; "S" ; "I"] 
      @ (f "R_" nr [])
      @ ["dS" ; "dI"]
      @ (f "dR_" nr [])
  end;;

module Default_Algp =
  struct
    let h0 = 1. /. (24. *. 60.);;
    let delta = 0.1;; (* an error of 0.1 person seems ok *)
    let min_step = 1. /. (24. *. 3600.);; (* more than one step a second looks overkill *)
    let max_step = 1.;; (* we always want at least a resolution of a day *)
  end;;


