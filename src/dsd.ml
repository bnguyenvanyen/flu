open Lacaml.D;;

let pi = 4. *. atan 1.;;

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
    val init_perturb : float
    val dilat_bound : float
  end;;

module Sys (Pars : PARS) : Dopri5.SYSTEM =
  struct
    open Pars;;

    assert (size >= 0. 
            && r0 >= 0. 
            && e >= 0. 
            && etaN1 >= 0. 
            && etaN2 >= 0. 
            && g1 >= 0. 
            && g2 >= 0. 
            && nu >= 0.
            && init_perturb > 0. 
            && dilat_bound > 0.);;

    let eta1 = etaN1 *. size;;
    let eta2 = etaN2 *. size;;
    let bet0 = r0 *. nu;;

    let g12 = ~-. g1 -. g2;;
    let nnu = ~-. nu;;
    let ng1 = ~-. nu -. g1;;
    let ng2 = ~-. nu -. g2;;
    let nq = ~-. q;;
    let qg1 = ~-. q -. g1;;
    let qg2 = ~-. q -. g2;;
    let a = Mat.of_array (* to simplify computations *)
      [| 
        [| 0. ; g1 ; g2 ; 0. ; 0. ; 0. ; 0. ; 0. ; q  ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; g2 ; 0. ; 0. ; 0. ; 0. ; 0. ; q  ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; g1 ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; q  ; 0. |] ;
        [| 0. ; 0. ; 0. ; g12; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; q  |] ;
        [| 0. ; 0. ; 0. ; 0. ; nnu; 0. ; g2 ; 0. ; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; nnu; 0. ; g1 ; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; ng2; 0. ; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; ng1; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; nq ; g1 ; g2 ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; nu ; 0. ; 0. ; 0. ; 0. ; qg1; 0. ; g2 |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; nu ; 0. ; 0. ; 0. ; 0. ; qg2; g1 |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; nu ; nu ; 0. ; 0. ; 0. ; g12 -. q |]
      |];;
    let j = Mat.of_array (* there's no copy function :( *)
      [| 
        [| 0. ; g1 ; g2 ; 0. ; 0. ; 0. ; 0. ; 0. ; q  ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; g2 ; 0. ; 0. ; 0. ; 0. ; 0. ; q  ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; g1 ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; q  ; 0. |] ;
        [| 0. ; 0. ; 0. ; g12; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; q  |] ;
        [| 0. ; 0. ; 0. ; 0. ; nnu; 0. ; g2 ; 0. ; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; nnu; 0. ; g1 ; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; ng2; 0. ; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; ng1; 0. ; 0. ; 0. ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; nq ; g1 ; g2 ; 0. |] ;
        [| 0. ; 0. ; 0. ; 0. ; nu ; 0. ; 0. ; 0. ; 0. ; qg1; 0. ; g2 |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; nu ; 0. ; 0. ; 0. ; 0. ; qg2; g1 |] ;
        [| 0. ; 0. ; 0. ; 0. ; 0. ; 0. ; nu ; nu ; 0. ; 0. ; 0. ; g12 -. q |]
      |];;

    let x = Vec.make0 12;; (* the s, i, r values *)
    let dx = Vec.make0 12;; (* the ds, di, dr values *)
    let tmp1 = Vec.make0 12;;
    let tmp2 = Vec.make0 12;;

    let n = 24

    let gi j = (* get i *)
      match j with
      | 10 -> 5
      | 20 -> 6
      | 12 -> 7
      | 21 -> 8
      | _ -> invalid_arg "expects 10, 20, 12 or 21"
    
    let gr j = (* get r *)
      match j with
      | 0 -> 1
      | 1 -> 2
      | 2 -> 3
      | 12 -> 4
      | _ -> invalid_arg "expects 0, 1, 2 or 12"

    let gq j = (* get q *)
      match j with
      | 0 -> 9
      | 1 -> 10
      | 2 -> 11
      | 12 -> 12
      | _ -> invalid_arg "expects 0, 1, 2, or 12"


    let f t y ~z =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      let beta_i1 = beta *. (y.{gi 10} +. y.{gi 12} +. eta1) /. size in
      let beta_i2 = beta *. (y.{gi 20} +. y.{gi 21} +. eta2) /. size in
      let beta_r0 = beta *. y.{gr 0} /. size in
      let beta_r1 = beta *. y.{gr 1} /. size in
      let beta_r2 = beta *. y.{gr 2} /. size in
      a.{gr 0, gr 0} <- ~-. beta_i1 -. beta_i2 ;
      a.{gr 1, gr 1} <- ~-. beta_i2 -. g1 ;
      a.{gr 2, gr 2} <- ~-. beta_i1 -. g2 ;
      a.{gi 10, gr 0} <- beta_i1 ;
      a.{gi 20, gr 0} <- beta_i2 ;
      a.{gi 12, gr 2} <- beta_i1 ;
      a.{gi 21, gr 1} <- beta_i2 ;
      
      j.{gr 0, gr 0} <- ~-. beta_i1 -. beta_i2 ;
      j.{gr 0, gi 10} <- ~-. beta_r0 ;
      j.{gr 0, gi 20} <- ~-. beta_r0 ;
      j.{gr 0, gi 12} <- ~-. beta_r0 ;
      j.{gr 0, gi 21} <- ~-. beta_r0 ;
      j.{gr 1, gr 1} <- ~-. beta_i2 -. g2 ;
      j.{gr 1, gi 20} <- ~-. beta_r1 ;
      j.{gr 1, gi 21} <- ~-. beta_r1 ;
      j.{gr 2, gr 2} <- ~-. beta_i1 -. g1 ;
      j.{gr 2, gi 10} <- ~-. beta_r2 ;
      j.{gr 2, gi 12} <- ~-. beta_r2 ;
      
      j.{gi 10, gr 0} <- beta_i1 ;
      j.{gi 10, gi 10} <- beta_r0 -. nu ;
      j.{gi 10, gi 12} <- beta_r0 +. g2 ;
      j.{gi 20, gr 0} <- beta_i2 ;
      j.{gi 20, gi 20} <- beta_r0 -. nu ;
      j.{gi 20, gi 21} <- beta_r0 +. g1 ;
      j.{gi 12, gr 2} <- beta_i1 ;
      j.{gi 12, gi 10} <- beta_r2 ;
      j.{gi 12, gi 12} <- beta_r2 -. nu -. g2 ;
      j.{gi 21, gr 1} <- beta_i2 ;
      j.{gi 21, gi 20} <- beta_r1 ;
      j.{gi 21, gi 21} <- beta_r1 -. nu -. g1 ;
      
      let x = copy ~y:x ~n:12 ~ofsx:1 y in
      let dx = copy ~y:dx ~n:12 ~ofsx:13 y in
      let tmp1 = gemv ~alpha:1. ~beta:0. ~y:tmp1 a x in
      let tmp2 = gemv ~alpha:1. ~beta:0. ~y:tmp2 j dx in
      let z = copy ~y:z ~ofsy:1 tmp1 in
      let z = copy ~y:z ~ofsy:13 tmp2 in
      z

    let norm1_var y =
      let dx = copy ~y:dx ~n:12 ~ofsx:13 y in
      amax dx
    
    let norm2_var y =
      let dx = copy ~y:dx ~n:(n/2) ~ofsx:(1 + n/2) y in
      sqrt (dot dx dx)

    let check_in_domain y =
      if Vec.min ~n:12 y < 0. then false else 
      if norm2_var y > init_perturb *. dilat_bound then false else
      true

    let shift_in_domain y ~z =   
      for i = 1 to 12 do
        if y.{i} < 0. then 
          let a = y.{i} /. 11. in
          for j = 1 to (i - 1) do
            z.{j} <- y.{j} +. a
          done;
          for j = (i + 1) to 12 do
            z.{j} <- y.{j} +. a
          done;
          z.{i} <- 0.
      done; 
      let nrm = norm2_var y in
        if nrm > init_perturb *. dilat_bound then
          for i = (1 + n/2) to n do
            z.{i} <- y.{i} *. init_perturb /. nrm
          done ;
      z

    let csv_init () =
      ["t" ; "h" ; 
       "R0" ; "R1" ; "R2" ; "R12" ;
       "I10" ; "I20" ; "I12" ; "I21" ;
       "Q0" ; "Q1" ; "Q2" ; "Q12" ;
       "dR0" ; "dR1" ; "dR2" ; "dR12" ;
       "dI10" ; "dI20" ; "dI12" ; "dI21" ;
       "dQ0" ; "dQ1" ; "dQ2" ; "dQ12"]
  end


module Default_Algp =
  struct
    let h0 = 1. /. (24. *. 60.)
    let delta = 0.1
    let min_step = 1. /. (24. *. 3600.)
    let max_step = 1.
  end;;

