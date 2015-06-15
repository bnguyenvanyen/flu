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
    open Pars
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

    let f t y ~z =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      let beta_i1 = beta *. (y.{5} +. y.{7} +. eta1) /. size in
      let beta_i2 = beta *. (y.{6} +. y.{8} +. eta2) /. size in
      let beta_r0 = beta *. y.{1} /. size in
      let beta_r1 = beta *. y.{2} /. size in
      let beta_r2 = beta *. y.{3} /. size in
      a.{1,1} <- ~-. beta_i1 -. beta_i2 ;
      a.{2,2} <- ~-. beta_i2 ;
      a.{3,3} <- ~-. beta_i1 ;
      a.{5,1} <- beta_i1 ;
      a.{6,1} <- beta_i2 ;
      a.{7,3} <- beta_i1 ;
      a.{8,2} <- beta_i2 ;
      
      j.{1,1} <- ~-. beta_i1 -. beta_i2 ;
      j.{1,5} <- ~-. beta_r0 ;
      j.{1,6} <- ~-. beta_r0 ;
      j.{1,7} <- ~-. beta_r0 ;
      j.{1,8} <- ~-. beta_r0 ;
      j.{2,2} <- ~-. beta_i2 -. g2 ;
      j.{2,6} <- ~-. beta_r1 ;
      j.{2,8} <- ~-. beta_r1 ;
      j.{3,3} <- ~-. beta_i1 -. g1 ;
      j.{3,5} <- ~-. beta_r2 ;
      j.{3,7} <- ~-. beta_r2 ;
      
      j.{5,1} <- beta_i1 ;
      j.{5,5} <- beta_r0 -. nu ;
      j.{5,7} <- beta_r0 +. g2 ;
      j.{6,1} <- beta_i2 ;
      j.{6,6} <- beta_r0 -. nu ;
      j.{6,8} <- beta_r0 +. g1 ;
      j.{7,3} <- beta_i1 ;
      j.{7,5} <- beta_r2 ;
      j.{7,7} <- beta_r2 -. nu -. g2 ;
      j.{8,2} <- beta_i2 ;
      j.{8,6} <- beta_r1 ;
      j.{8,8} <- beta_r1 -. nu -. g1 ;
      
      let x = copy ~y:x ~n:12 ~ofsx:1 y in
      let dx = copy ~y:dx ~n:12 ~ofsx:13 y in
      let tmp1 = gemv ~alpha:1. ~beta:0. ~y:tmp1 a x in
      let tmp2 = gemv ~alpha:1. ~beta:0. ~y:tmp2 j dx in
      let z = copy ~y:z ~ofsy:1 tmp1 in
      let z = copy ~y:z ~ofsy:13 tmp2 in
      z

    let norm_var y =
      let dx = copy ~y:dx ~n:12 ~ofsx:13 y in
      amax dx
    
    (* FIXME breaks if two values are < 0 at once ... *)
    (* FIXME also, useless manipulations of y ... *)
    let check_in_domain y =
      if Vec.min ~n:12 y < 0. then false else 
      if norm_var y > init_perturb *. dilat_bound then false else
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
      let ampl =
        let n = norm_var y in
        if n > init_perturb *. dilat_bound then
          n /. init_perturb
        else 
          1.
      in 
      for i = 13 to 24 do
        z.{i} <- y.{i} /. ampl
      done;
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

