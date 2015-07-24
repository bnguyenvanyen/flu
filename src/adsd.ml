(* Simulate the double-strain age-structured ODE system using Dopri5 *)
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
    val sensi_v : Lacaml_float64.vec
    val prop_v : Lacaml_float64.vec
    val cont_m : Lacaml_float64.mat
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
            && q >= 0.
            && init_perturb > 0. 
            && dilat_bound > 0.);;

    assert (Mat.dim1 cont_m = Mat.dim2 cont_m
            && Vec.dim sensi_v = Vec.dim prop_v
            && Vec.dim prop_v = Mat.dim1 cont_m);;
    let m = Vec.dim prop_v
    let n = 12 * m * 2

    let eta1 = etaN1 *. size;;
    let eta2 = etaN2 *. size;;
    let bet0 = r0 *. nu

    let g12 = ~-. g1 -. g2;;
    let nnu = ~-. nu;;
    let ng1 = ~-. nu -. g1;;
    let ng2 = ~-. nu -. g2;;
    let nq = ~-. q;;
    let qg1 = ~-. q -. g1;;
    let qg2 = ~-. q -. g2;;
    let a_base = Mat.of_array (* to simplify computations *)
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
  
    let a = Array.make m a_base

    let j_base = Mat.of_array (* there's no copy function :( *)
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

    let j = Mat.make0 (n / 2) (n / 2);;
    (* the real "base j" is this on the diagonal of the age-clases *)
    for k = 1 to m do
      ignore (lacpy ~b:j ~br:(1 + (k - 1) * 12) ~bc:(1 + (k - 1) * 12) j_base)
    done

    let gi j k = (* get i *)
      let ind = 
        match j with
        | 10 -> 5
        | 20 -> 6
        | 12 -> 7
        | 21 -> 8
        | _ -> invalid_arg "expects 10, 20, 12 or 21"
      in 12 * (k - 1) + ind

    
    let gr j k = (* get r *)
      let ind = 
        match j with
        | 0 -> 1
        | 1 -> 2
        | 2 -> 3
        | 12 -> 4
        | _ -> invalid_arg "expects 0, 1, 2 or 12"
      in 12 * (k - 1) + ind

    let gq j k = (* get q *)
      let ind = 
        match j with
        | 0 -> 9
        | 1 -> 10
        | 2 -> 11
        | 12 -> 12
        | _ -> invalid_arg "expects 0, 1, 2, or 12"
      in 12 * (k - 1) + ind

    let eta1_v = Vec.make m eta1;;
    let eta2_v = Vec.make m eta2;;
    (* prop_v will now contain the number of hosts in each age class 
     * (and not the proportion) *)
    scal size prop_v 

    (* allocate stuff *)
    (* initial values (by age-class) *)
    let r0_v = Vec.make0 m
    let r1_v = Vec.make0 m
    let r2_v = Vec.make0 m
    let i1_v = Vec.make0 m
    let i2_v = Vec.make0 m  
    (* infectivities by susceptible *)
    let beta_i = Vec.make0 m
    let beta_s = Vec.make0 m
    (* end values *)
    let s_dot_v = Vec.make0 m
    let i_dot_v = Vec.make0 m
    let r_dot_v = Vec.make0 m

    let x = Vec.make0 (n / 2)
    let dx = Vec.make0 (n / 2) (* the dr, di, dq values *)

    (* same. FIXME can we drop some ? *)
    let beta_i1 = Vec.make0 m
    let beta_i2 = Vec.make0 m
    let beta_r0 = Vec.make0 m
    let beta_r1 = Vec.make0 m
    let beta_r2 = Vec.make0 m
    let tmp_m1 = Vec.make0 m
    let tmp_m2 = Vec.make0 m

    let f t y ~z =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in

      (* Age-class infections *)
      let tmp_m2 = copy ~y:tmp_m2 sensi_v in
      scal beta tmp_m2 ; 
      (* tmp_m2 now contains infectivities for each age class *)
      let tmp_m2 = Vec.div ~z:tmp_m2 tmp_m2 prop_v in
      (* tmp_m2 now contains the per infectious per susceptible infectivity *)

      (* For the first strain *)
      let i1_v = Vec.add ~n:m ~z:i1_v ~ofsx:5 ~incx:12 ~ofsy:7 ~incy:12 y y in
      (* i1_v contains the added infected individuals (by history) by age class *)
      let tmp_m1 = copy ~y:tmp_m1 eta1_v in (* fill tmp_m1 with eta1 *)
      (* compute the number of contacts : 
       * matrix of contacts * each number of infected *)
      let tmp_m1 = gemv ~y:tmp_m1 ~beta:1. ~m:m cont_m i1_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i1 = Vec.mul ~z:beta_i1 tmp_m1 tmp_m2 in

      (* For the second strain *)
      let i2_v = Vec.add ~n:m ~z:i1_v ~ofsx:6 ~incx:12 ~ofsy:8 ~incy:12 y y in
      (* i2_v contains the added infected individuals (by history) by age class *)
      let tmp_m1 = copy ~y:tmp_m1 eta2_v in (* fill tmp_m1 with eta2 *)
      (* compute the number of contacts : matrix of contacts * each number of infected *)
      let tmp_m1 = gemv ~y:tmp_m1 ~beta:1. ~m:m cont_m i2_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i2 = Vec.mul ~z:beta_i2 tmp_m1 tmp_m2 in

      (* these vectors we need for the jacobian *)
      let r0_v = copy ~n:m ~y:r0_v ~ofsx:1 ~incx:12 y in
      let r1_v = copy ~n:m ~y:r1_v ~ofsx:2 ~incx:12 y in
      let r2_v = copy ~n:m ~y:r2_v ~ofsx:3 ~incx:12 y in
      let beta_r0 = Vec.mul ~z:beta_r0 tmp_m2 r0_v in
      let beta_r1 = Vec.mul ~z:beta_r1 tmp_m2 r1_v in
      let beta_r2 = Vec.mul ~z:beta_r2 tmp_m2 r2_v in

      for k = 1 to m do
        (* the infection terms are not linear and vary in time *)
        a.(k - 1).{gr 0 1, gr 0 1} <- ~-. (beta_i1.{k}) -. beta_i2.{k} ;
        a.(k - 1).{gr 1 1, gr 1 1} <- ~-. (beta_i2.{k}) -. g1 ;
        a.(k - 1).{gr 2 1, gr 2 1} <- ~-. (beta_i1.{k}) -. g2 ;
        a.(k - 1).{gi 10 1, gr 0 1} <- beta_i1.{k} ;
        a.(k - 1).{gi 20 1, gr 0 1} <- beta_i2.{k} ;
        a.(k - 1).{gi 12 1, gr 2 1} <- beta_i1.{k} ;
        a.(k - 1).{gi 21 1, gr 1 1} <- beta_i2.{k} ;
        (* in the jacobian too *)
        (* but in the jacobian we keep all terms around ... *)
        (* r_k against itself *)
        j.{gr 0 k, gr 0 k} <- ~-. (beta_i1.{k}) -. beta_i2.{k} ;
        j.{gr 1 k, gr 1 k} <- ~-. (beta_i2.{k}) -. g2 ;
        j.{gr 2 k, gr 2 k} <- ~-. (beta_i1.{k}) -. g1 ;
        (* i_k against r_k *)
        j.{gi 10 k, gr 0 k} <- beta_i1.{k} ;
        j.{gi 20 k, gr 0 k} <- beta_i2.{k} ;
        j.{gi 12 k, gr 2 k} <- beta_i1.{k} ;
        j.{gi 21 k, gr 1 k} <- beta_i2.{k} ;

        for l = 1 to m do
          let beta_r0_kl = beta_r0.{k} *. cont_m.{k, l} (* I guess *)
          and beta_r1_kl = beta_r1.{k} *. cont_m.{k, l}
          and beta_r2_kl = beta_r2.{k} *. cont_m.{k, l}
          in
          (* r_k against i_l *)
          j.{gr 0 k, gi 10 l} <- ~-. beta_r0_kl ;
          j.{gr 0 k, gi 20 l} <- ~-. beta_r0_kl ;
          j.{gr 0 k, gi 12 l} <- ~-. beta_r0_kl ;
          j.{gr 0 k, gi 21 l} <- ~-. beta_r0_kl ;
          j.{gr 1 k, gi 20 l} <- ~-. beta_r1_kl ;
          j.{gr 1 k, gi 21 l} <- ~-. beta_r1_kl ;
          j.{gr 2 k, gi 10 l} <- ~-. beta_r2_kl ;
          j.{gr 2 k, gi 12 l} <- ~-. beta_r2_kl ;
          (* i_k against i_l *)
          j.{gi 12 k, gi 10 l} <- beta_r2_kl ;
          j.{gi 21 k, gi 20 l} <- beta_r1_kl ;
          if k = l then 
            (j.{gi 10 k, gi 10 l} <- beta_r0_kl -. nu ;
             j.{gi 10 k, gi 12 l} <- beta_r0_kl +. g2 ;
             j.{gi 20 k, gi 20 l} <- beta_r0_kl -. nu ;
             j.{gi 20 k, gi 21 l} <- beta_r0_kl +. g1 ;
             j.{gi 12 k, gi 12 l} <- beta_r2_kl -. nu -. g2 ;
             j.{gi 21 k, gi 21 l} <- beta_r1_kl -. nu -. g1)
          else
            (j.{gi 10 k, gi 10 l} <- beta_r0_kl ;
             j.{gi 10 k, gi 12 l} <- beta_r0_kl ;
             j.{gi 20 k, gi 20 l} <- beta_r0_kl ;
             j.{gi 20 k, gi 21 l} <- beta_r0_kl ;
             j.{gi 12 k, gi 12 l} <- beta_r2_kl ;
             j.{gi 21 k, gi 21 l} <- beta_r1_kl)
        done ;
      done ; 

      let x = copy ~y:x ~n:(n / 2) ~ofsx:1 y in
      for k = 1 to m do 
       ignore (gemv 
           ~n:12 
           ~alpha:1. 
           ~beta:0. 
           ~y:z
           ~ofsy:(1 + (k - 1) * 12)
           ~ofsx:(1 + (k - 1) * 12) 
           a.(k - 1) x)
      done ;
      let dx = copy ~y:dx ~n:(n / 2) ~ofsx:(n / 2 + 1) y in
      let z = gemv ~n:(n / 2) ~alpha:1. ~beta:0. ~y:z ~ofsy:(n / 2 + 1) j dx in
      (*
      print_string "start\n" ;
      Vec.iteri (fun i -> fun x -> if i > 36 then (print_float x ; print_string ", ")) z ;
      print_string "\nend\n" ;
      *)
      z

    let norm1_var y =
      let dx = copy ~y:dx ~n:(n/2) ~ofsx:(1 + n/2) y in
      amax dx

    let norm2_var y =
      let dx = copy ~y:dx ~n:(n/2) ~ofsx:(1 + n/2) y in
      sqrt (dot dx dx)
    
    let check_in_domain y =
      if Vec.min ~n:(n/2) y < 0. then false else 
      if norm2_var y > init_perturb *. dilat_bound then false else
      true

    let shift_in_domain y ~z =   
      for i = 1 to (n/2) do
        if y.{i} < 0. then
          let a = y.{i} /. (float_of_int (n/2 - 1)) in
          for j = 1 to (i - 1) do
            z.{j} <- y.{j} +. a
          done;
          for j = (i + 1) to (n/2) do
            z.{j} <- y.{j} +.a
          done;
          z.{i} <- 0.;
      done;  
      let nrm = norm2_var y in
      if nrm > init_perturb *. dilat_bound then
        for i = (1 + n/2) to n do
          z.{i} <- y.{i} *. init_perturb /. nrm
        done ;
      z

    let csv_init () =
      let rec f s n s_l =
        match n > 0 with
        | true -> 
            f s (n - 1) ((s ^ string_of_int n) :: s_l)
        | false -> 
            s_l
      in
      ["t" ; "h"] 
      @ List.concat
      (List.map (fun s -> f s m []) ["R0_" ; "R1_" ; "R2_" ; "R12_" ; 
                                     "I10_" ; "I20_" ; "I12_" ; "I21_" ; 
                                     "Q0_" ; "Q1_" ; "Q2_" ; "Q12_" ; 
                                     "dR0_" ; "dR1_" ; "dR2_" ; "dR12_" ; 
                                     "dI10_" ; "dI20_" ; "dI12_" ; "dI21_" ;
                                     "dQ0_" ; "dQ1_" ; "dQ2_" ; "dQ12_"])
  end;;

module Default_Algp =
  struct
    let h0 = 1. /. (24. *. 60.);;
    let delta = 0.1;; (* an error of 0.1 person seems ok *)
    let min_step = 1. /. (24. *. 3600.);; (* more than one step a second looks overkill *)
    let max_step = 1.;; (* we always want at least a resolution of a day *)
  end;;


