(* Simulate the double-strain age-structured city-structured ODE system using Dopri5 *)
(* expected structure of y : [
 *                            [
 *                             [
 *                              [strains] 
 *                             age] 
 *                            city] 
 *                           sys / var ]
 *)

open Lacaml.D;;

let pi = 4. *. atan 1.;;

module type PARS =
  sig
    val a : int
    val c : int
    val size : float
    val r0 : float
    val e : float
    val g1 : float
    val g2 : float
    val nu : float
    val q : float
    val sensi_base_v : Lacaml_float64.vec
    val age_prop_v : Lacaml_float64.vec
    val cont_base_m : Lacaml_float64.mat
    val eta_base_m : Lacaml_float64.mat
    val city_prop_v : Lacaml_float64.vec
    val init_perturb : float
    val dilat_bound : float
  end;;

module Sys (Pars : PARS) : Dopri5.SYSTEM =
  struct
    open Pars;;
    
    (* check sign of parameters *)
    assert (a > 0
            && c > 0
            && size >= 0. 
            && r0 >= 0. 
            && e >= 0. 
            && g1 >= 0. 
            && g2 >= 0. 
            && nu >= 0.
            && q >= 0.
            && init_perturb > 0. 
            && dilat_bound > 0.);;

    (* check dimensions for age structure *)
    assert (a = Vec.dim age_prop_v
            && a = Mat.dim1 cont_base_m 
            && a = Mat.dim2 cont_base_m
            && a = Vec.dim sensi_base_v);;

    (* check the dimensions for city structure *)
    assert (c = Mat.dim1 eta_base_m
            && c = Mat.dim2 eta_base_m
            && c = Vec.dim city_prop_v);;

    (* there shouldn't be any flow from a city to itself *)
    assert (asum (Mat.copy_diag eta_base_m) = 0.)

    let n = 12 * a * c * 2
    let m = 2 (* for now ... *)

    let bet0 = r0 *. nu

    let g12 = ~-. g1 -. g2;;
    let nnu = ~-. nu;;
    let ng1 = ~-. nu -. g1;;
    let ng2 = ~-. nu -. g2;;
    let nq = ~-. q;;
    let qg1 = ~-. q -. g1;;
    let qg2 = ~-. q -. g2;;
    let a_base_m = Mat.of_array (* to simplify computations *)
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
  
    let lin_f_m = Array.make (a * c) a_base_m

    let jac_base_m = Mat.of_array (* there's no copy function :( *)
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

    let jac_m = Mat.make0 (n / 2) (n / 2);;
    (* the real "base j" is this on the diagonal of the age-clases *)
    for i = 1 to c do
      for k = 1 to a do
        ignore (lacpy 
                  ~b:jac_m
                  ~br:(1 + (i - 1) * 12 * a + (k - 1) * 12) 
                  ~bc:(1 + (i - 1) * 12 * a + (k - 1) * 12) 
                  jac_base_m
               )
      done
    done

    let sensi_v = Vec.make0 (a * c);;
    for i = 1 to c do
      for k = 1 to a do
        sensi_v.{(i - 1) * a + k} <- sensi_base_v.{k}
      done
    done

    let cont_m = Mat.make0 (a * c) (a * c);;
    for i = 1 to c do
      for k = 1 to a do
        for l = 1 to a do
          cont_m.{(i - 1) * a + k, (i - 1) * a + l} <- cont_base_m.{k, l}
        done
      done
    done

    let etlin_f_m = Mat.make0 (a * c) (a * c);;
    for i = 1 to c do
      for j = 1 to c do
        for k = 1 to a do
          for l = 1 to a do
            (* we suppose that airflow is independent of age-class (wrong) *)
            etlin_f_m.{(i - 1) * a + k, (j - 1) * a + l} 
              <- eta_base_m.{i, j} *. cont_base_m.{k, l}
          done
        done
      done
    done

    let prop_v = Vec.make0 (a * c);;
    for i = 1 to c do
      for k = 1 to a do
        prop_v.{(i - 1) * a + k} <- city_prop_v.{i} *. age_prop_v.{k}
      done
    done

    let gi j k i = (* get i *)
      let ind = 
        match j with
        | 10 -> 5
        | 20 -> 6
        | 12 -> 7
        | 21 -> 8
        | _ -> invalid_arg "expects 10, 20, 12 or 21"
      in 12 * a * (i - 1) + 12 * (k - 1) + ind

    let gr j k i = (* get r *)
      let ind = 
        match j with
        | 0 -> 1
        | 1 -> 2
        | 2 -> 3
        | 12 -> 4
        | _ -> invalid_arg "expects 0, 1, 2 or 12"
      in 12 * a * (i - 1) + 12 * (k - 1) + ind

    let gq j k i = (* get q *)
      let ind = 
        match j with
        | 0 -> 9
        | 1 -> 10
        | 2 -> 11
        | 12 -> 12
        | _ -> invalid_arg "expects 0, 1, 2, or 12"
      in 12 * a * (i - 1) + 12 * (k - 1) + ind;;

    (* prop_v will now contain the number of hosts in each age class 
     * (and not the proportion) *)
    scal size prop_v

    (* allocate stuff *)
    (* initial values (by age-class) *)
    let r0_v = Vec.make0 (a * c)
    let r1_v = Vec.make0 (a * c)
    let r2_v = Vec.make0 (a * c)
    let i1_v = Vec.make0 (a * c)
    let i2_v = Vec.make0 (a * c) 
    (* infectivities by susceptible *)
    let beta_i = Vec.make0 (a * c)
    let beta_s = Vec.make0 (a * c)
    (* end values *)
    let s_dot_v = Vec.make0 (a * c)
    let i_dot_v = Vec.make0 (a * c)
    let r_dot_v = Vec.make0 (a * c)

    let x = Vec.make0 (n / 2)
    let dx = Vec.make0 (n / 2) (* the dr, di, dq values *)

    (* same. FIXME can we drop some ? *)
    let eta1_v = Vec.make0 (a * c)
    let eta2_v = Vec.make0 (a * c)
    let beta_i1_v = Vec.make0 (a * c)
    let beta_i2_v = Vec.make0 (a * c)
    let beta_r0_v = Vec.make0 (a * c)
    let beta_r1_v = Vec.make0 (a * c)
    let beta_r2_v = Vec.make0 (a * c)
    let tmp_mc1 = Vec.make0 (a * c)
    let tmp_mc2 = Vec.make0 (a * c)

    let f ?(z=Vec.make0 n) t y =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in

      (* Age-class infections *)
      let tmp_mc2 = copy ~y:tmp_mc2 sensi_v in
      scal beta tmp_mc2 ; 
      (* tmp_mc2 now contains infectivities for each age class *)
      let tmp_mc2 = Vec.div ~z:tmp_mc2 tmp_mc2 prop_v in
      (* tmp_mc2 now contains the per infectious per susceptible infectivity *)

      (* For the first strain *)
      let i1_v = Vec.add ~n:(a * c) ~z:i1_v ~ofsx:5 ~incx:12 ~ofsy:7 ~incy:12 y y in
      let eta1_v = gemv ~y:eta1_v ~beta:0. ~m:(a * c) etlin_f_m i1_v in
      (* i1_v contains the added infected individuals (by history) by age class *)
      let tmp_mc1 = copy ~y:tmp_mc1 eta1_v in (* fill tmp_mc1 with eta1 *)
      (* compute the number of contacts : 
       * matrix of contacts * each number of infected *)
      let tmp_mc1 = gemv ~y:tmp_mc1 ~beta:1. ~m:a cont_m i1_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i1_v = Vec.mul ~z:beta_i1_v tmp_mc1 tmp_mc2 in

      (* For the second strain *)
      let i2_v = Vec.add ~n:(a * c) ~z:i2_v ~ofsx:6 ~incx:12 ~ofsy:8 ~incy:12 y y in
      (* i2_v contains the added infected individuals (by history) by age class *)
      let eta2_v = gemv ~y:eta2_v ~beta:0. ~m:(a * c) etlin_f_m i2_v in
      let tmp_mc1 = copy ~y:tmp_mc1 eta2_v in (* fill tmp_mc1 with eta2 *)
      (* compute the number of contacts : matrix of contacts * each number of infected *)
      let tmp_mc1 = gemv ~y:tmp_mc1 ~beta:1. ~m:a cont_m i2_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i2_v = Vec.mul ~z:beta_i2_v tmp_mc1 tmp_mc2 in

      (* these vectors we need for the jacobian *)
      let r0_v = copy ~n:(a * c) ~y:r0_v ~ofsx:1 ~incx:12 y in
      let r1_v = copy ~n:(a * c) ~y:r1_v ~ofsx:2 ~incx:12 y in
      let r2_v = copy ~n:(a * c) ~y:r2_v ~ofsx:3 ~incx:12 y in
      let beta_r0_v = Vec.mul ~z:beta_r0_v tmp_mc2 r0_v in
      let beta_r1_v = Vec.mul ~z:beta_r1_v tmp_mc2 r1_v in
      let beta_r2_v = Vec.mul ~z:beta_r2_v tmp_mc2 r2_v in

      for i = 1 to c do
        for k = 1 to a do
          let ind = (i - 1) * a + k in
          let beta_i1 = beta_i1_v.{ind} in
          let beta_i2 = beta_i2_v.{ind} in
          (* the infection terms are not linear and vary in time *)
          lin_f_m.(ind - 1).{gr 0 1 1, gr 0 1 1} <- ~-. beta_i1 -. beta_i2 ;
          lin_f_m.(ind - 1).{gr 1 1 1, gr 1 1 1} <- ~-. beta_i2 -. g1 ;
          lin_f_m.(ind - 1).{gr 2 1 1, gr 2 1 1} <- ~-. beta_i1 -. g2 ;
          lin_f_m.(ind - 1).{gi 10 1 1, gr 0 1 1} <- beta_i1 ;
          lin_f_m.(ind - 1).{gi 20 1 1, gr 0 1 1} <- beta_i2 ;
          lin_f_m.(ind - 1).{gi 12 1 1, gr 2 1 1} <- beta_i1 ;
          lin_f_m.(ind - 1).{gi 21 1 1, gr 1 1 1} <- beta_i2 ;
          (* in the jacobian too *)
          (* but in the jacobian we keep all terms around ... *)
          (* r_k against itself *)
          jac_m.{gr 0 k i, gr 0 k i} <- ~-. beta_i1 -. beta_i2 ;
          jac_m.{gr 1 k i, gr 1 k i} <- ~-. beta_i2 -. g2 ;
          jac_m.{gr 2 k i, gr 2 k i} <- ~-. beta_i1 -. g1 ;
          (* i_k against r_k *)
          jac_m.{gi 10 k i, gr 0 k i} <- beta_i1 ;
          jac_m.{gi 20 k i, gr 0 k i} <- beta_i2 ;
          jac_m.{gi 12 k i, gr 2 k i} <- beta_i1 ;
          jac_m.{gi 21 k i, gr 1 k i} <- beta_i2 ;
          for j = 1 to c do
            for l = 1 to a do
              (* I guess *)
              let beta_r0 =
                if i = j then beta_r0_v.{k} *. cont_m.{k, l}
                else beta_r0_v.{k} *. etlin_f_m.{i, j}
              and beta_r1 = 
                if i = j then beta_r1_v.{k} *. cont_m.{k, l}
                else beta_r1_v.{k} *. etlin_f_m.{i, j}
              and beta_r2 = 
                if i = j then beta_r2_v.{k} *. cont_m.{k, l} 
                else beta_r2_v.{k} *. etlin_f_m.{i, j}
              in
              (* r_k against i_l *)
              jac_m.{gr 0 k i, gi 10 l j} <- ~-. beta_r0 ;
              jac_m.{gr 0 k i, gi 20 l j} <- ~-. beta_r0 ;
              jac_m.{gr 0 k i, gi 12 l j} <- ~-. beta_r0 ;
              jac_m.{gr 0 k i, gi 21 l j} <- ~-. beta_r0 ;
              jac_m.{gr 1 k i, gi 20 l j} <- ~-. beta_r1 ;
              jac_m.{gr 1 k i, gi 21 l j} <- ~-. beta_r1 ;
              jac_m.{gr 2 k i, gi 10 l j} <- ~-. beta_r2 ;
              jac_m.{gr 2 k i, gi 12 l j} <- ~-. beta_r2 ;
              (* i_k against i_l *)
              jac_m.{gi 12 k i, gi 10 l j} <- beta_r2 ;
              jac_m.{gi 21 k i, gi 20 l j} <- beta_r1 ;
              if k = l then 
                (jac_m.{gi 10 k i, gi 10 l j} <- beta_r0 -. nu ;
                 jac_m.{gi 10 k i, gi 12 l j} <- beta_r0 +. g2 ;
                 jac_m.{gi 20 k i, gi 20 l j} <- beta_r0 -. nu ;
                 jac_m.{gi 20 k i, gi 21 l j} <- beta_r0 +. g1 ;
                 jac_m.{gi 12 k i, gi 12 l j} <- beta_r2 -. nu -. g2 ;
                 jac_m.{gi 21 k i, gi 21 l j} <- beta_r1 -. nu -. g1)
              else
                (jac_m.{gi 10 k i, gi 10 l j} <- beta_r0 ;
                 jac_m.{gi 10 k i, gi 12 l j} <- beta_r0 ;
                 jac_m.{gi 20 k i, gi 20 l j} <- beta_r0 ;
                 jac_m.{gi 20 k i, gi 21 l j} <- beta_r0 ;
                 jac_m.{gi 12 k i, gi 12 l j} <- beta_r2 ;
                 jac_m.{gi 21 k i, gi 21 l j} <- beta_r1)
            done ;
          done ;
        done ;
      done ; 

      let x = copy ~y:x ~n:(n / 2) ~ofsx:1 y in
      for k = 1 to a do
       ignore (gemv 
           ~n:12 
           ~alpha:1. 
           ~beta:0. 
           ~y:z
           ~ofsy:(1 + (k - 1) * 12)
           ~ofsx:(1 + (k - 1) * 12) 
           lin_f_m.(k - 1) x)
      done ;
      let dx = copy ~y:dx ~n:(n / 2) ~ofsx:(n / 2 + 1) y in
      let z = gemv ~n:(n / 2) ~alpha:1. ~beta:0. ~y:z ~ofsy:(n / 2 + 1) jac_m dx in
      (*
      print_string "start\n" ;
      Vec.iteri (fun i -> fun x -> if i > 36 then (print_float x ; print_string ", ")) z ;
      print_string "\nend\n" ;
      *)
      z

    (* FIXME we'd like two values per city *)
    let aux ?(z=Vec.make0 2) t y =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      (* Age-class infections *)
      let tmp_mc2 = copy ~y:tmp_mc2 sensi_v in
      scal beta tmp_mc2 ; 
      (* tmp_mc2 now contains infectivities for each age class *)
      let tmp_mc2 = Vec.div ~z:tmp_mc2 tmp_mc2 prop_v in
      (* tmp_mc2 now contains the per infectious per susceptible infectivity *)
      (* For the first strain *)
      let i1_v = Vec.add ~n:(a * c) ~z:i1_v ~ofsx:5 ~incx:12 ~ofsy:7 ~incy:12 y y in
      let eta1_v = gemv ~y:eta1_v ~beta:0. ~m:(a * c) etlin_f_m i1_v in
      (* i1_v contains the 1infected individuals (by history) by age class *)
      let tmp_mc1 = copy ~y:tmp_mc1 eta1_v in (* fill tmp_mc1 with eta1 *)
      (* compute the number of contacts : 
       * matrix of contacts * each number of infected *)
      let tmp_mc1 = gemv ~y:tmp_mc1 ~beta:1. ~m:a cont_m i1_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i1_v = Vec.mul ~z:beta_i1_v tmp_mc1 tmp_mc2 in

      (* For the second strain *)
      let i2_v = Vec.add ~n:(a * c) ~z:i2_v ~ofsx:6 ~incx:12 ~ofsy:8 ~incy:12 y y in
      (* i2_v contains the 2infected individuals (by history) by age class *)
      let eta2_v = gemv ~y:eta2_v ~beta:0. ~m:(a * c) etlin_f_m i2_v in
      let tmp_mc1 = copy ~y:tmp_mc1 eta2_v in (* fill tmp_mc1 with eta2 *)
      (* compute the number of contacts : matrix of contacts * each number of infected *)
      let tmp_mc1 = gemv ~y:tmp_mc1 ~beta:1. ~m:a cont_m i2_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i2_v = Vec.mul ~z:beta_i2_v tmp_mc1 tmp_mc2 in

      (* these vectors we need for the jacobian *)
      let r0_v = copy ~n:(a * c) ~y:r0_v ~ofsx:1 ~incx:12 y 
      and r1_v = copy ~n:(a * c) ~y:r1_v ~ofsx:2 ~incx:12 y
      and r2_v = copy ~n:(a * c) ~y:r2_v ~ofsx:3 ~incx:12 y
      in
      let infct1 = dot (Vec.add ~z:tmp_mc1 r0_v r2_v) beta_i1_v
      and infct2 = dot (Vec.add ~z:tmp_mc2 r0_v r1_v) beta_i2_v
      in
      z.{1} <- infct1 *. 7. *. 100000. /. size ;
      z.{2} <- infct2 *. 7. *. 100000. /. size ;
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

    let shift_in_domain ?(z=Vec.make0 n) y =   
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
      @ ["inc1" ; "inc2"]
      @ List.concat
      (List.map (fun s -> f s a []) ["R0_" ; "R1_" ; "R2_" ; "R12_" ; 
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


