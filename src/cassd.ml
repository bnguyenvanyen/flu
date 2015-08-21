
(* Simulate the single-strain age-structured city-structured ODE system using Dopri5 *)
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
    val g : float
    val nu : float
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
            && g >= 0. 
            && nu >= 0.
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

    let n = 3 * a * c * 2
    let m = 1 (* for now ... *)

    let bet0 = r0 *. nu

    let ng = ~-. g
    let nnu = ~-. nu
    let a_base_m = Mat.of_array (* to simplify computations *)
      [| 
        [| 0. ; 0.  ; g  |] ;
        [| 0. ; nnu ; 0. |] ;
        [| 0. ; nu  ; ng |]
      |];;
  
    let lin_f_m = Array.make (a * c) a_base_m

    let jac_base_m = Mat.of_array (* there's no copy function :( *)
      [| 
        [| 0. ; 0.  ; g  |] ;
        [| 0. ; nnu ; 0. |] ;
        [| 0. ; nu  ; ng |]
      |];;

    let jac_m = Mat.make0 (n / 2) (n / 2);;
    (* the real "base j" is this on the diagonal of the age-clases *)
    for i = 1 to c do
      for k = 1 to a do
        ignore (lacpy 
                  ~b:jac_m
                  ~br:(1 + (i - 1) * 3 * a + (k - 1) * 3) 
                  ~bc:(1 + (i - 1) * 3 * a + (k - 1) * 3) 
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
    done ;

    (* prop_v will now contain the number of hosts in each age class 
     * (and not the proportion) *)
    scal size prop_v

    (* allocate stuff *)
    (* initial values (by age-class) *)
    let s_v = Vec.make0 (a * c)
    let i_v = Vec.make0 (a * c)
    let r_v = Vec.make0 (a * c)
    (* infectivities by susceptible *)
    let beta_i_v = Vec.make0 (a * c)
    let beta_s_v = Vec.make0 (a * c)
    (* end values *)
    let s_dot_v = Vec.make0 (a * c)
    let i_dot_v = Vec.make0 (a * c)
    let r_dot_v = Vec.make0 (a * c)

    let x = Vec.make0 (n / 2)
    let dx = Vec.make0 (n / 2) (* the dr, di, dq values *)

    (* same. FIXME can we drop some ? *)
    let eta_v = Vec.make0 (a * c)
    let tmp_ac1 = Vec.make0 (a * c)
    let tmp_ac2 = Vec.make0 (a * c) ;;

    let f ?(z=Vec.make0 n) t y =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in

      (* Age-class infections *)
      let tmp_ac2 = copy ~y:tmp_ac2 sensi_v in
      scal beta tmp_ac2 ; 
      (* tmp_ac2 now contains infectivities for each age class *)
      let tmp_ac2 = Vec.div ~z:tmp_ac2 tmp_ac2 prop_v in
      (* tmp_ac2 now contains the per infectious per susceptible infectivity *)

      let i_v = copy ~n:(a * c) ~y:i_v ~ofsx:2 ~incx:3 y in
      let eta_v = gemv ~y:eta_v ~beta:0. ~m:(a * c) etlin_f_m i_v in
      (* i_v contains the infected individuals (by history) by age class *)
      let tmp_ac1 = copy ~y:tmp_ac1 eta_v in (* fill tmp_ac1 with eta *)
      (* compute the number of contacts : 
       * matrix of contacts * each number of infected *)
      let tmp_ac1 = gemv ~y:tmp_ac1 ~beta:1. ~m:a cont_m i_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i_v = Vec.mul ~z:beta_i_v tmp_ac1 tmp_ac2 in

      for i = 1 to c do
        for k = 1 to a do
          let ind = (i - 1) * a + k in
          let beta_i = beta_i_v.{ind} in
          (* the infection terms are not linear and vary in time *)
          lin_f_m.(ind - 1).{1, 1} <- ~-. beta_i ;
          lin_f_m.(ind - 1).{2, 1} <- beta_i ;
          (* in the jacobian too *)
          (* but in the jacobian we keep all terms around ... *)
          (* r_k against itself *)
          jac_m.{(ind - 1) * 3 + 1, (ind - 1) * 3 + 1} <- ~-. beta_i ;
          (* i_k against r_k *)
          jac_m.{(ind - 1) * 3 + 2, (ind - 1) * 3 + 1} <- beta_i ;
          for j = 1 to c do
            for l = 1 to a do
              let ind_jl = (j - 1) * a + l in
              (* I guess *)
              let beta_s =
                if i = j then beta_s_v.{k} *. cont_m.{k, l}
                else beta_s_v.{k} *. etlin_f_m.{i, j}
              in
              (* s_k against i_l *)
              jac_m.{(ind - 1) * 3 + 1, (ind_jl - 1) * 3 + 2} <- ~-. beta_s ;
              (* i_k against i_l *)
              jac_m.{(ind - 1) * 3 + 2, (ind_jl - 1) * 3 + 2} <- beta_s ;
              if k = l then 
                (jac_m.{(ind - 1) * 3 + 2, (ind_jl - 1) * 3 + 2} <- beta_s -. nu)
              else
                (jac_m.{(ind - 1) * 3 + 2, (ind_jl - 1) * 3 + 2} <- beta_s)
            done ;
          done ;
        done ;
      done ; 

      let x = copy ~y:x ~n:(n / 2) ~ofsx:1 y in
      for k = 1 to a do
       ignore (gemv 
           ~n:3
           ~alpha:1. 
           ~beta:0. 
           ~y:z
           ~ofsy:(1 + (k - 1) * 3)
           ~ofsx:(1 + (k - 1) * 3) 
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
    let aux ?(z=Vec.make0 1) t y =
      let beta = bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.)) in
      (* Age-class infections *)
      let tmp_ac2 = copy ~y:tmp_ac2 sensi_v in
      scal beta tmp_ac2 ; 
      (* tmp_ac2 now contains infectivities for each age class *)
      let tmp_ac2 = Vec.div ~z:tmp_ac2 tmp_ac2 prop_v in
      (* tmp_ac2 now contains the per infectious per susceptible infectivity *)
      let i_v = copy ~n:(a * c) ~y:i_v ~ofsx:2 ~incx:3 y in
      let eta_v = gemv ~y:eta_v ~beta:0. ~m:(a * c) etlin_f_m i_v in
      (* i1_v contains the 1infected individuals (by history) by age class *)
      let tmp_ac1 = copy ~y:tmp_ac1 eta_v in (* fill tmp_mc1 with eta1 *)
      (* compute the number of contacts : 
       * matrix of contacts * each number of infected *)
      let tmp_ac1 = gemv ~y:tmp_ac1 ~beta:1. ~m:a cont_m i_v in
      (* compute the "per susceptible" number of infections *)
      let beta_i_v = Vec.mul ~z:beta_i_v tmp_ac1 tmp_ac2 in

      (* these vectors we need for the jacobian *)
      let s_v = copy ~n:(a * c) ~y:s_v ~ofsx:1 ~incx:3 y 
      in
      let infct = dot s_v beta_i_v
      in
      z.{1} <- infct *. 7. *. 100000. /. size ;
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
      let rec strange n =
        match n > 0 with
        | true ->
          (strange (n - 1)) @ [string_of_int n]
        | false ->
          []
      in
      let rec f s_l_l s =
        match s_l_l with
        | [] ->
          [s]
        | [] :: tl ->
          []
        | (s_start :: s_l_red) :: tl ->
          (f tl (s_start ^ s)) @ (f (s_l_red :: tl) s)
      in
      ["t" ; "h"] 
      @ ["inc"]
      @ (f [strange c ; ["_"] ; strange a ; ["_"] ; ["S" ; "I" ; "R"]] "")
      @ (f [strange c ; ["_"] ; strange a ; ["_"] ; ["dS" ; "dI" ; "dR"]] "")
  end;;

module Default_Algp =
  struct
    let h0 = 1. /. (24. *. 60.);;
    let delta = 0.1;; (* an error of 0.1 person seems ok *)
    let min_step = 1. /. (24. *. 3600.);; (* more than one step a second looks overkill *)
    let max_step = 1.;; (* we always want at least a resolution of a day *)
  end;;


