(* A bit of a simplified port of dopri5.f : less options *)
(* Implements the Dormand-Prince method of order 5 (4) *)

open Lacaml.D;;

module type SYSTEM =
  sig
    val n : int
    val m : int
    val f : 
            ?z:Lacaml_float64.vec -> 
            float -> 
            Lacaml_float64.vec -> 
            Lacaml_float64.vec
    val aux : 
            ?z:Lacaml_float64.vec ->
            float ->
            Lacaml_float64.vec ->
            Lacaml_float64.vec  
    val check_in_domain :
            Lacaml_float64.vec ->
            bool
    val shift_in_domain : 
            ?z:Lacaml_float64.vec ->
            Lacaml_float64.vec ->
            Lacaml_float64.vec
    val csv_init :
            unit ->
            string list
  end;;

module type ALGPARAMS =
  sig
    val h0 : float
    val delta : float
    val min_step : float
    val max_step : float
  end;;

module type INTEGR =
  sig
   (* 
    val update_k : 
          float -> 
          float -> 
          Lacaml_float64.mat -> 
          Lacaml_float64.vec -> 
          int -> 
          Lacaml_float64.mat
    val compute_k :
          float ->
          float ->
          Lacaml_float64.mat ->
          Lacaml_float64.vec ->
          Lacaml_float64.mat
    val runge_kutta : 
          float ->
          Lacaml_float64.mat -> 
          Lacaml_float64.vec -> 
          Lacaml_float64.vec -> 
          Lacaml_float64.vec ->
          Lacaml_float64.vec
    val loop : 
          float -> 
          float -> 
          Lacaml_float64.mat ->
          Lacaml_float64.vec ->
          Lacaml_float64.vec ->
          Lacaml_float64.vec ->
          (float * float * Lacaml_float64.vec)
   *)
    val simulate : 
          out_channel ->
          float -> 
          Lacaml_float64.vec -> 
          (float * Lacaml_float64.vec)
  end;;

(* declare the coefficients in the top namespace *)
let c = Vec.of_array 
        [| 
          0.; 1. /. 5.; 3. /. 10.; 4. /. 5.; 8. /. 9.; 1.; 1. 
        |];;
let b5 = Vec.of_array 
         [| 
           35. /. 384. ; 
           0. ; 
           500. /. 1113. ; 
           125. /. 192. ; 
           ~-. 2187. /. 6784. ;
           11. /. 84. ; 
           0. 
         |];;

let b4 = Vec.of_array 
          [|
            5179. /. 57600. ;
            0. ;
            7571. /. 16695. ;
            393. /. 640. ;
            ~-. 92097. /. 339200. ;
            187. /. 2100. ;
            1. /. 40.
          |]

let a = Mat.of_array 
        [|
          [| 0.; 0.; 0.; 0.; 0.; 0.; 0. |] ;
          [| 1. /. 5.; 0.; 0.; 0.; 0.; 0.; 0. |] ;
          [| 3. /. 40.; 9. /. 40.; 0.; 0.; 0.; 0.; 0. |] ;
          [| 44. /. 45.; ~-. 56. /. 15.; 32. /. 9.; 0.; 0.; 0.; 0. |] ;
          [| 19372. /.6561.; ~-. 25360. /. 2187.; 64448. /. 6561.; 
            ~-. 212. /. 729.; 0.; 0.; 0. |] ; 
          [| 9017. /. 3168.; ~-. 355. /. 33.; 46732. /. 5247.; 
            49. /. 176.; ~-. 5103. /. 18656.; 0.; 0. |] ;
          [| 35. /. 384.; 0.; 500. /. 1113.; 125. /. 192.; 
            ~-. 2187. /. 6784.; 11. /. 84.; 0. |] 
        |]	

(* create a unit test suite somewhere *)
module Integrator (Sys : SYSTEM) (Algp : ALGPARAMS) : INTEGR =
  struct
    open Algp

    let tmp1 = Vec.make0 Sys.n (* a temp storage vector of size n *)
    let tmp2 = Vec.make0 Sys.n (* a temp storage vector of size n *)

    let update_k t h k y i =
      let tmp1 = copy ~y:tmp1 y in
      let tmp1 = gemv ~trans:`N ~alpha:h ~beta:1. ~y:tmp1 k (Mat.copy_row a i) in
      let tmp2 = Sys.f ~z:tmp2 (t +. h *. c.{i}) tmp1 in
      for j = 1 to Sys.n do
        k.{j,i} <- tmp2.{j}
      done;
      k
    
    let compute_k t h k y =
      (* we start with a clean k again *)
      Mat.fill k 0.;
      for i = 1 to 7 do
        ignore (update_k t h k y i) ;
      done;
      k
    
    (* fourth-order evaluation *)
    let runge_kutta h k ny y b =
      let tmp1 = copy ~y:tmp1 y in
      let tmp1 = gemv ~trans:`N ~alpha:h ~beta:1. ~y:tmp1 k b in
      copy ~y:ny tmp1
    
    (* compute dt *)
    
    let loop t h k y4 y5 y z =
      let k = compute_k t h k y in
      let y4 = runge_kutta h k y4 y b4 in
      let y5 = runge_kutta h k y5 y b5 in
      let tmp1 = Vec.sub ~z:tmp1 y5 y4 in (* tmp gets y5 - y4 *)
      let e = abs_float (amax tmp1) in      (* or 1-norm ? *)
      let t = t +. h in
      (*
      let h = 
        if e = 0. then max_step else 0.9 *. h *. (delta /. e) ** (1. /. 5.) in
      *)
      let h = 0.9 *. h *. (delta /. e) ** (1. /. 5.) in
      let h = 
        if h < min_step then min_step else 
        if (h > max_step) || ((compare h nan) = 0) then max_step else 
        h in
      let y = copy ~y y5 in
      let y =
        (* can we manage to have this for free with an exception ? *)
        if not (Sys.check_in_domain y) then 
          Sys.shift_in_domain ~z:y y
        else
          y
      in
      let z = Sys.aux ~z t y
      in 
      (t, h, y, z)

    let csv_info () =
        ["n=" ^ string_of_int Sys.n ; "m=" ^ string_of_int Sys.m]

    let csv_line t h y z =
      let l = ref [] in
      Vec.iter (fun x -> l := Printf.sprintf "%.3e" x :: !l) z ;
      Vec.iter (fun x -> l := Printf.sprintf "%.3e" x :: !l) y ;
      [Printf.sprintf "%f" t; Printf.sprintf "%.2e" h] @ (List.rev !l)

    let simulate st_chan tf y0 =
      (* outputs the data as csv to chan *)
      let chan = Csv.to_channel st_chan in
      Csv.output_record chan (csv_info ()) ;
      Csv.output_record chan (Sys.csv_init ()) ;
      let t = ref 0. in
      let print_t = ref 0. in
      let h = ref Algp.h0 in
      let y4 = copy y0 in
      let y5 = copy y0 in
      let y = copy y0  in
      let z = Sys.aux 0. y0 in
      let k = Mat.make0 Sys.n 7 in
      while !t < tf do
        let nt, nh, y, z = loop !t !h k y4 y5 y z in
        (* We want to print only every max_step *)
        if nt -. !print_t > Algp.max_step then 
          (Csv.output_record chan (csv_line nt nh y z) ;
           print_t := nt) ;
        t := nt ;
        h := nh ;
      done;
      (!t, y) 
  end;;


