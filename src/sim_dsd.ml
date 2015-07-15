open Lacaml.D

(* system parameter values *)

(* Empirical parameter set *)

let size_r = ref (10. ** 5.)
let r0_r = ref (2.)
let e_r = ref (0.15)
let etaN1_r = ref (10. ** (-7.6))
let etaN2_r = ref (10. ** (-7.6))
let g1_r = ref (1. /. (20. *. 365.))
let g2_r = ref (1. /. (20. *. 365.))
let nu_r = ref (1. /. 2.77) 
let q_r = ref (2. /. 365.)


(* Theoretical parameter set *)
(*
let size_r = ref (10. ** 5.)
let r0_r = ref (5.)
let e_r = ref (0.35)
let etaN1_r = ref (10. ** (-.7.))
let etaN2_r = ref (10. ** (-.7.))
let g1_r = ref (1. /. (20. *. 365.))
let g2_r = ref (1. /. (20. *. 365.))
let nu_r = ref (1. /. 8.) 
let q_r = ref (2. /. 365.) 
*)

(* variational system behaviour *)
let init_perturb_r = ref (10. ** (~-. 8.))
let dilat_bound_r = ref 10.

(* simulation arguments *)
let tf_r = ref (365. *. 200.)

let f = fun n -> Random.float 2. -. 1. 
let dx_0 = Array.init 12 f 
let s = Array.fold_left (+.) 0. dx_0 
let dy_0 = Array.map (fun x -> (x -. s)) dx_0
let y0 = Vec.of_array 
           (Array.append 
             [| 0.2 ; 0.2 ; 0.2 ; 0.2 ;
                0.001 ; 0.001 ; 0.001 ; 0.001 ;
                0.049 ; 0.049 ; 0.049 ; 0.049 |]
            dy_0)
  

(* Algorithm parameters *)
let h0_r = ref (1. /. (24. *. 60.))
let delta_r = ref 0.1
let min_step_r = ref (1. /. (24. *. 3600.))
let max_step_r = ref 1.

let verbose_r = ref false

let main () =
  let change_chan_to_file co_r s =
    co_r := open_out s
  in
  let chan_r = ref stdout in
  let specy0 = 
        [Arg.Float (fun x -> y0.{1} <- x);
         Arg.Float (fun x -> y0.{2} <- x);
         Arg.Float (fun x -> y0.{3} <- x);
         Arg.Float (fun x -> y0.{4} <- x);
         Arg.Float (fun x -> y0.{5} <- x);
         Arg.Float (fun x -> y0.{6} <- x);
         Arg.Float (fun x -> y0.{7} <- x);
         Arg.Float (fun x -> y0.{8} <- x);
         Arg.Float (fun x -> y0.{9} <- x);
         Arg.Float (fun x -> y0.{10} <- x);
         Arg.Float (fun x -> y0.{11} <- x);
         Arg.Float (fun x -> y0.{12} <- x)] in
  let specl = 
        [("-v", Arg.Set verbose_r,
                ": verbose output (head and tail) to stdout.");
         ("-dest", Arg.String (change_chan_to_file chan_r),
                ": location of the destination CSV file.\n" ^ 
                "      If not given, outputs to standard output.");
         ("-tf", Arg.Set_float tf_r,
                ": Simulate until (in days)");
         ("-y0", Arg.Tuple specy0,
                ": Initial proportions in each compartment (Should sum to 1)");
         ("-N", Arg.Set_float size_r, 
                ": Total number of hosts in the population");
         ("-R0", Arg.Set_float r0_r, 
                ": Basic reproductive ratio");
         ("-e", Arg.Set_float e_r, 
                ": Strength of the seasonal forcing");
         ("-etaN1", Arg.Set_float etaN1_r, 
                ": Intensity of immigration for strain 1 (per host)");
         ("-etaN2", Arg.Set_float etaN2_r, 
                ": Intensity of immigration for strain 2 (per host)");
         (* it doesn't make sense to use both -etaN and -etaN1 or -etaN2 *)
         ("-etaN", Arg.Float (fun x -> etaN1_r := x; etaN2_r := x),
                ": Intensity of immigration for strains 1 and 2 (per host).");
         ("-g1", Arg.Set_float g1_r, 
                ": Frequency of immunity loss for strain 1 (1/days)");
         ("-g2", Arg.Set_float g2_r, 
                ": Frequency of immunity loss for strain 2 (1/days)");
         (* it doesn't make sense to use both -g and -g1 or -g2 *)
         ("-g", Arg.Float (fun x -> g1_r := x; g2_r := x),
                ": Frequency of immunity loss for strains 1 and 2 (per host).");
         ("-nu", Arg.Set_float nu_r, 
                ": Frequency of recovery from infection (1/days)");
         ("-q", Arg.Set_float q_r, 
                ": Frequency of loss of cross immunity (1/days)");
         ("-init_perturb", Arg.Set_float init_perturb_r,
                ": Initial norm (1) of the perturbation");
         ("-dil", Arg.Set_float dilat_bound_r, 
                ": Dilatation factor before rescaling the variational system");
         ("-h0", Arg.Set_float h0_r,
                ": Initial step size");
         ("-delta", Arg.Set_float delta_r,
                ": Maximum tolerated local error. Must be > 0.");
         ("-min_step", Arg.Set_float min_step_r,
                ": Minimum resolution in time. Must be > 0.");
         ("-max_step", Arg.Set_float max_step_r,
                ": Maximum resolution in time. Must be > 0.")]
  in
  (* simply ignore anonymous arguments *)
  let anon_print s = print_endline ("Ignored anonymous argument : " ^ s) in
  (* printed before the help message : *)
  let usage_msg = "  Simulate using Dopri5(.ml) a double strain " ^
                  "seasonally forced SIR model approximating (for example) " ^
                  "influenza dynamics." ^
                  "\nFor more info, look into dsd.mli and dopri5.mli.\n" ^
                  "Available options :" in
  (* parse the command line and update the parameter values *)
  Arg.parse specl anon_print usage_msg ;
  (* sanity check *)
  if (1. -. 10. ** (~-. !size_r -. 1.)  < Vec.sum ~n:12 y0) 
  && (Vec.sum ~n:12 y0 < 1. -. 10. ** (~-. !size_r -. 1.))  then 
    failwith "The user did not pass a proportion tuple (sums to 1) as y0 : \n" ;

  (* We scale the population size appropriately *)
  scal ~n:12 ~ofsx:1 !size_r y0 ;
  (* We scale the perturbation appropriately *)
  scal ~n:12 ~ofsx:13 !init_perturb_r y0 ;
  let module Pars = 
    struct
      let size = !size_r
      let r0 = !r0_r
      let e = !e_r
      let etaN1 = !etaN1_r
      let etaN2 = !etaN2_r
      let g1 = !g1_r
      let g2 = !g2_r
      let nu = !nu_r
      let q = !q_r
      let init_perturb = !init_perturb_r
      let dilat_bound = !dilat_bound_r
    end
  in
  let module DsdSys = Dsd.Sys (Pars) in
  let module Algp =
    struct
      let h0 = !h0_r
      let delta = !delta_r
      let min_step = !min_step_r
      let max_step = !max_step_r
    end
  in
  if !verbose_r = true then
    (* change to a formatted print *)
    print_string ("The parameter values used are\n"
      ^ "N :" ^ (string_of_float !size_r) ^ "\n"
      ^ "R0 :" ^ (string_of_float !r0_r) ^ "\n"
      ^ "e :" ^ (string_of_float !e_r) ^ "\n"
      ^ "etaN1 :" ^ (string_of_float !etaN1_r) ^ "\n"
      ^ "etaN2 :" ^ (string_of_float !etaN2_r) ^ "\n"
      ^ "g1 :" ^ (string_of_float !g1_r) ^ "\n"
      ^ "g2 :" ^ (string_of_float !g2_r) ^ "\n"
      ^ "nu :" ^ (string_of_float !nu_r) ^ "\n"
      ^ "q :" ^ (string_of_float !q_r) ^ "\n") ;
  let module Gen = Dopri5.Integrator (DsdSys) (Algp) in
  Gen.simulate !chan_r !tf_r y0

let () = ignore (main ())
