(* system parameter values *)
let n_r = ref (10. ** 6.);;
let r0_r = ref (2.);;
let e_r = ref (0.15);;
let etaN1_r = ref (10. ** (-7.6));;
let etaN2_r = ref (10. ** (-7.6));;
let g1_r = ref (1. /. (10. *. 365.));;
let g2_r = ref (1. /. (10. *. 365.));;
let nu_r = ref (1. /. 2.77);;
let q_r = ref (2. /. 365.);;

(* simulation arguments *)
let dest_r = ref "./sim_sss_default_dest.csv"
let tf_r = ref (365. *. 100.);;
(* FIXME need to recompute y0 later (if size_r has been changed) *)
let y0 =   [| 0.2 ; 0.2 ; 0.2 ; 0.2 ;
              0.001 ; 0.001 ; 0.001 ; 0.001 ;
              0.049 ; 0.049 ; 0.049 ; 0.049 |]
(* Algorithm parameters *)
let min_step_r = ref 1.;;

let main () =
  let change_chan_to_file co_r s =
    co_r := open_out s
  in
  let chan_r = ref stdout in
  let specy0 =
        [Arg.Float (fun x -> y0.(1) <- x);
         Arg.Float (fun x -> y0.(2) <- x);
         Arg.Float (fun x -> y0.(3) <- x);
         Arg.Float (fun x -> y0.(4) <- x);
         Arg.Float (fun x -> y0.(5) <- x);
         Arg.Float (fun x -> y0.(6) <- x);
         Arg.Float (fun x -> y0.(7) <- x);
         Arg.Float (fun x -> y0.(8) <- x);
         Arg.Float (fun x -> y0.(9) <- x);
         Arg.Float (fun x -> y0.(10) <- x);
         Arg.Float (fun x -> y0.(11) <- x);
         Arg.Float (fun x -> y0.(12) <- x)] in
  let specl = 
        [("-dest", Arg.String (change_chan_to_file chan_r),
                ": location of the destination CSV file");
         ("-tf", Arg.Set_float tf_r,
                ": Simulate until (in days)");
         ("-y0", Arg.Tuple specy0,
                ": Initial conditions (Should sum to 1)");
         ("-N", Arg.Set_float n_r, 
                ": Total number of hosts in the population");
         ("-R0", Arg.Set_float r0_r, 
                ": Basic reproductive ratio");
         ("-e", Arg.Set_float e_r, 
                ": Strength of the seasonal forcing");
         ("-etaN1", Arg.Set_float etaN1_r, 
                ": Intensity of immigration for strain 1 (per host)");
         ("-etaN2", Arg.Set_float etaN1_r, 
                ": Intensity of immigration for strain 2 (per host)");
         ("-g1", Arg.Set_float g1_r, 
                ": Frequency of immunity loss for strain 1 (1/days)");
         ("-g2", Arg.Set_float g2_r, 
                ": Frequency of immunity loss for strain 2 (1/days)");
         ("-nu", Arg.Set_float nu_r, 
                ": Frequency of recovery from infection (1/days)");
         ("-q", Arg.Set_float q_r, 
                ": Frequency of cross-immunity loss (1/days)");
         ("-min_step", Arg.Set_float min_step_r,
                ": Minimum resolution in time. Must be > 0.")]
  in
  (* simply ignore anonymous arguments *)
  let anon_print s = print_endline ("Ignored anonymous argument : " ^ s) in
  (* printed before the help message : *)
  let usage_msg = "  Simulate using Gillespie(.ml) a single strain " ^
                  "seasonally forced SIR model approximating (for example) " ^
                  "influenza dynamics." ^
                  "\nFor more info, look into sss.mli and gill.mli" in
  (* parse the command line and update the parameter values *)
  Arg.parse specl anon_print usage_msg ;
  (* sanity check *)
  if (1. +. 1. /. (!n_r *. 10.)  < Array.fold_left (+.) 0. y0 ) 
  || (Array.fold_left (+.) 0. y0 < 1. -. 1. /. (!n_r *. 10.))  then 
    failwith "The initial compartment proportions don't sum to 1 \n" ;
  let y0 = Array.map (fun x -> int_of_float (x *. !n_r)) y0 in
  let module Pars = 
    struct
      let n = !n_r
      let r0 = !r0_r
      let e = !e_r
      let etaN1 = !etaN1_r
      let etaN2 = !etaN2_r
      let g1 = !g1_r
      let g2 = !g2_r
      let nu = !nu_r
      let q = !q_r
    end
  in
  let module DssSys = Dss.Sys (Pars) in
  let module Algp =
    struct
      let min_step = !min_step_r
    end
  in
  let module Gen = Gill.Integrator (DssSys) (Algp) in
  ignore(Gen.simulate !chan_r !tf_r y0)

let () = main ()
