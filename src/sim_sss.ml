(* system parameter values *)
let n_r = ref (10. ** 6.)
let rr0_r = ref (2.)
let e_r = ref (0.15)
let etaN_r = ref (10. ** (-7.1))
let g_r = ref (1. /. (14. *. 365.))
let nu_r = ref (1. /. 2.77)

(* simulation arguments *)
let dest_r = ref "./sim_sss_default_dest.csv"
let tf_r = ref (365. *. 10.)
(* FIXME need to recompute y0 later (if size_r has been changed) *)

let s0_r = ref 0.5
let i0_r = ref 0.001
let r0_r = ref 0.499

(* Algorithm parameters *)
let min_step_r = ref 1.

let main () =
  let change_chan_to_file co_r s =
    co_r := open_out s
  in
  let chan_r = ref stdout in
  let specy0 =
        [Arg.Float (fun x -> s0_r := x);
         Arg.Float (fun x -> i0_r := x);
         Arg.Float (fun x -> r0_r := x);] in
  let specl = 
        [("-dest", Arg.String (change_chan_to_file chan_r),
                ": location of the destination CSV file");
         ("-tf", Arg.Set_float tf_r,
                ": Simulate until (in days)");
         ("-y0", Arg.Tuple specy0,
                ": Initial conditions (Should sum to 1)");
         ("-N", Arg.Set_float n_r, 
                ": Total number of hosts in the population");
         ("-R0", Arg.Set_float rr0_r, 
                ": Basic reproductive ratio");
         ("-e", Arg.Set_float e_r, 
                ": Strength of the seasonal forcing");
         ("-etaN", Arg.Set_float etaN_r, 
                ": Intensity of immigration (per host)");
         ("-g", Arg.Set_float g_r, 
                ": Frequency of immunity loss (1/days)");
         ("-nu", Arg.Set_float nu_r, 
                ": Frequency of recovery from infection (1/days)");
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
  if (1. -. 1. /. (!n_r *. 10.)  < !s0_r +. !i0_r +. !r0_r ) 
  && (!s0_r +. !i0_r +. !r0_r < 1. -. 1. /. (!n_r *. 10.))  then 
    failwith "The initial compartment proportions don't sum to 1 \n" ;

  let f x = int_of_float (!x *. !n_r) in
  let (s0, i0, r0) = (f s0_r, f i0_r, f r0_r) in
  let module Pars = 
    struct
      let n = !n_r
      let r0 = !rr0_r
      let e = !e_r
      let etaN = !etaN_r
      let g = !g_r
      let nu = !nu_r
    end
  in
  let module SssSys = Sss.Sys (Pars) in
  let module Algp =
    struct
      let min_step = !min_step_r
    end
  in
  let module Gen = Gill.Integrator (SssSys) (Algp) in
  ignore(Gen.simulate !chan_r !tf_r (s0, i0, r0))

let () = main ()
