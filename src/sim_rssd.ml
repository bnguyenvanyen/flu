open Lacaml.D

let load_vec_from_file v_r fname =
  let data = Csv.load fname in
  match data with
  | [] -> 
      failwith "no data read from file"
  | l1 :: l2 :: _ -> 
      failwith "file has more than one line"
  | l1 :: [] ->
      let l = List.map (fun s -> float_of_string s) l1 in
      let v = Vec.of_list l in
      v_r := v

(* system parameter values *)
let size_r = ref (10. ** 5.)
let r0_r = ref (2.)
let e_r = ref (0.15)
let b_r = ref (0.)
let etaN_r = ref (10. ** (-7.1))
let phi_r = ref (0.)
let g_r = ref (1. /. (10. *. 365.))
let nu_r = ref (1. /. 2.77)
let nr_r = ref 3
let r_prop_v_r = ref (Vec.of_array [| 1. /. 3. ; 1. /. 3. ; 1. /. 3. |])

(* variational system behaviour *)
let init_perturb_r = ref (10. ** (~-. 8.))
let dilat_bound_r = ref 10.

(* simulation arguments *)
let tf_r = ref (365. *. 200.)


(* Algorithm parameters *)
let h0_r = ref (1. /. (24. *. 60.))
let delta_r = ref 0.1
let min_step_r = ref (1. /. (24. *. 3600.))
let max_step_r = ref 1.

let sir_v = Vec.of_array 
             [| 0.5 ; 0.001 ; 0.499 |]

let main () =
  let change_chan_to_file co_r s =
    co_r := open_out s
  in
  let chan_r = ref stdout in
  let spec_sir = 
        [Arg.Float (fun x -> sir_v.{1} <- x);
         Arg.Float (fun x -> sir_v.{2} <- x);
         Arg.Float (fun x -> sir_v.{3} <- x);] in
  let specl = 
        [("-dest", Arg.String (change_chan_to_file chan_r),
                ": location of the destination CSV file.\n" ^ 
                "      If not given, outputs to standard output.");
         ("-tf", Arg.Set_float tf_r,
                ": Simulate until (in days)");
         ("-N", Arg.Set_float size_r, 
                ": Total number of hosts in the population");
         ("-R0", Arg.Set_float r0_r, 
                ": Basic reproductive ratio");
         ("-e", Arg.Set_float e_r, 
                ": Strength of the seasonal forcing");
         ("-b", Arg.Set_float b_r,
                ": Strength of the immigration seasonal forcing");
         ("-etaN", Arg.Set_float etaN_r, 
                ": Intensity of immigration (per host)");
         ("-phi", Arg.Set_float phi_r,
                ": Phase shift between infectivity and immigration");
         ("-g", Arg.Set_float g_r, 
                ": Frequency of immunity loss (1/days)");
         ("-nu", Arg.Set_float nu_r, 
                ": Frequency of recovery from infection (1/days)");
         ("-nr", Arg.Set_int nr_r,
                ": Number of removed classes");
         ("-sir_prop", Arg.Tuple spec_sir,
                ": Initial proportions of S, I and R");
         ("-r_prop", Arg.String (load_vec_from_file r_prop_v_r),
                ": Initial proportions of the removed compartments");
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
  let usage_msg = "  Simulate using Dopri5(.ml) a single strain " ^
                  "seasonally forced SIR model approximating (for example) " ^
                  "influenza dynamics, with multiple removed compartments." ^
                  "\nFor more info, look into rssd.mli and dopri5.mli.\n" ^
                  "Available options :" in
  (* parse the command line and update the parameter values *)
  Arg.parse specl anon_print usage_msg ;
  (* sanity check *)
  if (1. +. 1. /. (!size_r *. 10.) < Vec.sum ~n:3 sir_v) 
  || (Vec.sum ~n:3 sir_v < 1. -. 1. /. (!size_r *. 10.))  then 
    failwith "The user did not pass a proportion tuple (sums to 1) as y0 \n" ;
  let f = fun n -> Random.float 2. -. 1. in
  let rdu = Array.init (2 + !nr_r) f in
  let s = Array.fold_left (+.) 0. rdu in
  let dx0 = Vec.of_array (Array.map (fun x -> (x -. s)) rdu) in
  scal ~n:(2 + !nr_r) ~ofsx:1 !init_perturb_r dx0 ;
  let y0 = Vec.make0 (2 * (2 + !nr_r)) in
  y0.{1} <- sir_v.{1} *. !size_r ;
  y0.{2} <- sir_v.{2} *. !size_r ;
  for i = 1 to !nr_r do
    y0.{2 + i} <- sir_v.{3} *. (!r_prop_v_r).{i} *. !size_r
  done ;
  ignore (copy ~n:(2 + !nr_r) 
               ~ofsy:(3 + !nr_r) 
               ~y:y0
               dx0) ;
  let module Pars = 
    struct
      let size = !size_r
      let r0 = !r0_r
      let e = !e_r
      let b = !b_r
      let etaN = !etaN_r
      let phi = !phi_r
      let g = !g_r
      let nu = !nu_r
      let nr = !nr_r
      let init_perturb = !init_perturb_r
      let dilat_bound = !dilat_bound_r
    end
  in
  let module RssdSys = Rssd.Sys (Pars) in
  let module Algp =
    struct
      let h0 = !h0_r
      let delta = !delta_r
      let min_step = !min_step_r
      let max_step = !max_step_r
    end
  in
  let module Gen = Dopri5.Integrator (RssdSys) (Algp) in
  Gen.simulate !chan_r !tf_r y0

let () = ignore (main ())
