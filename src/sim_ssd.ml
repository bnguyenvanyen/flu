open Lacaml.D;;

(* system parameter values *)
let size_r = ref (10. ** 5.);;
let r0_r = ref (2.);;
let e_r = ref (0.15);;
let etaN_r = ref (10. ** (-7.1));;
let g_r = ref (1. /. (14. *. 365.));;
let nu_r = ref (1. /. 2.77);;

(* variational system behaviour *)
let init_perturb_r = ref (10. ** (~-. 6.));;
let dilat_bound_r = ref 100.;;
let contract_bound_r = ref 100.;;

(* simulation arguments *)
let dest_r = ref "./sim_ssd_default_dest.csv"
let tf_r = ref (365. *. 10.);;

let ds_0 = Random.float 2. -. 1.;;
let di_0 = Random.float 2. -. 1.;;
let dr_0 = ~-. ds_0 -. di_0;;
let y0 = Vec.of_array 
          [| 
            0.5 *. !size_r ; 0.001 *. !size_r ; 0.499 *. !size_r ;
            ds_0 *. !init_perturb_r ; di_0 *. !init_perturb_r ; dr_0 *. !init_perturb_r
          |];;

(* Algorithm parameters *)
let h0_r = ref (1. /. (24. *. 60.));;
let delta_r = ref 0.1;;
let min_step_r = ref (1. /. (24. *. 3600.));;
let max_step_r = ref 7.;;

let main () =
  let specy0 = 
        [Arg.Float (fun x -> y0.{1} <- x);
         Arg.Float (fun x -> y0.{2} <- x);
         Arg.Float (fun x -> y0.{3} <- x);] in
  let specl = 
        [("-dest", Arg.Set_string dest_r,
                ": location of the destination CSV file");
         ("-tf", Arg.Set_float tf_r,
                ": Simulate until (in days)");
         ("-y0", Arg.Tuple specy0,
                ": Initial conditions (Should sum to 1)");
         ("-N", Arg.Set_float size_r, 
                ": Total number of hosts in the population");
         ("-R0", Arg.Set_float r0_r, 
                ": Basic reproductive ratio");
         ("-e", Arg.Set_float e_r, 
                ": Strength of the seasonal forcing");
         ("-etaN", Arg.Set_float etaN_r, 
                ": Intensity of immigration (per host)");
         ("-g", Arg.Set_float g_r, 
                ": Frequency of immunity loss (1/days)");
         ("-nu", Arg.Set_float nu_r, 
                ": Frequency of recovery from infection (1/days)");
         ("-init_perturb", Arg.Set_float init_perturb_r,
                ": Initial norm (1) of the perturbation");
         ("-dil", Arg.Set_float dilat_bound_r, 
                ": Dilatation factor before rescaling the variational system");
         ("-contr", Arg.Set_float contract_bound_r,
                ": Contraction factor before rescaling the variational system");

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
                  "influenza dynamics." ^
                  "\nFor more info, look into ssd.mli and dopri5.mli" in
  (* parse the command line and update the parameter values *)
  Arg.parse specl anon_print usage_msg ;
  (* sanity check *)
  let init_sz = y0.{1} +. y0.{2} +. y0.{3} in
  if not (init_sz = !size_r) then
    failwith ("The announced population size is not equal to the initial population size : \n" 
              ^ (string_of_float !size_r) ^ " != " ^ (string_of_float init_sz));
  let module Pars = 
    struct
      let size = !size_r
      let r0 = !r0_r
      let e = !e_r
      let etaN = !etaN_r
      let g = !g_r
      let nu = !nu_r
      let init_perturb = !init_perturb_r
      let dilat_bound = !dilat_bound_r
      let contract_bound = !contract_bound_r
    end
  in
  let module SsdSys = Ssd.Sys (Pars) in
  let module Algp =
    struct
      let h0 = !h0_r
      let delta = !delta_r
      let min_step = !min_step_r
      let max_step = !max_step_r
    end
  in
  let module Gen = Dopri5.Integrator (SsdSys) (Algp) in
  Gen.simulate !dest_r !tf_r y0

let () = ignore (main ())
