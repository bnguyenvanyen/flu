type time = float;;

let rec rand_float top =
  let u = Random.float top in
  match u with
  | 0. -> rand_float top
  | _ -> u;;

let rand_exp lambda =
  let u = rand_float 1. in
  -1. *. (log u) /. lambda;;

let rec successive_sum l =
  match l with
  | a :: b :: tl -> a :: (successive_sum ((a +. b) :: tl))   
  | hd :: [] -> l
  | [] -> l (* this case should never happen *);;

let choose c p =
  let ss = successive_sum p in
  let u = rand_float 1. in
  let rec f l1 l2 =
    match l1, l2 with
    | (_ :: _, []) | ([], _ :: _) -> 
        failwith "the two lists have different lengths"
    | [], [] -> 
        failwith "empty arguments"
    | hd1 :: tl1, hd2 :: tl2 ->
        (match u < hd2 with
         | true -> hd1 
         | false -> f tl1 tl2)
  in
  f c ss

let multinomial n c p =
  let rec f n l =
    match n with (* problem if n < 0 *)
    | 0 -> l
    | _ -> f (n - 1) (choose c p :: l)
  in f n []


module type SYSTEM =
  sig
    type state
    val csv_init : unit -> string list
    val csv_line : time -> state -> string list
    val fl : (time -> state -> float) list
    val ml : (state -> state) list
  end;;

module type ALGPARAMS =
  sig
    val step_size : float
  end;;

module type INTEGR =
  sig
    type state
    type process = {past : (time * state) list ;
                    present : time * state ; 
                    future : (time * state) list}
    val backward : process -> process
    val forward : process -> process
    val csv_of_proc : process -> string list list
    val simulate :
          out_channel ->
          float ->
          state ->
          (time * state)
  end;;

module Integrator (Sys : SYSTEM) (Algp : ALGPARAMS) : 
  (INTEGR with type state = Sys.state) =
  struct
    type state = Sys.state
    type process = {past : (time * state) list ;
                    present : time * state ;
                    future : (time * state) list}
  
    (* Base functions *)
    (* Random generator *)
    
    let backward p =
      match p.past, p.future with
      | [], _ -> invalid_arg "Initial state has been reached"
      | hd :: tl, l -> {past = tl ; present = hd ; future = p.present :: l}
    
    let forward p =
      match p.past, p.future with
      | _, [] -> invalid_arg "Final state has been reached"
      | l, hd :: tl -> {past = p.present :: l ; present = hd ; future = tl}
    
    (* tail_recursive*)
    let rec backward_to_start p =
      let np = backward p in
      match np.past with
      | [] -> np
      | _  -> backward_to_start np
    
    (* tail_recursive *)
    let rec forward_to_end p =
      let np = forward p in
      match np.future with
      | [] -> np
      | _ -> forward_to_end np

    let rec choose_modif w srl ml =
      match srl, ml with
      | [], [] -> 
          failwith "w is out of bounds"
      | ([], _ :: _) | (_ :: _, []) -> 
          failwith "rate list and modif list are incompatible"
      | sr :: tl1, m :: tl2 -> 
          (match w < sr with
           | true -> m
           | false -> choose_modif w tl1 tl2)

    let draw_future prst =
      let t, st = prst in
      let rl = List.map (fun f -> f t st) Sys.fl in
      let srl = successive_sum rl in
      let lambd = List.hd (List.rev srl) in
      let next_t = t +. Algp.step_size in
      let n_events = int_of_float (Algp.step_size *. lambd) in
      let modifs = multinomial n_events Sys.ml rl in
      let next_st = List.fold_left (fun st -> fun m -> m st) st modifs  in
      (next_t, next_st)
      
    let csv_of_proc proc =
      let rec add_line proc csv =
        let t, st = proc.present in  
        match proc.past with
        | [] -> csv
        | _ -> add_line (backward proc) ((Sys.csv_line t st) :: csv)
      in let csv_data = [] in
      (Sys.csv_init ()) :: (add_line proc csv_data)

    let rec to_the_end tf (proc : process) =
      (* are we just changind the next future state, or lenghtening the series ? *)
      let t, st = proc.present in
      match t < tf with
      | false -> (* final time reached *)
        proc
      | true -> (* computations continue *)
        let (nnt, nnst), tl = 
          match proc.future with
          | [] -> 
            draw_future (t, st), [] 
          | (nt, nst) :: tl -> 
            failwith "Trying to draw future in the middle of a process"
        in (to_the_end tf) @@ forward @@
           {past = proc.past ;
            present = proc.present ;
            future = (nnt, nnst) :: tl}

    let simulate st_chan tf x0 =
      let chan = Csv.to_channel st_chan in
      Csv.output_record chan (Sys.csv_init ());
      let prst = (0., x0) in
      let proc0 = {past = []; present = prst; future = []} in
      let full_proc = to_the_end tf proc0 in
      let csv_data = csv_of_proc full_proc in
      Csv.output_all chan csv_data ;
      full_proc.present


  end;;

