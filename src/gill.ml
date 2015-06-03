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
    val min_step : float
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
          string ->
          float ->
          state ->
          process
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
      let elapsed_t = rand_exp lambd in
      let next_t = t +. elapsed_t in
      let w = rand_float lambd in
      let m = choose_modif w srl Sys.ml in
      let next_st = m st in
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
      let rec action t nnt =
        match nnt -. t < Algp.min_step with
        | false -> fun p -> (to_the_end tf) @@ forward @@ p
        | true -> to_the_end tf
      in
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
            draw_future (nt, nst), tl
        in (action t nnt) @@
           {past = proc.past ;
            present = proc.present ;
            future = (nnt, nnst) :: tl}

    let simulate fname tf x0 =
      let prst = (0., x0) in
      let proc0 = {past = []; present = prst; future = []} in
      let full_proc = to_the_end tf proc0 in
      let csv_data = csv_of_proc full_proc in
      Csv.save fname csv_data ;
      full_proc

          

  end;;

