let pi = 4. *. atan 1.;;


let a_incr a i =
  a.(i) <- a.(i) + 1

let a_decr a i =
  a.(i) <- a.(i) - 1

let foi = float_of_int

module type PARS =
  sig
    val n : float
    val r0 : float
    val e : float
    val etaN1 : float
    val etaN2 : float
    val g1 : float
    val g2 : float
    val nu : float
    val q : float
  end

module Sys (Pars : PARS) : (Gill.SYSTEM with type state = int array
                                         and type aux = float * float) =
  struct
    open Pars
    type state = int array;;
    type aux = float * float;;
    let eta1 = etaN1 *. n
    let eta2 = etaN2 *. n
    let bet0 = r0 *. nu
    
    let cnt_incr = ref 0
    let cnt_decr = ref 0

    let beta t =
      bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.))

    let gi j = (* get i *)
      match j with
      | 10 -> 4
      | 20 -> 5
      | 12 -> 6
      | 21 -> 7
      | _ -> invalid_arg "expects 10, 20, 12 or 21"
    
    let gr j =
      match j with
      | 0 -> 0
      | 1 -> 1
      | 2 -> 2
      | 12 -> 3
      | _ -> invalid_arg "expects 0, 1, 2 or 12"

    let gq j =
      match j with
      | 0 -> 8
      | 1 -> 9
      | 2 -> 10
      | 12 -> 11
      | _ -> invalid_arg "expects 0, 1, 2, or 12"

    let fi snk st =
      match snk with
      | 10 -> foi (st.(gi 10) + st.(gi 12)) +. eta1
      | 20 -> foi (st.(gi 20) + st.(gi 21)) +. eta2
      | 12 -> foi (st.(gi 10) + st.(gi 12)) +. eta1
      | 21 -> foi (st.(gi 20) + st.(gi 21)) +. eta2
      | _ -> invalid_arg "expects 10, 20, 12 or 21" 

    let imloss1_rate src snk t st =
      ignore (gr snk) ;
      g1 *. foi st.(gr src)

    let imloss2_rate src snk t st =
      ignore (gr snk) ;
      g2 *. foi st.(gr src)

    let imloss_modif src snk st =
      let nst = Array.copy st in
      a_decr nst (gr src) ;
      a_incr nst (gr snk) ;
      nst

    let imlossi1_rate src snk t st =
      g1 *. foi st.(gi src)

    let imlossi2_rate src snk t st =
      g2 *. foi st.(gi src)

    let imlossi_modif src snk st =
      let nst = Array.copy st in
      a_decr nst (gi src) ;
      a_incr nst (gi snk) ;
      nst

 
    let qloss_rate src snk t st =
      ignore (gr snk) ;
      q *. foi st.(gq src)
  
    let qloss_modif src snk st =
      let nst = Array.copy st in
      a_decr nst (gq src) ;
      a_incr nst (gr snk) ;
      nst

    let inf_rate src snk t st =
      beta t *. foi st.(gr src) /. n *. (fi snk st)

    let inf_modif src snk st =
      let nst = Array.copy st in
      a_decr nst (gr src) ;
      a_incr nst (gi snk) ;
      nst

    let recov_rate src snk t st =
      ignore (gq snk) ; (* in case bad sink *)
      nu *. foi st.(gi src)

    let recov_modif src snk st =
      let nst = Array.copy st in
      a_decr nst (gi src) ;
      a_incr nst (gq snk) ;
      nst
 
    let min_step = 1.
    let aux_fun t st =
      let infct1 = (inf_rate 0 10 t st) +. (inf_rate 2 12 t st)
      and infct2 = (inf_rate 0 20 t st) +. (inf_rate 1 21 t st)
      in
      (infct1 *. 7. *. 100000. /. n, infct2 *. 7. *. 100000. /. n)

    let fl = [imloss1_rate 1 0; imloss2_rate 2 0; imloss1_rate 12 2; 
              imloss2_rate 12 1;
              imlossi1_rate 21 20 ; imlossi2_rate 12 10 ; 
              qloss_rate 0 0; qloss_rate 1 1 ; qloss_rate 2 2 ; 
              qloss_rate 12 12 ;
              inf_rate 0 10 ; inf_rate 0 20; inf_rate 2 12 ; inf_rate 1 21 ; 
              recov_rate 10 1 ; recov_rate 20 2 ; recov_rate 21 12 ; 
              recov_rate 12 12]
    let ml = [imloss_modif 1 0 ; imloss_modif 2 0 ; imloss_modif 12 2; 
              imloss_modif 12 1 ; 
              imlossi_modif 21 20 ; imlossi_modif 12 10 ; 
              qloss_modif 0 0 ; qloss_modif 1 1 ; qloss_modif 2 2 ; 
              qloss_modif 12 12 ;
              inf_modif 0 10 ; inf_modif 0 20 ; inf_modif 2 12 ; inf_modif 1 21 ; 
              recov_modif 10 1 ; recov_modif 20 2 ; recov_modif 21 12 ; 
              recov_modif 12 12]
    
    let csv_init () =
      [
        ["n=12" ; "m=2"] ;
        ["t" ; 
         "inc1" ; "inc2" ;
         "R0" ; "R1" ; "R2" ; "R12" ; 
         "I10" ; "I20" ; "I12" ; "I21" ; 
         "Q0" ; "Q1" ; "Q2" ; "Q12"]
      ]
    
    let csv_line t au st =
      let au1, au2 = au in
      (string_of_float t) :: 
      (string_of_float au1) :: (string_of_float au2) ::
      Array.to_list (Array.map string_of_int st) 
  end;;

module Default_Algp =
  struct
    let min_step = 1.
  end
