let pi = 4. *. atan 1.;;

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

module Sys (Pars : PARS) : (Gill.SYSTEM with type state = int array) =
  struct
    open Pars
    type state = int array;;
    let eta1 = etaN1 *. n;;
    let eta2 = etaN2 *. n;;
    let bet0 = r0 *. nu;;
    
    let foi = float_of_int;;   
 
    let beta t =
      bet0 *. (1. +. e *. cos (2. *. pi *. t /. 365.));;

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

    let fi1 st =
      foi (st.(gi 10) + st.(gi 12))

    let fi2 st =
      foi (st.(gi 20) + st.(gi 21))
    
    let imloss10_rate t st =
      g1 *. foi st.(gr 1);;

    let imloss10_modif st =
      st.(gr 1) <- st.(gr 1) - 1;
      st.(gr 0) <- st.(gr 0) + 1;
      st
  
    let imloss20_rate t st =
      g2 *. foi st.(gr 2);;
    
    let imloss20_modif st =
      st.(gr 2) <- st.(gr 2) - 1;
      st.(gr 0) <- st.(gr 0) + 1;
      st
  
    let imloss122_rate t st =
      g1 *. foi st.(gr 12);;

    let imloss122_modif st =
      st.(gr 12) <- st.(gr 12) - 1;
      st.(gr 2) <- st.(gr 2) + 1;
      st
  
    let imloss121_rate t st =
      g2 *. foi st.(gr 12);;
    
    let imloss121_modif st =
      st.(gr 12) <- st.(gr 12) - 1;
      st.(gr 1) <- st.(gr 1) + 1;
      st
  
    let imlossi2120_rate t st =
      g1 *. foi st.(gi 21);;

    let imlossi2120_modif st =
      st.(gi 21) <- st.(gi 21) - 1;
      st.(gi 20) <- st.(gi 20) + 1;
      st
  
    let imlossi1210_rate t st =
      g2 *. foi st.(gi 12);;

    let imlossi1210_modif st =
      st.(gi 12) <- st.(gi 12) - 1;
      st.(gi 10) <- st.(gi 10) + 1;
      st
  
    let qloss0_rate t st =
      q *. foi (gq 0)

    let qloss0_modif st =
      st.(gq 0) <- st.(gq 0) - 1;
      st.(gr 0) <- st.(gr 0) + 1;
      st
  

    let qloss1_rate t st =
      q *. foi st.(gq 1)

    let qloss1_modif st =
      st.(gq 1) <- st.(gq 1) - 1;
      st.(gr 1) <- st.(gr 1) + 1;
      st
  
    let qloss2_rate t st =
      q *. foi st.(gq 2)

    let qloss2_modif st =
      st.(gq 2) <- st.(gq 2) - 1;
      st.(gr 2) <- st.(gr 2) + 1;
      st
  
    let qloss12_rate t st =
      q *. foi st.(gq 12)

    let qloss12_modif st =
      st.(gq 12) <- st.(gq 12) - 1;
      st.(gr 12) <- st.(gr 12) + 1;
      st

    let inf010_rate t st =
      beta t *. foi st.(gr 0) /. n *. (fi1 st +. eta1);;

    let inf010_modif st =
      st.(gr 0) <- st.(gr 0) - 1;
      st.(gi 10) <- st.(gi 10) + 1;
      st
       
    let inf020_rate t st =
      beta t *. foi st.(gr 0) /. n *. (fi2 st +. eta2);;

    let inf020_modif st =
      st.(gr 0) <- st.(gr 0) - 1;
      st.(gi 20) <- st.(gi 20) + 1;
      st
       
    let inf212_rate t st =
      beta t *. foi st.(gr 2) /. n *. (fi1 st +. eta1);;

    let inf212_modif st =
      st.(gr 2) <- st.(gr 2) - 1;
      st.(gi 12) <- st.(gi 12) + 1;
      st
       
    let inf121_rate t st =
      beta t *. foi st.(gr 1) /. n *. (fi2 st +. eta2);;

    let inf121_modif st =
      st.(gr 1) <- st.(gr 1) - 1;
      st.(gi 21) <- st.(gi 21) + 1;
      st
       
    let recov101_rate t st =
      nu *. foi st.(gi 10);;

    let recov101_modif st =
      st.(gi 10) <- st.(gi 10) - 1;
      st.(gq 1) <- st.(gq 1) + 1;
      st
    
    let recov202_rate t st =
      nu *. foi st.(gi 20);;

    let recov202_modif st =
      st.(gi 20) <- st.(gi 20) - 1;
      st.(gq 2) <- st.(gq 2) + 1;
      st
    
    let recov2112_rate t st =
      nu *. foi st.(gi 21);;

    let recov2112_modif st =
      st.(gi 21) <- st.(gi 21) - 1;
      st.(gq 12) <- st.(gq 12) + 1;
      st
    
    let recov1212_rate t st =
      nu *. foi st.(gi 12);;

    let recov1212_modif st =
      st.(gi 12) <- st.(gi 12) - 1;
      st.(gq 12) <- st.(gq 12) + 1;
      st
    
    let min_step = 1.
    let fl = [imloss10_rate ; imloss20_rate ; imloss122_rate ; imloss121_rate ;
              imlossi2120_rate ; imlossi1210_rate ; 
              qloss0_rate ; qloss1_rate ; qloss2_rate ; qloss12_rate ;
              inf010_rate ; inf020_rate ; inf212_rate ; inf121_rate ; 
              recov101_rate ; recov202_rate ; recov2112_rate ; recov1212_rate]
    let ml = [imloss10_modif ; imloss20_modif ; imloss122_modif ; 
              imloss121_modif ; imlossi2120_modif ; imlossi1210_modif ; 
              qloss0_modif ; qloss1_modif ; qloss2_modif ; qloss12_modif ;
              inf010_modif ; inf020_modif ; inf212_modif ; inf121_modif ; 
              recov101_modif ; recov202_modif ; recov2112_modif ; recov1212_modif]
    
    let csv_init () =
      ["t" ; 
       "R0" ; "R1" ; "R2" ; "R12" ; 
       "I10" ; "I20" ; "I12" ; "I21" ; 
       "Q0" ; "Q1" ; "Q2" ; "Q12"]
    
    let csv_line t st =
      (string_of_float t) :: 
      Array.to_list (Array.map string_of_int st) 
  end;;

module Default_Algp =
  struct
    let min_step = 1.
  end
