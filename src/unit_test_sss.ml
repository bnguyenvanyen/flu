(* A small executable to test and debug *)
#require "flu";;

let tf = 365. *. 100.;;
let s0 = int_of_float (5. *. (10. ** 5.));;
let i0 = int_of_float (10. ** 3.);;
let r0 = int_of_float (4.99 *. (10. ** 5.));;
let x0 = (s0, i0, r0);;
let fname = "/home/queue/Documents/2015stage/code/sss_106.csv";;

Sss.Gen.simulate fname tf x0;;

