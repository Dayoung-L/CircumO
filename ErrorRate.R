kk = fn(100)
bb = list()

temp = which(PC_res1 == "BEEP")
A = 1
for (i in temp){
  bb[[A]] = kk[[i]] 
  A = A + 1
}

library(purrr)
possibly_some_function = possibly(CircumOLS,
                                  otherwise = "BEEP")
PCSB_res1 = map(bb, possibly_some_function)