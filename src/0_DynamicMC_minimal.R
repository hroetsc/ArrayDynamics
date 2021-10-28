### ARRAY DYNAMICS SIMULATION ###
# input:        Julia script for array dynamics simulation
# output:       simulation results
# description:  minimal framework for the execution of the Julia simulation script 
# author:       HR


library(doParallel)
library(JuliaCall)
library(rhdf5)
library(stringr)


# ----- initialisation -----
julia <- julia_setup()


# ----- execution -----
tm = proc.time()
julia_source("src/0_DynamicMC.jl")
elap = proc.time() - tm


# ----- load and parse results -----
source('src/argparse.R')

print("LOAD SIMULATION RESULTS AND PARSE THEM INTO RDATA")

for (i in 1:length(allresults)) {
  cnt = allresults[[i]]
  
  SIMresults = list(A = split(cnt$A, rep(1:nrow(cnt$A), each=ncol(cnt$A))),
                    M = split(cnt$M, rep(1:nrow(cnt$M), each=ncol(cnt$M))))
  names(SIMresults[["A"]]) = cnt$Ts
  names(SIMresults[["M"]]) = cnt$Ts
  
  outfile = paste0(outfol,"_c",c,"_met-",met,"_rep",i,".RData")
  save(SIMresults, file = outfile)
}


# ----- output -----
write(paste0("Julia: ",
             round(elap[3], 4), " - r0 ", r_0, ", J ", paste(J, collapse = '-') , ", c ", c, ", met ", met),
      file = paste0(outfol, 'time'), append = T)

write(paste0('simulation finished sucessfully at ', Sys.time()),
      file = paste0('logs/simulation_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'.txt'))

