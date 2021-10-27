### ARRAY DYNAMICS SIMULATION ###
# input:        Julia script for array dynamics simulation
# description:  minimal framework for the execution of the Julia simulation script 
# author:       HR

library(doParallel)
library(JuliaCall)

# ----- initialisation -----

julia <- julia_setup()

# ----- execution -----

tm = proc.time()
julia_source("src/0_DynamicMC.jl")
elap = proc.time() - tm


# ----- output -----

write(paste0("Julia: ",
             round(elap[3], 4), " - r0 ", r_0, ", J ", paste(J, collapse = '-') , ", c ", c, ", met ", met),
      file = paste0(outfol, 'time'), append = T)

write(paste0('simulation finished sucessfully at ', Sys.time()),
      file = paste0('logs/simulation_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'.txt'))

