source("../run_simulation.R")

run_simulation(generate_het_gaussian,"../simulations/het.RDS", n=1000,d=500,n_test=5000,B=100)
