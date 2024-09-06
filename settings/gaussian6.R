source("../run_simulation.R")

run_simulation(simulator=generate_hom_gaussian_3,"../simulations/gaussian6.RDS", n=1000,d=800,n_test=5000,B=100)
