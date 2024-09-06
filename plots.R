source("libraries.R")
source("distributions.R")
source("methods.R")

results <- readRDS("simulations/gaussian1.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian1.png",width = 10,height = 6,bg = "white")

results <- readRDS("simulations/gaussian2.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian2.png",width = 10,height = 6,bg = "white")

results <- readRDS("simulations/gaussian3.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian3.png",width = 10,height = 6,bg = "white")

results <- readRDS("simulations/gaussian4.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian4.png",width = 10,height = 6,bg = "white")


results <- readRDS("simulations/gaussian5.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian5.png",width = 10,height = 6,bg = "white")

results <- readRDS("simulations/gaussian6.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian6.png",width = 10,height = 6,bg = "white")

results <- readRDS("simulations/gaussian7.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/gaussian7.png",width = 10,height = 6,bg = "white")



results <- readRDS("simulations/het.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/het.png",width = 10,height = 6,bg = "white")

results <- readRDS("simulations/het2.RDS")
results <- results %>% filter(X2!="Elastic (conformal)")
create_plot(results)
ggsave("figures/het2.png",width = 10,height = 6,bg = "white")
