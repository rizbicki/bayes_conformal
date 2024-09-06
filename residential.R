source("libraries.R")
source("distributions.R")
source("methods.R")
library(readxl)
library(readr)

#dados <- readxl::read_xlsx("real_data/Residential-Building-Data-Set.xlsx", range = cell_rows(1:2))
dados <- readxl::read_xlsx("real_data/ResidentialClean.xlsx")


ids <- sample(1:nrow(dados))
nTrain <- round(length(ids)*0.7)

dados_train <- dados[ids[1:nTrain],]
dados_test <- dados[-ids[1:nTrain],]

set.seed(0)
ids <- sample(c("train","val"),size = nrow(dados_train),prob = c(0.7,0.3),replace = TRUE)
xTrain <- dados_train[ids=="train",] %>% select(-c("V-9","V-10")) %>% as.matrix()
yTrain <- dados_train$`V-10`[ids=="train",drop=FALSE]
xVal <- dados_train[ids=="val",] %>% select(-c("V-9","V-10"))%>% as.matrix()
yVal <- dados_train$`V-10`[ids=="val",drop=FALSE]
xTrain_full <- dados_train %>% select(-c("V-9","V-10")) %>% as.matrix()
yTrain_full <- dados_train$`V-10`[drop=FALSE]

xTest <- dados_test %>% select(-c("V-9","V-10"))%>% as.matrix()
yTest <- dados_test$`V-10`

# Elastic Net
fitted_elastic_initial <- cv.glmnet(x=xTrain,y=yTrain)
l <- sort(c(fitted_elastic_initial$lambda,exp(seq(log(1e-13),log(max(lambda_path(xTrain,yTrain))^2),length.out=5000))))
fitted_elastic <- cv.glmnet(x=xTrain,y=yTrain,lambda = l)
fitted_elastic_conformal <- conformal(fitted_elastic,xVal,yVal)
pred_elastic_conformal <- predict(fitted_elastic_conformal,
                                  xTest)

# Ridge 
fitted_ridge_initial <- cv.glmnet(x=xTrain,y=yTrain,alpha = 0)
l <- sort(c(fitted_ridge_initial$lambda,exp(seq(log(1e-5),log(max(lambda_path(xTrain,yTrain))^2),length.out=5000))))
fitted_ridge <- cv.glmnet(x=xTrain,y=yTrain,alpha = 0,lambda = l)
fitted_ridge_conformal <- conformal(fitted_ridge,xVal,yVal)
pred_ridge_conformal <- predict(fitted_ridge_conformal,
                                xTest)


# LM Conformal
fitted_lm <- glmnet(x=xTrain,y=yTrain,alpha = 0,lambda = 0)
fitted_lm_conformal <- conformal(fitted_lm,xVal,yVal)
pred_lm_conformal <- predict(fitted_lm_conformal,
                             xTest)

# LM Exact predictive
data <- data.frame(y=yTrain_full,x=xTrain_full)
fitted_lm_exact <- lm(y~.,data = data)
data_full_test <- data.frame(y=yTest,x=xTest)
pred_lm_exact <- predict(fitted_lm_exact,data_full_test,
                         interval = "prediction",level=0.9) %>% 
  as.data.frame()
colnames(pred_lm_exact) <- c("pred","lower","upper")

best_l_se <- fitted_ridge_initial$lambda.1se
best_l <- fitted_ridge_initial$lambda.min

# Bayes
gamma_0_grid=sort(c(best_l_se,best_l,seq(0.1,100,length.out=20),2*round(0.9*nrow(xTrain_full))*best_l,
                    2*round(0.9*nrow(xTrain_full))*l[seq(1,length(l),length.out=40)]))
gamma_0_grid <- gamma_0_grid[-c(1:7)]
fitted_bayes <- bayesian_regression_cv(xTrain_full,yTrain_full,
                                       gamma_0_grid=gamma_0_grid)
pred_bayes <- predict(fitted_bayes,xTest)

fitted_bayes_conformal <- bayesian_regression_cv(xTrain,yTrain,gamma_0_grid=gamma_0_grid)
fitted_bayes_conformal <- conformal(fitted_bayes_conformal,xVal,yVal)
pred_bayes_conformal <- predict(fitted_bayes_conformal,xTest)


results <- rbind(evaluate2(pred_ridge_conformal,yTest,"Ridge (conformal)"),
      evaluate2(pred_lm_conformal,yTest,"lm (conformal)"),
      evaluate2(pred_lm_exact,yTest,"lm (exact)"),
      evaluate2(pred_bayes,yTest,"Bayes"),
      evaluate2(pred_bayes_conformal,yTest,"Bayes (conformal)"))
#,
 #     evaluate2(pred_elastic_conformal,yTest,"Elastic (conformal)"))

ggplot(results,aes(x=method, y=value)) +
  geom_point( fill="blue", pch=21, size=3)+
  geom_errorbar(aes(ymin=value-2*se, ymax=value+2*se), width=0.03, colour="blue", alpha=0.7)+
  facet_wrap(vars(metric), scales = "free")+
  theme_minimal(base_size = 16)+
  xlab("Method")+ylab("Statistic")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.spacing = unit(2, "lines"))
ggsave("figures/residential.png",width = 10,height = 6,bg = "white")
