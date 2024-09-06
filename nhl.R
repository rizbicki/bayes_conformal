source("libraries.R")
source("distributions.R")
source("methods.R")

dados_train <- read_csv("real_data/nhl/train.csv")
dados_testX <- read_csv("real_data/nhl/archive/test.csv")
dados_testY <- read_csv("real_data/nhl/archive/test_salaries.csv")
dados_test <- bind_cols(dados_testY,dados_testX)

dados_train$Hand <- recode(dados_train$Hand,L=0,R=1)
dados_test$Hand <- recode(dados_test$Hand,L=0,R=1)
dados_train$Salary <- log(dados_train$Salary)
dados_test$Salary <- log(dados_test$Salary)

cov_to_remove <- c("Born", "City", "Pr/St", "Cntry", "Nat","Last Name", "First Name", "Position",
                   "Team")
dados_train <- dados_train %>% select(-one_of(cov_to_remove)) %>% 
  filter(complete.cases(.))
dados_test <- dados_test %>% select(-one_of(cov_to_remove)) %>% 
  filter(complete.cases(.))

set.seed(0)
ids <- sample(c("train","val"),size = nrow(dados_train),prob = c(0.7,0.3),replace = TRUE)
xTrain <- dados_train[ids=="train",] %>% select(-Salary) %>% as.matrix()
yTrain <- dados_train$Salary[ids=="train",drop=FALSE]
xVal <- dados_train[ids=="val",] %>% select(-Salary)%>% as.matrix()
yVal <- dados_train$Salary[ids=="val",drop=FALSE]
xTrain_full <- dados_train %>% select(-Salary) %>% as.matrix()
yTrain_full <- dados_train$Salary[drop=FALSE]

xTest <- dados_test %>% select(-Salary)%>% as.matrix()
yTest <- dados_test$Salary

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
fitted_bayes <- bayesian_regression_cv(xTrain_full,yTrain_full,
                                       gamma_0_grid=gamma_0_grid)
pred_bayes <- predict(fitted_bayes,xTest)

fitted_bayes_conformal <- bayesian_regression_cv(xTrain,yTrain,gamma_0_grid=gamma_0_grid)
fitted_bayes_conformal <- conformal(fitted_bayes_conformal,xVal,yVal)
pred_bayes_conformal <- predict(fitted_bayes_conformal,xTest)

results <- as.data.frame(matrix(NA,nrow=6,ncol=7))
results[1,] <- cbind(method="Ridge (conformal)",
                                      evaluate(pred_ridge_conformal,yTest))
results[2,] <- cbind(method="lm (conformal)",
                                      evaluate(pred_lm_conformal,yTest))
results[3,] <- cbind(method="lm (exact)",
                                      evaluate(pred_lm_exact,yTest))
results[4,] <- cbind(method="Bayes",
                                      evaluate(pred_bayes,yTest))
results[5,] <- cbind(method="Bayes (conformal)",
                                      evaluate(pred_bayes_conformal,yTest))
results[6,] <- cbind(method="Elastic (conformal)",
                                      evaluate(pred_elastic_conformal,yTest))
colnames(results) <- c("Method",names(evaluate(pred_elastic_conformal,yTest)))

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
ggsave("figures/nhl.png",width = 10,height = 6,bg = "white")
