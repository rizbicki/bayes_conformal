run_simulation <- function(simulator,file_name, n,d,n_test=5000,B=100) {
  source("../libraries.R")
  source("../distributions.R")
  source("../methods.R")
  
  n_methods <- 6
  if(file.exists(file_name))
  {
    results <- readRDS(file_name)
  } else {
    results <- data.frame(matrix(NA,nrow = B*n_methods,ncol=8))  
  }
  for(ii in 1:B)
  {
    print(ii/B)
    if(!is.na(results[n_methods*(ii-1)+1,1]))
      next;
    data_test <- simulator(n_test,d)
    data_train_and_val <- simulator(n,d)
    ids <- sample(c("train","val"),size = n,prob = c(0.7,0.3),replace = TRUE)
    xTrain <- data_train_and_val$x[ids=="train",]
    yTrain <- data_train_and_val$y[ids=="train",drop=FALSE]
    xVal <- data_train_and_val$x[ids=="val",]
    yVal <- data_train_and_val$y[ids=="val",drop=FALSE]
    
    
    # Elastic Net
    fitted_elastic_initial <- cv.glmnet(x=xTrain,y=yTrain)
    l <- sort(c(fitted_elastic_initial$lambda,exp(seq(log(1e-13),log(max(lambda_path(xTrain,yTrain))^2),length.out=5000))))
    fitted_elastic <- cv.glmnet(x=xTrain,y=yTrain,lambda = l)
    fitted_elastic_conformal <- conformal(fitted_elastic,xVal,yVal)
    pred_elastic_conformal <- predict(fitted_elastic_conformal,
                                      data_test$x)
    
    # Ridge 
    fitted_ridge_initial <- cv.glmnet(x=xTrain,y=yTrain,alpha = 0)
    l <- sort(c(fitted_ridge_initial$lambda,exp(seq(log(1e-5),log(max(lambda_path(xTrain,yTrain))^2),length.out=5000))))
    fitted_ridge <- cv.glmnet(x=xTrain,y=yTrain,alpha = 0,lambda = l)
    fitted_ridge_conformal <- conformal(fitted_ridge,xVal,yVal)
    pred_ridge_conformal <- predict(fitted_ridge_conformal,
                                    data_test$x)
    
    # LM Conformal
    fitted_lm <- glmnet(x=xTrain,y=yTrain,alpha = 0,lambda = 0)
    fitted_lm_conformal <- conformal(fitted_lm,xVal,yVal)
    pred_lm_conformal <- predict(fitted_lm_conformal,
                                 data_test$x)
    
    # LM Exact predictive
    data <- data.frame(y=data_train_and_val$y,x=data_train_and_val$x)
    fitted_lm_exact <- lm(y~.,data = data)
    data_full_test <- data.frame(y=data_test$y,x=data_test$x)
    pred_lm_exact <- predict(fitted_lm_exact,data_full_test,
                             interval = "prediction",level=0.9) %>% 
      as.data.frame()
    colnames(pred_lm_exact) <- c("pred","lower","upper")
    
    best_l_se <- fitted_ridge_initial$lambda.1se
    best_l <- fitted_ridge_initial$lambda.min
    
    # Bayes
    gamma_0_grid=sort(c(best_l_se,best_l,seq(0.1,100,length.out=20),2*round(0.9*nrow(data_train_and_val$x))*best_l,
                        2*round(0.9*nrow(data_train_and_val$x))*l[seq(1,length(l),length.out=40)]))
    fitted_bayes <- bayesian_regression_cv(data_train_and_val$x,data_train_and_val$y,
                                           gamma_0_grid=gamma_0_grid)
    pred_bayes <- predict(fitted_bayes,data_test$x)
    
    fitted_bayes_conformal <- bayesian_regression_cv(xTrain,yTrain,gamma_0_grid=gamma_0_grid)
    fitted_bayes_conformal <- conformal(fitted_bayes_conformal,xVal,yVal)
    pred_bayes_conformal <- predict(fitted_bayes_conformal,data_test$x)
    
    results[n_methods*(ii-1)+1,] <- cbind(it=ii,method="Ridge (conformal)",
                                          evaluate(pred_ridge_conformal,data_test$y))
    results[n_methods*(ii-1)+2,] <- cbind(it=ii,method="lm (conformal)",
                                          evaluate(pred_lm_conformal,data_test$y))
    results[n_methods*(ii-1)+3,] <- cbind(it=ii,method="lm (exact)",
                                          evaluate(pred_lm_exact,data_test$y))
    results[n_methods*(ii-1)+4,] <- cbind(it=ii,method="Bayes",
                                          evaluate(pred_bayes,data_test$y))
    results[n_methods*(ii-1)+5,] <- cbind(it=ii,method="Bayes (conformal)",
                                          evaluate(pred_bayes_conformal,data_test$y))
    results[n_methods*(ii-1)+6,] <- cbind(it=ii,method="Elastic (conformal)",
                                          evaluate(pred_elastic_conformal,data_test$y))
    gc()
    saveRDS(results,file_name)
  }
}