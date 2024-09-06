conformal <- function(fit,xCal,yCal,alpha=0.1)
{
  predictions <- predict(fit,xCal)
  if("pred"%in%names(predictions))
  {
    predictions=predictions$pred 
  }
  residuals <- abs(predictions-yCal)
  output <- list(cutoff=quantile(residuals,probs = 1-alpha),fit=fit)
  class(output) <- "conformal"
  return(output)
}

predict.conformal <- function(fit,xNew)
{
  pred <- predict(fit$fit,xNew)
  if("pred"%in%names(pred))
  {
    pred=pred$pred 
  }
  lower <- pred-fit$cutoff
  upper <- pred+fit$cutoff
  return(list(pred=pred,
              lower=lower,
              upper=upper)
  )
}

lambda_path <- function(x,y)
{
  mysd <- function(z) sqrt(sum((z-mean(z))^2)/length(z))
  sx <- scale(x, scale = apply(x, 2, mysd))
  sx <- as.matrix(sx, ncol = 20, nrow = 100)
  
  ## Calculate lambda path (first get lambda_max):
  lambda_max <- max(abs(colSums(sx*y)))/length(y)
  epsilon <- .0001
  K <- 100
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                              length.out = K)), digits = 10)
  return(lambdapath)
}

fit_bayes <- function(xTrain,yTrain,a,b,gamma_0)
{
  mu_0 <- rep(0,ncol(xTrain))
  reescale_factor <- mean(yTrain)
  yTrain <- yTrain-reescale_factor
  Gamma_0 <- diag(gamma_0,ncol(xTrain))
  a_star <- a+nrow(xTrain)/2
  Gamma_star <- t(xTrain)%*%xTrain+Gamma_0
  mu_star <- solve(Gamma_star)%*%(Gamma_0%*%mu_0+t(xTrain)%*%yTrain)
  b_star <- b+(t(yTrain)%*%yTrain + t(mu_0)%*% Gamma_0%*%mu_0- t(mu_star)%*%Gamma_star%*%mu_star )/2
  output <- list(a_star=a_star,
                 Gamma_star=Gamma_star,
                 mu_star=mu_star,
                 b_star=b_star,
                 reescale_factor=reescale_factor)
  class(output) <- "fitbayes"
  return(output)
}


bayesian_regression <- function(xTrain,yTrain,xCal,yCal,
                                a_grid=seq(0.00001,5,length.out=20),
                                b_grid=seq(0.00001,5,length.out=20),
                                gamma_0_grid=seq(0.1,100,length.out=20))
{
  grid <- expand.grid(a=a_grid,b=b_grid,gamma_0=gamma_0_grid)
  mse <- rep(NA,nrow(grid))
  for(gg in 1:nrow(grid))
  {
    fit <- fit_bayes(xTrain,yTrain,grid$a[gg],grid$b[gg],grid$gamma_0[gg])
    prediction_calibration <- predict(fit,xCal)
    
    mse[gg] <- mean((prediction_calibration$pred-yCal)^2)
  }
  best_par <- grid[which.min(mse)[1],]
  fit <- fit_bayes(xTrain,yTrain,best_par$a,best_par$b,best_par$gamma_0)
  
  output <- list(a_star=fit$a_star,
                 Gamma_star=fit$Gamma_star,
                 mu_star=fit$mu_star,
                 b_star=fit$b_star,
                 reescale_factor=fit$reescale_factor,
                 a_0=best_par$a,
                 b_0=best_par$b,
                 gamma_0=best_par$gamma_0,
                 mse=mse)
  class(output) <- "fitbayes"
  return(output)
}


bayesian_regression_cv <- function(xTrain,yTrain,
                                a_0=0.00001,
                                b_0=1,nFolders=10,gamma_0_grid=seq(0.001,100,length.out=100))
{
  mse <- rep(NA,length(gamma_0_grid))
  cv <- vfold_cv(xTrain, v = nFolders)
  predictions <- rep(NA,length(yTrain))
  for(gg in 1:length(mse))
  {
    for(folder in 1:nFolders)
    {
      fit <- fit_bayes(xTrain[cv$splits[[folder]]$in_id,],yTrain[cv$splits[[folder]]$in_id],
                       a_0,b_0,gamma_0_grid[gg])
      predictions[-cv$splits[[folder]]$in_id] <- predict(fit,xTrain[-cv$splits[[folder]]$in_id,])$pred
    }
    mse[gg] <- mean((predictions-yTrain)^2)  
  }
  gamma_0 <- gamma_0_grid[which.min(mse)[1]]
  fit <- fit_bayes(xTrain,yTrain,a_0,b_0,gamma_0)
  
  output <- list(a_star=fit$a_star,
                 Gamma_star=fit$Gamma_star,
                 mu_star=fit$mu_star,
                 b_star=fit$b_star,
                 reescale_factor=fit$reescale_factor,
                 a_0=a_0,
                 b_0=b_0,
                 gamma_0=gamma_0,
                 mse=mse)
  class(output) <- "fitbayes"
  return(output)
}


predict.fitbayes <- function(fit,xNew,alpha=0.1)
{
  pred=fit$reescale_factor+xNew%*%fit$mu_star
  cov=as.numeric(fit$b_star/fit$a_star)*
    (diag(1,nrow(xNew))+xNew%*%solve(fit$Gamma_star)%*%t(xNew))
  lower <- qlst(alpha/2, df=2*fit$a_star, 
                mu = as.numeric(pred), 
                sigma = sqrt(diag(cov)))
  upper <- qlst(1-alpha/2, df=2*fit$a_star, 
                mu = as.numeric(pred), 
                sigma = sqrt(diag(cov)))
  return(list(pred=pred,
              cov=cov,
              lower=lower,
              upper=upper)
  )
}

evaluate <- function(pred,yTest)
{
  mse <- mean((pred$pred-yTest)^2)
  mse_se <- sd((pred$pred-yTest)^2)/sqrt(length(pred$upper))
  size <- mean(pred$upper-pred$lower)
  size_se <- sd(pred$upper-pred$lower)/sqrt(length(pred$upper))
  coverage <- mean((pred$lower <= yTest)&(yTest <= pred$upper))
  coverage_se <- sd((pred$lower <= yTest)&(yTest <= pred$upper))/sqrt(length(pred$upper))
  return(data.frame(mse=mse,mse_se=mse_se,
                    size=size,size_se=size_se,
                    coverage=coverage,coverage_se=coverage_se))
}

evaluate2 <- function(pred,yTest,method_name)
{
  mse <- mean((pred$pred-yTest)^2)
  mse_se <- sd((pred$pred-yTest)^2)/sqrt(length(pred$upper))
  size <- mean(pred$upper-pred$lower)
  size_se <- sd(pred$upper-pred$lower)/sqrt(length(pred$upper))
  coverage <- mean((pred$lower <= yTest)&(yTest <= pred$upper))
  coverage_se <- sd((pred$lower <= yTest)&(yTest <= pred$upper))/sqrt(length(pred$upper))
  d <- data.frame(cbind(method_name,rbind(c("MSE",mse,mse_se),
                        c("Size",size,size_se),
                        c("Coverage",coverage,coverage_se))))
  colnames(d) <- c("method","metric","value","se")
  d$value <- as.numeric(d$value)
  d$se <- as.numeric(d$se)
  return(d)
}

results_tidy <- function(results)
{
  results <- results[complete.cases(results),]
  colnames(results) <- c("it","method",
                         "mse","mse_se",
                         "size","size_se",
                         "coverage","coverage_se")
  results <-  results %>% 
    rename(MSE = mse,Coverage=coverage,Size=size)
  
  results_tidy <- results %>% 
    pivot_longer(
      cols = c("MSE","Size","Coverage"), 
      names_to = "measurement", 
      values_to = "value"
    )
  results_tidy$measurement <- factor(results_tidy$measurement,
                                     levels = c("Coverage","Size","MSE"))
  
  return(results_tidy)
}

results_tidy_noit <- function(results)
{
  results <- results[complete.cases(results),]
  colnames(results) <- c("method",
                         "mse","mse_se",
                         "size","size_se",
                         "coverage","coverage_se")
  results <-  results %>% 
    rename(MSE = mse,Coverage=coverage,Size=size)
  
  results_tidy <- results %>% 
    pivot_longer(
      cols = c("MSE","Size","Coverage"), 
      names_to = "measurement", 
      values_to = "value"
    )
  results_tidy$measurement <- factor(results_tidy$measurement,
                                     levels = c("Coverage","Size","MSE"))
  
  return(results_tidy)
}


create_plot <- function(results)
{
  results <- results_tidy(results)
  return(ggplot(results,aes(x=method, y=value)) +
           stat_summary(fun.data=function(x){mean_cl_normal(x, conf.int=.95)}, 
                        geom="errorbar", 
                        width=0.03, colour="blue", alpha=0.7) +
           stat_summary(fun.y=mean, geom="point", fill="blue", pch=21, size=3) +
           facet_wrap(vars(measurement), scales = "free")+
           theme_minimal(base_size = 16)+
           xlab("Method")+ylab("Statistic")+
           theme(axis.text.x = element_text(angle = 45, hjust=1),
                 panel.spacing = unit(2, "lines")))
}


create_plot_noit <- function(results)
{
  results <- results_tidy_noit(results)
  return(ggplot(results,aes(x=method, y=value)) +
           stat_summary(fun.data=function(x){mean_cl_normal(x, conf.int=.95)}, 
                        geom="errorbar", 
                        width=0.03, colour="blue", alpha=0.7) +
           stat_summary(fun.y=mean, geom="point", fill="blue", pch=21, size=3) +
           facet_wrap(vars(measurement), scales = "free")+
           theme_minimal(base_size = 16)+
           xlab("Method")+ylab("Statistic")+
           theme(axis.text.x = element_text(angle = 45, hjust=1),
                 panel.spacing = unit(2, "lines")))
}