generate_bimodal <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-1.5,1.5),n,d)
  }
  f=(x[,1]-1)^2*(x[,1]+1)
  g=rep(0,n)
  g[x[,1]> -0.5]=2*sqrt(x[x[,1]> -0.5,1]+0.5)
  s=1/4+abs(x[,1])
  # response
  y=ifelse(runif(n)>0.5,f-g,f+g)+rnorm(n,0,sqrt(s))
  return(list(x=x,y=y))
}

generate_het_gaussian <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-5,5),n,d)
  }
  # response
  y=rnorm(nrow(x),x[,1],1+abs(x[,1]))
  return(list(x=x,y=y))
}

generate_hom_gaussian <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-5,5),n,d)
  }
  # response
  y=x[,1]+rnorm(nrow(x),0,1)
  return(list(x=x,y=y))
}


generate_hom_gaussian_2 <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-5,5),n,d)
  }
  # response
  y=4*x[,1]+2*x[,2]+rnorm(nrow(x),0,1)
  return(list(x=x,y=y))
}

generate_hom_gaussian_3 <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-5,5),n,d)
  }
  # response
  y=rowSums(0.5*x[,1:20])+rnorm(nrow(x),0,1)
  return(list(x=x,y=y))
}

generate_gamma <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-5,5),n,d)
  }
  # response
  y=5*x[,1]+rgamma(nrow(x),1+2*abs(x[,1]),1+2*abs(x[,1]))
  return(list(x=x,y=y))
}


