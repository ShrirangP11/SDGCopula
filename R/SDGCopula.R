# Prerequisite packages
.packages = c('MASS', 'mvtnorm', 'copula', 'fitdistrplus', 'plotrix', 'psych', 'ggcorrplot', 'gridExtra', 'reshape2', 'ggplot2', 'lsa')
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session 
# lapply(.packages, require, character.only=TRUE)
lapply(.packages, FUN = function(X) {
  do.call("require", list(X)) 
})


#Function to for parametric density estimation of marginals
get.marginal = function(x){
  if(length(x)>1){
    aic.vec = c()
    if (any(x < 0)){
      distributions <- c("norm", "unif","cauchy","logis")
      for(i in 1:length(distributions)){
        fit = fitdistrplus::fitdist(x,distributions[i],method = "mle")
        aic.vec[i] = fit$aic
      }
    }else{
      if(all(x >= 0 & x <= 1)){
        distributions <- c("beta","norm", "exp", "gamma", "weibull", "unif","cauchy","lnorm","logis")
        for(i in 1:length(distributions)){
          fit = fitdistrplus::fitdist(x,distributions[i],method = "mle")
          aic.vec[i] = fit$aic
        }}else{
          distributions <- c("norm", "exp", "gamma", "weibull", "unif","cauchy","lnorm","logis")
          for(i in 1:length(distributions)){
            fit = fitdistrplus::fitdist(x,distributions[i],method = "mle")
            aic.vec[i] = fit$aic
          }
        }
    }
    fit = fitdistrplus::fitdist(x,distributions[which.min(aic.vec)],method = "mle")
    return(list(distribution=distributions[which.min(aic.vec)],par = fit$estimate,AIC=fit$aic))
  }else(return("Error: Invalid Vector"))
}


#Function to return pairwise correlations between distinct variables
cor_struct <- function(mat) {
  # Check if the matrix is symmetric
  if (!is.matrix(mat) || !all(mat == t(mat))) {
    stop("Input must be a symmetric matrix.")
  }
  non_diag_entries <- c()
  for (i in 1:(nrow(mat) - 1)) {
    non_diag_entries <- c(non_diag_entries, mat[i, (i+1):ncol(mat)])
  }
  return(non_diag_entries)
}



#Define a function to obtain KDE, sample from it, then return quantiles.
qdist <- function(u, input, type){
  if(type==0){ # 0:numeric
    quantile <- c()
    sampled <- sample(density(input)$x,prob=density(input)$y, 99999, replace=T)
    my_ecdf <- ecdf(input)
    dataframe <- data.frame(sort(sampled),sort(my_ecdf(sampled)))
    for(i in u){
      quantile <- append(quantile, dataframe[min(which(dataframe[,2]>=i)),1])
    }
    return(quantile)
  }else if(type==1){ # 1:categorical
    quantile <- c()
    prob <- c()
    support <- c(1:length(input))
    for(i in support){
      prob <- append(prob, sum(ifelse(input==i,1,0))/length(input))
    }
    cdf <- cumsum(prob)
    dataframe <- data.frame(support, cdf)
    for(j in u){
      quantile <- append(quantile, dataframe[min(which(dataframe[,2]>=j)),1])
    }
    return(quantile)
  }else if(type==2){ # 2:integer
    quantile <- c()
    prob <- c() 
    support <- sort(unique(input)) 
    for(i in support){
      prob <- append(prob, sum(ifelse(input==i,1,0))/length(input))
    }
    cdf <- cumsum(prob)
    dataframe <- data.frame(support, cdf)
    for(j in u){
      quantile <- append(quantile, dataframe[min(which(dataframe[,2]>=j)),1])
    }
    return(quantile)
  }
}



#Function to get parametric estimate of marginals
marginals <- function(dataframe){
  dist <- c()
  param <- list()
  for(i in 1:ncol(dataframe)){
    dist[i] <- get.marginal(dataframe[,i])$distribution
    param[i] <- list(as.list(get.marginal(dataframe[,i])$par))
  }
  return(list(dist=dist,par=param))
}


#Function to convert Categorical variables to numeric.
tonumeric <- function(x){
  isfactor <- !sapply(x, is.numeric)
  codes <- lapply(x[isfactor], \(x)levels(factor(x)))
  x[isfactor] <- data.matrix(x[isfactor])
  structure(x, codes = codes)
}

#Function to convert numeric back to original categorical variables.
tocategory <- function(x){
  codes <- attr(x, 'codes')
  x[names(codes)] <- Map(\(i,j)j[i], x[names(codes)], codes)
  x
} 


#Function to fit Normal/t-copula for dim>=2
fitCop <- function(dataframe, copula, parametric, dof){
  
  dataframe1 <- tonumeric(dataframe)
  
  if(missing(dof)){ #dof=4 by default
    dof=4
  }
  if(missing(parametric)){ #parametric=FALSE by default
    parametric=FALSE
  }
  if(missing(copula)){ #copula='normal' by default
    copula='normal'
  }
  #For normal-Copula
  if(copula=='normal'){
    if(parametric==T){
      myCop <- copula::normalCopula(param=cor_struct(cor(dataframe1, method=c('spearman'))), dim = ncol(dataframe1), dispstr = "un")
      myMvd <- copula::mvdc(myCop, margins=marginals(dataframe1)$dist, paramMargins = marginals(dataframe1)$par)
      Z <- copula::rMvdc(nrow(dataframe1), myMvd)
      colnames(Z) <- colnames(dataframe1)
      return(tocategory(Z))
    }else if(parametric==F){
      u <- data.frame(matrix(ncol=ncol(dataframe),nrow=nrow(dataframe)))
      sample <- MASS::mvrnorm(nrow(dataframe1),mu=rep(0,ncol(dataframe1)),Sigma=cor(dataframe1,method=c('spearman')),empirical=T)
      for(i in 1:ncol(u)){
        u[,i] <- pnorm(sample[,i])
      }
      for(j in 1:ncol(dataframe1)){
        if(class(dataframe[,j])=='numeric'){
          dataframe1[,j] <- qdist(u[,j],dataframe1[,j],type = 0)
        }else if(class(dataframe[,j])=='factor'){
          dataframe1[,j] <- qdist(u[,j],dataframe1[,j],type = 1)
        }else if(class(dataframe[,j])=='integer'){
          dataframe1[,j] <- qdist(u[,j],dataframe1[,j],type = 2)
        }
      }
      dataframe1 <- tocategory(dataframe1)
      dataframe1[sapply(dataframe1, is.character)] <- lapply(dataframe1[sapply(dataframe1, is.character)], as.factor)
      return(dataframe1)
    }
    
  #For t-Copula
  }else if(copula=='t'){ 
    if(parametric==T){
      myCop <- copula::tCopula(df=dof, param=cor_struct(cor(dataframe1, method=c('spearman'))), dim = ncol(dataframe1), dispstr = "un")
      myMvd <- copula::mvdc(myCop, margins=marginals(dataframe1)$dist, paramMargins = marginals(dataframe1)$par)
      Z <- copula::rMvdc(nrow(dataframe1), myMvd)
      colnames(Z) <- colnames(dataframe1)
      return(tocategory(Z))
    }else if(parametric==F){
      u <- data.frame(matrix(ncol=ncol(dataframe1),nrow=nrow(dataframe1)))
      sample <- mvtnorm::rmvt(nrow(dataframe1), df=dof, sigma=cor(dataframe1,method=c('spearman')))
      for(i in 1:ncol(u)){
        u[,i] <- pt(sample[,i], df=dof)
      }
      for(j in 1:ncol(dataframe1)){
        if(class(dataframe[,j])=='numeric'){
          dataframe1[,j] <- qdist(u[,j],dataframe1[,j],type = 0)
        }else if(class(dataframe[,j])=='factor'){
          dataframe1[,j] <- qdist(u[,j],dataframe1[,j],type = 1)
        }else if(class(dataframe[,j])=='integer'){
          dataframe1[,j] <- qdist(u[,j],dataframe1[,j],type = 2)
        }
      }
      dataframe1 <- tocategory(dataframe1)
      dataframe1[sapply(dataframe1, is.character)] <- lapply(dataframe1[sapply(dataframe1, is.character)], as.factor)
      return(dataframe1)
    }
  
  }else{
    stop('Error: Undefined Copula')
  }
}



#Function to return diagnostic histograms(2 histograms side-by-side)
histfit <- function(df,position,ylim,xlab){
  l <- list(df[,1], df[,2])
  if(missing(position)){
    position='topright'
  }
   if(missing(xlab)){
    xlab='Value'
  }
  if(missing(ylim)){
    plotrix::multhist(l, xlab=xlab, ylab='Frequency', col=c('red', 'blue'))
    legend(position, legend=colnames(df), fill=c('red', 'blue'))
  }else{
    plotrix::multhist(l, xlab=xlab, ylab='Frequency', col=c('red', 'blue'),ylim=ylim)
    legend(position, legend=colnames(df), fill=c('red', 'blue'))
  }
}



#Function to compare correlation heat maps.
corrfit <- function(orig, syn){
  orig <- tonumeric(orig)
  syn <- tonumeric(syn)
  p1 <- ggcorrplot::ggcorrplot(cor(orig, method=c('spearman')), title = 'Original Dataset', ggtheme = ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5)))
  p2 <- ggcorrplot::ggcorrplot(cor(syn, method=c('spearman')), title='Synthetic Dataset', ggtheme = ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5)))
  gridExtra::grid.arrange(grobs=list(p1, p2), ncol=2)
}


#Function to plot (PC1 vs. PC2)  for Original & Synthetic datasets.
pcafit <- function(orig, syn){
  orig <- tonumeric(orig)
  syn <- tonumeric(syn)
  foo <- function(x) return(x/sqrt(sum(x^2)))
  cor_orig <- round(cor(orig, method='spearman'),2)
  cor_syn <- round(cor(syn, method='spearman'),2)
  eigen1 <- eigen(cor_orig)$vectors
  apply(eigen1,2,foo)
  eigen2 <- eigen(cor_syn)$vectors
  apply(eigen2,2,foo)
  PC_orig <- as.matrix(orig)%*%eigen1
  PC_syn <- as.matrix(syn)%*%eigen2
  par(mfrow=c(1,2))
  plot(PC_orig[,1], PC_orig[,2], xlab='PC1', ylab='PC2', main='Original Dataset', col='red')
  plot(PC_syn[,1], PC_syn[,2], xlab='PC1', ylab='PC2', main='Synthetic Dataset', col='blue')
  par(mfrow=c(1,1))
  
}

#Function to return pair plots
pairs <- function(df) psych::pairs.panels(df, method='spearman')

#Return RMSE of synthetic data
RMSE <- function(orig, syn) return(sqrt(sum((tonumeric(orig)-tonumeric(syn))^2)))

#Function to return cosine similarity between two datasets
cosineSim <- function(orig, syn){
  csim <- c()
  for(i in 1:ncol(orig)){
    c1 <- as.matrix(tonumeric(data.frame(orig[,i])))
    c2 <- as.matrix(tonumeric(data.frame(syn[,i])))
    csim[i] <- lsa::cosine(cbind(c1,c2))[1,2]
  }
  return(sum(csim)/ncol(orig))
}
