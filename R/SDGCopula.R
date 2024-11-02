# Prerequisite packages
.packages = c('MASS', 'mvtnorm', 'copula', 'fitdistrplus', 'plotrix', 'psych', 'ggcorrplot', 'gridExtra', 'reshape2', 'ggplot2', 'lsa')
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(.packages, require, character.only=TRUE)


#Function to for parametric density estimation of marginals
#' Function to obtain the marginal distribution from the search space.
#' Search space = c('normal', 'cauchy', 'uniform', 'logistic', 'beta', 'exponential', 'gamma', 'Weibull', 'lognormal'). The method used for parameter estimation is maximum likelihood estimation.
#'
#' @param x numeric vector
#'
#' @return distribution name, parametr estimates and AIC value
#' @export
#'
#' @examples get.marginal(x)
get.marginal = function(x){
  if(length(x)>1){
    aic.vec = c()
    if (any(x < 0)){
      distributions <- c("norm", "unif","cauchy","logis")
      for(i in 1:length(distributions)){
        fit = fitdist(x,distributions[i],method = "mle")
        aic.vec[i] = fit$aic
      }
    }else{
      if(all(x >= 0 & x <= 1)){
        distributions <- c("beta","norm", "exp", "gamma", "weibull", "unif","cauchy","lnorm","logis")
        for(i in 1:length(distributions)){
          fit = fitdist(x,distributions[i],method = "mle")
          aic.vec[i] = fit$aic
        }}else{
          distributions <- c("norm", "exp", "gamma", "weibull", "unif","cauchy","lnorm","logis")
          for(i in 1:length(distributions)){
            fit = fitdist(x,distributions[i],method = "mle")
            aic.vec[i] = fit$aic
          }
        }
    }
    fit = fitdist(x,distributions[which.min(aic.vec)],method = "mle")
    return(list(distribution=distributions[which.min(aic.vec)],par = fit$estimate,AIC=fit$aic))
  }else(return("Error: Invalid Vector"))
}


#Function to return pairwise correlations between distinct variables
#' "A function that returns pairwise correlations between distinct variables.
#'
#' @param mat data matrix
#'
#' @return
#' @export
#'
#' @examples
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
#' Function to obtain the kernel density estimate. From that density, draw a sample and return the quantiles.
#'
#' @param u
#' @param input
#' @param type
#'
#' @return quantiles
#' @export
#'
#' @examples
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
#' Function to obtain the parametric estimates of marginals.
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
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
#' Function to convert categorical variables into numeric values.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
tonumeric <- function(x){
  isfactor <- !sapply(x, is.numeric)
  codes <- lapply(x[isfactor], \(x)levels(factor(x)))
  x[isfactor] <- data.matrix(x[isfactor])
  structure(x, codes = codes)
}

#Function to convert numeric back to original categorical variables.
#' Function to convert numeric values to the original categorical variables.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
tocategory <- function(x){
  codes <- attr(x, 'codes')
  x[names(codes)] <- Map(\(i,j)j[i], x[names(codes)], codes)
  x
}


#Function to fit Normal/t-copula for dim>=2
#' This function returns a synthetic data frame.
#'
#' @param dataframe from which synthetic data is to be generated.
#' @param copula normal or t, default = 'normal'
#' @param parametric TRUE or FALSE, default = FALSE
#' @param dof degree of freedom for t copula, default dof = 4
#'
#' @return A data frame containing synthetic data.
#' @export
#'
#' @examples fitCop(dataframe = df,copula = normal, parametric = F)
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
      myCop <- normalCopula(param=cor_struct(cor(dataframe1, method=c('spearman'))), dim = ncol(dataframe1), dispstr = "un")
      myMvd <- mvdc(myCop, margins=marginals(dataframe1)$dist, paramMargins = marginals(dataframe1)$par)
      Z <- rMvdc(nrow(dataframe1), myMvd)
      colnames(Z) <- colnames(dataframe1)
      return(tocategory(Z))
    }else if(parametric==F){
      u <- data.frame(matrix(ncol=ncol(dataframe),nrow=nrow(dataframe)))
      sample <- mvrnorm(nrow(dataframe1),mu=rep(0,ncol(dataframe1)),Sigma=cor(dataframe1,method=c('spearman')),empirical=T)
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
      myCop <- tCopula(df=dof, param=cor_struct(cor(dataframe1, method=c('spearman'))), dim = ncol(dataframe1), dispstr = "un")
      myMvd <- mvdc(myCop, margins=marginals(dataframe1)$dist, paramMargins = marginals(dataframe1)$par)
      Z <- rMvdc(nrow(dataframe1), myMvd)
      colnames(Z) <- colnames(dataframe1)
      return(tocategory(Z))
    }else if(parametric==F){
      u <- data.frame(matrix(ncol=ncol(dataframe1),nrow=nrow(dataframe1)))
      sample <- rmvt(nrow(dataframe1), df=dof, sigma=cor(dataframe1,method=c('spearman')))
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
#' A function to compare the histograms of original and synthetic data.
#'
#' @param df Data frame containing both original and synthetic data.
#' @param position Position of the legend.
#' @param ylim Set y-axis limits if necessary.
#'
#' @return Histogram plot comparing original and synthetic data.
#' @export
#'
#' @examples histfit(df = dataframe, position = 'topright')
histfit <- function(df,position,ylim){
  l <- list(df[,1], df[,2])
  if(missing(position)){
    position='topright'
  }
  if(missing(ylim)){
    multhist(l, xlab='Value', ylab='Frequency', col=c('red', 'blue'))
    legend(position, legend=colnames(df), fill=c('red', 'blue'))
  }else{
    multhist(l, xlab='Value', ylab='Frequency', col=c('red', 'blue'),ylim=ylim)
    legend(position, legend=colnames(df), fill=c('red', 'blue'))
  }
}



#Function to compare correlation heat maps.
#' Function to create a side-by-side heatmap of original and synthetic data for comparison.
#'
#' @param orig Original data.
#' @param syn Synthetic data.
#'
#' @return Side-by-side heatmap plots of both data sets.
#' @export
#'
#' @examples corrfit(orig=data,syn=synthetic_data)
corrfit <- function(orig, syn){
  orig <- tonumeric(orig)
  syn <- tonumeric(syn)
  p1 <- ggcorrplot(cor(orig, method=c('spearman')), title = 'Original Dataset')+theme(plot.title = element_text(hjust=0.5))
  p2 <- ggcorrplot(cor(syn, method=c('spearman')), title='Synthetic Dataset')+theme(plot.title = element_text(hjust=0.5))
  grid.arrange(grobs=list(p1, p2), ncol=2)
}


#Function to plot (PC1 vs. PC2)  for Original & Synthetic datasets.
#' Function to create a plot of principal component PC1 vs. PC2 for comparing original and synthetic data.
#'
#' @param orig original data.
#' @param syn synthetic data.
#'
#' @return Side-by-side PC1 vs. PC2 plots of both data sets.
#' @export
#'
#' @examples pcafit(orig=data,syn=synthetic_data)
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
#' Function to create a scatter plot of data.
#'
#' @param df data.
#'
#' @return Pairs plot.
#' @export
#'
#' @examples pairs(data)
pairs <- function(df) pairs.panels(df, method='spearman')

#Return RMSE of synthetic data
#' Function to calculate root mean square error.
#'
#' @param orig original data.
#' @param syn synthetic data.
#'
#' @return RMSE value.
#' @export
#'
#' @examples RMSE(orig=data,syn=synthetic_data)
RMSE <- function(orig, syn) return(sqrt(sum((tonumeric(orig)-tonumeric(syn))^2)))


#Function to return cosine similarity between two datasets
#' Compare angles
#'
#' @param orig original data.
#' @param syn synthetic data
#'
#' @return value between -1 to 1 which indicate similarity between two data frames.
#' @export
#'
#' @examples cosineSim(orig, syn)
cosineSim <- function(orig, syn){
  csim <- c()
  for(i in 1:ncol(orig)){
    c1 <- as.matrix(tonumeric(data.frame(orig[,i])))
    c2 <- as.matrix(tonumeric(data.frame(syn[,i])))
    csim[i] <- cosine(cbind(c1,c2))[1,2]
  }
  return(sum(csim)/ncol(orig))
}
