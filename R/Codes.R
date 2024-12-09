#' Multivariate FAR estimation
#'
#' This function estimates the functional coefficients of a multivariate FAR model at specific reference signal value u0.
#'
#' @param y A matrix containing the multivariate time series data
#' @param u A vector containing the reference signal series
#' @param p The order of FAR(p) model
#' @param d The lag of the reference signal
#' @param bwp The bandwidth proportion from 0 to 1 (with respect to the range of observations)
#' @param numpoints The number of discretized points in the range of the reference signal 'u' that will be used to evaluate the local linear approximation
#' @return A list of outputs including the functional coefficient estimates, the evaluation points and model residuals.
#' @export
FAR.est <- function(y,u,p,d,bwp = 0.1,numpoints = 50){

  lagrm <- max(p,d)
  k <- ncol(y)
  Tlength <- nrow(y)

  Mybounds <- as.numeric(stats::quantile(u,probs = c(0.05,0.95)))
  fhat_int <- seq(Mybounds[1],Mybounds[2],length = numpoints)
  fhat_points <- as.numeric(stats::filter(fhat_int, c(0.5,0.5)))[-length(fhat_int)]
  comdiff <- fhat_points[2] - fhat_points[1]
  fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)
  fhat_mean <- array(rep(NA,(k^2)*p*length(fhat_points)),dim = c(k,k*p,length(fhat_points)))

  for(i in 1:length(fhat_points)){
    try({
      fhat_mean[,,i] <- t(FARest(y,u,u0 = fhat_points[i],p,d,bwp))[,1:(k*p)]
    },silent = TRUE)
  }

  mat_Y <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
  mat_X <- matrix(NA,nrow = Tlength-lagrm,ncol = k*p)
  mat_U <- matrix(NA,nrow = Tlength-lagrm,ncol = 1)

  mat_Y[1:(Tlength-lagrm),] <- y[(lagrm+1):Tlength,]
  mat_U[1:(Tlength-lagrm),] <- u[(lagrm+1-d):(Tlength-d)]
  for(l in 1:p){
    for(j in 1:k){
      mat_X[,k*(l-1)+j] <- y[(lagrm+1-l):(Tlength-l),j]
    }
  }
  cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int,Inf)))
  mat_Yhat <- matrix(NA,nrow = Tlength-lagrm,ncol = k)

  for(t in 1:(Tlength-lagrm)){
    fhat_t <- fhat_mean[,,cat_U[t]]
    mat_Yhat[t,] <- t(fhat_t%*%matrix(mat_X[t,],ncol=1))
  }

  mat_resid <- mat_Y - mat_Yhat

  output <- list(mat_Y = mat_Y,mat_Yhat = mat_Yhat,mat_U = mat_U,
                 mat_resid = mat_resid, fhat_int = fhat_int,fhat_points = fhat_points,
                 fhat_mean = fhat_mean)

  return(output)

}

#' Multivariate MXFAR estimation
#'
#' This function estimates the functional coefficients of a multivariate mixed-effects FAR model at specific reference signal value u0.
#'
#' @param SeriesNum A numeric vector stating the number of series in a 'group' data
#' @param Tlength A numeric scalar stating the number of time points per series
#' @param y A matrix containing the 'stacked' multivariate time series data from multiple series
#' @param u A vector containing the 'stacked' reference signal series
#' @param p The order of MXFAR(p) model
#' @param d The lag of the reference signal
#' @param bwp The bandwidth proportion from 0 to 1 (with respect to the range of observations)
#' @param numpoints The number of discretized points in the range of the reference signal 'u' that will be used to evaluate the local linear approximation
#' @return A list of outputs including the functional coefficient (mean and series-specific) estimates, the evaluation points and model residuals.
#' @export
MXFAR.est <- function(SeriesNum,Tlength,y,u,p,d,bwp = 0.1,numpoints = 50){

  lagrm <- max(p,d)
  k <- ncol(y)
  g_num <- length(SeriesNum)
  SeriesNum_Total <- sum(SeriesNum)
  SeriesNum_cum <- cumsum(SeriesNum)
  SeriesNum_cumlag <- SeriesNum_cum - SeriesNum
  ts_desc <- list(SeriesNum = SeriesNum,Tlength = Tlength,p = p,
                  d = d, u = u, k = k)

  Mybounds <- as.numeric(stats::quantile(u,probs = c(0.05,0.95)))
  fhat_int <- seq(Mybounds[1],Mybounds[2],length = numpoints)
  fhat_points <- as.numeric(stats::filter(fhat_int, c(0.5,0.5)))[-length(fhat_int)]
  comdiff <- fhat_points[2] - fhat_points[1]
  fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)
  fhat_all <- array(rep(NA,2*(k^2)*p*(SeriesNum_Total+g_num)*length(fhat_points)),dim = c(k,2*k*p*(SeriesNum_Total+g_num),length(fhat_points)))
  fhat_mean <- array(rep(NA,(k^2)*p*length(fhat_points)*g_num),dim = c(k,k*p,length(fhat_points),g_num))
  fhat_indv <- array(rep(NA,(k^2)*p*length(fhat_points)*SeriesNum_Total),dim = c(k,k*p,length(fhat_points),SeriesNum_Total))

  for(i in 1:length(fhat_points)){
    try({
      fhat_all[,,i] <- MXFARest(SeriesNum,Tlength,y,u,u0 = fhat_points[i],p,d,bwp)
      for(s in 1:g_num){
        fhat_mean[,,i,s] <- fhat_all[,((s-1)*2*k*p+1):((s-1)*2*k*p+k*p),i]
        for(j in (SeriesNum_cumlag[s]+1):SeriesNum_cum[s]){
          fhat_indv[,,i,j] <- fhat_all[,((s-1)*2*k*p+1):((s-1)*2*k*p+k*p),i] + fhat_all[,((g_num+j-1)*2*k*p+1):((g_num+j-1)*2*k*p+k*p),i]
        }
      }
    },silent = TRUE)
  }

  mat_Y <- matrix(NA,nrow = SeriesNum_Total*(Tlength-lagrm),ncol = k)
  mat_X <- matrix(NA,nrow = SeriesNum_Total*(Tlength-lagrm),ncol = k*p)
  mat_U <- matrix(NA,nrow = SeriesNum_Total*(Tlength-lagrm),ncol = 1)

  for(i in 1:SeriesNum_Total){
    mat_Y[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),] <- y[((i-1)*Tlength+lagrm+1):(i*Tlength),]
    mat_U[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),] <- u[((i-1)*Tlength+lagrm+1-d):(i*Tlength-d)]
    for(l in 1:p){
      for(j in 1:k){
        mat_X[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),k*(l-1)+j] <- y[((i-1)*Tlength+lagrm+1-l):(i*Tlength-l),j]
      }
    }
  }

  cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int,Inf)))
  mat_Yhat <- matrix(NA,nrow = SeriesNum_Total*(Tlength-lagrm),ncol = k)

  for(s in 1:g_num){
    for(i in (SeriesNum_cumlag[s]+1):SeriesNum_cum[s]){
      for(t in 1:(Tlength-lagrm)){
        fhat_u <- fhat_all[,,cat_U[(i-1)*(Tlength-lagrm)+t]]
        fhat_t <- fhat_u[,((s-1)*2*k*p+1):((s-1)*2*k*p+k*p)] + fhat_u[,((g_num+i-1)*2*k*p+1):((g_num+i-1)*2*k*p+k*p)]
        mat_Yhat[(i-1)*(Tlength-lagrm)+t,] <- t(fhat_t%*%matrix(mat_X[(i-1)*(Tlength-lagrm)+t,],ncol=1))
      }
    }
  }

  mat_resid <- mat_Y - mat_Yhat

  output <- list(ts_desc = ts_desc,mat_Y = mat_Y,mat_Yhat = mat_Yhat,mat_U = mat_U,
                 mat_resid = mat_resid, fhat_int = fhat_int,fhat_points = fhat_points,
                 fhat_all = fhat_all,fhat_mean = fhat_mean,fhat_indv = fhat_indv)

  return(output)

}


#' Multi-fold cross validation for the FAR model
#'
#' This function calculates the accumulated prediction error (APE) from a fitted FAR model
#'
#' @param y A matrix containing the 'stacked' multivariate time series data from multiple series
#' @param u A vector containing the 'stacked' reference signal series
#' @param p The order of MXFAR(p) model
#' @param d The lag of the reference signal
#' @param r The number of time points for predictions
#' @param Q The number of segments considered for the time series
#' @param bwp The bandwidth proportion from 0 to 1 (with respect to the range of observations)
#' @param numpoints The number of discretized points in the range of the reference signal 'u' that will be used to evaluate the local linear approximation
#' @return The calculated APE metric.
aux_APE_indv <- function(y,u,p,d,r,Q,bwp,numpoints = 50){

  lagrm <- max(p,d)
  k <- ncol(y)
  Tlength <- nrow(y)
  ape_tot <- 0

  for(q in 1:Q){
    y_q <- y[1:(Tlength-q*r),]
    u_q <- u[1:(Tlength-q*r)]

    Mybounds <- as.numeric(stats::quantile(u_q,probs = c(0.05,0.95)))
    fhat_int_q <- seq(Mybounds[1],Mybounds[2],length = numpoints)
    fhat_points_q <- as.numeric(stats::filter(fhat_int_q, c(0.5,0.5)))[-length(fhat_int_q)]
    comdiff <- fhat_points_q[2] - fhat_points_q[1]
    fhat_points_q <- c(fhat_points_q[1]-comdiff,fhat_points_q,fhat_points_q[length(fhat_points_q)]+comdiff)
    fhat_mean_q <- array(rep(NA,(k^2)*p*length(fhat_points_q)),dim = c(k,k*p,length(fhat_points_q)))

    for(i in 1:length(fhat_points_q)){
      try({
        fhat_mean_q[,,i] <- t(FARest(y_q,u_q,u0 = fhat_points_q[i],p,d,bwp=bwp*(Tlength/(Tlength-q*r))^(1/5)))[,1:(k*p)]
      },silent = TRUE)
    }

    y_q_n <- y[(Tlength-q*r+1-lagrm):(Tlength-q*r+r),]
    u_q_n <- u[(Tlength-q*r+1-lagrm):(Tlength-q*r+r)]

    mat_Y <- matrix(NA,nrow = r,ncol = k)
    mat_X <- matrix(NA,nrow = r,ncol = k*p)
    mat_U <- matrix(NA,nrow = r,ncol = 1)

    mat_Y[1:r,] <- y_q_n[lagrm + 1:r,]
    mat_U[1:r,] <- u_q_n[lagrm + 1:r - d]
    for(l in 1:p){
      for(j in 1:k){
        mat_X[,k*(l-1)+j] <- y_q_n[lagrm + 1:r - l,j]
      }
    }

    cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int_q,Inf)))
    mat_Yhat <- matrix(NA,nrow = r,ncol = k)

    for(t in 1:r){
      fhat_t_q <- fhat_mean_q[,,cat_U[t]]
      mat_Yhat[t,] <- t(fhat_t_q%*%matrix(mat_X[t,],ncol=1))
    }

    mat_resid <- mat_Y - mat_Yhat

    ape_tot <- ape_tot + sum(mat_resid^2,na.rm = TRUE)

  }
  return(ape_tot)
}

#' Multi-fold cross validation for the MXFAR model
#'
#' This function calculates the accumulated prediction error (APE) from a fitted MXFAR model
#'
#' @param SeriesNum A numeric vector stating the number of series in a 'group' data
#' @param Tlength A numeric scalar stating the number of time points per series
#' @param y A matrix containing the 'stacked' multivariate time series data from multiple series
#' @param u A vector containing the 'stacked' reference signal series
#' @param p The order of MXFAR(p) model
#' @param d The lag of the reference signal
#' @param r The number of time points for predictions
#' @param Q The number of segments considered for the time series
#' @param bwp The bandwidth proportion from 0 to 1 (with respect to the range of observations)
#' @param numpoints The number of discretized points in the range of the reference signal 'u' that will be used to evaluate the local linear approximation
#' @return The calculated APE metric.
#' @export
APE <- function(SeriesNum,Tlength,y,u,p,d,r,Q,bwp=0.1,numpoints=50){

  SeriesNum_Total <- sum(SeriesNum)

  ape_tot <- 0
  for(i in 1:SeriesNum_Total){
    y_i <- y[((i-1)*Tlength+1):(i*Tlength),]
    u_i <- u[((i-1)*Tlength+1):(i*Tlength)]
    ape_tot <- ape_tot + aux_APE_indv(y_i,u_i,p,d,r,Q,bwp,numpoints)
  }

  return(ape_tot/SeriesNum_Total)

}


#' Bootstrap-based Non-linearity Test
#'
#' A nonparametric test for the presence of non-linearity in the multiple vector time series setting.
#'
#' @param SeriesNum A numeric vector stating the number of series in a 'group' data
#' @param Tlength A numeric scalar stating the number of time points per series
#' @param y A matrix containing the 'stacked' multivariate time series data from multiple series
#' @param u A vector containing the 'stacked' reference signal series
#' @param p The order of MXFAR(p) model
#' @param d The lag of the reference signal
#' @param bwp The bandwidth proportion from 0 to 1 (with respect to the range of observations)
#' @param numpoints The number of discretized points in the range of the reference signal 'u' that will be used to evaluate the local linear approximation
#' @param maxboot The number of bootstrap replicates to be considered.
#' @return A list containing the calculated test statistic and its p-value.
#' @export
NLTest <- function(SeriesNum,Tlength,y,u,p,d,bwp=0.1,numpoints = 50,maxboot=200){

  lagrm <- max(c(p,d))
  k <- ncol(y)
  SeriesNum_Total <- sum(SeriesNum)

  MXFAR_res <- array(NA,dim = c(SeriesNum_Total*(Tlength-lagrm),k))
  VAR_res <- array(NA,dim = c(SeriesNum_Total*(Tlength-lagrm),k))

  MXFAR_res_test <- array(NA,dim = c(SeriesNum_Total*(Tlength-2*lagrm),k))
  VAR_res_test <- array(NA,dim = c(SeriesNum_Total*(Tlength-2*lagrm),k))

  VAR_Yhat <- array(NA,dim = c(SeriesNum_Total*(Tlength-lagrm),k))
  MXFAR_u <- vector("numeric",length = SeriesNum_Total*(Tlength-lagrm))

  for(i in 1:SeriesNum_Total){
    y_i <- y[((i-1)*Tlength+1):(i*Tlength),]
    u_i <- u[((i-1)*Tlength+1):(i*Tlength)]
    VAR_est <- vars::VAR(y_i, p, type = "none")
    VAR_res[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),] <- stats::residuals(VAR_est)[-c(1:abs(p-lagrm)),]
    VAR_Yhat[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),] <- stats::fitted.values(VAR_est)[-c(1:abs(p-lagrm)),]
    FAR_est <- FAR.est(y_i,u_i,p,d,bwp,numpoints)
    MXFAR_res[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),] <- FAR_est$mat_resid-colMeans(FAR_est$mat_resid)
    MXFAR_u[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm))] <- FAR_est$mat_U[,1]

    MXFAR_res_test[((i-1)*(Tlength-2*lagrm)+1):(i*(Tlength-2*lagrm)),] <- FAR_est$mat_resid[-c(1:lagrm),]
    VAR_res_test[((i-1)*(Tlength-2*lagrm)+1):(i*(Tlength-2*lagrm)),] <- stats::residuals(VAR_est)[-c(1:abs(p-lagrm)),][-c(1:lagrm),]
  }

  teststat <- sum(diag(crossprod(VAR_res_test)),na.rm = TRUE)/sum(diag(crossprod(MXFAR_res_test)),na.rm = TRUE) - 1
  teststat_comp <- rep(NA,maxboot)

  for(b in 1:maxboot){

    samp_ind <- sample(1:SeriesNum_Total,SeriesNum_Total,replace = TRUE)
    y_b <- VAR_Yhat + MXFAR_res[(Tlength-lagrm)*(rep(samp_ind,each=Tlength-lagrm)-1)+rep(1:(Tlength-lagrm),times=SeriesNum_Total),]
    u_b <- MXFAR_u[(Tlength-lagrm)*(rep(samp_ind,each=Tlength-lagrm)-1)+rep(1:(Tlength-lagrm),times=SeriesNum_Total)]

    MXFAR_res_b <- array(NA,dim = c(SeriesNum_Total*(Tlength-2*lagrm),k))
    VAR_res_b <- array(NA,dim = c(SeriesNum_Total*(Tlength-2*lagrm),k))
    VAR_Yhat_b <- array(NA,dim = c(SeriesNum_Total*(Tlength-2*lagrm),k))

    for(i in 1:SeriesNum_Total){
      y_i_b <- y_b[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm)),]
      u_i_b <- u_b[((i-1)*(Tlength-lagrm)+1):(i*(Tlength-lagrm))]
      VAR_est_b <- vars::VAR(y_i_b, p, type = "none")
      VAR_res_b[((i-1)*(Tlength-2*lagrm)+1):(i*(Tlength-2*lagrm)),] <- stats::residuals(VAR_est_b)[-c(1:abs(p-lagrm)),]
      FAR_est_b <- FAR.est(y_i_b,u_i_b,p,d,bwp,numpoints)
      MXFAR_res_b[((i-1)*(Tlength-2*lagrm)+1):(i*(Tlength-2*lagrm)),] <- FAR_est_b$mat_resid
    }

    teststat_comp[b] <- sum(diag(crossprod(VAR_res_b)),na.rm = TRUE)/sum(diag(crossprod(MXFAR_res_b)),na.rm = TRUE) - 1

  }

  output <- list(teststat=teststat,teststat_comp=teststat_comp)
  output$pval <- sum(output$teststat_comp >= output$teststat,na.rm = TRUE)/(maxboot - sum(is.na(output$teststat_comp),na.rm = TRUE))

  return(output)

}


#' Multivariate FAR simulation
#'
#' This function simulates a multivariate time series from a multivariate FAR model with a chosen reference signal.
#'
#' @param Tlength A numeric scalar stating the number of time points
#' @param d The lag of the reference signal
#' @param Y_d The index of the time series to be used as endogenous reference signal. Put '0' for an exogenous standard normal reference signal
#' @param Fmat A list containing all functional coefficient matrices. The length of 'Fmat' dictates the order of the FAR model and the dimension of each matrix should match the dimension of the multivariate time series.
#' @param reff A numeric vector of random effects in 'Fmat'.
#' @param eps.Sigma A numeric vector of random error variances for the individual series.
#' @param nburn The number of time points that will be discarded after simulation to ensure convergence of the Markov Chain.
#' @return A list including the simulated time series and reference signal series.
#' @export
FAR.sim <- function(Tlength,d,Y_d,Fmat,reff,eps.Sigma,nburn = 500){

  P <- length(Fmat(0,reff))
  K <- ncol(Fmat(0,reff)[[1]])
  y <- matrix(NA,nrow = nburn+Tlength,ncol = K)
  ypast <- array(1,dim = c(P,K))
  ylag <- rep(1,d)

  e <- MASS::mvrnorm(n = nburn+Tlength,mu = rep(0,K+1),Sigma = diag(c(eps.Sigma,1)))

  for(t in 1:(nburn+Tlength)){

    Fmat_u <- Fmat(ylag[d],reff)
    FX <- matrix(0,nrow = K,ncol = 1)
    for(p in 1:P){
      FX <- FX + Fmat_u[[p]]%*%matrix(ypast[p,],nrow = K)
    }
    y[t,] <- t(FX + matrix(e[t,1:K],nrow = K))

    if(P>1){
      for(p in P:2){
        ypast[p,] <- ypast[p-1,]
      }
    }
    ypast[1,] <- y[t,]

    for(i in d:2){
      ylag[i] <- ylag[i-1]
    }
    if(Y_d == 0){
      ylag[1] <- e[t,3]
    } else {
      ylag[1] <- y[t,Y_d]
    }
  }

  y <- y[-c(1:500),]
  e <- e[-c(1:500),]

  if(Y_d == 0){
    output <- list(y = y, ref = e[,3])
  } else {
    output <- list(y = y, ref = y[,Y_d])
  }
  return(output)
}




#' Multivariate MXFAR simulation
#'
#' This function simulates multiple multivariate time series from a multivariate MXFAR model with a chosen reference signal.
#'
#' @param SeriesNum A numeric scalar stating the number of individual series to be simulated
#' @param Tlength A numeric scalar stating the number of time points per series
#' @param d The lag of the reference signal
#' @param Y_d The index of the time series to be used as endogenous reference signal. Put '0' for an exogenous standard normal reference signal
#' @param Fmat A list containing all functional coefficient matrices. The length of 'Fmat' dictates the order of the FAR model and the dimension of each matrix should match the dimension of the multivariate time series.
#' @param reff.Sigma A numeric vector of variances for the random effects in 'Fmat'.
#' @param eps.Sigma A numeric vector of random error variances for the individual series.
#' @param nburn The number of time points that will be discarded after simulation to ensure convergence of the Markov Chain.
#' @return A list including the simulated time series and reference signal series.
#' @export
MXFAR.sim <- function(SeriesNum,Tlength,d,Y_d,Fmat,reff.Sigma,eps.Sigma,nburn=500){

  SeriesNum_Total <- sum(SeriesNum)
  reff_s <- MASS::mvrnorm(n = SeriesNum_Total,mu = rep(0,length(reff.Sigma)),Sigma = diag(reff.Sigma))

  P <- length(Fmat(0,reff_s[1,]))
  K <- ncol(Fmat(0,reff_s[1,])[[1]])

  y <- matrix(NA,nrow = SeriesNum_Total*Tlength,ncol = K)
  ref <- matrix(NA,nrow = SeriesNum_Total*Tlength,ncol = 1)

  for(i in 1:SeriesNum_Total){
    entry <- FAR.sim(Tlength,d,Y_d,Fmat,reff_s[i,],eps.Sigma,nburn = nburn)
    y[((i-1)*Tlength+1):(i*Tlength),] <- entry$y
    ref[((i-1)*Tlength+1):(i*Tlength),] <- entry$ref
  }
  output <- list(y = y, ref = ref)
  return(output)
}


