#' CirucmOLS
#'
#' @title Circumplex Models Ordinary Least Square Estimation
#' @name CircumOLS
#' @description The function is to compute ordinary least square estimates of circumplex model parameters and their standard errors. It allows both continuous data and ordinal data.
#' @usage CircumOLS(rawdt, m=1, mcsc="unconstrained", type="ordinal", simulation_based = F, icorrection = 0, ncore = 8, maxit = 1000, factr = 1e7 ,pgtol = 0, lmm = NULL, N_star = 10000)
#' @param rawdt indicates the raw data which is a n by p matrix where n is the number of participants and p is the number of manifest variables.
#' @param m indicates the number of cosine function coefficients. Defaults to 1.
#' @param mcsc minimum common score correlation value: \code{"unconstrained"} (default), \code{-1}
#' @param type indicates the type of data: \code{"ordinal"} (default), \code{"continuous"}
#' @param simulation_based determines whether an asymptotic covariance matrix estimate is a sample-based ACM estimate or a Monte Carlo ACM estimate: \code{False} (default, sample-based ACM), \code{True} (Monte Carlo ACM)
#' @param icorrection Methods to adjust for empty cells: a scalar where 0 is no adjustment is done (default), 1  adds 1/(nc*nr) to all cells where nc and nr is the number of columns and rows of the contingency table respectively, 2 adds 0.1 to all cells, 3 adds 0.5 to all cells, 11 adds 1/(nc*nr) to only zero cells, 12 adds 0.1 to only zero cells, and 13 adds 0.5 to only zero cells
#' @param ncore indicates the number of cores for parallel computing. It defaults to 8.
#' @param maxit The maximum number of iterations. It defaults to 1000.
#' @param factr controls the convergence of the "L-BFGS-B" method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default is 1e7, that is a tolerance of about 1e-8.
#' @param pgtol helps control the convergence of the "L-BFGS-B" method. It is a tolerance on the projected gradient in the current search direction. This defaults to zero, when the check is suppressed.
#' @param lmm is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method. It defaults to number of free parameters.
#' @param N_star is the sample size of a simulated sample for Monte Carlo ACM. It defaults to 10000.
#' @author Dayoung Lee \email{dlee33@nd.edu}
#' @details Code modified from function CircE.BFGS obtained from https://cran.r-project.org/src/contrib/Archive/CircE/CircE_1.1.tar.gz
#' @references Browne, M. W. (1992). Circumplex models for correlation matrices. Psychometrika, 57, 469–497. doi: 10.1007/BF02294416
#' @references Monroe, S. (2018). Contributions to estimation of polychoric correlations. Multivariate Behavioral Research, 53, 247–266. doi: 10.1080/00273171.2017.1419851
#' @references Grassi, M., Luccio, R., & Di Blas, L. (2010). Circe: An r implementation of browne’s circular stochastic process model. Behavior Research Methods, 42(1), 55–73. doi:10.3758/BRM.42.1.55
#' @references Lee, D., & Zhang, G. (in preparation). Circumplex models with ordinal data.
#' @references Zhang, G., Trichtinger, L., Lee, D., & Jiang, G. (2021). Polychoricrm: A computationally efficient r function for estimating polychoric correlations and their asymptotic covariance matrix. Structural Equation Modeling: A Multidisciplinary Journal. doi:10.1080/10705511.2021.1929996
#' @examples #Examples using the data set included in the packages: data("emotions") CircumOLS(emotions, m = 1, type = "ordinal")    # Big-five inventory (N = 228)
#' @importFrom EFAutilities efa
#' @importFrom MASS mvrnorm
#' @importFrom Turbofuns PolychoricRM
#' @export CircumOLS
#' @method print CircumOLS
#' @export

CircumOLS <- function(rawdt, m=1, mcsc="unconstrained",
                      type="ordinal", simulation_based = F, icorrection = 0, ncore = 8,
                      maxit = 1000, factr = 1e7 ,pgtol = 0, lmm = NULL, N_star = 10000){

  m = m
  N = nrow(rawdt)
  p = ncol(rawdt)

  library(Turbofuns)
  library(EFAutilities)
  library(MASS)

  if (type == "ordinal"){
    Lpoly = PolychoricRM(iRaw = rawdt,IAdjust = 0, NCore = 8,estimate.acm=T)
    dt = Lpoly[[2]]
  } else {dt = cor(rawdt)}

  k=3
  K=pi/180
  r=1
  mcsc= mcsc
  title="Circumplex Estimation"
  print.level=1

  start.value = function(dt){
    r=r
    if (type == "ordinal"){
      if (sum(eigen(Lpoly[[4]])$values < 0) == 0){
        lambda = efa(covmat = dt, factors = 3, n.obs = N, mtest = F, se = "none")$unrotated
      } else {lambda = efa(rawdt, dist = "ordinal", factors = 3, n.obs = N, mtest = F, se = "none")$unrotated}
    } else {lambda = efa(covmat = dt, factors = 3, n.obs = N, mtest = F, se = "none")$unrotated}


    k=3
    dpsi = diag(diag(dt - lambda%*%t(lambda)))

    dzeta = sqrt(diag(diag(lambda%*%t(lambda))))
    dupsilon = solve(dzeta)%*%dpsi

    lambdastar = solve(dzeta)%*%lambda
    lambdastarbar = (1/p)*t(lambdastar)%*%c(rep(1,p))

    cmating = lambdastar - (c(rep(1,p))%*%t(lambdastarbar))
    cmat = (1/p)*t(cmating)%*%cmating

    U = eigen(cmat)$vectors
    lambdasim = lambdastar%*%U
    lambdadot = matrix(data=NA, nrow = p, ncol = 2)

    for(i in 1:p)
    {for(j in 1:2)
    {lambdadot[i,j] = lambdasim[i,j] / sqrt(lambdasim[i,1]^2 + lambdasim[i,2]^2)}
    }

    theta = function(lambdadot)
    {a = rep(0, nrow(lambdadot))

    for (i in seq(1:nrow(lambdadot))[-r])
    {if ((lambdadot[i,2]*lambdadot[r,1] - lambdadot[i,1]*lambdadot[r,2]) >= 0)
    {a[i] = acos(lambdadot[i,1]*lambdadot[r,1] + lambdadot[i,2]*lambdadot[r,2])*180/pi}
      else
      {a[i] = 360 - acos(lambdadot[i,1]*lambdadot[r,1] + lambdadot[i,2]*lambdadot[r,2])*180/pi}
    }

    return(a)
    }

    angle = theta(lambdadot)
    v = diag(dupsilon)

    if(mcsc=="unconstrained"){
      beta0 = ((1/p)*sum(lambdasim[,3]))^2
      beta1 = 1 - beta0

      betas=matrix(0,1,m+1);
      betas[1] = beta0
      betas[2] = beta1
      betas = betas/betas[2]
      par <- c(angle[-c(r)]*(pi/180),betas[-c(2)],v)
    }
    else {
      betas=matrix(0,1,m+1);
      betas[1] = (1+mcsc)/2
      betas[2] = (1-mcsc)/2
      betas = betas/betas[2]
      par <- c(angle[-c(r)]*(pi/180),betas[-c(2)],v)
    }
    return(par)
  }

  par = start.value(dt)

  ols<-
    function(par){
      R = dt
      p = ncol(dt)
      ang1=par[1:(p-1)]
      ang=append(ang1,0,r-1)   #ang=ang2*(pi/180)
      if (m==1) {alpha = c(par[(p-1+1)],1)
      } else {alpha=  c(par[((p-1)+1)],1,par[((p-1)+2):((p-1)+(m))])}
      if (mcsc=="unconstrained") {b=alpha/sum(alpha)
      } else if (mcsc==-1) {alpha[seq(1,(m+1),by=2)]=0; b=alpha/sum(alpha)
      if (m<=2) b[2] =1 else b[2] = 1 - sum(b[seq(4,(m+1),by=2)])
      } else {if (m==1){b=c((mcsc+1)/2,(1-mcsc)/2)
      }else if (m==2){b=alpha/sum(alpha);b[1]=(mcsc+1)/2 - b[3]; b[2]=(1-mcsc)/2
      } else {b = alpha/sum(alpha);b[1] = (mcsc+1)/2 - sum(b[seq(3,(m+1),by=2)])
      b[2] = (1-mcsc)/2 - sum(b[seq(4,(m+1),by=2)])}}

      v=par[(length(par)-p+1):length(par)]
      z2 = rep(0,p)
      for (i in 1:p)
      {z2[i] = (1/(1+v[i]))^(1/2)}
      M=matrix(c(0),p,p,byrow=TRUE)
      for(i in 1:p){
        for(j in 1:p){
          M[i,j]=c(b[-1])%*%cos(c(1:m)*(ang[j]-ang[i]))
        }}
      Pc=M+matrix(b[1],p,p)
      Dv=diag(v,p)
      Dz2=diag(z2);
      Residual = R - Dz2%*%(Pc+Dv)%*%Dz2
      f=sum(Residual^2)/2;
      f}

  gradient<-
    function(par){

      R = dt
      p = ncol(dt)
      ang1=par[1:(p-1)]
      ang=append(ang1,0,r-1)
      if (m==1) {alpha = c(par[(p-1+1)],1)
      } else {alpha=  c(par[((p-1)+1)],1,par[((p-1)+2):((p-1)+(m))])}

      if (mcsc=="unconstrained") {b=alpha/sum(alpha)
      } else if (mcsc==-1) {alpha[seq(1,(m+1),by=2)]=0; b=alpha/sum(alpha)
                            if (m<=2) b[2] =1 else b[2] = 1 - sum(b[seq(4,(m+1),by=2)])
      } else {if (m==1){b=c((mcsc+1)/2,(1-mcsc)/2)
              } else if (m==2){b=alpha/sum(alpha);b[1]=(mcsc+1)/2 - b[3]; b[2]=(1-mcsc)/2
              } else {b = alpha/sum(alpha);b[1] = (mcsc+1)/2 - sum(b[seq(3,(m+1),by=2)])
                      b[2] = (1-mcsc)/2 - sum(b[seq(4,(m+1),by=2)])}}

      v=par[(length(par)-p+1):length(par)]
      z2 = rep(0,p)
      for (i in 1:p)
      {z2[i] = (1/(1+v[i]))^(1/2)}
      M=matrix(c(0),p,p,byrow=TRUE)
      for(i in 1:p){
        for(j in 1:p){
          M[i,j]=c(b[-1])%*%cos(c(1:m)*(ang[j]-ang[i]))
        }}
      Pc=M+matrix(b[1],p,p)
      Dv=diag(v,p)
      Dz2=diag(z2)

      # theta

      angi <- function(ang,i,j,m){
        x = rep(0, m)

        for (k in 1:m)
        {x[k] = k*b[k+1]*sin(k*(ang[j]-ang[i]))}

        z2[i]*z2[j]*sum(x)}

      angj <- function(ang,i,j,m)
      {x = rep(0, m)

      for (k in 1:m)
      {x[k] = -k*b[k+1]*sin(k*(ang[j]-ang[i]))}

      z2[i]*z2[j]*sum(x)}

      dPxdthe <- list()

      for (ii in 1:(p-1)){dPxdthe[[ii]] <- matrix(0, (p-ii), (p))}

      for (ii in 1:(p-1))
      {
        for (l in 1:(p-ii))
        { dPxdthe[[ii]][l,ii] <- angi(ang,ii,(ii+l),m)
        dPxdthe[[ii]][l,(ii+l)] <- angj(ang,ii,(ii+l),m)
        }
      }

      dPxdthe <- do.call(rbind, dPxdthe)

      # beta

      if (mcsc == "unconstrained") {
        beta = list()
        for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), (m))}

        for (ii in 1:(p-1))
        {
          for (l in 1:(p-ii))
          {
            bt <- rep(0, m)
            for (mm in 1:m) {bt[(mm)] <- z2[ii]*z2[ii+l]*(cos(mm*(ang[ii+l]-ang[ii]))/sum(alpha) - (alpha[-1]%*%(cos((1:m)*(ang[ii+l]-ang[ii]))))/sum(alpha)^2 - alpha[1]/sum(alpha)^2)}
            bt[1] <- z2[ii]*z2[ii+l]*(1/sum(alpha)-alpha[1]/sum(alpha)^2 -1*sum( (alpha[-c(1)]/sum(alpha)^2)*cos(c(1:m)*(ang[ii+l]-ang[ii]) ) ) )

            beta[[ii]][l,] <- bt
          }
        }


        dPxdb <- do.call(rbind, beta)
      } else if (mcsc == -1){

        beta = list()
        for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

        for (ii in 1:(p-1))
        {
          for (l in 1:(p-ii))
          {
            bt <- rep(0, m)
            for (mm in 1:m) {bt[(mm)] <- z2[ii]*z2[ii+l]*(cos(mm*(ang[ii+l]-ang[ii]))/sum(alpha) - (alpha[-1]%*%(cos((1:m)*(ang[ii+l]-ang[ii]))))/sum(alpha)^2 - alpha[1]/sum(alpha)^2)}
            bt[1] <- z2[ii]*z2[ii+l]*(1/sum(alpha)-alpha[1]/sum(alpha)^2 -1*sum( (alpha[-c(1)]/sum(alpha)^2)*cos(c(1:m)*(ang[ii+l]-ang[ii]) ) ) )

            beta[[ii]][l,] <- bt
          }
        }

        dPxdb <- do.call(rbind, beta)

        if (m <= 2) {dPxdb[,1:(m)] <- 0} else {dPxdb <- do.call(rbind, beta)
        dPxdb[,seq(2,m,by=2)] = 0; dPxdb[,1] = 0
        }
      } else {

        beta = list()
        for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

        for (ii in 1:(p-1))
        {
          for (l in 1:(p-ii))
          {
            bt <- rep(0, m)
            if (m==1) {for (mm in seq(1,m,by=2)) {bt[mm] <- z2[ii]*z2[ii+l]*(cos(mm*(ang[ii+l]-ang[ii]))-cos(ang[ii+l]-ang[ii]))}
            }else {
              for (mm in seq(1,m,by=2)) {bt[mm] <- z2[ii]*z2[ii+l]*(((1/sum(alpha))-alpha[mm]/sum(alpha)^2)*(cos(mm*(ang[ii+l]-ang[ii]))-cos(ang[ii+l]-ang[ii])))}
              for (mm in seq(2,m,by=2)) {bt[mm] <- z2[ii]*z2[ii+l]*(((1/sum(alpha))-alpha[mm]/sum(alpha)^2)*(cos(mm*(ang[ii+l]-ang[ii]))-1))}}
            bt[1] <- z2[ii]*z2[ii+l]*(1/sum(alpha)-alpha[1]/sum(alpha)^2 -1*sum( (alpha[-c(1)]/sum(alpha)^2)*cos(c(1:m)*(ang[ii+l]-ang[ii]) ) ) )
            beta[[ii]][l,] <- bt}
        }


        dPxdb <- do.call(rbind, beta)
        dPxdb[,1] = 0
      }

      # nu

      dPxdv = matrix(0,p*(p-1)/2,p)

      for (l in 1:p) {

        kk <- matrix(0,p,p)
        diag(kk)[l] <- (-1/2)*((1+v[l])^(-3/2))

        jj <- matrix(0,p,p)
        diag(jj)[l] <- 1

        dPxdvl = (kk%*%(Pc+Dv)%*%Dz2) + t(kk%*%(Pc+Dv)%*%Dz2) + Dz2%*%jj%*%Dz2

        dPxdv[,l] <-  as.vector(dPxdvl[lower.tri(dPxdvl)])

      }


      dPxdr = cbind(dPxdthe[,-r], dPxdb, dPxdv)

      res = dt - Dz2%*%(Pc+Dv)%*%Dz2
      res = as.matrix(res[lower.tri(res)])

      dfdr = -2*t(dPxdr)%*%res

      dfdr}

  up = c(rep(Inf,(2*p + m -1)))
  low = c(rep(-Inf,(p-1)),rep(0,m),rep(0,p))
  control <- list(trace=print.level,
                  REPORT=20,
                  maxit = maxit,
                  factr = factr,
                  pgtol = pgtol,
                  lmm = if(is.null(lmm)){length(par)} else lmm)

  est <- optim(par, fn = ols, gr = gradient, method="L-BFGS-B", lower = low, upper = up, control = control)
  if (sum(est$par[(p+m):(2*p+m-1)]==0) != 0) warning("Heywood case occurs")
  if (est$convergence != 0) warning("Convergence error occurs")

  test.stat = function(par)
  {
    R = dt
    p = ncol(dt)
    ang1=par[1:(p-1)]
    ang=append(ang1,0,r-1)
    b = par[(p-1+1):(p-1+m)]/(1+sum(par[(p-1+1):(p-1+m)]))
    if (m==1) {b = c(b[1],1-sum(b[1]))
    } else {b = c(b[1],1-sum(b),b[-1])}

    v=par[(length(par)-p+1):length(par)]
    z2 = rep(0,p)
    for (i in 1:p)
    {z2[i] = (1/(1+v[i]))^(1/2)}
    M=matrix(c(0),p,p,byrow=TRUE)
    for(i in 1:p){
      for(j in 1:p){
        M[i,j]=c(b[-1])%*%cos(c(1:m)*(ang[j]-ang[i]))
      }}
    Pc=M+matrix(b[1],p,p)
    Dv=diag(v,p)
    Dz2=diag(z2)

    if (type == "ordinal"){
      if (simulation_based == T){
        set.seed(1)
        if (sum(eigen(dt)$values >= 0) == p){samp = mvrnorm(n = N_star, mu = rep(0,p), dt, tol = 1e-06, empirical = F)
        } else {samp = mvrnorm(n = N_star, mu = rep(0,p), Dz2%*%(Pc+Dv)%*%Dz2, tol = 1e-06, empirical = F)}

        samp2 = samp

        lmt = which(Lpoly[[1]][,1] == 1e+10)

        for (j in 1:p){
          for (i in 2:lmt){
            samp2[which(samp[,j]>=Lpoly[[1]][(i-1),j]
                        & samp[,j]<Lpoly[[1]][i,j]),j] = i - 2}
        }
        Y.Hat = PolychoricRM(iRaw=samp2, IAdjust=3, NCore=ncore, estimate.acm=T)[[4]]
      } else {
        Y.Hat = Lpoly[[4]]}
    } else {
      p2 = p * p
      p.star = p*(p-1)/2
      EliU <- function(MP, eta = 1) {

        ## Browne, M. W. & Shapiro, A. (1986). The asymptotic Covariance matrix of
        ## sample correlation coefficients under general conditions. Linear Algebra
        ## and its applications, 82, 169-176.
        ## Equations (4.1) and (4.3)

        ## EliU -> The asymptotic covariance matrix of sample correlations if manifest variables are
        ## of an elliptical distribution.

        # It does not require any external functions.

        p = dim(MP)[1]

        Ms= matrix(0,p*p,p*p)

        for (j in 1:p) {
          for (i in 1:p)  {

            if (j==i) {
              ii = (i-1)*p + i
              Ms[ii,ii] = 1
            } else
            {
              ij = (j-1)*p + i
              ji = (i-1)*p + j
              Ms[ij,ij] = 0.5
              Ms[ij,ji] = 0.5
            }

          } # i
        } # j


        Kd = matrix(0,p*p,p)
        for (i in 1:p) {
          ii = (i-1) * p + i
          Kd[ii,i] = 1
        }


        A = Ms %*% (MP %x% diag(p)) %*% Kd

        Gamma = 2 * Ms %*% (MP %x% MP)

        if (eta != 1 ) {
          MP.v = array(MP)
          Gamma = eta * Gamma + (eta -1) * outer(MP.v, MP.v)
        }

        B = Gamma %*% Kd
        G = t(Kd) %*% Gamma %*% Kd

        Cov.r = Gamma - A %*% t(B) - B %*% t(A) + A %*% G %*% t(A)

      } # EliU

      u.r = EliU(dt)

      M.Select = matrix(0,p2,p.star)
      ij=0
      ij.new = 0
      for (j in 1:p) {
        for (i in 1:p) {
          ij = ij + 1
          if (i<j) {
            ij.new = ij.new + 1
            M.Select[ij, ij.new] = 1
          }
        } # i
      } # j

      Y.Hat = t(M.Select) %*% u.r %*% M.Select}


    angi <- function(ang,i,j,m){
      x = rep(0, m)

      for (k in 1:m)
      {x[k] = k*b[k+1]*sin(k*(ang[j]-ang[i]))}

      z2[i]*z2[j]*sum(x)}

    angj <- function(ang,i,j,m)
    {x = rep(0, m)

    for (k in 1:m)
    {x[k] = -k*b[k+1]*sin(k*(ang[j]-ang[i]))}

    z2[i]*z2[j]*sum(x)}

    dPxdthe <- list()

    for (ii in 1:(p-1)){dPxdthe[[ii]] <- matrix(0, (p-ii), (p))}

    for (ii in 1:(p-1))
    {
      for (l in 1:(p-ii))
      { dPxdthe[[ii]][l,ii] <- angi(ang,ii,(ii+l),m)
      dPxdthe[[ii]][l,(ii+l)] <- angj(ang,ii,(ii+l),m)
      }
    }

    dPxdthe <- do.call(rbind, dPxdthe)

    # beta

    if (mcsc == "unconstrained") {
      beta = list()
      for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

      for (ii in 1:(p-1))
      {
        for (l in 1:(p-ii))
        {
          bt <- rep(0, m)
          for (mm in 1:(m)) {bt[mm] <- z2[ii]*z2[ii+l]*(cos(mm*(ang[ii+l]-ang[ii]))-1)}
          beta[[ii]][l,] <- bt
        }
      }


      dPxdb <- do.call(rbind, beta)
    } else {

      beta = list()
      for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

      for (ii in 1:(p-1))
      {
        for (l in 1:(p-ii))
        {
          bt <- rep(0, m)
          for (mm in 1:(m)) {bt[mm] <- z2[ii]*z2[ii+l]*(cos(mm*(ang[ii+l]-ang[ii]))-cos(ang[ii+l]-ang[ii]))}
          beta[[ii]][l,] <- bt
        }
      }

      dPxdb <- do.call(rbind, beta)

      if (m <= 2) {dPxdb[,1:m] <- 0} else {dPxdb <- do.call(rbind, beta)
      dPxdb[,seq(2,m,by=2)] = 0; dPxdb[,1] = 0}}

    # nu part

    dPxdv = matrix(0,p*(p-1)/2,p)

    for (l in 1:p) {

      kk <- matrix(0,p,p)
      diag(kk)[l] <- (-1/2)*((1+v[l])^(-3/2))

      jj <- matrix(0,p,p)
      diag(jj)[l] <- 1

      dPxdvl = (kk%*%(Pc+Dv)%*%Dz2) + t(kk%*%(Pc+Dv)%*%Dz2) + Dz2%*%jj%*%Dz2

      dPxdv[,l] <-  as.vector(dPxdvl[lower.tri(dPxdvl)])

    }


    dPxdr = cbind(dPxdthe[,-1], dPxdb, dPxdv)

    res = dt - Dz2%*%(Pc+Dv)%*%Dz2
    res = as.matrix(res[upper.tri(res)])

    ex = matrix(0,p,p)
    ex[lower.tri(ex)] = 1:(p*(p-1)/2)
    ex[upper.tri(ex)] = 1:(p*(p-1)/2)
    ex2 = t(ex)

    Delta = dPxdr[ex2[upper.tri(ex2)],]

    null <-function(M)
    {
      tmp <- qr(M)
      set <- if (tmp$rank == 0L)
        seq_len(ncol(M))
      else -seq_len(tmp$rank)
      qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
    }

    Del.c = null(Delta)

    F0 = t(res) %*% Del.c %*% solve(t(Del.c) %*% Y.Hat %*% Del.c) %*% t(Del.c) %*% res
    ts.st = N * F0
    pstar = dim(Del.c)[1]
    q = pstar - dim(Del.c)[2]
    df = pstar - q

    if (ts.st < df) {

      RMSEA = 0

    } else {

      RMSEA = sqrt((ts.st - df)/(df*N))

    }
    p.perfect = 1 - pchisq(ts.st,df,ncp=0)
    p.close = 1 - pchisq(ts.st,df,ncp=0.0025*N*df)

    # Hessian
    # row: dndn dndb dndt dbdb dbdt dtdt
    # col: correlation coeffiicient

    # dtdt
    # i =/ j case
    k = c(1:m)
    hh = list()

    for (i in 1:(p-1))
    {hh[[i]] = rep(0,p-i)
    for (j in (i+1):p)
    {hh[[i]][(j-i)] = (k^2%*%(b[2:(m+1)]*cos(k*(ang[j]-ang[i]))))*z2[j]*z2[i]}
    }

    xx = rep(0,p*(p-1)/2)
    xx[(i*p - i*(i+1)/2 + (j-i))] = hh[[i]][(j-i)]

    dPxdtdt.1 = matrix(0,(p*(p-1)/2),p*(p-1)/2)
    diag(dPxdtdt.1) <- unlist(hh)
    dPxdtdt.1 <- dPxdtdt.1[-c(1:p-1),]

    # i=j

    nn=list()

    for (i in (1:(p-1)))
    {nn[[i]] = matrix(0,p,(p-i))
    for (j in 1:(p-i))
    {nn[[i]][i,j] = (-k^2%*%(b[2:(m+1)]*cos(k*(ang[(i+j)]-ang[i]))))*z2[(i+j)]*z2[i]}
    }

    dPxdtdt.2 <- do.call(cbind, nn)
    dPxdtdt.2 <- dPxdtdt.2[-1,]

    dPxdtdt.2[1,1] = (-k^2%*%(b[2:(m+1)]*cos(k*(ang[(2)]-ang[1]))))*z2[(2)]*z2[1]

    for (i in 2:(p-1))
    { NN = c(i,(p-2):((p-2)-(i-2)))
    WW = rep(0,(i))
    for (MM in 1:i){WW[MM] = sum(NN[1:MM])}
    for (JJ in WW)
    { j = which(JJ == WW)
    dPxdtdt.2[i,JJ] = (-k^2%*%(b[2:(m+1)]*cos(k*(ang[(i+1)]-ang[j]))))*z2[(i+1)]*z2[j]}
    }

    d2Pxdtdt <- rbind(dPxdtdt.1, dPxdtdt.2)

    # dtdb

    dPxdtdbi=list()

    for (i in 1:(p-1))
     {dPxdtdbi[[i]] = matrix(0,(p-1)*m,(p-i))
      for (j in ((1+m*(i-2)):(1+m*(i-2)+m-1)))
        {q = j - m*(i-2)
         for (K in ((1+i):p))
          {dPxdtdbi[[i]][j,(K-i)] = (q*sin(q*(ang[K]-ang[i])))*z2[K]*z2[i]}
        }
     }

    d2Pxdbdti <- do.call(cbind, dPxdtdbi)
    d2Pxdbdti[,c(1:(p-1))] <- 0

    dPxdtdbj=list()

    for (i in 1:(p-1))
     {dPxdtdbj[[i]] = matrix(0,m*(p-1),(p-i))
      for (j in 1:(p-i))
       {dPxdtdbj[[i]][(1+m*(i-1)+m*(j-1)):(m+m*(i-1)+m*(j-1)),j] = (-c(1:m)*sin(c(1:m)*(ang[(i+j)]-ang[i])))*z2[(i+j)]*z2[i]}
     }

    d2Pxdbdtj <- do.call(cbind, dPxdtdbj)

    if (mcsc == "unconstrained"){
      d2Pxdbdt = d2Pxdbdti + d2Pxdbdtj
    } else if (mcsc == -1){
      if (m <=2) {d2Pxdbdt = matrix(0,m*(p-1),p*(p-1)/2)
      } else  {
        dPxdtdbi=list()

        for (i in 1:(p-1))
        {dPxdtdbi[[i]] = matrix(0,(p-1)*m,(p-i))
        for (j in ((1+m*(i-2)):(1+m*(i-2)+m-1)))
        {q = j - m*(i-2)
        for ( K in ((1+i):p))
        {dPxdtdbi[[i]][j,(K-i)] = (q*sin(q*(ang[K]-ang[i])))*z2[K]*z2[i]}
        dPxdtdbi[[i]][seq(2,((p-1)*m),by=2),] = 0
        dPxdtdbi[[i]][1,] =0
        } }

        d2Pxdbdti <- do.call(cbind, dPxdtdbi)
        d2Pxdbdti[,c(1:(p-1))] <- 0

        dPxdtdbj=list()

        for (i in 1:(p-1))
        {dPxdtdbj[[i]] = matrix(0,m*(p-1),(p-i))
        for (j in 1:(p-i))
        {dPxdtdbj[[i]][(1+m*(i-1)+m*(j-1)):(m+m*(i-1)+m*(j-1)),j] = (-c(1:m)*sin(c(1:m)*(ang[(i+j)]-ang[i])))*z2[(i+j)]*z2[i]}
        dPxdtdbj[[i]][seq(2,((p-1)*m),by=2),] = 0
        dPxdtdbj[[i]][1,] =0
        }

        d2Pxdbdtj <- do.call(cbind, dPxdtdbj)
        d2Pxdbdt = d2Pxdbdti + d2Pxdbdtj
      }}  else {

        dPxdtdbi=list()

        for (i in 1:(p-1))
        {dPxdtdbi[[i]] = matrix(0,(p-1)*m,(p-i))
        for (j in ((1+m*(i-2)):(1+m*(i-2)+m-1)))
        {q = j - m*(i-2)
        for ( K in ((1+i):p))
        {dPxdtdbi[[i]][j,(K-i)] = (q*sin(q*(ang[K]-ang[i])))*z2[K]*z2[i]}
        dPxdtdbi[[i]][1,] =0
        } }

        d2Pxdbdti <- do.call(cbind, dPxdtdbi)
        d2Pxdbdti[,c(1:(p-1))] <- 0

        dPxdtdbj=list()

        for (i in 1:(p-1))
        {dPxdtdbj[[i]] = matrix(0,m*(p-1),(p-i))
        for (j in 1:(p-i))
        {dPxdtdbj[[i]][(1+m*(i-1)+m*(j-1)):(m+m*(i-1)+m*(j-1)),j] = (-c(1:m)*sin(c(1:m)*(ang[(i+j)]-ang[i])))*z2[(i+j)]*z2[i]}
        dPxdtdbj[[i]][1,] =0
        }

        d2Pxdbdtj <- do.call(cbind, dPxdtdbj)
        d2Pxdbdt = d2Pxdbdti + d2Pxdbdtj
      }

    # dtdn

    angi <- function(ang,i,j,m){
      x = rep(0, m)

      for (K in 1:m)
      {x[K] = K*b[K+1]*sin(K*(ang[j]-ang[i]))}

      sum(x)}

    angj <- function(ang,i,j,m)
    {x = rep(0, m)

    for (K in 1:m)
    {x[K] = -K*b[K+1]*sin(K*(ang[j]-ang[i]))}

    sum(x)}

    dPxdthe <- list()

    for (ii in 1:(p-1)){dPxdthe[[ii]] <- matrix(0, (p-ii), (p))}

    for (ii in 1:(p-1))
    {
      for (l in 1:(p-ii))
      { dPxdthe[[ii]][l,ii] <- angi(ang,ii,(ii+l),m)
      dPxdthe[[ii]][l,(ii+l)] <- angj(ang,ii,(ii+l),m)
      }
    }

    dPxdthe <- do.call(rbind, dPxdthe)

    tt = list()

    for (i in 1:p)
    {
      tt[[i]] = matrix(0,p,p)

      tt[[i]][lower.tri(tt[[i]])] <- dPxdthe[,i]
      tt[[i]] = tt[[i]] + t(tt[[i]])
    }

    dPxdthedv = list()

    for (i in 1:p){
      dPxdthedv[[i]] =  matrix(0, p, p*(p-1)/2)

      for(l in 1:p){

        KK <- matrix(0,p,p)
        diag(KK)[l] <- (-1/2)*(1+v[l])^(-3/2)

        dPxdthedvl = KK%*%tt[[i]]%*%Dz2 + t(KK%*%tt[[i]]%*%Dz2)

        dPxdthedv[[i]][l,] = as.vector(dPxdthedvl[lower.tri(dPxdthedvl)])
      }
    }

    d2Pxdtdv <- do.call(rbind, dPxdthedv)
    d2Pxdtdv <- d2Pxdtdv[-c(1:p),]

    # dbdn
    # beta part

    if (mcsc == "unconstrained") {
      beta = list()
      for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

      for (ii in 1:(p-1))
      {
        for (l in 1:(p-ii))
        {
          bt <- rep(0, m)
          for (mm in 1:(m)) {bt[mm] <- (cos(mm*(ang[ii+l]-ang[ii]))-1)}
          beta[[ii]][l,] <- bt
        }
      }

      dPxdb <- do.call(rbind, beta)
    } else if (mcsc == -1){

      beta = list()
      for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

      for (ii in 1:(p-1))
      {
        for (l in 1:(p-ii))
        {
          bt <- rep(0, m)
          for (mm in 1:(m)) {bt[mm] <- (cos(mm*(ang[ii+l]-ang[ii]))-cos(ang[ii+l]-ang[ii]))}
          beta[[ii]][l,] <- bt
        }
      }

      dPxdb <- do.call(rbind, beta)

      if (m <= 2) {dPxdb[,1:m] <- 0} else {dPxdb <- do.call(rbind, beta)
      dPxdb[,seq(2,m,by=2)] = 0; dPxdb[,1] = 0
      }
    } else {

      beta = list()
      for (ii in 1:(p-1)){beta[[ii]] <- matrix(0, (p-ii), m)}

      for (ii in 1:(p-1))
      {
        for (l in 1:(p-ii))
        {
          bt <- rep(0, m)
          if (m==1) {for (mm in seq(1,m,by=2)) {bt[mm] <- (cos(mm*(ang[ii+l]-ang[ii]))-cos(ang[ii+l]-ang[ii]))}
          }else {
            for (mm in seq(1,m,by=2)) {bt[mm] <- (cos(mm*(ang[ii+l]-ang[ii]))-cos(ang[ii+l]-ang[ii]))}
            for (mm in seq(2,m,by=2)) {bt[mm] <- (cos(mm*(ang[ii+l]-ang[ii]))-1)}}
          beta[[ii]][l,] <- bt}
      }
      dPxdb <- do.call(rbind, beta)
      dPxdb[,1] = 0
    }

    bb = list()

    for (i in 1:m)
    {
      bb[[i]] = matrix(0,p,p)

      bb[[i]][lower.tri(bb[[i]])] <- dPxdb[,i]
      bb[[i]] = bb[[i]] + t(bb[[i]])
    }

    dPxdbdv = list()

    for (i in 1:m){

      dPxdbdv[[i]] =  matrix(0, p, p*(p-1)/2)

      for(l in 1:p){

        KK <- matrix(0,p,p)
        diag(KK)[l] <- (-1/2)*(1+v[l])^(-3/2)

        dPxdbdvl = KK%*%bb[[i]]%*%Dz2 + t(KK%*%bb[[i]]%*%Dz2)

        dPxdbdv[[i]][l,] = as.vector(dPxdbdvl[lower.tri(dPxdbdvl)])
      }
    }

    d2Pxdbdv <- do.call(rbind, dPxdbdv)

    # dndn

    d2Pxdvdv = list()

    for (l in 1:p) {

      d2Pxdvdv[[l]] = matrix(0, (p-l+1), p*(p-1)/2)

      KK1 <- matrix(0,p,p)
      diag(KK1)[l] <- (-1/2)*(1+v[l])^(-3/2)

      JJ1 <- matrix(0,p,p)
      diag(JJ1)[l] <- 1

      for (q in l:p) {
        KK2 <- matrix(0,p,p)
        diag(KK2)[q] <- (-1/2)*(1+v[q])^(-3/2)

        JJ2 <- matrix(0,p,p)
        diag(JJ2)[q] <- 1

        if(q == l)
        {dKdv = matrix(0,p,p)
        diag(dKdv)[l] <- (3/4)*(1+v[l])^(-5/2)
        d2Pxdvdvl = dKdv%*%(Pc+Dv)%*%Dz2 + KK1%*%JJ2%*%Dz2 + KK1%*%(Pc+Dv)%*%KK2 + KK2%*%JJ1%*%Dz2 +
          t(dKdv%*%(Pc+Dv)%*%Dz2 + KK1%*%JJ2%*%Dz2 + KK1%*%(Pc+Dv)%*%KK2 + KK2%*%JJ1%*%Dz2)

        d2Pxdvdv[[l]][(q-l+1),] <-  as.vector(d2Pxdvdvl[lower.tri(d2Pxdvdvl)])} ##change into lower.triangle

        else
        {d2Pxdvdvl = KK1%*%JJ2%*%Dz2 + KK1%*%(Pc+Dv)%*%KK2 + KK2%*%JJ1%*%Dz2 +
          t(KK1%*%JJ2%*%Dz2 + KK1%*%(Pc+Dv)%*%KK2 + KK2%*%JJ1%*%Dz2)

        d2Pxdvdv[[l]][(q-l+1),] <-  as.vector(d2Pxdvdvl[lower.tri(d2Pxdvdvl)])}
      }
    }

    d2Pxdvdv <- do.call(rbind, d2Pxdvdv)

    DD = (t(Delta)%*%Delta)
    AA = matrix(0,length(par),length(par))
    BB = t(Delta)%*%Y.Hat%*%Delta
    hh[[1]] <- NULL

    for (i in 1:(length(par)))
    {for (j in i:(length(par)))
    {if (i <= (p-1) & j <= (p-1)){
      if (i == j) {d2Pxdrdr = dPxdtdt.2[i,]
      } else {
        d2Pxdrdr = rep(0,p*(p-1)/2)
        d2Pxdrdr[(i*p - i*(i+1)/2 + (j-i))] = hh[[i]][(j-i)]}
    }else if (i <= (p-1) & j > (p-1) & j <= (p+m-1))
    {d2Pxdrdr = d2Pxdbdt[(m*(i-1)+(j-(p-1))),] # dtdb
    }else if (i <= (p-1) & j > (p+m-1))
    {d2Pxdrdr = d2Pxdtdv[(p*(i-1)+(j-(p+m-1))),]
    }else if (i > (p-1) & i <= (p+m-1) & j > (p-1) & j <= (p+m-1))
    {d2Pxdrdr = rep(0,p*(p-1)/2)
    }else if (i > (p-1) & i <= (p+m-1) & j > (p+m-1))
    {d2Pxdrdr = d2Pxdbdv[(p*(i-(p-1)-1)+(j-(p+m-1))),]
    }else if (i > (p+m-1) & j > (p+m-1))
    {ii = i - (p+m-1); jj = j - (p+m-1)
    d2Pxdrdr = d2Pxdvdv[((2*p-(ii-2))*(ii-1)/2 + j - i + 1),]
    }

      AA[i,j] = (DD[i,j] + d2Pxdrdr%*%res)}
    }

    AA.t = t(AA)
    diag(AA.t) = 0
    AA = AA + AA.t

    ##se of beta1 ~
    se = sqrt(diag(solve(AA)%*%BB%*%solve(AA)))
    se[1:(p-1)] = se[1:(p-1)]*180/pi

    resultA = list()
    resultA$F0 = F0
    resultA$RMSEA = RMSEA
    resultA$perfect.fit = p.perfect
    resultA$close.fit = p.close
    resultA$s.e = se/sqrt(N)
    resultA$df= df
    return(resultA)
  }

  estim <- est$par

  if (m==1){
    estim[(p-1+1):(p-1+m)] <- 1 - estim[(p-1+1):(p-1+m)]/(1+sum(estim[(p-1+1):(p-1+m)]))
  } else estim[(p-1+1):(p-1+m)] <- c(1, estim[(p-1+2):(p-1+m)])/(1+sum(estim[(p-1+1):(p-1+m)]))

  estim2 = estim
  op = which(estim[1:(p-1)]>2*pi)
  estim[op] = estim[op] - 2*pi
  op2 = which(estim[1:(p-1)]<0)
  estim[op2] = 2*pi + estim[op2]
  estim[1:(p-1)] = estim[1:(p-1)]*180/pi
  com = sqrt(1/(1+estim[(length(estim)-p+1) : length(estim)]))

  v.se = test.stat(est$par)$s.e[(length(estim)-p+1) : length(estim)]
  com.se = 1/2*v.se*(1+estim[(length(estim)-p+1) : length(estim)])^(-3/2)

  result = list()
  result$coefficients = data.frame(point.estimates = c(estim[1:(p+m-1)],com), SE = c(test.stat(est$par)$s.e[1:(p+m-1)],com.se))
  result$v = data.frame(v = estim[(p+m):(length(estim))], SE = v.se)
  result$test.stat = test.stat(est$par)
  result$radians = estim2
  result$optim = est

  if (sum(is.na(result$coef$SE)) != 0) result$SE_NA = 1

  class(result) = "CircumOLS"
  return(result)
}


print.CircumOLS = function(x, ...) print(x[c(1,3)])
