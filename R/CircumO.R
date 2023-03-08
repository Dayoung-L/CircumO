#' CircumO
#'
#' @importFrom EFAutilities efa
#' @importFrom MASS mvrnorm
#' @importFrom Turbofuns PolychoricRM
#' @export CircumO
#' @method print CircumO
#' @export

CircumO <- function(rawdt, m=1, mcsc="unconstrained", type="N",
                     simulation_based = F, icorrection = 0, ncore = 2,
                     maxit = 1000, factr = 1e9, pgtol = 0, lmm = NULL,
                     N_star = 10000, maxlag = 20, SE = 1){

  if (!type %in% c("N", "O", "TS")) stop("Specify the type of your data in the 'type' argument. N: Normal, O: Ordinal, TS: Time Series")

  SE = SE
  maxlag = maxlag
  m = m
  N = nrow(rawdt)
  p = ncol(rawdt)

  library(Turbofuns)
  library(EFAutilities)
  library(MASS)

  if (type == "O"){
    Lpoly = PolychoricRM(iRaw = rawdt,IAdjust = 0, NCore = ncore, estimate.acm=T)
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
    if (type == "O"){
      if (sum(eigen(Lpoly[[4]])$values < 0) == 0){
        lambda = efa(covmat = dt, factors = 3, n.obs = N, mtest = F)$unrotated
      } else {lambda = efa(rawdt, dist = "ordinal", factors = 3, n.obs = N, mtest = F)$unrotated}
    } else {lambda = efa(covmat = dt, factors = 3, n.obs = N, mtest = F)$unrotated}


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

  test.stat = function(par, maxlag = maxlag)
  {
    maxlag = maxlag
    R = dt
    p = ncol(dt)
    ang1=par[1:(p-1)]
    ang=append(ang1,0,r-1)
    b = par[(p-1+1):(p-1+m)]/(1+sum(par[(p-1+1):(p-1+m)]))
    if (m==1) {b = c(b[1],1-sum(b[1]))
    } else {b = c(b[1],1-sum(b),b[-1])}

    v=par[(length(par)-p+1):length(par)]

    z2 = rep(0,p)
    for (i in 1:p) {z2[i] = (1/(1+v[i]))^(1/2)}

    M=matrix(c(0),p,p,byrow=TRUE)

    for(i in 1:p){
      for(j in 1:p){
        M[i,j]=c(b[-1])%*%cos(c(1:m)*(ang[j]-ang[i]))
      }
    }

    Pc=M+matrix(b[1],p,p)
    Dv=diag(v,p)
    Dz2=diag(z2)

    if (type == "O"){
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

    } else if (type == "TS") {

      p2 = p * p
      p.star = p*(p-1)/2

      TSACM = function(X,nlagin=maxlag) {
        find.r <- function(dumlag, dumrow, dumcol) {

          if (dumlag < 0) {
            truelag = dumlag * (-1) + 1
            dumtemp = dumrow
            dumrow = dumcol
            dumcol = dumtemp
          } else {
            truelag = dumlag + 1
          }

          Gamma[truelag,dumrow,dumcol]

        } # find.r

        # declare variables

        np = ncol(X)
        nlagout = 0
        nlagout1 = nlagout + 1

        # Calculate lagged correlations
        Gamma = acf(X,lag.max=nlagin,plot = FALSE)$acf

        CovCorr = array( array(0, (nlagout1*np*np)**2), dim=c(nlagout1,np,np,nlagout1,np,np))

        # Declare the logical array LogCovCorr and set its initial values
        LogCovCorr <- CovCorr > 5
        for (i in 1:np) {
          LogCovCorr[1,i,i,,,] = TRUE
          LogCovCorr[,,,1,i,i] = TRUE
        }

        nsumterm = nlagin - nlagout

        m = 0; n = 0; m1 = 1; n1 = 1

        find.temp = array(NA,dim=c((2*nlagin+1),np,np))

        for (a in -nlagin:nlagin){
          for (b in 1:np){
            for (d in 1:np){
              find.temp[a+(nlagin+1),b,d] = find.r(a,b,d)
            }
          }
        }

        for (j in 1:np) {
          for (i in 1:np) {
            for (l in 1:np) {
              for (k in 1:np) {

                if (! LogCovCorr[m1,i,j,n1,k,l]) {

                  for (u.actual in -nsumterm:nsumterm) {
                    r.m.ij = Gamma[m1,i,j]
                    r.n.kl = Gamma[n1,k,l]

                    u = u.actual + (nlagin+1)

                    CovCorr[m1,i,j,n1,k,l] = CovCorr[m1,i,j,n1,k,l] +
                      0.5 * r.m.ij * r.n.kl * ((find.temp[u,i,k])**2 + (find.temp[u,j,k])**2 + (find.temp[u,i,l])**2 + (find.temp[u,j,l])**2) - # terms 1, 2, 3, 4
                      r.n.kl * (find.temp[u,j,k] * find.temp[u+m,i,k]+find.temp[u,j,l]*find.temp[u+m,i,l]) - # terms 5 and 6
                      r.m.ij * (find.temp[u,i,l] * find.temp[u-n,i,k]+find.temp[u,j,l] * find.temp[u-n,j,k]) + # terms 7 and 8
                      find.temp[u,j,l] * find.temp[u-n+m,i,k] + find.temp[u-n,j,k] * find.temp[u+m,i,l]
                  } # u

                  CovCorr[n1,k,l,m1,i,j] = CovCorr[m1,i,j,n1,k,l]
                  LogCovCorr[m1,i,j,n1,k,l] = TRUE
                  LogCovCorr[n1,k,l,m1,i,j] = TRUE

                } # (! LogCovCorr[m,i,j,n,k,l])

              } # k
            } # l
          } # i
        } # j

        temp.arr = array(CovCorr[1,,,1,,])
        asc = array(temp.arr,dim=c(np*np, np*np))

        list(corr = Gamma[1,,], asc = asc)
      }

      u.r = TSACM(rawdt)$asc

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

      Y.Hat = t(M.Select) %*% u.r %*% M.Select
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

      Y.Hat = t(M.Select) %*% u.r %*% M.Select
    }

    # Delta

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

    if (SE == 2){

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
    } else {
      cov.est = solve( t(Delta) %*% Delta ) %*% ( t(Delta) %*% Y.Hat %*% Delta) %*% solve( t(Delta) %*% Delta )

      ## se of beta1 ~
      se = sqrt(diag(cov.est))
      se[1:(p-1)] = se[1:(p-1)]*180/pi
    }

    resultA = list(test.stat = ts.st, RMSEA = RMSEA, perferct.fit = p.perfect,
                   close.fit = p.close, df = df, s.e = se/sqrt(N))

    return(resultA)
  }

  estim <- est$par

  # alpha to beta
  if (m==1){
    estim[(p-1+1):(p-1+m)] <- 1 - estim[(p-1+1):(p-1+m)]/(1+sum(estim[(p-1+1):(p-1+m)]))
  } else estim[(p-1+1):(p-1+m)] <- c(1, estim[(p-1+2):(p-1+m)])/(1+sum(estim[(p-1+1):(p-1+m)]))

  op = which(estim[1:(p-1)]>2*pi)
  estim[op] = estim[op] - 2*pi

  op2 = which(estim[1:(p-1)]<0)
  estim[op2] = 2*pi + estim[op2]
  estim[1:(p-1)] = estim[1:(p-1)]*180/pi

  resultT = test.stat(est$par, maxlag)

  # v: unique variance
  v.se = resultT$s.e[(length(estim)-p+1) : length(estim)]

  # com: communality index
  com = sqrt(1/(1+estim[(length(estim)-p+1) : length(estim)]))
  com.se = 1/2*v.se*(1+estim[(length(estim)-p+1) : length(estim)])^(-3/2)

  coef = rbind(data.frame(Estimates = c(estim[1:(p+m-1)],com), SE = c(resultT$s.e[1:(p+m-1)],com.se)),
               data.frame(Estimates = estim[(p+m):(length(estim))], SE = v.se))
  rownames(coef) = c(paste("a", 2:p, sep = ""), paste("b", 1:m, sep = ""), paste("z", 1:p, sep = ""), paste("u", 1:p, sep = ""))

  if (any(is.na(coef$SE))) warning("NA in the standard error estimates")

  result = list(coef = coef, test = resultT[1:5], optim = est)

  class(result) = "CircumO"
  return(result)
}


print.CircumO = function(x, ...) print(list(coef = round(x[[1]],2), test = x[[2]][1:5]))
