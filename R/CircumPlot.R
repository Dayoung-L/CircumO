#' CircumPlot
#'
#' @export

CircumPlot <-
  function(C.object,m=1,r=1, labels = TRUE, pchar=NULL,bg.points="red",ef=0.4,
           big.points=15,big.labels=15,bg.plot="white",
           col.axis="black",color="black",col.text="white",
           twodim=TRUE,bound=TRUE,label=NULL,reverse=FALSE, title=NULL){

    ## Grassi, M., Luccio, R., & Di Blas, L. (2010). Circe: An r implementation of browne’s circular stochastic process model.
    ## Behavior Research Methods,42(1), 55–73.

    words <- c(rownames(C.object$coefficients))
    starts_with <- "a"
    matching_words <- grep(paste0("^", starts_with), words, value = TRUE)
    p = length(matching_words) + 1

    par = C.object$optim$par
    ang1=par[1:(p-1)]
    ang=append(ang1,0,r-1)

    b1 = C.object$coef[c(p+(m-1)),1]
    b = c(1-sum(b1), b1)

    v = par[(length(par)-p+1):length(par)]
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
    Px = Dz2%*%(Pc+Dv)%*%Dz2

    object = list()
    object$m = m
    object$equal.ang = F
    object$equal.com = F
    object$b = b
    object$v.names = if (is.null(label)) {c(1:p)} else {label}
    object$upper = NULL
    object$lower = NULL
    object$R = Px

    par2 = c(0,C.object$optim$par)
    par2[1:p] = par2[1:p]*180/pi + 90
    object$coeff<-data.frame(c(round(par2,5)))

    k = object$m
    K = pi/180

    equal.ang=object$equal.ang
    equal.com=object$equal.com
    #------------------  Plot of correlation function  ----------------
    b.iter=c(object$b)
    #beta=res$parres$parres$par[(p+1):(p+(k+1))]
    theta=K*c(0:360)
    rho=rep(0,length(theta))
    for(i in 1:length(rho)){

      rho[i]=b.iter[1]+c(rep(1,length(b.iter[-c(1)])))%*%(b.iter[-c(1)]*cos(c(seq(1,k))*theta[i] ))

    }
    #plot(theta,rho,xlim=c(0*K,360*K),ylim=c(0,1),bty="n")
    mincorr180=2*(sum(b.iter[c(seq(1,(k+1),by=2))]))-1;summary(rho);
    mincorr180=round(mincorr180,3)

    grad=c(0,45,90,135,180,225,270,315,360)
    grad*K
    pos=grad*K

    v.names=object$v.names
    #--------------------------- circular plot ----------------------

    # par(mar = c(0,0,0,0), bg = bg.plot)
    if(equal.com==TRUE) com.ind=sqrt(1/(1+object$coeff[(p+m+1):(p+m+p)])) else com.ind=sqrt(1/(1+object$coeff[(p+m+1):(p+m+p),1]))
    plot(c(-1,1),c(-1,1),type="n",bty="n",axes=FALSE, main = title, ylab = " ", xlab = " ");
    x1=seq(0,1,by=0.00001);x2=seq(-1,0,by=0.001);x=c(x2,x1)
    points(c(x,x),c(sqrt(1^2-x^2),c(-1)*sqrt(1^2-x^2)),cex=0.1,col=color,type="l",lwd=2,lty="solid");if(twodim==TRUE){abline(h=0,v=0,col=color,lty="solid",lwd=2)};

    if(bound==TRUE){
      if(!is.null(object$upper)){
        up<-unique(object$upper)
        segments(rep(0,length(up)),rep(0,length(up)),cos(up*K)*1,sin(up*K)*1)
      }
      if(!is.null(object$lower)){
        low<-unique(object$lower)
        segments(rep(0,length(low)),rep(0,length(low)),cos(low*K)*1,sin(low*K)*1)
      }
    }
    angular.points=matrix(0,dim(object$R)[1],2)
    for(i in 1:dim(object$R)[1]){
      if(reverse==FALSE){angular.points[i,]=c(cos(K*object$coeff[i,1]),sin(K*object$coeff[i,1]))}
      if(reverse==TRUE){angular.points[i,]=c(cos(K*(360-object$coeff[i,1])),sin(K*(360-object$coeff[i,1])))}
    }

    row.names(angular.points)=object$v.names
    if(labels==TRUE){
      text(angular.points[,1]*(ef*com.ind),angular.points[,2]*(ef*com.ind),col=color,pch=17,cex=(big.labels)/nrow(object$R),labels=v.names)}
    segments(c(rep(0,dim(object$R)[1])),c(rep(0,dim(object$R)[1])),angular.points[,1]*com.ind,angular.points[,2]*com.ind,col=color,lty="solid")
    if(labels==TRUE){
      segments(angular.points[,1]*com.ind,angular.points[,2]*com.ind,angular.points[,1],angular.points[,2],col=color,lty="dotted")}
    points(angular.points[,1]*com.ind,angular.points[,2]*com.ind,col="black",pch=if(is.null(pchar))2 else pchar,bg=bg.points,cex=(big.points)/nrow(object$R))
    text(-1.09,0.92,labels=substitute(list(rho[180])==list(r),list(r=mincorr180)),col=col.text,pos=4);
    text(-1.09,0.86,labels=substitute(list("max "*h)==list(maxcom),list(maxcom=round(1,2))),col=col.text,pos=4)
    text(-1.09,0.80,labels=substitute(list("max "*h^2)==list(maxcom),list(maxcom=round(max(com.ind^2),2))),col=col.text,pos=4)
  }
