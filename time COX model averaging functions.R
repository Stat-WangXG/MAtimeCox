MAtdcCOX<-function(data,K.set=c(5:10),test.plot=F,compute.S=T)
{
  #remove the missing data
  data=na.omit(data)

  #retain the subjects whose survival times are greater than zero
  data=data[data[,ncol(data)-1]>0,]

  #deal with tied events
  set.seed(-10)
  t_ties=as.numeric(names(table(data[,ncol(data)-1]))[table(data[,ncol(data)-1])>1])
  if (length(t_ties)>0)
  {
    row_ties=which(data[,ncol(data)-1] %in% t_ties)
    time_after=data[,ncol(data)-1]
    while (1)
    {
      time_after[row_ties]=data[row_ties,ncol(data)-1]+data[row_ties,ncol(data)-1]*runif(length(row_ties),0.0001,0.001)
      model_c_ties=survfit(Surv(time_after,1-delta)~1,data = data.frame(cbind(time_after,delta=data[,ncol(data)] )))
      if ( length(model_c_ties$surv) == length(unique(time_after)) )
        break
    }
  }
  data[,ncol(data)-1]=time_after

  n=nrow(data) # sample size
  p=ncol(data)-2 # number of covariates
  submodel_all=as.list(rep(0,p)) #candidate models
  names(submodel_all)=paste("submodel",1:p,sep="")
  z=as.matrix(data[,1:p])
  time=data[,ncol(data)-1] # observed time
  delta=data[,ncol(data)] # censoring indicator
  censoring_rate=1-mean(delta)
  delta=delta[order(time)]
  z=z[order(time),]
  time=time[order(time)]
  data_original=data.frame(cbind(z,time,delta))

  vfit = coxph(Surv(time,delta) ~.,data=data_original,iter.max = 40)
  if ( test.plot==T )
  {
    plot(cox.zph(vfit),col = "red")
  }

  model_c=survfit(Surv(time,1-delta)~1,data = data_original)

  #determine the number of basis functions based on AIC
  K_set=unique(sort(K.set))
  AIC=1e+100
  K=K_set[1]
  for ( i in 1:length(K_set) )
  {
    Bspline=bs(time,df=K_set[i],knots =quantile(time[delta==1],probs = seq(1:(K_set[i]-4))/(K_set[i]-3)),intercept=T,degree = 3)
    data_new<-data.expand(delta, time, z, Bspline,ncol(Bspline))
    data_new<-data.frame(data_new)
    index=which(data_new$time_start==data_new$time_stop)
    if (length(index)!=0)
      data_new<-data_new[-index,]
    for (m in 1:p)
    {
      formula_m=paste("Surv(time_start,time_stop,delta_new)~",
                      paste(colnames(data_new)[1:p][-m],sep="",collapse = "+"),
                      "+",paste(rep(colnames(data_new)[m],ncol(Bspline)),":",
                                colnames(data_new)[(p+1):(p+ncol(Bspline))],sep="",collapse = "+"),sep = "")
      formula_m=as.formula(formula_m)
      r_temp=suppressMessages(try(coxph(formula_m,data = data_new,control = coxph.control(timefix = F),singular.ok=F,iter.max = 40),silent = T))
      if (class(r_temp)=="try-error")
        break
      else
      {
        submodel_m=coxph(formula_m,data = data_new,control = coxph.control(timefix = F),singular.ok=F,iter.max = 40)
        submodel_all[[m]]=submodel_m
      }
    }
    if (class(r_temp)=="try-error")
    {
      break
    }
    AIC_temp=sum(vapply(1:p,function(x) extractAIC(submodel_all[[x]])[2],as.numeric(1)))
    if ( AIC_temp < AIC )
    {
      K=K_set[i]
      AIC=AIC_temp
    }
  }

  Bspline=bs(time,df=K,knots =quantile(time[delta==1],probs = seq(1:(K-4))/(K-3)),intercept=T,degree = 3)
  data_new=data.expand(delta, time, z, Bspline,ncol(Bspline))
  data_new=data.frame(data_new)
  index=which(data_new$time_start==data_new$time_stop)
  if (length(index)!=0)
    data_new<-data_new[-index,]
  for (m in 1:p)
  {
    formula_m=paste("Surv(time_start,time_stop,delta_new)~",
                    paste(colnames(data_new)[1:p][-m],sep="",collapse = "+"),
                    "+",paste(rep(colnames(data_new)[m],ncol(Bspline)),":",
                              colnames(data_new)[(p+1):(p+ncol(Bspline))],sep="",collapse = "+"),sep = "")
    formula_m=as.formula(formula_m)
    r_temp=suppressWarnings(try(coxph(formula_m,data = data_new,control = coxph.control(timefix = F),singular.ok=F,iter.max = 40),silent = T))
    if (class(r_temp)=="try-error")
    {
      break
    }
    else
    {
      submodel_m=coxph(formula_m,data = data_new,control = coxph.control(timefix = F),singular.ok=F,iter.max = 40)
      submodel_all[[m]]=submodel_m
    }
  }
  if (class(r_temp)=="try-error")
  {
    cat("Submodels can not convergence!","\n")
    return(NULL)
  }

  G=matrix(rep(model_c$surv,n),byrow = F,ncol = n)
  if (n-nrow(G)!=0)
  {
    G=rbind(G,matrix(rep(G[nrow(G),],n-nrow(G)),nrow = n-nrow(G),byrow = T))
  }
  Dmat=0
  dvec=0
  term_C=0
  for (i in 1:n)
  {
    S_i=vapply(1:p
               ,FUN = function(x)   Breslow.S.m(n,delta,z,Bspline,x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K-1)],z[i,])
               ,FUN.VALUE = as.double(1:n))
    W_i=diag(vapply(1:n,FUN =function(x) w.t.i.G(time[x],i,time,delta,G),FUN.VALUE = as.numeric(1)))
    I_i=time<time[i]
    Dmat=Dmat+2*t(S_i)%*%W_i%*%S_i
    dvec=dvec+2*t(I_i)%*%W_i%*%S_i
    term_C=term_C+t(I_i)%*%W_i%*%I_i
  }
  Amat=as.matrix(rbind(rep(1,p),diag(p)))
  bvec=c(1,rep(0,p))
  r_temp=suppressWarnings(try(solve.QP(Dmat,dvec,t(Amat),bvec=bvec,meq=1),silent = T))
  if (class(r_temp)=="try-error")
  {
    cat("There is no solution to the quadratic programming problem!","\n")
    return(NULL)
  }
  answer=solve.QP(Dmat,dvec,t(Amat),bvec=bvec,meq=1)
  w_best=answer$solution # the weight vector of candidate models
  S_MA_train=NULL # each column denotes a MA estimator of the conditional survival function
  if (compute.S==T)
  {
    S_MA_train=vapply(1:n,function(i) vapply(1:p,FUN = function(x)
      Breslow.S.m(n,delta,z,Bspline,x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K-1)],z[i,]),FUN.VALUE = as.double(1:n))%*%w_best,as.double(1:n))
  }
  return( list(data_train=data,z.order=z,Bspline=Bspline,time.order=time,delta.order=delta,n=n,p=p,K_n=K,test=cox.zph(vfit),candidate_models=submodel_all,MA_weights=w_best,S_MA_train=S_MA_train) )

}





MApredict<-function(MAtdcCOX.object,newdata,t_star) #newdata can be a p-dimensional vector, n*p dimensional matrix or data frame
{
  if (sum(is.na(newdata))>0)
  {
    cat("Warning:There is NA in newdata!","\n")
    return(NULL)
  }

  t_star=t_star[order(t_star)]
  time=MAtdcCOX.object$time.order
  location_t_star=vapply(1:length(t_star),function(x) sum(t_star[x]>=time),FUN.VALUE = as.numeric(1) )
  n=MAtdcCOX.object$n
  p=MAtdcCOX.object$p
  delta=MAtdcCOX.object$delta.order
  z=MAtdcCOX.object$z.order
  Bspline=MAtdcCOX.object$Bspline
  K=MAtdcCOX.object$K_n
  submodel_all=MAtdcCOX.object$candidate_models
  w_best=MAtdcCOX.object$MA_weights

  if ( is.vector(newdata) & !is.list(newdata) )
  {
    S_MA_test=vapply(1:p,FUN = function(x)
      Breslow.S.m(n,delta,z,Bspline,x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K-1)],newdata),FUN.VALUE = as.double(1:n))%*%w_best
    if (max(location_t_star)==0)
    {
      S_MA_t_star=rep(1,length(t_star))
    }
    else
    {
      S_MA_t_star=c(rep(1,sum(location_t_star==0)),S_MA_test[location_t_star[ (sum(location_t_star==0)+1):length(t_star) ]])
    }
    S_MA_t_star=as.matrix(S_MA_t_star)
    rownames(S_MA_t_star)=t_star

    return( list(n_test=1,t_star.order=t_star,S_MA_t_star=S_MA_t_star,time.order=time,delta.order=delta,S_MA_test=as.vector(S_MA_test)) )
  }

  else if ( is.matrix(newdata) | is.data.frame(newdata)  )
  {
    newdata=as.matrix(newdata)
    n_test=nrow(newdata)
    S_MA_test=vapply(1:n_test,function(i) vapply(1:p,FUN = function(x)
      Breslow.S.m(n,delta,z,Bspline,x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K-1)],newdata[i,]),FUN.VALUE = as.double(1:n))%*%w_best,as.double(1:n))
    if(max(location_t_star)==0)
    {
      S_MA_t_star=matrix(1,nrow = length(t_star),ncol=n_test)
    }
    else
    {
      S_MA_t_star=rbind( matrix(1,nrow=sum(location_t_star==0),ncol=n_test) , S_MA_test[location_t_star[ (sum(location_t_star==0)+1):length(t_star) ],] )
    }
    rownames(S_MA_t_star)=t_star

    return( list(n_test=n_test,t_star.order=t_star,S_MA_t_star=S_MA_t_star,time.order=time,delta.order=delta,S_MA_test=S_MA_test) )
  }

  else
  {
    cat("Warning:The newdata is not a vector or matrix or data frame!","\n")
    return(NULL)
  }
}





Breslow.S.m<-function(n,delta,z, B_spline,m,beta_m,theta_m,covariate) # the Breslow estimator of S_m(t|covariate)
{
  Surv.m.covariate=rep(0,n)
  if (ncol(z)==2)
  {
    hazard0<-vapply(c(1:n),FUN =function(x) sum( exp( z[x:n,-m]*beta_m + z[x:n,m]*as.numeric(theta_m%*%B_spline[x,]) )),FUN.VALUE =numeric(1) )
  }
  else
  {
    hazard0<-vapply(c(1:n),FUN =function(x) sum( exp( z[x:n,-m]%*%beta_m + z[x:n,m]*as.numeric(theta_m%*%B_spline[x,]) )),FUN.VALUE =numeric(1) )
  }
  hazard0<-delta/hazard0
  hazard0[delta==0]=0
  hazard0[hazard0==Inf]=1e+300
  Surv.m.covariate<-exp(-cumsum( as.vector(hazard0) *  as.vector( exp( as.numeric(covariate[-m]%*%beta_m) + covariate[m]*B_spline[,]%*%theta_m ) ) ))
  Surv.m.covariate
}





w.t.i.G<-function(t,i,time,delta,G) #Inverse of Probability of Censoring Weighting
{
  if ( time[i]<=min(time) & t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(1)+(time[i]>t)/(1)
  if ( time[i]<=min(time) & !t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(1)+ifelse(time[i]<=t,0,(time[i]>t)/(G[rev(which(t>=time))[1],i]))
  if ( !time[i]<=min(time) & t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(G[rev(which(time[i]>time))[1],i])+(time[i]>t)/(1)
  if ( !time[i]<=min(time) & !t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(G[rev(which(time[i]>time))[1],i])+ifelse(time[i]<=t,0,(time[i]>t)/(G[rev(which(t>=time))[1],i]))
  return(weight.t.i)
}





data.expand<- function(delta2, time2, z2, bs7_2,K) {

  n=length(time2)
  delta_new = numeric(n*(n+1)/2)
  time_start = numeric(n*(n+1)/2)
  time_stop = numeric(n*(n+1)/2)
  z_new = matrix(0, n*(n+1)/2,ncol(z2))
  bs_new = matrix(0, n*(n+1)/2,K)
  index = numeric(n*(n+1)/2)
  key=0
  i=0
  repeat {
    i=i+1

    j=0
    repeat{
      j=j+1

      delta_new[key+j] = delta2[i]*(i==j)
      if (j==1) {
        time_start[key+j]=0
      }

      if (j>1) {
        time_start[key+j]=time2[j-1]
      }

      time_stop[key+j]=time2[j]
      z_new[key+j,]=z2[i,]
      bs_new[key+j,]=bs7_2[j,]
      index[key+j]=j
      if(j==i) break
    }

    key=key+j
    if(i==n) break
  }
  list(z_new=z_new, bs_new=bs_new, delta_new=delta_new, time_start=time_start, time_stop=time_stop, index=index)
}

