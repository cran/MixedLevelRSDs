
#################################
#' Measure of rotatability Q based on a second order model
#'
#' @param design Second order design matrix without intercept
#'
#' @return Rotatability measure Q
#' @description Calculates the measure of rotatability (measure Q, 0 <= Q <= 1) given by Draper and Pukelsheim(1990) for given design based on a second order model.
#' @export
#'
#' @examples
#' \dontrun{
#' library(MixedLevelRSDs)
#'RotatabilityQ(design)
#'}
#'@references
#'1)	Norman R. Draper and Friedrich Pukelsheim(1990), <doi: 10.1080/00401706.1990.10484635>. "Another look at rotatability".
#'
#'2)	M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<doi: 10.1080/03610926.2021.1944213>." Theoretical developments in response surface designs:    an informative review and further thoughts.".


RotatabilityQ<-function(design){
  #####taking all the combination
  m=as.matrix(design)
  v=ncol(m)
  mm=m
  for(i in 1:ncol(m)){
    for(j in 1:ncol(m)){
      x=mm[,i]*mm[,j]
      mm=cbind(mm,x)
    }
  }
  colnames(mm)=NULL
  mm=cbind(1,mm)
  x_mat=mm
  ######x prime x
  x_prime_x=t(x_mat)%*% x_mat
  ###################A matrix
  A=(x_prime_x)/nrow(x_mat)
  ############ v0 matrix
  v_0=matrix(0,nrow=nrow(A),ncol=ncol(A))
  v_0[1,1]=1
  v_0
  ########v2 matrix
  totalcol=1+v+(v*v)
  sq=c()
  for(i in 1:(v)){
    if(i==1){
      sq=c(sq,(v+2))
    }
    if(i>1){
      x=(sq[length(sq)]+(v+1))
      sq=c(sq,x)
    }
  }
  v_2=matrix(0,nrow=nrow(A),ncol=ncol(A))
  v_2[1,c(sq)]<-(3*v)^(-0.5)
  v_2[c(sq),1]<-(3*v)^(-0.5)
  #######
  for(j in 2:(v+1)){
    v_2[j,j]=(3*v)^(-0.5)
  }
  v_2
  ##################v4 matrix
  v_4=matrix(0,nrow=nrow(A),ncol=ncol(A))
  ########################
  totalcol=1+v+(v*v)
  sqtrm=c()
  for(i in 1:(v)){
    if(i==1){
      sqtrm=c(sqtrm,(v+2))
    }
    if(i>1){
      x=(sqtrm[length(sqtrm)]+(v+1))
      sqtrm=c(sqtrm,x)
    }
  }
  intterm=c(setdiff(seq(min(sqtrm),max(sqtrm)),sqtrm))
  #####swap
  swap<-function(a,b){
    x<-c(b,a)
    return(x)
  }
  ##########for square terms
  matint=matrix(min(sqtrm):max(sqtrm),nrow=v,byrow=T)
  matintdia=c(diag(matint))
  for(i in matintdia){
    for(j in matintdia){
      if(i==j){
        v_4[i,j]<-3*((3*v*(v+2))^(-0.5))
      }else{
        v_4[i,j]<-((3*v*(v+2))^(-0.5))
      }
    }
  }
  ################for interaction terms
  matint=matrix(min(sqtrm):max(sqtrm),nrow=v,byrow=T)
  same=c()
  for(k in intterm){
    pos=c(which(matint==k,arr.ind = T))
    pos1=c(swap(pos[1],pos[2]))
    x1=matint[pos[1],pos[2]]
    x2=matint[pos1[1],pos1[2]]
    same=c(x1,x2)
    for(i in same){
      for(j in same){
        v_4[i,j]<-((3*v*(v+2))^(-0.5))
      }
    }
  }
  v_4

  #################A bar
  A_bar=v_0+v_2*(sum(diag(A%*%v_2)))+v_4*(sum(diag(A%*%v_4)))
  ##########
  A_bar_v_0_sq=(A_bar-v_0)%*%(A_bar-v_0)
  A_v_0_sq=(A-v_0)%*%(A-v_0)
  #########Q star
  Q=(sum(diag(A_bar_v_0_sq))/sum(diag(A_v_0_sq)))
  message("Rotatability Measure Q*")
  return(round(Q,digits = 5))
}
