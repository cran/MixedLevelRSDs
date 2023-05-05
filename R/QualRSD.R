
#' Second Order RSDs with qualitative factors
#'
#' @param v Total number of input factors
#' @param k Number of qualitative factors, 1<=k<=v-2
#' @param Interaction To specify whether to generate a design for fitting second order model which include Interaction term between qualitative and quantitative factors. The interaction = T means the generated designs will be suitable to fit a second order model which include interaction between qualitative factors and the linear terms of quantitative factors
#'
#' @return Second Order RSDs with qualitative factors along with D-efficiency and G-efficiency
#' @description Generate a Second Order Design where first v-k column represent v-k quantitative factors and the last k column represents the k qualitative factors (1<=k<=v-2). It also gives D-efficiency and G-efficiency for the generated design.
#' @export
#'
#' @examples
#'library(MixedLevelRSDs)
#'QualRSD(5,2, Interaction = FALSE )
#'@references
#'  1)	M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<doi: 10.1080/03610926.2021.1944213>."Theoretical developments in response surface designs: an informative review and further thoughts".
#'
#'  2)	Jyoti Divecha and Bharat Tarapara (2017). < doi:10.1080/08982112.2016.1217338>. "Small, balanced, efficient, optimal, and near rotatable response surface designs for factorial experiments asymmetrical in some quantitative, qualitative factors".

QualRSD=function(v,k,Interaction=FALSE){
  ql=k
  qt=v-ql
  if(v==3){
    fmatrix=suppressMessages(FrF2::FrF2(2^(v),v))
  }
  if(v==4){
    fmatrix=suppressMessages(FrF2::FrF2(2^(v),v))
  }
  if(v==5){
    fmatrix=suppressMessages(FrF2::FrF2(2^(v),v))
  }
  if(v==6){
    fmatrix=FrF2::FrF2(2^(v-1),v)
  }
  if(v==7){
    fmatrix=FrF2::FrF2(2^(v-1),v)
  }
  if(v==8){
    fmatrix=FrF2::FrF2(2^(v-2),v)
  }
  if(v==9){

    fmatrix=FrF2::FrF2(2^(v-2),v)
  }
  if(v==10){
    fmatrix=FrF2::FrF2(2^(v-3),v)
  }
  row=nrow(fmatrix)
  col=ncol(fmatrix)
  fmatrix=as.matrix(fmatrix)
  fmatrix=c(fmatrix)
  fmatrix=as.integer(fmatrix)
  fmatrix=matrix((fmatrix),ncol=col,nrow=row)
  qtmat=fmatrix[,1:qt]
  qlmat=fmatrix[,(qt+1):v]
  a=2
  axial=qt
  axial1<-diag(a,nrow=axial)
  axial2<-diag(-a,nrow=axial)
  axmat=NULL
  for(i in 1:(qt)){
    axmat=rbind(axmat,axial1[i,],axial2[i,])
  }
  #################
  #axmat=cbind(axmat,matrix(0,nrow=2*(qt-ql),ncol=ncol(qtmat)-ncol(axmat)))
  qtmat=rbind(as.matrix(qtmat),axmat,matrix(0,nrow=4,ncol=ncol(qtmat)))
  ###########
  qtmatsq=qtmat
  for(l in 1:qt){
    qtmatsq=cbind(qtmatsq,qtmatsq[,l]^2)
  }
  qtfulmat<-NULL
  for(i in 1:(qt-1)){
    for(j in (i+1):qt){
      qtfulmat=cbind(qtfulmat,qtmat[,i]*qtmat[,j])
    }
  }
  qtfulmat=cbind(qtmatsq,qtfulmat)
  ###############
  addzero=matrix(0,nrow=(nrow(qtmat)-nrow(fmatrix)),ncol=ql)
  qlmat=rbind(t(t(qlmat)),addzero)
  ############

  qltfulmat=NULL
  for(i in 1:ql){
    qltfulmat=cbind(qltfulmat,qlmat[,i])
  }
  if(Interaction==T){
    for(k in 1:ql){
      for(j in 1:qt){
        qltfulmat<-cbind(qltfulmat,qtmat[,j]*qlmat[,k])
      }
    }
  }
  ########################################
  x_mat=cbind(1,qtfulmat,qltfulmat)

  ############################################
  x_matrix=x_mat
  x_prime_x=t(x_mat)%*%x_mat
  k1=1
  var<-c()
  while(k1<=nrow(x_matrix))
  {
    V=t(x_matrix[k1,])
    b<-t(V)
    v_y_hat<-V %*%MASS::ginv(x_prime_x) %*% b
    var<-c(var,v_y_hat)
    k1<-k1+1
  }
  variance_of_esitmated_response<-var

  ##############
  ########


  #############
  maxpredictedvar<- max(variance_of_esitmated_response)
  #############
  N=nrow(x_mat)
  #######No of parameters
  p=ncol(x_mat)
  #return(x_mat)
  ###########
  message("Second Order Design")
  cat("\n")
  fmat=cbind(qtmat,qlmat)
  #fmat=round(fmat,digits = 5)
  fmat=format(fmat,nsmall=5)
  fmat=noquote(fmat)
  print(fmat)
  cat("\n")
  ################G efficiency
  message("G-efficiency")
  cat("\n")
  G_eff=p/(maxpredictedvar*N)
  G_eff=round(G_eff,digits = 5)
  print(G_eff)
  cat("\n")
  ###########D efficiency
  message("D-efficiency")
  cat("\n")
  D_eff=((det(x_prime_x))^(1/p))/N
  D_eff=round(D_eff,digits = 5)
  print(D_eff)
}


