treekde<- function(dados,alpha,MinPts){
  dados<-data.frame(dados)
  fr <- function(x,X){
    
    sum(1/(m*h*sqrt(2*pi))*exp(-((x-X)^2)/(2*h^2)))
  } #function kde
  fr2 <- function(x,X){
    -sum(1/(m*h^3*sqrt(2*pi))*(x-X)*exp(-((x-X)^2)/(2*h^2)))
  }#derivative kde
  
  n<-dim(dados)[2]
  m<-dim(dados)[1]
  cluster<-rep(1,m)
  Ctudo<-cbind(dados,cluster)
  parar<-0#number of leaf clusters
  completo<-{}#set of leaf clusters
  k<-1#number of parent clusters
  int<-0#iterations
  
  while (parar<k) {
    rep<-k
    for (j in 1:rep) {
      t<-completo==j
      w<-any(t)
      corte<-0
      variavel<-0
      fmin<-10^8
      
      if(w==FALSE){
        C<-Ctudo[Ctudo[,n+1]==j,-(n+1)]
        if(is.null(dim(C)[1])==FALSE){
          m<-dim(C)[1]
          for (i in 1:n) {
            X<-C[,i]
            sigma<-sd(X)
            if(is.na(sigma)==FALSE){
              h<-alpha*((4/3)^(1/5)*sigma*length(X)^(-1/5))
              l<-min(X)#lower bound
             u<-max(X)#upper bound
            r<-optimise(fr,c(l,u),tol = 10^-6,X)#optimization
            verifica<-abs(fr2(r$minimum,X))#second order derivative
              if(is.na(verifica)==TRUE){
               verifica<-10^4
              }
              
            }else{verifica<-10^4}
            if(verifica<10^-6){
              if(r$objective<fmin){
                corte<-r$minimum
                variavel<-i
                fmin<-r$objective
                grad_corte<-verifica
                lm<-l
                um<-u
              }
            }
          }
          Cauxiliar<-Ctudo
          
          if(variavel>0){
            Ctudo[Ctudo[,variavel]>corte&Ctudo[,n+1]==j,n+1]<-k+1
            Cj<-Ctudo[Ctudo[,n+1]==j,n+1]
            Cj1<-Ctudo[Ctudo[,n+1]==k+1,n+1]
            if(min(length(Cj),length(Cj1))<MinPts){
              Ctudo<-Cauxiliar
              completo<-cbind(completo,j)
            }else{k<-k+1
            }
          }else{completo<-cbind(completo,j)
          }
          parar<-length(completo)
        }
      }
    }
    int<-int+1
  }
  return(list(cluster=Ctudo[,n+1]))
}

