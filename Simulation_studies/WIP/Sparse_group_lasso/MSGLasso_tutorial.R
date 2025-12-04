install.packages('MSGLasso')

library('MSGLasso')
set.seed(sample(1:100,1))
G.arr <- c(0,20,20,20,20,20,20,20,20,20,20)

data("Beta.m")

######## generate data set for model fitting

simDataGen<-function(N, Beta, rho, s, G.arr, seed=1){
  
  P<-nrow(Beta)
  Q<-ncol(Beta)
  gsum<-0
  X.m<-NULL
  
  set.seed(seed)
  
  Sig<-matrix(0,P,P)
  jstart <-1
  
  for(g in 1:length(G.arr)-1){
    X.m<-cbind(X.m, matrix(rnorm(N*G.arr[g+1]),N,G.arr[g+1], byrow=TRUE))
    
    for(i in 2:P){ for(j in jstart: (i-1)){
      
      Sig[i,j]<-rho^(abs(i-j))
      
      Sig[j,i]<-Sig[i,j]
      
    }}
    jstart <- jstart + G.arr[g+1]
  }
  
  
  diag(Sig)<-1
  R<-chol(Sig)
  
  X.m<-X.m%*%R
  
  SVsum <-0
  
  for (q in 1:Q){SVsum <-SVsum+var(X.m %*% Beta[,q])}
  sdr =sqrt(s*SVsum/Q)
  
  E.m <- matrix(rnorm(N*Q,0,sdr),N, Q, byrow=TRUE)
  
  Y.m<-X.m%*%Beta+E.m
  
  return(list(X=X.m, Y=Y.m, E=E.m))
}

N <-150

rho=0.5; 
s=4;

Data <- simDataGen(N, Beta.m,rho, s, G.arr, seed=sample(1:100,1))
X.m<-Data$X
Y.m<-Data$Y

dim(X.m)
dim(Y.m)


############################################################
## fit model for one set of (lam1, lam.G) using example data
############################################################

P <- dim(Beta.m)[1]
Q <- dim(Beta.m)[2]
G <- 10
R <- 10

gmax <- 1
cmax <- 20
GarrStarts <- c(0,20,40,60,80,100,120,140,160,180)
GarrEnds <- c(19,39,59,79,99,119,139,159,179,199)
RarrStarts <- c(0,20,40,60,80,100,120,140,160,180)
RarrEnds <- c(19,39,59,79,99,119,139,159,179,199)

tmp <- FindingPQGrps(P, Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
PQgrps <- tmp$PQgrps

tmp1 <- Cal_grpWTs(P, Q, G, R, gmax, PQgrps)
grpWTs <- tmp1$grpWTs

tmp2 <- FindingGRGrps(P, Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
GRgrps <- tmp2$GRgrps

Pen_L <- matrix(rep(1,P*Q),P,Q, byrow=TRUE)
Pen_G <- matrix(rep(1,G*R),G,R, byrow=TRUE)
grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)

MSGLassolam1 <- 1.6
MSGLassolamG <- 0.26
MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)

system.time(try <-MSGLasso(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQgrps, GRgrps, grp_Norm0, 
                           MSGLassolam1, MSGLassolamG.m, Beta0=NULL))
