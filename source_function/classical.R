#####################################################
#         PCA, SIR, SAVE, DR    4/4/10              #
#####################################################
#####################################################
#                 Subroutine 1:                     #
#           Dimension Reduction Function            #
#         PCA,  SIR, SAVE, DR 6/8/09                #
#  yunit:      the name of the classes, eg c(0,6,9) #
#  method: 0 = pca, 1 = sir, 2 = save, 3 = dr       #  
#  This function gives the full dimension reduction # 
#  matrix                                           #
#  NB: xc must be centered                          #
#####################################################
reduction = function(xc,y,yunit,d,sig,signrt,method){
if(method==0) return(eigen(var(xc))$vectors[,1:d]) else
if(method==1) return(sirmat(xc,y,yunit,sig,signrt)[,1:d]) else
if(method==2) return(savmat(xc,y,yunit,sig,signrt)[,1:d]) else
return(drmat(xc,y,yunit,sig,signrt)[,1:d])
}
#####################################################
#          Subroutine 2:                            #
#       sliced inverse regression                   #
#  sig = variance of xc; signrt = sig^(-1/2)        #
#  xc must be centered                              #
#####################################################
sirmat = function(xc,y,yunit,sig,signrt){
xstand = xc%*%signrt
xgy = slav(xstand,y,yunit)
prob = slprob(y,yunit)
p = ncol(xc)
out = matrix(0,p,p)
nslice = length(yunit)
for(i in 1:nslice){
out = out + prob[i]*xgy[i,]%*%t(xgy[i,])}
out = eigen(out)$vectors 
out = signrt%*%out
return(out)
}
#####################################################
#            Subroutine 3:                          #
#            slice average                          #
#####################################################
slav = function(x,y,yunit){
n = nrow(x)
p = ncol(x)
nslice = length(yunit)
xgy = matrix(0,nslice,p)
for(i in 1:nslice){
xgy[i,] = apply(x[y==yunit[i],],2,mean)}
return(xgy)
}
#####################################################
#          Subroutine 4:                            #
#       sliced average variance estimate            #
#  sig = variance of xc; signrt = sig^(-1/2)        #
#  xc must be centered                              #
#####################################################
savmat = function(xc,y,yunit,sig,signrt){
xstand = xc%*%signrt
vxgy = slco(xstand,xstand,y,yunit)
prob = slprob(y,yunit)
p = ncol(xc)
out = matrix(0,p,p)
nslice = length(yunit)
for(i in 1:nslice){
out = out + prob[i]*(vxgy[,,i]-diag(p))%*%
                    (vxgy[,,i]-diag(p))}
out = eigen(out)$vectors 
out = signrt%*%out
return(out)
}
#####################################################
#           Subroutine 5:                           #
#           slice covariance                        #
#####################################################
slco = function(x1,x2,y,yunit){
n = nrow(x1)
p = ncol(x1)
nslice = length(yunit)
cx1x2y = array(0,c(p,p,nslice))
for(i in 1:nslice){
cx1x2y[,,i] = cov(x1[y==yunit[i],],x2[y==yunit[i],])}
return(cx1x2y)
}
#####################################################
#      Subroutine 6:                                #
#     directional regression                        #
#     requirements same as before                   #
#####################################################
drmat = function(xc,y,yunit,sig,signrt){
xstand = xc%*%signrt
exy = slav(xstand,y,yunit)
vxy = slco(xstand,xstand,y,yunit)
p = ncol(xc)
nslice = length(yunit)
prob = slprob(y,yunit)
mat1 = matrix(0,p,p)
mat2 = matrix(0,p,p)
for(i in 1:nslice){
mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
          (vxy[,,i]+exy[i,]%*%t(exy[i,]))
mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-
      2*diag(p)
out = eigen(out)$vectors
out = signrt%*%out
return(out)
}
###################################################
#           Subroutine 7:                         #
#          Compute slice proportions              #
###################################################
slprob = function(y,yunit){
n = length(y)
nslice = length(yunit)
out = rep(0,nslice)
for(i in 1:nslice){
out[i] = length(y[y==yunit[i]])/n}
return(out)
}
###################################################
#            Subroutine 8:  LZC distance          #
#  dism if v1 v2 are matrices; disv if they are   #
#                    vectors                      #
###################################################
dism = function(v1,v2){
p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
d <- sqrt(sum((p1-p2)*(p1-p2)))
      return(d)
}
disv = function(v1,v2){
p1 = v1%*%t(v1)/c(t(v1)%*%v1)
p2 = v2%*%t(v2)/c(t(v2)%*%v2)
d = sqrt(sum((p1-p2)*(p1-p2)))
return(d)
}
###################################################
#           Subroutine 9:                         #
#          trace multiple correlation             #
###################################################
dis1 = function(v1,v2,x){
u1 = x%*%v1
u2 = x%*%v2
sig1 = var(u1)
sig2 = var(u2)
sig12 = cov(u1,u2)
sig21 = t(sig12)
out = matpower(sig2,-1/2)%*%sig21%*%solve(sig1)%*%
      sig12%*%matpower(sig2,-1/2)
out = sum(diag(out))/ncol(v1)
return(out)}
###################################################
#  SUBROUTINE 10: power of a matrix               #
###################################################
matpower = function(a,alpha){
a = (a + t(a))/2
tmp = eigen(a)
return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
t(tmp$vectors))}
###################################################
#  SUBROUTINE 12: discretize Y                    #
###################################################
discretize = function(y,yunit){
nsli=length(yunit)
yord = y[order(y)]
n = length(y)
nwith = round(n/nsli)
divpt = rep(0,nsli-1)
for(i in 1:(nsli-1)){
divpt[i] = yord[i*nwith+1]}
y1 = rep(0,n)
y1[y>=divpt[nsli-1]]=nsli
y1[y<divpt[1]]=1
for(i in 2:(nsli-1)){
y1[(y>=divpt[i-1])&(y<divpt[i])]=i}
return(y1)
}
###################################################
#  SUBROUTINE: Moore-Penrose type power           #
#  Taking power ignoring 0 eigenvalues;           #
#    ignoring criterion=ignore                    #
###################################################
mppower = function(matrix,power,ignore){
eig = eigen(matrix)
eval = eig$values
evec = eig$vectors
m = length(eval[abs(eval)>ignore])
tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
t(evec[,1:m])
return(tmp)
}
###################################################
#   SUBROUTINE 13: center X (n by p matrix        #
###################################################
center = function(x){
return(t(t(x)-apply(x,2,mean)))}





