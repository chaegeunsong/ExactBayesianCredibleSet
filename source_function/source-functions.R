#####################################################-
##  Source functions                             ####
#####################################################-
#####################################################
#          function 1:                               
#       sliced inverse regression                    
#  sig = variance of xc; signrt = sig^(-1/2)         
#  y must be discretized according to yunit for  
#  continuous y; for discrete y, use original y
#  with yunit being the distinct values in 
#  discrete y                                        
#####################################################
sir=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  yless=ydis;ylabel=numeric()
  for(i in 1:n) {if(var(yless)!=0) {ylabel=c(ylabel,yless[1]);yless=yless[yless!=yless[1]]}}
  ylabel=c(ylabel,yless[1])
  prob=numeric();exy=numeric()
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  sirmat=t(exy)%*%diag(prob)%*%exy
  return(signrt%*%eigen(sirmat)$vectors[,1:r])}
#plot(c(x%*%sir(x,y,6,1,"continuous")),y)
####################################################
#     function 7: discretize Y                    
#     note that i added a small perturbation to 
#     y; when y is discrete it will create wierd
#     slices if you treat it as continuous
####################################################
discretize = function(y,h){
  n=length(y);m=round(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1]) 
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i 
  return(y1)
}
#####################################################
#    function: choose a method
#####################################################
method=function(mtype){
  if(mtype=="sir") beta=sir(x,y,3,3,"categorical")
  if(mtype=="save") beta=save(x,y,3,3,"categorical")
  if(mtype=="sirii") beta=sirii(x,y,3,3,"categorical")
}



#####################################################################
#   2/20/2014: This code create 3-d scatterplot                       
#####################################################################
#####################################################################
#    function: connect two points by a line                      
#####################################################################
connect = function(a,b){
  lines(c(a[1],b[1]),c(a[2],b[2]))}
#####################################################################
#    function: connect two points by a dotted line                      
#####################################################################
dotconnect = function(a,b){
  points(seq(from=a[1],to=b[1],length=50),
         seq(from=a[2],to=b[2],length=50),pch=".",cex=1.2)}
#####################################################################
#    function: convert a 3d position to a 2d position                     
#####################################################################
onepos=function(v,th){
  theta=th*(pi/180)
  v1=c(v[1]*cos(theta),v[1]*sin(theta))
  theta1=pi-theta
  v2=c(v[2]*cos(theta1),v[2]*sin(theta1))
  v3=c(0,v[3])
  w1=v1
  w2=v1+v2
  w3=v1+v2+v3
  return(w3)
}
#####################################################################
#        function: plot one point                               
#        input v: 3-d vector representing a single point           
#        th: angle of the cube. typical angle is 30                
#        color: 1 = red, 2 = green, 3 = blue
#####################################################################
onept=function(v,th,color){
  w = onepos(v,th)
  if(color=="red")
    points(w[1],w[2],col="red") else
      if(color=="green")
        points(w[1],w[2],col="green") else
          points(w[1],w[2],col="blue")
}
####################################################################
#       function:             create frame                            
#                  h1, h2, h3 are lengths in x, y, z axis  
#                  (typical values are 1.5,1,2)             
#                  th is angle (typical value is 30)
####################################################################
frame3d = function(h1=1.5,h2=1,h3=1.3,th=30,mtype){
  point = numeric()
  for(i in c(0,h1)){
    for(j in c(0,h2)){
      for(k in c(0,h3)){
        point = rbind(point,onepos(c(i,j,k),th))}}}
  x1 = min(point[,1]-0.05*h1)
  x2 = max(point[,1])
  y2=max(point[,2]);y1=min(point[,2])
  scale=(x2-x1)*0.9/(y2-y1)
  plot(0,0,ylim=c(min(point[,2]),max(point[,2])),
       xlim=c(x1,x2),
       axes=F,xlab=" ",ylab=" ",pch=" ")
  connect(point[1,],point[2,])
  connect(point[2,],point[4,])
  connect(point[4,],point[3,])
  connect(point[3,],point[1,])
  connect(point[1,],point[5,])
  connect(point[5,],point[6,])
  connect(point[6,],point[8,])
  dotconnect(point[8,],point[7,])
  dotconnect(point[7,],point[5,])
  dotconnect(point[7,],point[3,])
  connect(point[2,],point[6,])
  connect(point[4,],point[8,])
  pt1=(point[1,]+point[5,])/2
  pt2=(point[2,]+point[6,])/2
  pt3=pt2+(pt1-pt2)*1.08
  point=rbind(point,pt1,pt2,pt3) 
  angle1=atan(scale*(point[5,2]-point[1,2])/(point[5,1]-point[1,1]));angle1=angle1*180/pi
  pt1=(point[1,]+point[2,])/2;pt2=(point[3,]+point[4,])/2;pt3=pt1+(pt2-pt1)*1.08
  point=rbind(point,pt1,pt2,pt3)
  angle2=atan(scale*(point[4,2]-point[2,2])/(point[4,1]-point[2,1]));angle2=angle2*180/pi
  pt1=(point[2,]+point[4,])/2;pt2=(point[1,]+point[3,])/2;pt3=pt1+(pt2-pt1)*1.06
  point=rbind(point,pt1,pt2,pt3)
  if(mtype=="sir") {
    text(point[11,1],point[11,2], "first SIR direction",srt=angle1)
    text(point[14,1],point[14,2],"third SIR direction",srt=90)
    text(point[17,1],point[17,2],"second SIR direction",srt=angle2)}
  if(mtype=="save") {
    text(point[11,1],point[11,2], "first SAVE direction",srt=angle1)
    text(point[14,1],point[14,2],"third SAVE direction",srt=90)
    text(point[17,1],point[17,2],"second SAVE direction",srt=angle2)}
  if(mtype=="sirii"){
    text(point[11,1],point[11,2], "first SIR-II direction",srt=angle1)
    text(point[14,1],point[14,2],"third SIR-II direction",srt=90)
    text(point[17,1],point[17,2],"second SIR-II direction",srt=angle2)}
  if(mtype=="cr"){
    text(point[11,1],point[11,2], "first CR direction",srt=angle1)
    text(point[14,1],point[14,2],"third CR direction",srt=90)
    text(point[17,1],point[17,2],"second CR direction",srt=angle2)}
  if(mtype=="dr"){
    text(point[11,1],point[11,2], "first DR direction",srt=angle1)
    text(point[14,1],point[14,2],"third DR direction",srt=90)
    text(point[17,1],point[17,2],"second DR direction",srt=angle2)}
  #for(i in 1:nrow(point)) text(point[i,1],point[i,2],i)
}

#####################################################################
#             function:       plot data                             
#  x is the 3-d data; h1,h2,h3: lenghts of 3 axes
#  th: angle; index: the subsets of x plotted - this is useful
#  when plotting several groups; color: the color of that subgroup
#####################################################################
data3d = function(x,h1=1.5,h2=1,h3=1.3,th=30,index,color){
  gcr=x
  gcrstd=gcr
  lim2=max(gcr[,1])
  lim1=min(gcr[,1])
  del=lim2-lim1
  lim0=lim1-del/10
  lim3=lim2+del/10
  gcrstd[,1]=(gcr[,1]-lim0)/(lim3-lim0)*h1
  lim2=max(gcr[,2])
  lim1=min(gcr[,2])
  del=lim2-lim1
  lim0=lim1-del/10
  lim3=lim2+del/10
  gcrstd[,2]=(gcr[,2]-lim0)/(lim3-lim0)*h2
  lim2=max(gcr[,3])
  lim1=min(gcr[,3])
  del=lim2-lim1
  lim0=lim1-del/10
  lim3=lim2+del/10
  gcrstd[,3]=(gcr[,3]-lim0)/(lim3-lim0)*h3
  n=nrow(gcrstd)
  for(i in index){
    v=c(gcrstd[i,])
    onept(v,th,color)
  }
}
#####################################################################
#   axis label
#####################################################################
axislab=function(mtype){
  pt1=c(.7,.25);pt2=c(-0.5,0.15);pt3=c(-.91,1.4)
  text(pt2[1],pt2[2],"first mtype predictor",srt=-17)
  text(pt1[1],pt1[2],"second SIR predictor",srt=19)
  text(pt3[1],pt3[2],"third SIR predictor",srt=270)}

###############################################-
#### matrix power
###############################################-
matpower=function(a,alpha){
  a=(a+t(a))/2;tmp=eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))
}

###############################################-
## discretize continuous response// h=round(n/50)
###############################################-
discretize=function(y,h){
  n=length(y);m=floor(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)
}

###############################################-
## Directional regression
###############################################-
dr <- function(x,y,h,r,ytype){
  ## x: n*p predictor matrix
  ## y: n*1 response
  ## h: number of slices// h=round(n/50)
  ## r: dimension of sc
  p <- ncol(x);n <- nrow(x)
  signrt <-  matpower(var(x),-1/2)
  xc <-  t(t(x)-apply(x,2,mean))
  xst <-  xc%*%signrt
  if(ytype=="continuous") ydis <- discretize(y,h)
  if(ytype=="categorical") ydis <- y
  yless <- ydis; ylabel <- numeric()
  for(i in 1:n) {
    if(var(yless)!=0) {
      ylabel <- c(ylabel,yless[1])
      yless <- yless[yless!=yless[1]]
    }
  }
  ylabel <- c(ylabel,yless[1])
  prob=numeric()
  for(i in 1:h) prob <- c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h)); exy=numeric()
  for(i in 1:h) {
    vxy[,,i] <- var(xst[ydis==ylabel[i],])
    exy <- rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  }
  mat1 <- matrix(0,p,p)
  mat2 <- matrix(0,p,p)
  for(i in 1:h){
    mat1 <- mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 <- mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out <- 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  return(signrt%*%eigen(out)$vectors[,1:r])
}


###############################################-
## evaluation
###############################################-
standmat <- function(x) {
  x <- as.matrix(x)
  return(x%*%(crossprod(x)%^%(-1/2)))
}
###############################################-
## frobenius norm diff btw projection matrix
###############################################-
eva <- function(x,y){
  x <- as.matrix(x); y <- as.matrix(y)
  x <- standmat(x); y <- standmat(y)
  error <- tcrossprod(x)-tcrossprod(y)
  return(norm(error, type = 'F'))
}


##################################################################### 
#          Contour Regression 5/15/2017            
#####################################################################
##################################################################### 
#         contour regression
#         percent: percent of small abs(y[i]-y[j]) used
#                  .1 to .2 usually works well         
#####################################################################
crmat = function(x,y,percent){
  z = standmat(x)
  n = dim(x)[1]
  p = dim(x)[2]
  deltax = numeric()
  deltay = numeric()
  for(i in 2:n){
    for(j in 1:(i-1)){
      deltax = rbind(deltax,z[i,]-z[j,])
      deltay = c(deltay,y[i]-y[j])}}
  nn = length(deltay)
  index = order(abs(deltay))[1:round(percent*nn)]
  cont = deltax[index,]
  crmat = matpower(var(x),-1/2)%*%eigen(t(cont)%*%cont)$vectors
  return(crmat)
}
##################################################################### 
#                      transform index 
#####################################################################
tradeindex12 = function(k,n){
  j = ceiling(k/n)
  i = k - (j-1)*n
  return(c(i,j))}
tradeindex21 = function(i,j,n){
  return((j-1)*n+i)}
##################################################################### 
#         contour regression: fast version
#         percent: percent of small abs(y[i]-y[j]) used
#                  .1 to .2 usually works well         
#####################################################################
crfastmat = function(x,y,percent){
  z = standmat(x)
  n = dim(x)[1]
  p = dim(x)[2]
  ymat = matrix(y,n,n)
  deltay  = c(abs(ymat - t(ymat)))
  singleindex = (1:n^2)[deltay < percent*mean(deltay)]
  contourmat = matrix(0,p,p)
  for(k in singleindex){
    doubleindex = tradeindex12(k,n) 
    deltaz = z[doubleindex[1],]-z[doubleindex[2],]
    contourmat = contourmat + deltaz %*% t(deltaz)}
  signrt = matpower(var(x),-1/2)
  return(signrt%*%eigen(contourmat)$vectors)
}


#########################################
cr = function(x,y,percent,r){                                  
  tradeindex12 = function(k,n){                                
    j = ceiling(k/n)                                             
    i = k - (j-1)*n                                              
    return(c(i,j))}                                                                                    
  mu=apply(x,2,mean);signrt=matpower(var(x),-1/2)              
  z=t(t(x)-mu)%*%signrt                                        
  n=dim(x)[1];p = dim(x)[2]                                    
  ymat=matrix(y,n,n)                                           
  deltay=c(abs(ymat - t(ymat)))                                
  singleindex=(1:n^2)[deltay < percent*mean(deltay)]           
  contourmat=matrix(0,p,p)                                     
  for(k in singleindex){                                       
    doubleindex=tradeindex12(k,n)                                
    deltaz=z[doubleindex[1],]-z[doubleindex[2],]                 
    contourmat=contourmat+deltaz %*% t(deltaz)}                  
  signrt=matpower(var(x),-1/2)                                 
  return(signrt%*%eigen(contourmat)$vectors[,p:(p-r+1)])                   
}     


#####################################################
#    function: choose a method
#####################################################
method=function(mtype){
  if(mtype=="save") beta=save(x.tra,y.tra,3,3,"categorical")
  if(mtype=="sirii") beta=sirii(x.tra,y.tra,3,3,"categorical")
  if(mtype=="cr") beta=cr(x.tra,y.tra,0.01,3)
  if(mtype=="dr") beta=dr(x.tra,y.tra,3,3,"categorical")
  return(beta)
}

