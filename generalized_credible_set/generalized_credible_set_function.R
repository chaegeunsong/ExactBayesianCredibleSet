#####################################################-
##  Read R code from other R files               ####
#####################################################-
source('../source_function/source-functions.R')     # source functions
source('../source_function/classical.R')            # get classical functions
source('../source_function/common-functions.R')     # get common functions

#####################################################-
##  Packages                                     ####
#####################################################-
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(grDevices)


#####################################################-
##   function: choose a method and return beta   ####
#####################################################-
method=function(mtype){
  if(mtype=="save") beta=save(x.tra,y.tra,3,3,"categorical")
  if(mtype=="sirii") beta=sirii(x.tra,y.tra,3,3,"categorical")
  if(mtype=="cr") beta=cr(x.tra,y.tra,0.01,3)
  if(mtype=="dr") beta=dr(x.tra,y.tra,3,3,"categorical")
  return(beta)
}


#####################################################-
##   QDA Bayes rule function                     ####
#####################################################-
## Bayes rule to predict class and get posterior  ---

## input 
#' @param x : data/predictors. n x p matrix
#' @param prior : prior dist'n of y
#' @param class : class/parameter space (ex) c(0,6,9)
#' @param mu : mean vectors list (ex) mu=list(m0, mu6, mu9)
#' @param sig: sigma list (ex) sig=list(sig0, sig6, sig9)
## output 
#' @return predv class predictions
#' @return posterior probaility

br <- function(x, prior, class, mu, sig){
  n <- nrow(x)
  k <- length(class)
  postp <- matrix(0, n, k) # posterior expected loss
  predv <- matrix(NA, n, 1) # prediction vector
  
  for(l in 1:n){ # number of nrow(x)
    y <- x[l,]
    for(j in 1:k){ # number of action space
      postp[l,j] <- 
        dmvnorm(y, mean = mu[[j]], sigma = sig[[j]]) *
        prior[j]
    }
    postp[l,] <- postp[l,]/sum(postp[l,])
    postp[l,] <- postp[l,] %>% round(.,3)
    predv[l] <- class[which.max(postp[l,])]
  }
  return(list('predv'=predv, 'postp'=postp))
}  


## .-------------------------------------------  ####
#####################################################-
##  Generalized credible set                     ####
#####################################################-
## .-------------------------------------------  ####

rho <- function(k, posterior){
  # input 
  #' @param k : P(posterior > k | x )
  #' @param posterior : posterior dist'n as table
  # output
  #' @return posteriror probability that posteriror density is greater than k 
  value <- sum(posterior[posterior > k])
  return(value)
}


phi <- function(posterior, alpha = 0.05){
  # input
  #' @param posteriror : posterior dist'n given x. n x p matrix. 
  #' @param alpha : (1-alpha) credible level
  # output
  #' @return phi value: generalized CS. 
  
  ## get ka(kappa alpah) first
  grid <- seq(0,1, by=.001) # k grid
  rho_vec <- rep(NA, length(grid)) # rho vector
  for(i in 1:length(rho_vec)){
      rho_vec[i] <- rho(grid[i], posterior)  
  }
  
  ka <- c(NA)
  ka_idx <- which.max(rho_vec <= 1-alpha)
  ka <- grid[ka_idx]
  ka <- round(ka, 3) # to avoid equality error
  
  ## get gamma
  gamma <- c(NA)
  denom <- sum(posterior[posterior==ka])
    if(denom==0) {gamma<- 0
    } else {
      gamma <- ( 1-alpha - rho(ka, posterior)) / denom
    }
  
  ## get phi
  value <- rep(NA, length(posterior))
  value[posterior>ka] <- 1
  value[posterior==ka] <- gamma
  value[posterior<ka] <- 0
  value <- value %>% round(.,3)
  ## return value
  return(value)
}

#####################################################-
##   plot one bar plot                           ####
#####################################################-
## input 
#' @param v : 3-d vector representing a single point     
#' @param th : angle of the cube. typical angle is 30  
#' @param phi : generalized credible set value

onebarplot <- function(v, th, phi){
  z <- v
  zz <- z
  
  ## horizonal line : +- 0.06 
  zz <- z
  zz[1] <- zz[1]-0.06
  zz[3] <- zz[3]+0.05
  ha <- onepos(zz, th) # 2d coordinates
  
  zz <- z
  zz[1] <- zz[1]+0.06
  zz[3] <- zz[3]+0.05
  hb <- onepos(zz, th) # 2d coordinates
  
  connect(ha, hb)
  
  ## 1,2,3 prob (0.2 weight on prob)
  zz <- z
  
  # 1
  zz <- z
  zz[1] <- zz[1]-0.03
  zz[3] <- zz[3]+0.05
  h1 <- onepos(zz, th)
  
  zz <- z
  zz[1] <- zz[1]-0.03
  zz[3] <- zz[3]+0.05+0.2*phi[1]
  v1 <- onepos(zz, th)
  
  connect(h1, v1)
  
  # 2
  zz <- z
  zz[1] <- zz[1]
  zz[3] <- zz[3]+0.05
  h2 <- onepos(zz, th)
  
  zz <- z
  zz[1] <- zz[1]
  zz[3] <- zz[3]+0.05+0.2*phi[2]
  v2 <- onepos(zz, th)
  
  connect(h2, v2)
  
  
  # 3
  zz <- z
  zz[1] <- zz[1]+0.03
  zz[3] <- zz[3]+0.05
  h3 <- onepos(zz, th)
  
  zz <- z
  zz[1] <- zz[1]+0.03
  zz[3] <- zz[3]+0.05+0.2*phi[3]
  v3 <- onepos(zz, th)
  
  connect(h3, v3)
}
#onebarplot(gcrstd[p1,], 30, phi_vec[p1,])


#####################################################-
##   plot bar plots                              ####
#####################################################-
## input 
#' @param x : 3-d data as a matrix     
#' @param h1,h2,h3 : lengths of 3 axes
#' @param th : angle
#' @param index : the subsets of x plotted - this is useful
#' @param phi : generalized credible set values as a matrix

barplot3d = function(x,h1=1.5,h2=1,h3=1.3,th=30,index,phi){
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
    p <- phi[i,]
    onebarplot(v,th,p)
  }
}

#####################################################-
## function: connect two points by a bold line ## ###           
#####################################################-
## input 
#' @param a : 2-d vector
#' @param b : 2-d vector
#' @param color : color

bconnect = function(a,b, color='black'){
  lines(c(a[1],b[1]),c(a[2],b[2]), lwd=3, col=color)}


#####################################################-
## function: plot one color point                ####           
#####################################################-
## input 
#' @param v: 3-d vector representing a single point 
#' @param th: angle of the cube. typical angle is 30  
#' @param color : color

conept=function(v,th,color){
  w = onepos(v,th)
  if(color=="red")
    points(w[1],w[2],col="red", pch=16) else
      if(color=="green")
        points(w[1],w[2],col="green", pch=16) else
          points(w[1],w[2],col="blue", pch=16)
  points(w[1], w[2])
}

#####################################################-
## function: rotation matrix                     ####           
#####################################################-
## input 
#' @param w: 2-d vector representing a single point 
#' @param theta: angle to rotate in circular measure. 2*pi=360
#' 
#' output
#' @return r: rotated 2-d vector
rotate <- function(w, theta){
  r <- c(NA, NA)
  r[1] <- w[1]*cos(theta) - w[2]*sin(theta)
  r[2] <- w[1]*sin(theta) + w[2]*cos(theta)
  return(r)
}

#####################################################-
##   flower plot                                 ####
#####################################################-
## input 
#' @param v: 3-d vector representing a single point 
#' @param th: angle of the cube. typical angle is 30  
#' @param phi : phi values= generalized credible set values
#' @param color : color
#' 

oneflower3dplot <- function(v, th, phi, color){
  z <- v
  zz <- z
  w <- onepos(zz, th) # 2d coordinates / original center
  
  ## 1,2,3 prob (0.2 weight on prob)
  # 1
  zz <- z
  zz[3] <- zz[3]+0.2*phi[1]
  h1 <- onepos(zz, th)
  
  bconnect(w, h1, 'red')
  
  # 2
  zz <- z
  zz[3] <- zz[3]+0.2*phi[2]
  h2 <- onepos(zz, th)
  v2 <- h2-w
  
  # rotate
  theta <- (2*pi) * (1/3)
  r2 <- rotate(v2, theta)
  
  # translate
  t2 <- w+r2
  bconnect(w, t2,'green')
  
  
  
  # 3
  zz <- z
  zz[3] <- zz[3]+0.2*phi[3]
  h3 <- onepos(zz, th)
  v3 <- h3-w
  
  # rotate
  theta <- (2*pi) * (2/3)
  r3 <- rotate(v3, theta)
  
  # translate
  t3 <- w+r3
  bconnect(w, t3, 'blue')
  
  # color point
  conept(v, th, color=color)
}

#####################################################-
##   plot flower plots                           ####
#####################################################-
## input 
#' @param x : 3-d data as a matrix     
#' @param h1,h2,h3 : lengths of 3 axes
#' @param th : angle
#' @param index : the subsets of x plotted - this is useful
#' @param phi : generalized credible set values as a matrix
#' @param color : color

flowerplot3d = function(x,h1=1.5,h2=1,h3=1.3,th=30,index,phi, color){
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
    p <- phi[i,]
    oneflower3dplot(v,th,p,color)
  }
}


#####################################################-
##   flower plot in 2d                           ####
#####################################################-
## input 
#' @param v: 2-d vector representing a single point 
#' @param phi : phi values= generalized credible set values
#' @param r : exterior radius
#' @param color : color for center point
#' @param colorset : colors for blades
#' @param cex : inner hub size
#' @param ir : innerhub radius

library(plotrix)


oneflower2dplot <- function(v, phi, r=0.5, color,
                            colorset = c('red','green','blue','black','orange'),
                            cex=1, ir=0.07){
  z <- v[c(1,2)] # 2d coordinates / original center
  zz <- z
  n <- length(phi)
  
  # the rest
  for(i in 1:n){
    zz <- z
    #zz[2] <- zz[2]+r*phi[i]
    zz[2] <- zz[2]+r*phi[i]+ir*cex*(phi[i]!=0)
    vi <- zz-z
    
    # rotate
    theta <- (2*pi) * ((i-1)/n)
    ri <- rotate(vi, theta)
    
    # translate
    ti <- z+ri
    
    col <- colorset[i]
    #col <- adjustcolor(colorset[i], alpha.f = 0.6)
    bconnect(z, ti, col)
  }
  
  # color point
  draw.circle(z[1],z[2],col=color,radius = ir*cex)
  # points(z[1],z[2], col=color, pch=19, cex=cex)
  # points(z[1],z[2], cex=cex)
  
  # draw exterior circle using plotrix package
  draw.circle(z[1],z[2], radius = r*1.05+ir*cex)
  # draw.circle(z[1],z[2], radius = r*1.05, 
  #             border = adjustcolor('black', alpha.f = 0.6))
}

#####################################################-
##   plot flower plots in 2-d                    ####
#####################################################-
## input 
#' @param w : 2-d data as a matrix augmented with label     
#' @param index : the subsets of x plotted - this is useful
#' @param phi : generalized credible set values as a matrix
#' @param r : radius
#' @param color : color
#' @param colorset : colors for blades
#' @param cex : inner hub size
#' @param ir : innerhub radius
flowerplot2d = function(w, index, phi, r=0.5, color,
                        colorset = c('red','green','blue','black','orange'),
                        cex=1, ir=0.07){
  for(i in index){
    v <- w[i,]
    p <- phi[i,]
    oneflower2dplot(v,p,r,color,colorset,cex,ir)
  }
}

#####################################################-
##   Draw multiple circles                       ####
#####################################################-
draw.circles <- function(x,y,ir=0.07,cex=1, border=NULL, col=NA){
  n <- length(x)
  for(i in 1:n){
    draw.circle(x[i],y[i],radius = ir*cex, border = border, col=col)  
  }
}

