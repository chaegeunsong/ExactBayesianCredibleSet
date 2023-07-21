## .-------------------------------------------  ####
#####################################################-
##   simulation 1: 3 classes                     ####
#####################################################-
## .-------------------------------------------  ####

## generate simulations ----------------------------
## define mean vector and covariance mat
m1 <- c(5,6); m2 <- c(4,5); m3 <- c(6,4)
s <- diag(1, 2, 2)

## generate data
nsim <- 10
set.seed(716)
x1 <- rmvnorm(nsim, m1, s)
x2 <- rmvnorm(nsim, m2, s)
x3 <- rmvnorm(nsim, m3, s)

x <- rbind(x1,x2,x3)
y <- rep(1:3, each=nsim)

data1 <- cbind(x,y)

## grid
# make squaregrid for pretty plot. +-1 to make a room for outer circle
squaregrid <- c(min(range(x[,1])-0.5, range(x[,2])-0.5), 
                max(range(x[,1])+0.5, range(x[,2])+0.5))

## plot data
par(family = "Times New Roman", mar = c(5, 5, 4, 2) + 0.1)
plot(squaregrid,squaregrid,pch=" ", asp = 1, 
     xlab=expression(x[1]), ylab=expression(x[2]),
     cex=2, cex.lab=2, cex.axis=2)
points(x[y==1,1], x[y==1,2], col="red", cex=1.25)
points(x[y==2,1], x[y==2,2], col="green", cex=1.25)
points(x[y==3,1], x[y==3,2], col="blue", cex=1.25)

## classification ----------------------------------------
p <- ncol(x)
n1 <- nrow(x[y==1,])
n2 <- nrow(x[y==2,])
n3 <- nrow(x[y==3,])

x1 <- x[y==1,]
x2 <- x[y==2,]
x3 <- x[y==3,]

## estimate mu
mu1 <- (1/n1) * ( x1 %>% colSums() ) # apply(pred1, 2, mean)
mu2 <- (1/n2) * ( x2 %>% colSums() ) # apply(pred2, 2, mean)
mu3 <- (1/n3) * ( x3 %>% colSums() ) # apply(pred3, 2, mean)
mu <- list(mu1, mu2, mu3)

## estimate sigma
sig1 <- matrix(1, p, p)
for(i in 1:n1){
  sig1 <- sig1+(x1[i,]-mu1) %*% t((x1[i,]-mu1))
}
sig1 <- (1/(n1-1)) * sig1

sig2 <- matrix(1, p, p)
for(i in 1:n2){
  sig2 <- sig2+(x2[i,]-mu2) %*% t((x2[i,]-mu2))
}
sig2 <- (1/(n2-1)) * sig2

sig3 <- matrix(1, p, p)
for(i in 1:n3){
  sig3 <- sig3+(x3[i,]-mu3) %*% t((x3[i,]-mu3))
}
sig3 <- (1/(n3-1)) * sig3
sig <- list(sig1, sig2, sig3)


## prediction
## prior distribution 
(prior <- table(y) / length(y))
prior <- prior %>% round(.,3)

qda <- br(x = x, prior = prior, class = c(1,2,3), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv,y) # confusion matrix

#hist(qda$predv) # histogram of qda

post <- table(qda$predv) / length(qda$predv) # posterior dist'n
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y)
wrong_idx


## check rho function -------------------------------
qda$postp %>% dim()

grid <- seq(0,1, by=.001) # k grid
rho_list <- vector(mode = 'list', length = nrow(qda$postp)) # rho list
rho_vec <- rep(NA, length(grid)) # rho vector
for (j in 1:length(rho_list)){
  for(i in 1:length(grid)){
    rho_vec[i] <- rho(grid[i], qda$postp[j,])  
  }
  rho_list[[j]] <- rho_vec
}

## get ka --------------------------------------
alpha <- 0.05
ka_vec <- rep(NA, nrow(qda$postp))
for(i in 1:length(ka_vec) ){
  ka_idx <- which.max(rho_list[[i]] <= 1-alpha)
  ka_vec[i] <- grid[ka_idx]
}


phi_vec <- matrix(NA, nrow(qda$postp), ncol(qda$postp))
for(i in 1:nrow(qda$postp)){
  phi_vec[i,] <- phi(qda$postp[i,])
}

cairo_ps("sim1original.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 2) + 0.1)
plot(c(1.5,8),c(1.5,8),pch=" ", #asp = 1, 
     xlab=expression(x[1]), ylab=expression(x[2]),
     #xlim=c(1.5,7.5), ylim=c(1.5,8),
     cex=2, cex.lab=2, cex.axis=2)
# points(x[y==1,1], x[y==1,2], col="red", cex=1.25) # 2:(4.61, 3.25), 4:(6.07, 7.40)
# points(x[y==2,1], x[y==2,2], col="green", cex=1.25)
# points(x[y==3,1], x[y==3,2], col="blue", cex=1.25)
draw.circles(x[y==1,1], x[y==1,2], border="red", cex=1.25, ir=0.07)
draw.circles(x[y==2,1], x[y==2,2], border="green", cex=1.25, ir=0.07)
draw.circles(x[y==3,1], x[y==3,2], border="blue", cex=1.25, ir=0.07)

text(6.07, 7.40+0.3, 'A', cex=1.5)
text(4.61, 3.25+0.3, 'B', cex=1.5)
dev.off()

## A and B
# 4th : A=(6.07, 7.40), phi = 1, 0.49, 1
# 2nd : B=(4.61, 3.25), phi = 0.967, 0, 0
phi_vec

## steering wheel plot ------------------------------------------

cairo_ps("sim1steer.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 2) + 0.1)
plot(c(1.5,8),c(1.5,8),pch=" ", #asp = 1,
     xlab=expression(x[1]), ylab=expression(x[2]),
     #xlim=c(2,8), ylim=c(2,8),
     cex=2, cex.lab=2, cex.axis=2)
draw.circles(x[qda$predv==1,1], x[qda$predv==1,2], border="red", cex=1.25, ir=0.07)
draw.circles(x[qda$predv==2,1], x[qda$predv==2,2], border="green", cex=1.25, ir=0.07)
draw.circles(x[qda$predv==3,1], x[qda$predv==3,2], border="blue", cex=1.25, ir=0.07)

n <- length(qda$predv)
flowerplot2d(cbind(x,y),(1:n)[qda$predv==1],phi_vec, r=0.3, 'red', cex=1.25, ir=0.07)
flowerplot2d(cbind(x,y),(1:n)[qda$predv==2],phi_vec, r=0.3,'green', cex=1.25, ir=0.07)
flowerplot2d(cbind(x,y),(1:n)[qda$predv==3],phi_vec, r=0.3,'blue', cex=1.25, ir=0.07)
dev.off()
