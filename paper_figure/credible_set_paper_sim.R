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
#par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
#par(family = "Times New Roman")
#plot(squaregrid,squaregrid,pch=" ", asp = 1)

## plot data
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x[,1],x[,2],pch=" ", asp = 1) 
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

# # most cases
# qda$postp[100,]
# rho_list[[100]]
# plot(grid, rho_list[[100]], type = 'l')
# # -> clustering is clear -> posterior 0 or 1 
# # -> ka is 1 -> phi is 0.95
# 
# # when there is ambiguity
# qda$postp[69,]
# rho_list[[69]]
# plot(grid, rho_list[[69]], type = 'l')

## get ka --------------------------------------
alpha <- 0.05
ka_vec <- rep(NA, nrow(qda$postp))
for(i in 1:length(ka_vec) ){
  ka_idx <- which.max(rho_list[[i]] <= 1-alpha)
  ka_vec[i] <- grid[ka_idx]
}

## check phi function -------------------------------
# # most cases
# qda$postp[100,]
# phi(qda$postp[100,])
# qda$postp[100,] %*% phi(qda$postp[100,]) # exp(phi)
# # -> clustering is clear -> posterior 0 or 1 
# # -> ka is 1 -> phi is 0.95
# 
# # when there is an ambiguity
# qda$postp[69,]
# phi(qda$postp[69,])
# qda$postp[69,] %*% phi(qda$postp[69,]) # exp(phi)

phi_vec <- matrix(NA, nrow(qda$postp), ncol(qda$postp))
for(i in 1:nrow(qda$postp)){
  phi_vec[i,] <- phi(qda$postp[i,])
}

## plot data
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x[,1],x[,2],pch=" ", asp = 1) 

# cairo_ps("sim1original.eps", family = "Times New Roman")
# par(family = "Times New Roman", mar = c(5, 5, 4, 2) + 0.1)
# plot(squaregrid,squaregrid,pch=" ", asp = 1, 
#      xlab=expression(x[1]), ylab=expression(x[2]),
#      xlim=c(1.5,7.5), ylim=c(1.5,8),
#      cex=2, cex.lab=2, cex.axis=2)
# points(x[y==1,1], x[y==1,2], col="red", cex=1.25) # (4.61, 3.25)
# points(x[y==2,1], x[y==2,2], col="green", cex=1.25)
# points(x[y==3,1], x[y==3,2], col="blue", cex=1.25)
# dev.off()

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
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x[,1],x[,2],pch=" ", asp = 1) 

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
# ## .-------------------------------------------  ####
# #####################################################-
# ##   simulation 2: 5 classes                     ####
# #####################################################-
# ## .-------------------------------------------  ####
# 
# ## generate simulations ----------------------------
# ## define mean vector and covariance mat
# m1 <- c(5,6); m2 <- c(4,5); m3 <- c(4,3); m4 <- c(6,3); m5 <- c(6,5)
# s <- diag(1/2, 2, 2)
# 
# ## generate data
# nsim <- 10
# set.seed(716)
# x1 <- rmvnorm(nsim, m1, s)
# x2 <- rmvnorm(nsim, m2, s)
# x3 <- rmvnorm(nsim, m3, s)
# x4 <- rmvnorm(nsim, m4, s)
# x5 <- rmvnorm(nsim, m5, s)
# 
# x <- rbind(x1,x2,x3,x4,x5)
# y <- rep(1:5, each=nsim)
# 
# data2 <- cbind(x,y)
# 
# ## grid
# # make squaregrid for pretty plot. +-1 to make a room for outer circle
# squaregrid <- c(min(range(x[,1])-0.5, range(x[,2])-0.5), 
#                 max(range(x[,1])+0.5, range(x[,2])+0.5))
# plot(squaregrid,squaregrid,pch=" ", asp = 1)
# ## plot train
# # par(pty="s") # to make the plot region square (independent of size and limits)
# # plot(x[,1],x[,2],pch=" ", asp = 1) 
# plot(squaregrid,squaregrid,pch=" ", asp = 1,
#      xlab=expression(x[1]), ylab=expression(x[2]),
#      cex=1.5, cex.lab=1.25, cex.axis=1.25)
# points(x[y==1,1], x[y==1,2], col="red")
# points(x[y==2,1], x[y==2,2], col="green")
# points(x[y==3,1], x[y==3,2], col="blue")
# points(x[y==4,1], x[y==4,2], col="black")
# points(x[y==5,1], x[y==5,2], col="orange")
# 
# 
# ## classification ----------------------------------------
# p <- ncol(x)
# n1 <- nrow(x[y==1,])
# n2 <- nrow(x[y==2,])
# n3 <- nrow(x[y==3,])
# n4 <- nrow(x[y==4,])
# n5 <- nrow(x[y==5,])
# 
# x1 <- x[y==1,]
# x2 <- x[y==2,]
# x3 <- x[y==3,]
# x4 <- x[y==4,]
# x5 <- x[y==5,]
# 
# ## estimate mu
# mu1 <- (1/n1) * ( x1 %>% colSums() ) # apply(pred1, 2, mean)
# mu2 <- (1/n2) * ( x2 %>% colSums() ) # apply(pred2, 2, mean)
# mu3 <- (1/n3) * ( x3 %>% colSums() ) # apply(pred3, 2, mean)
# mu4 <- (1/n4) * ( x4 %>% colSums() ) # apply(pred4, 2, mean)
# mu5 <- (1/n5) * ( x5 %>% colSums() ) # apply(pred5, 2, mean)
# 
# mu <- list(mu1, mu2, mu3, mu4, mu5)
# 
# ## estimate sigma
# sig1 <- matrix(1, p, p)
# for(i in 1:n1){
#   sig1 <- sig1+(x1[i,]-mu1) %*% t((x1[i,]-mu1))
# }
# sig1 <- (1/(n1-1)) * sig1
# 
# sig2 <- matrix(1, p, p)
# for(i in 1:n2){
#   sig2 <- sig2+(x2[i,]-mu2) %*% t((x2[i,]-mu2))
# }
# sig2 <- (1/(n2-1)) * sig2
# 
# sig3 <- matrix(1, p, p)
# for(i in 1:n3){
#   sig3 <- sig3+(x3[i,]-mu3) %*% t((x3[i,]-mu3))
# }
# sig3 <- (1/(n3-1)) * sig3
# 
# sig4 <- matrix(1, p, p)
# for(i in 1:n4){
#   sig4 <- sig4+(x4[i,]-mu4) %*% t((x4[i,]-mu4))
# }
# sig4 <- (1/(n4-1)) * sig4
# 
# sig5 <- matrix(1, p, p)
# for(i in 1:n5){
#   sig5 <- sig5+(x5[i,]-mu5) %*% t((x5[i,]-mu5))
# }
# sig5 <- (1/(n5-1)) * sig5
# 
# sig <- list(sig1, sig2, sig3, sig4, sig5)
# 
# 
# ## prediction
# ## prior distribution 
# (prior <- table(y) / length(y))
# prior <- prior %>% round(.,3)
# 
# qda <- br(x = x, prior = prior, class = c(1,2,3,4,5), mu = mu, sig = sig)
# head(qda$predv)
# 
# table(qda$predv)
# table(qda$predv, y) # confusion matrix
# 
# #hist(qda$predv) # histogram of qda
# 
# post <- table(qda$predv) / length(qda$predv) # posterior dist'n
# post %>% round(.,3)
# 
# ## posterior probability
# qda$postp %>% head()
# # Most of them cosists of 0 and 1 because it is easy to 
# # do clustering. 
# 
# ## wrong clustering
# wrong_idx <- which(qda$predv != y)
# wrong_idx
# 
# 
# ## check rho function -------------------------------
# qda$postp %>% dim()
# 
# grid <- seq(0,1, by=.001) # k grid
# rho_list <- vector(mode = 'list', length = nrow(qda$postp)) # rho list
# rho_vec <- rep(NA, length(grid)) # rho vector
# for (j in 1:length(rho_list)){
#   for(i in 1:length(grid)){
#     rho_vec[i] <- rho(grid[i], qda$postp[j,])  
#   }
#   rho_list[[j]] <- rho_vec
# }
# 
# # # most cases
# # qda$postp[100,]
# # rho_list[[100]]
# # plot(grid, rho_list[[100]], type = 'l')
# # # -> clustering is clear -> posterior 0 or 1 
# # # -> ka is 1 -> phi is 0.95
# # 
# # # when there is ambiguity
# # qda$postp[69,]
# # rho_list[[69]]
# # plot(grid, rho_list[[69]], type = 'l')
# 
# ## get ka --------------------------------------
# alpha <- 0.05
# ka_vec <- rep(NA, nrow(qda$postp))
# for(i in 1:length(ka_vec) ){
#   ka_idx <- which.max(rho_list[[i]] <= 1-alpha)
#   ka_vec[i] <- grid[ka_idx]
# }
# 
# 
# ## check phi function -------------------------------
# # # most cases
# # qda$postp[100,]
# # phi(qda$postp[100,])
# # qda$postp[100,] %*% phi(qda$postp[100,]) # exp(phi)
# # # -> clustering is clear -> posterior 0 or 1 
# # # -> ka is 1 -> phi is 0.95
# # 
# # # when there is an ambiguity
# # qda$postp[69,]
# # phi(qda$postp[69,])
# # qda$postp[69,] %*% phi(qda$postp[69,]) # exp(phi)
# 
# phi_vec <- matrix(NA, nrow(qda$postp), ncol(qda$postp))
# for(i in 1:nrow(qda$postp)){
#   phi_vec[i,] <- phi(qda$postp[i,])
# }
# 
# ## steering wheel plot ------------------------------------------
# #par(pty="s") # to make the plot region square (independent of size and limits)
# #plot(x[,1],x[,2],pch=" ", asp = 1) 
# 
# plot(squaregrid,squaregrid,pch=" ", asp = 1,
#      xlab=expression(x[1]), ylab=expression(x[2]),
#      cex=1.5, cex.lab=1.25, cex.axis=1.25)
# points(x[qda$predv==1,1], x[qda$predv==1,2], col="red")
# points(x[qda$predv==2,1], x[qda$predv==2,2], col="green")
# points(x[qda$predv==3,1], x[qda$predv==3,2], col="blue")
# points(x[qda$predv==4,1], x[qda$predv==4,2], col="black")
# points(x[qda$predv==5,1], x[qda$predv==5,2], col="orange")
# 
# n <- length(qda$predv)
# flowerplot2d(cbind(x,y),(1:n)[qda$predv==1],phi_vec, 'red',r=0.3)
# flowerplot2d(cbind(x,y),(1:n)[qda$predv==2],phi_vec, 'green',r=0.3)
# flowerplot2d(cbind(x,y),(1:n)[qda$predv==3],phi_vec, 'blue',r=0.3)
# flowerplot2d(cbind(x,y),(1:n)[qda$predv==4],phi_vec, 'black',r=0.3)
# flowerplot2d(cbind(x,y),(1:n)[qda$predv==5],phi_vec, 'orange',r=0.3)
# 

