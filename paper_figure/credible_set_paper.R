## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##  Read R code from other R files               ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
source('credible-set.R')
source('pendigit.R')

## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   Simulation                                  ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
library(mvtnorm)
library(dplyr)

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
nsim <- 50
set.seed(716)
x1 <- rmvnorm(nsim, m1, s)
x2 <- rmvnorm(nsim, m2, s)
x3 <- rmvnorm(nsim, m3, s)

x <- rbind(x1,x2,x3)
y <- rep(1:3, each=nsim)

data1 <- cbind(x,y)

## split data into train and test
set.seed(1113)
index.tra <- sample(ntotal, size = ntra)

train1 <- data1[index.tra,]
test1 <- data1[-index.tra,]

ntotal <- nsim*length(unique(y))
ntra <- ntotal*0.8
ntes <- ntotal*0.2

x_train1 <- train1[,1:2]
y_train1 <- train1[,3]

x_test1 <- test1[,1:2]
y_test1 <- test1[,3]

x.tra1 <- as.matrix(x_train1)
y.tra1 <- as.matrix(y_train1)

x.tes1 <- as.matrix(x_test1)
y.tes1 <- as.matrix(y_test1)

## grid
# make squaregrid for pretty plot. +-1 to make a room for outer circle
squaregrid <- c(min(range(x[,1])-0.5, range(x[,2])-0.5), 
                max(range(x[,1])+0.5, range(x[,2])+0.5))
plot(squaregrid,squaregrid,pch=" ", asp = 1)

## plot train
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x.tra1[,1],x.tra1[,2],pch=" ", asp = 1) 
plot(squaregrid,squaregrid,pch=" ", asp = 1, 
     xlab=expression(x[1]), ylab=expression(x[2]))
points(x.tra1[y.tra1==1,1], x.tra1[y.tra1==1,2], col="red")
points(x.tra1[y.tra1==2,1], x.tra1[y.tra1==2,2], col="green")
points(x.tra1[y.tra1==3,1], x.tra1[y.tra1==3,2], col="blue")

## plot test
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x.tes1[,1],x.tes1[,2],pch=" ", asp = 1) 
plot(squaregrid,squaregrid,pch=" ", asp = 1,
     xlab=expression(x[1]), ylab=expression(x[2]))
points(x.tes1[y.tes1==1,1], x.tes1[y.tes1==1,2], col="red")
points(x.tes1[y.tes1==2,1], x.tes1[y.tes1==2,2], col="green")
points(x.tes1[y.tes1==3,1], x.tes1[y.tes1==3,2], col="blue")

## classification ----------------------------------------
p <- ncol(x.tra1)
n1 <- nrow(x.tra1[y.tra1==1,])
n2 <- nrow(x.tra1[y.tra1==2,])
n3 <- nrow(x.tra1[y.tra1==3,])

x1 <- x.tra1[y.tra1==1,]
x2 <- x.tra1[y.tra1==2,]
x3 <- x.tra1[y.tra1==3,]

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
(prior <- table(y.tra) / length(y.tra))
prior <- prior %>% round(.,3)

qda <- br(x = x.tes1, prior = prior, class = c(1,2,3), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv,y.tes1) # confusion matrix

#hist(qda$predv) # histogram of qda

post <- table(qda$predv) / length(qda$predv) # posterior dist'n
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y.tes1)
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

## steering wheel plot ------------------------------------------
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x.tes1[,1],x.tes1[,2],pch=" ", asp = 1) 

plot(squaregrid,squaregrid,pch=" ", asp = 1,
     xlab=expression(x[1]), ylab=expression(x[2]))
points(x.tes1[qda$predv==1,1], x.tes1[qda$predv==1,2], col="red")
points(x.tes1[qda$predv==2,1], x.tes1[qda$predv==2,2], col="green")
points(x.tes1[qda$predv==3,1], x.tes1[qda$predv==3,2], col="blue")

n <- length(qda$predv)
flowerplot2d(cbind(x.tes1,y.tes1),(1:n)[qda$predv==1],phi_vec, r=0.33, 'red')
flowerplot2d(cbind(x.tes1,y.tes1),(1:n)[qda$predv==2],phi_vec, r=0.33,'green')
flowerplot2d(cbind(x.tes1,y.tes1),(1:n)[qda$predv==3],phi_vec, r=0.33,'blue')

## .-------------------------------------------  ####
#####################################################-
##   simulation 2: 5 classes                     ####
#####################################################-
## .-------------------------------------------  ####

## generate simulations ----------------------------
## define mean vector and covariance mat
m1 <- c(5,6); m2 <- c(4,5); m3 <- c(4,3); m4 <- c(6,3); m5 <- c(6,5)
s <- diag(1/2, 2, 2)

## generate data
nsim <- 30
set.seed(716)
x1 <- rmvnorm(nsim, m1, s)
x2 <- rmvnorm(nsim, m2, s)
x3 <- rmvnorm(nsim, m3, s)
x4 <- rmvnorm(nsim, m4, s)
x5 <- rmvnorm(nsim, m5, s)

x <- rbind(x1,x2,x3,x4,x5)
y <- rep(1:5, each=nsim)

data2 <- cbind(x,y)

## split data into train and test
set.seed(78)
index.tra <- sample(ntotal, size = ntra)

train2 <- data2[index.tra,]
test2 <- data2[-index.tra,]

ntotal <- nsim*length(unique(y))
ntra <- ntotal*0.8
ntes <- ntotal*0.2

x_train2 <- train2[,1:2]
y_train2 <- train2[,3]

x_test2 <- test2[,1:2]
y_test2 <- test2[,3]

x.tra2 <- as.matrix(x_train2)
y.tra2 <- as.matrix(y_train2)

x.tes2 <- as.matrix(x_test2)
y.tes2 <- as.matrix(y_test2)

## grid
# make squaregrid for pretty plot. +-1 to make a room for outer circle
squaregrid <- c(min(range(x[,1])-0.5, range(x[,2])-0.5), 
                max(range(x[,1])+0.5, range(x[,2])+0.5))
plot(squaregrid,squaregrid,pch=" ", asp = 1)
## plot train
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x.tra2[,1],x.tra2[,2],pch=" ", asp = 1) 
plot(squaregrid,squaregrid,pch=" ", asp = 1,
     xlab=expression(x[1]), ylab=expression(x[2]))
points(x.tra2[y.tra2==1,1], x.tra2[y.tra2==1,2], col="red")
points(x.tra2[y.tra2==2,1], x.tra2[y.tra2==2,2], col="green")
points(x.tra2[y.tra2==3,1], x.tra2[y.tra2==3,2], col="blue")
points(x.tra2[y.tra2==4,1], x.tra2[y.tra2==4,2], col="black")
points(x.tra2[y.tra2==5,1], x.tra2[y.tra2==5,2], col="orange")

## plot test
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(x.tes2[,1],x.tes2[,2],pch=" ", asp = 1) 
plot(squaregrid,squaregrid,pch=" ", asp = 1,
     xlab=expression(x[1]), ylab=expression(x[2]))
points(x.tes2[y.tes2==1,1], x.tes2[y.tes2==1,2], col="red")
points(x.tes2[y.tes2==2,1], x.tes2[y.tes2==2,2], col="green")
points(x.tes2[y.tes2==3,1], x.tes2[y.tes2==3,2], col="blue")
points(x.tes2[y.tes2==4,1], x.tes2[y.tes2==4,2], col="black")
points(x.tes2[y.tes2==5,1], x.tes2[y.tes2==5,2], col="orange")


## classification ----------------------------------------
p <- ncol(x.tra2)
n1 <- nrow(x.tra2[y.tra2==1,])
n2 <- nrow(x.tra2[y.tra2==2,])
n3 <- nrow(x.tra2[y.tra2==3,])
n4 <- nrow(x.tra2[y.tra2==4,])
n5 <- nrow(x.tra2[y.tra2==5,])

x1 <- x.tra2[y.tra2==1,]
x2 <- x.tra2[y.tra2==2,]
x3 <- x.tra2[y.tra2==3,]
x4 <- x.tra2[y.tra2==4,]
x5 <- x.tra2[y.tra2==5,]

## estimate mu
mu1 <- (1/n1) * ( x1 %>% colSums() ) # apply(pred1, 2, mean)
mu2 <- (1/n2) * ( x2 %>% colSums() ) # apply(pred2, 2, mean)
mu3 <- (1/n3) * ( x3 %>% colSums() ) # apply(pred3, 2, mean)
mu4 <- (1/n4) * ( x4 %>% colSums() ) # apply(pred4, 2, mean)
mu5 <- (1/n5) * ( x5 %>% colSums() ) # apply(pred5, 2, mean)

mu <- list(mu1, mu2, mu3, mu4, mu5)

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

sig4 <- matrix(1, p, p)
for(i in 1:n4){
  sig4 <- sig4+(x4[i,]-mu4) %*% t((x4[i,]-mu4))
}
sig4 <- (1/(n4-1)) * sig4

sig5 <- matrix(1, p, p)
for(i in 1:n5){
  sig5 <- sig5+(x5[i,]-mu5) %*% t((x5[i,]-mu5))
}
sig5 <- (1/(n5-1)) * sig5

sig <- list(sig1, sig2, sig3, sig4, sig5)


## prediction
## prior distribution 
(prior <- table(y.tra2) / length(y.tra2))
prior <- prior %>% round(.,3)

qda <- br(x = x.tes2, prior = prior, class = c(1,2,3,4,5), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv, y.tes2) # confusion matrix

#hist(qda$predv) # histogram of qda

post <- table(qda$predv) / length(qda$predv) # posterior dist'n
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y.tes2)
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

## steering wheel plot ------------------------------------------
#par(pty="s") # to make the plot region square (independent of size and limits)
#plot(x.tes2[,1],x.tes2[,2],pch=" ", asp = 1) 

plot(squaregrid,squaregrid,pch=" ", asp = 1,
     xlab=expression(x[1]), ylab=expression(x[2]))
points(x.tes2[qda$predv==1,1], x.tes2[qda$predv==1,2], col="red")
points(x.tes2[qda$predv==2,1], x.tes2[qda$predv==2,2], col="green")
points(x.tes2[qda$predv==3,1], x.tes2[qda$predv==3,2], col="blue")
points(x.tes2[qda$predv==4,1], x.tes2[qda$predv==4,2], col="black")
points(x.tes2[qda$predv==5,1], x.tes2[qda$predv==5,2], col="orange")

n <- length(qda$predv)
flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==1],phi_vec, 'red',r=0.3)
flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==2],phi_vec, 'green',r=0.3)
flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==3],phi_vec, 'blue',r=0.3)
flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==4],phi_vec, 'black',r=0.3)
flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==5],phi_vec, 'orange',r=0.3)


## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   Real data: Pen digit data                   ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####

## .-------------------------------------------  ####
#####################################################-
##   1.target 0,6,7,8,9                          ####
#####################################################-
## .-------------------------------------------  ####

#################################################-
## Upload the full pen digit data from online####
#################################################-
#install.packages("jsonlite", repos="https://cran.rstudio.com/")
library("jsonlite")

json_file <- 'https://datahub.io/machine-learning/pendigits/datapackage.json'
json_data <- fromJSON(paste(readLines(json_file), collapse=""))

# get list of all resources:
print(json_data$resources$name)

# print all tabular data(if exists any)
for(i in 1:length(json_data$resources$datahub$type)){
  if(json_data$resources$datahub$type[i]=='derived/csv'){
    path_to_file = json_data$resources$path[i]
    data <- read.csv(url(path_to_file))
    print(data)
  }
}

# load the data
path_to_file = json_data$resources$path[7]
data <- read.csv(url(path_to_file))
dim(data) # 10992 x 17
head(data)
str(data)

###############################################-
##  training and test                      ####
###############################################-

train_idx <- 1:7494
train <- data[train_idx,] # m=7494
test <- data[-train_idx,] # n=3498

###############################################-
##  Focus on 0,6,7,8,9 -> reduce the data  ####
###############################################-

library(dplyr)
target <- c(0,6,7,8,9)
train <- filter(train, class %in% target) # m=2219
test <- filter(test, class %in% target)   # n=1035

###############################################-
##  Make x and y                           ####
###############################################-
x_train <- train[,1:16]
y_train <- train[,17]

x_test <- test[,1:16]
y_test <- test[,17]

x.tra <- as.matrix(x_train)
y.tra <- as.matrix(y_train)

x.tes <- as.matrix(x_test)
y.tes <- as.matrix(y_test)

dim(x.tra)
dim(y.tra)

dim(x.tes)
dim(y.tes)

###############################################-
##  DR on training to get 2d               ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,2,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tra) %*% beta
dim(pred)


###############################################-
##  plot training in 2d                    ####
###############################################-
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction") 
# points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
# points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
# points(pred[y.tra==7,1],pred[y.tra==7,2],col="blue")
# points(pred[y.tra==8,1],pred[y.tra==8,2],col="black")
# points(pred[y.tra==9,1],pred[y.tra==9,2],col="orange")
# title(main='wine training in 2d')

# make squaregrid for pretty plot. +-1 to make a room for outer circle
predtemp <- pred
squaregrid <- c(min(range(predtemp[,1])-0.5, range(predtemp[,2])-0.5), 
                max(range(predtemp[,1])+0.5, range(predtemp[,2])+0.5))
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)

## plot
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction") 
plot(squaregrid,squaregrid,pch=" ", asp = 1)
points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
points(pred[y.tra==7,1],pred[y.tra==7,2],col="blue")
points(pred[y.tra==8,1],pred[y.tra==8,2],col="black")
points(pred[y.tra==9,1],pred[y.tra==9,2],col="orange")
#title(main='wine training in 2d')

###############################################-
##  DR on test to get 2d               ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,2,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tes) %*% beta
dim(pred)

###############################################-
##  plot test in 2d                        ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,2,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tes) %*% beta
dim(pred)

## true test data
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction") 
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(pred[y.tes==0,1],pred[y.tes==0,2],col="red")
points(pred[y.tes==6,1],pred[y.tes==6,2],col="green")
points(pred[y.tes==7,1],pred[y.tes==7,2],col="blue")
points(pred[y.tes==8,1],pred[y.tes==8,2],col="black")
points(pred[y.tes==9,1],pred[y.tes==9,2],col="orange")
#title(main='wine test in 2d')

#####################################################-
##  QDA                                          ####
#####################################################-
## prior distribution -------------------------------
(prior <- table(y.tes) / length(y.tes))
prior <- prior %>% round(.,3)


## mu, sigma estimation------------------------------
## assume x follows multivariate normal(MVN) dist'n like QDA
## estimate MVN parameters: mu and sigma

## split the data conditioned on 0,6,7,8,9
pred0 <- pred[y.tes==0,];dim(pred0)
pred6 <- pred[y.tes==6,];dim(pred6)
pred7 <- pred[y.tes==7,];dim(pred7)
pred8 <- pred[y.tes==8,];dim(pred8)
pred9 <- pred[y.tes==9,];dim(pred9)

## define dimension objects
p <- ncol(pred)
n0 <- nrow(pred0)
n6 <- nrow(pred6)
n7 <- nrow(pred7)
n8 <- nrow(pred8)
n9 <- nrow(pred9)

## estimate mu
mu0 <- (1/n0) * ( pred0 %>% colSums() ) # apply(pred0, 2, mean)
mu6 <- (1/n6) * ( pred6 %>% colSums() ) # apply(pred6, 2, mean)
mu7 <- (1/n7) * ( pred7 %>% colSums() ) # apply(pred7, 2, mean)
mu8 <- (1/n8) * ( pred8 %>% colSums() ) # apply(pred8, 2, mean)
mu9 <- (1/n9) * ( pred9 %>% colSums() ) # apply(pred9, 2, mean)
mu <- list(mu0, mu6, mu7, mu8, mu9)


## estimate sigma
sig0 <- matrix(1, p, p)
for(i in 1:n0){
  sig0 <- sig0+(pred0[i,]-mu0) %*% t((pred0[i,]-mu0))
}
sig0 <- (1/(n0-1)) * sig0

sig6 <- matrix(1, p, p)
for(i in 1:n6){
  sig6 <- sig6+(pred6[i,]-mu6) %*% t((pred6[i,]-mu6))
}
sig6 <- (1/(n6-1)) * sig6

sig7 <- matrix(1, p, p)
for(i in 1:n7){
  sig7 <- sig7+(pred7[i,]-mu7) %*% t((pred7[i,]-mu7))
}
sig7 <- (1/(n7-1)) * sig7

sig8 <- matrix(1, p, p)
for(i in 1:n8){
  sig8 <- sig8+(pred8[i,]-mu8) %*% t((pred8[i,]-mu8))
}
sig8 <- (1/(n8-1)) * sig8

sig9 <- matrix(1, p, p)
for(i in 1:n9){
  sig9 <- sig9+(pred9[i,]-mu9) %*% t((pred9[i,]-mu9))
}
sig9 <- (1/(n9-1)) * sig9

sig <- list(sig0, sig6, sig7, sig8, sig9)

## prediction ------------------------------
qda <- br(x = pred, prior = prior, class = c(0,6,7,8,9), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv,y.tes) # confusion matrix

#hist(qda$predv) # histogram of qda

# posterior dist'n
post <- table(qda$predv) / length(qda$predv) 
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y.tes)
wrong_idx

#####################################################-
##  Phi: SGCS                                    ####
#####################################################-
## rho function -------------------------------
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

## ka --------------------------------------
alpha <- 0.05
ka_vec <- rep(NA, nrow(qda$postp))
for(i in 1:length(ka_vec) ){
  ka_idx <- which.max(rho_list[[i]] <= 1-alpha)
  ka_vec[i] <- grid[ka_idx]
}

## phi --------------------------------------
phi_vec <- matrix(NA, nrow(qda$postp), ncol(qda$postp))
for(i in 1:nrow(qda$postp)){
  phi_vec[i,] <- phi(qda$postp[i,])
}



#####################################################-
##  Plot the predctoin                           ####
#####################################################-
## test data prediction
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction", asp = 1) 
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(pred[qda$predv==0,1],pred[qda$predv==0,2],col="red")
points(pred[qda$predv==6,1],pred[qda$predv==6,2],col="green")
points(pred[qda$predv==7,1],pred[qda$predv==7,2],col="blue")
points(pred[qda$predv==8,1],pred[qda$predv==8,2],col="black")
points(pred[qda$predv==9,1],pred[qda$predv==9,2],col="orange")
#title(main='wine predicion in 2d')

dim(y.tes)
dim(qda$predv)
###############################################-
##  steering wheel plot                    ####
###############################################-
## plot with true data
# n <- length(y.tes)
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==0],phi_vec, 'red')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==6],phi_vec, 'green')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==7],phi_vec, 'blue')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==8],phi_vec, 'black')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==9],phi_vec, 'orange')

## plot with test data prediction
n=length(y.tes)

colororder <-  c('red','blue','black','green','orange')
phi_vec_order <- phi_vec[,c(1,3,4,2,5)]

flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==0],
             phi_vec_order, r=0.25,'red',colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==6],
             phi_vec_order, r=0.25,'green',colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==7],
             phi_vec_order, r=0.25,'blue',colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==8],
             phi_vec_order, r=0.25,'black',colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==9],
             phi_vec_order, r=0.25,'orange',colororder)


###############################################-
##  selected prediction case               ####
###############################################-
x0 <- cbind(pred[qda$predv==0,1], pred[qda$predv==0,2]) 
x6 <- cbind(pred[qda$predv==6,1], pred[qda$predv==6,2])
x7 <- cbind(pred[qda$predv==7,1], pred[qda$predv==7,2])
x8 <- cbind(pred[qda$predv==8,1], pred[qda$predv==8,2])
x9 <- cbind(pred[qda$predv==9,1], pred[qda$predv==9,2])

set.seed(434)
n0 <- sample(nrow(x0), 10)
n6 <- sample(nrow(x6), 10)
n7 <- sample(nrow(x7), 10)
n8 <- sample(nrow(x8), 10)
n9 <- sample(nrow(x9), 10)
       
x0 <- x0[n0, ]
x6 <- x6[n6, ]
x7 <- x7[n7, ]
x8 <- x8[n8, ]
x9 <- x9[n9,]
## test data prediction
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction", asp = 1) 
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(x0[,1],x0[,2],col="red")
points(x6[,1],x6[,2],col="green")
points(x7[,1],x7[,2],col="blue")
points(x8[,1],x8[,2],col="black")
points(x9[,1],x9[,2],col="orange")
#title(main='wine selected predicion in 2d')


## plot with true data
# n <- length(y.tes)
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==0],phi_vec, 'red')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==6],phi_vec, 'green')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==7],phi_vec, 'blue')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==8],phi_vec, 'black')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==9],phi_vec, 'orange')

## plot with test data prediction
colororder <-  c('red','blue','black','green','orange')
phi_vec_order <- phi_vec[,c(1,3,4,2,5)]

flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==0][n0],phi_vec_order, 
             r=0.25,'red', colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==6][n6],phi_vec_order, 
             r=0.25,'green', colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==7][n7],phi_vec_order, 
             r=0.25,'blue', colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==8][n8],phi_vec_order, 
             r=0.25,'black', colororder)
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==9][n9],phi_vec_order, 
             r=0.25,'orange', colororder)

## .-------------------------------------------  ####
#####################################################-
##   2.target 0,6,8,9                          ####
#####################################################-
## .-------------------------------------------  ####

#################################################-
## Upload the full pen digit data from online####
#################################################-
#install.packages("jsonlite", repos="https://cran.rstudio.com/")
library("jsonlite")

json_file <- 'https://datahub.io/machine-learning/pendigits/datapackage.json'
json_data <- fromJSON(paste(readLines(json_file), collapse=""))

# get list of all resources:
print(json_data$resources$name)

# print all tabular data(if exists any)
for(i in 1:length(json_data$resources$datahub$type)){
  if(json_data$resources$datahub$type[i]=='derived/csv'){
    path_to_file = json_data$resources$path[i]
    data <- read.csv(url(path_to_file))
    print(data)
  }
}

# load the data
path_to_file = json_data$resources$path[7]
data <- read.csv(url(path_to_file))
dim(data) # 10992 x 17
head(data)
str(data)

###############################################-
##  training and test                      ####
###############################################-

train_idx <- 1:7494
train <- data[train_idx,] # m=7494
test <- data[-train_idx,] # n=3498

###############################################-
##  Focus on 0,6,7,8,9 -> reduce the data  ####
###############################################-

library(dplyr)
target <- c(0,6,8,9)
train <- filter(train, class %in% target) # m=2219
test <- filter(test, class %in% target)   # n=1035

###############################################-
##  Make x and y                           ####
###############################################-
x_train <- train[,1:16]
y_train <- train[,17]

x_test <- test[,1:16]
y_test <- test[,17]

x.tra <- as.matrix(x_train)
y.tra <- as.matrix(y_train)

x.tes <- as.matrix(x_test)
y.tes <- as.matrix(y_test)

dim(x.tra)
dim(y.tra)

dim(x.tes)
dim(y.tes)

###############################################-
##  DR on training to get 2d               ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,2,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tra) %*% beta
dim(pred)


###############################################-
##  plot training in 2d                    ####
###############################################-
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction") 
# points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
# points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
# points(pred[y.tra==7,1],pred[y.tra==7,2],col="blue")
# points(pred[y.tra==8,1],pred[y.tra==8,2],col="blue")
# points(pred[y.tra==9,1],pred[y.tra==9,2],col="black")
# title(main='wine training in 2d')

# make squaregrid for pretty plot. +-1 to make a room for outer circle
predtemp <- pred
squaregrid <- c(min(range(predtemp[,1])-0.5, range(predtemp[,2])-0.5), 
                max(range(predtemp[,1])+0.5, range(predtemp[,2])+0.5))
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)

## plot
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction") 
plot(squaregrid,squaregrid,pch=" ", asp = 1)
points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
points(pred[y.tra==8,1],pred[y.tra==8,2],col="blue")
points(pred[y.tra==9,1],pred[y.tra==9,2],col="black")
#title(main='wine training in 2d')

###############################################-
##  DR on test to get 2d               ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,2,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tes) %*% beta
dim(pred)

###############################################-
##  plot test in 2d                        ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,2,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tes) %*% beta
dim(pred)

## true test data
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction") 
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(pred[y.tes==0,1],pred[y.tes==0,2],col="red")
points(pred[y.tes==6,1],pred[y.tes==6,2],col="green")
points(pred[y.tes==8,1],pred[y.tes==8,2],col="blue")
points(pred[y.tes==9,1],pred[y.tes==9,2],col="black")
#title(main='wine test in 2d')

#####################################################-
##  QDA                                          ####
#####################################################-
## prior distribution -------------------------------
(prior <- table(y.tes) / length(y.tes))
prior <- prior %>% round(.,3)


## mu, sigma estimation------------------------------
## assume x follows multivariate normal(MVN) dist'n like QDA
## estimate MVN parameters: mu and sigma

## split the data conditioned on 0,6,7,8,9
pred0 <- pred[y.tes==0,];dim(pred0)
pred6 <- pred[y.tes==6,];dim(pred6)
pred8 <- pred[y.tes==8,];dim(pred8)
pred9 <- pred[y.tes==9,];dim(pred9)

## define dimension objects
p <- ncol(pred)
n0 <- nrow(pred0)
n6 <- nrow(pred6)
n8 <- nrow(pred8)
n9 <- nrow(pred9)

## estimate mu
mu0 <- (1/n0) * ( pred0 %>% colSums() ) # apply(pred0, 2, mean)
mu6 <- (1/n6) * ( pred6 %>% colSums() ) # apply(pred6, 2, mean)
mu8 <- (1/n8) * ( pred8 %>% colSums() ) # apply(pred8, 2, mean)
mu9 <- (1/n9) * ( pred9 %>% colSums() ) # apply(pred9, 2, mean)
mu <- list(mu0, mu6, mu8, mu9)


## estimate sigma
sig0 <- matrix(1, p, p)
for(i in 1:n0){
  sig0 <- sig0+(pred0[i,]-mu0) %*% t((pred0[i,]-mu0))
}
sig0 <- (1/(n0-1)) * sig0

sig6 <- matrix(1, p, p)
for(i in 1:n6){
  sig6 <- sig6+(pred6[i,]-mu6) %*% t((pred6[i,]-mu6))
}
sig6 <- (1/(n6-1)) * sig6


sig8 <- matrix(1, p, p)
for(i in 1:n8){
  sig8 <- sig8+(pred8[i,]-mu8) %*% t((pred8[i,]-mu8))
}
sig8 <- (1/(n8-1)) * sig8

sig9 <- matrix(1, p, p)
for(i in 1:n9){
  sig9 <- sig9+(pred9[i,]-mu9) %*% t((pred9[i,]-mu9))
}
sig9 <- (1/(n9-1)) * sig9

sig <- list(sig0, sig6, sig8, sig9)

## prediction ------------------------------
qda <- br(x = pred, prior = prior, class = c(0,6,8,9), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv,y.tes) # confusion matrix

#hist(qda$predv) # histogram of qda

# posterior dist'n
post <- table(qda$predv) / length(qda$predv) 
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y.tes)
wrong_idx

#####################################################-
##  Phi: SGCS                                    ####
#####################################################-
## rho function -------------------------------
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

## ka --------------------------------------
alpha <- 0.05
ka_vec <- rep(NA, nrow(qda$postp))
for(i in 1:length(ka_vec) ){
  ka_idx <- which.max(rho_list[[i]] <= 1-alpha)
  ka_vec[i] <- grid[ka_idx]
}

## phi --------------------------------------
phi_vec <- matrix(NA, nrow(qda$postp), ncol(qda$postp))
for(i in 1:nrow(qda$postp)){
  phi_vec[i,] <- phi(qda$postp[i,])
}



#####################################################-
##  Plot the predctoin                           ####
#####################################################-
## test data prediction
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction", asp = 1) 
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(pred[qda$predv==0,1],pred[qda$predv==0,2],col="red")
points(pred[qda$predv==6,1],pred[qda$predv==6,2],col="green")
points(pred[qda$predv==8,1],pred[qda$predv==8,2],col="blue")
points(pred[qda$predv==9,1],pred[qda$predv==9,2],col="black")
#title(main='wine predicion in 2d')

dim(y.tes)
dim(qda$predv)
###############################################-
##  steering wheel plot                    ####
###############################################-
## plot with true data
# n <- length(y.tes)
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==0],phi_vec, 'red')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==6],phi_vec, 'green')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==7],phi_vec, 'blue')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==8],phi_vec, 'blue')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==9],phi_vec, 'black')

## plot with test data prediction
n=length(y.tes)

# colororder <-  c('red','green','blue','black')
# phi_vec_order <- phi_vec[,c(1,3,4,2,5)]

flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==0],phi_vec,r=0.25,'red')
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==6],phi_vec,r=0.25,'green')
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==8],phi_vec,r=0.25,'blue')
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==9],phi_vec,r=0.25,'black')


###############################################-
##  selected prediction case               ####
###############################################-
x0 <- cbind(pred[qda$predv==0,1], pred[qda$predv==0,2]) 
x6 <- cbind(pred[qda$predv==6,1], pred[qda$predv==6,2])
x8 <- cbind(pred[qda$predv==8,1], pred[qda$predv==8,2])
x9 <- cbind(pred[qda$predv==9,1], pred[qda$predv==9,2])

set.seed(434)
n0 <- sample(nrow(x0), 10)
n6 <- sample(nrow(x6), 10)
n8 <- sample(nrow(x8), 10)
n9 <- sample(nrow(x9), 10)

x0 <- x0[n0, ]
x6 <- x6[n6, ]
x8 <- x8[n8, ]
x9 <- x9[n9,]
## test data prediction
# par(pty="s") # to make the plot region square (independent of size and limits)
# plot(pred[,1],pred[,2],pch=" ",xlab="first DR direction",
#      ylab="second DR direction", asp = 1) 
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(x0[,1],x0[,2],col="red")
points(x6[,1],x6[,2],col="green")
points(x8[,1],x8[,2],col="blue")
points(x9[,1],x9[,2],col="black")
#title(main='wine selected predicion in 2d')


## plot with true data
# n <- length(y.tes)
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==0],phi_vec, 'red')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==6],phi_vec, 'green')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==7],phi_vec, 'blue')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==8],phi_vec, 'blue')
# flowerplot2d(cbind(pred,y.tes),(1:n)[y.tes==9],phi_vec, 'black')

## plot with test data prediction
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==0][n0],phi_vec,r=0.25,'red')
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==6][n6],phi_vec,r=0.25,'green')
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==8][n8],phi_vec,r=0.25,'blue')
flowerplot2d(cbind(pred,y.tes),(1:n)[qda$predv==9][n9],phi_vec,r=0.25,'black')

## .-------------------------------------------  ####
#####################################################-
##   3. 3D/ target 0,6,7,8,9                     ####
#####################################################-
## .-------------------------------------------  ####

#################################################-
## Edit perspective plot functions           ####
#################################################-

#####################################################################
#        function: plot one point                               
#        input v: 3-d vector representing a single point           
#        th: angle of the cube. typical angle is 30                
#        color: 1 = red, 2 = green, 3 = blue
#####################################################################
onept=function(v,th,color){
  w = onepos(v,th)
    points(w[1],w[2],col=color)
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


#####################################################-
##   flower plot                                 ####
#####################################################-
# #onebarplot(gcrstd[p1,], 30, phi_vec[p1,])
# p1 <- phi1_idx[1]
# v=c(gcrstd[p1,])
# phi <- phi_vec[p1,]
# th <- 30

## input 
#' @param v: 3-d vector representing a single point 
#' @param th: angle of the cube. typical angle is 30  
#' @param phi : phi values= generalized credible set values
#' @param r: radius
#' @param color : color for center point
#' @param colorset : colors for blades

oneflower3dplot <- function(v, th, phi, r=0.5,color, 
                            colorset = c('red','green','blue','black','orange')){
  z <- v
  zz <- z
  w <- onepos(zz, th) # 2d coordinates / original center
  
  ## 1,2,3 prob (0.2 weight on prob)
  # 1
  zz <- z
  zz[3] <- zz[3]+r*phi[1]
  n <- length(phi)
  # h1 <- onepos(zz, th)
  # 
  # bconnect(w, h1, 'red')
  # 
  # # 2
  # zz <- z
  # zz[3] <- zz[3]+r*phi[2]
  # h2 <- onepos(zz, th)
  # v2 <- h2-w
  # 
  # # rotate
  # theta <- (2*pi) * (1/3)
  # r2 <- rotate(v2, theta)
  # 
  # # translate
  # t2 <- w+r2
  # bconnect(w, t2,'green')
  # 
  # 
  # 
  # # 3
  # zz <- z
  # zz[3] <- zz[3]+r*phi[3]
  # h3 <- onepos(zz, th)
  # v3 <- h3-w
  # 
  # # rotate
  # theta <- (2*pi) * (2/3)
  # r3 <- rotate(v3, theta)
  # 
  # # translate
  # t3 <- w+r3
  # bconnect(w, t3, 'blue')
  # 
  # the rest
  for(i in 1:n){
    zz <- z
    zz[3] <- zz[3]+r*phi[i]
    hi <- onepos(zz,th)
    vi <- hi-w
    
    # rotate
    theta <- (2*pi) * ((i-1)/n)
    ri <- rotate(vi, theta)
    
    # translate
    ti <- w+ri
    
    col <- colorset[i]
    bconnect(w, ti, col)
  }
  
  # color point
  #conept(v, th, color=color)
  points(w[1],w[2],col=color, pch=16)
  points(w[1],w[2])
  
  # draw exterior circle using plotrix package
  draw.circle(w[1],w[2], radius = r*1.05)
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
#' @param r: radius
#' @param color : color for center point
#' @param colorset : colors for blades

flowerplot3d = function(x,h1=1.5,h2=1,h3=1.3,th=30,index,phi,r=0.5,color,
                        colorset = c('red','green','blue','black','orange')){
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
    oneflower3dplot(v,th,p,r,color,colorset)
  }
}



#################################################-
## Upload the full pen digit data from online####
#################################################-
#install.packages("jsonlite", repos="https://cran.rstudio.com/")
library("jsonlite")

json_file <- 'https://datahub.io/machine-learning/pendigits/datapackage.json'
json_data <- fromJSON(paste(readLines(json_file), collapse=""))

# get list of all resources:
print(json_data$resources$name)

# print all tabular data(if exists any)
for(i in 1:length(json_data$resources$datahub$type)){
  if(json_data$resources$datahub$type[i]=='derived/csv'){
    path_to_file = json_data$resources$path[i]
    data <- read.csv(url(path_to_file))
    print(data)
  }
}

# load the data
path_to_file = json_data$resources$path[7]
data <- read.csv(url(path_to_file))
dim(data) # 10992 x 17
head(data)
str(data)

###############################################-
##  training and test                      ####
###############################################-

train_idx <- 1:7494
train <- data[train_idx,] # m=7494
test <- data[-train_idx,] # n=3498

###############################################-
##  Focus on 0,6,7,8,9 -> reduce the data  ####
###############################################-

library(dplyr)
target <- c(0,6,7,8,9)
train <- filter(train, class %in% target) # m=3716
test <- filter(test, class %in% target)   # n=1735

###############################################-
##  Make x and y                           ####
###############################################-
x_train <- train[,1:16]
y_train <- train[,17]

x_test <- test[,1:16]
y_test <- test[,17]

x.tra <- as.matrix(x_train)
y.tra <- as.matrix(y_train)

x.tes <- as.matrix(x_test)
y.tes <- as.matrix(y_test)

dim(x.tra)
dim(y.tra)

dim(x.tes)
dim(y.tes)

###############################################-
##  DR on training to get 3d               ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,3,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tra) %*% beta
dim(pred)


###############################################-
##  plot training in 3d                    ####
###############################################-
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tra),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tra==0],color = "red")
data3d(cbind(pred,y.tra),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tra==6],color = "green")
data3d(cbind(pred,y.tra),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tra==7],color = "blue")
data3d(cbind(pred,y.tra),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tra==8],color = "black")
data3d(cbind(pred,y.tra),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tra==9],color = "orange")

###############################################-
##  DR on test to get 3d               ####
###############################################-
## beta
beta <- dr(x.tra,y.tra,3,3,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x.tes) %*% beta
dim(pred)

###############################################-
##  plot test in 3d                        ####
###############################################-
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==0],color = "red")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==6],color = "green")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==7],color = "blue")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==8],color = "black")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==9],color = "orange")


#####################################################-
##  QDA                                          ####
#####################################################-
## prior distribution -------------------------------
(prior <- table(y.tes) / length(y.tes))
prior <- prior %>% round(.,3)


## mu, sigma estimation------------------------------
## assume x follows multivariate normal(MVN) dist'n like QDA
## estimate MVN parameters: mu and sigma

## split the data conditioned on 0,6,7,8,9
pred0 <- pred[y.tes==0,];dim(pred0)
pred6 <- pred[y.tes==6,];dim(pred6)
pred7 <- pred[y.tes==7,];dim(pred7)
pred8 <- pred[y.tes==8,];dim(pred8)
pred9 <- pred[y.tes==9,];dim(pred9)

## define dimension objects
p <- ncol(pred)
n0 <- nrow(pred0)
n6 <- nrow(pred6)
n7 <- nrow(pred7)
n8 <- nrow(pred8)
n9 <- nrow(pred9)

## estimate mu
mu0 <- (1/n0) * ( pred0 %>% colSums() ) # apply(pred0, 2, mean)
mu6 <- (1/n6) * ( pred6 %>% colSums() ) # apply(pred6, 2, mean)
mu7 <- (1/n7) * ( pred7 %>% colSums() ) # apply(pred7, 2, mean)
mu8 <- (1/n8) * ( pred8 %>% colSums() ) # apply(pred8, 2, mean)
mu9 <- (1/n9) * ( pred9 %>% colSums() ) # apply(pred9, 2, mean)
mu <- list(mu0, mu6, mu7, mu8, mu9)


## estimate sigma
sig0 <- matrix(1, p, p)
for(i in 1:n0){
  sig0 <- sig0+(pred0[i,]-mu0) %*% t((pred0[i,]-mu0))
}
sig0 <- (1/(n0-1)) * sig0

sig6 <- matrix(1, p, p)
for(i in 1:n6){
  sig6 <- sig6+(pred6[i,]-mu6) %*% t((pred6[i,]-mu6))
}
sig6 <- (1/(n6-1)) * sig6

sig7 <- matrix(1, p, p)
for(i in 1:n7){
  sig7 <- sig7+(pred7[i,]-mu7) %*% t((pred7[i,]-mu7))
}
sig7 <- (1/(n7-1)) * sig7

sig8 <- matrix(1, p, p)
for(i in 1:n8){
  sig8 <- sig8+(pred8[i,]-mu8) %*% t((pred8[i,]-mu8))
}
sig8 <- (1/(n8-1)) * sig8

sig9 <- matrix(1, p, p)
for(i in 1:n9){
  sig9 <- sig9+(pred9[i,]-mu9) %*% t((pred9[i,]-mu9))
}
sig9 <- (1/(n9-1)) * sig9

sig <- list(sig0, sig6, sig7, sig8, sig9)

## prediction ------------------------------
qda <- br(x = pred, prior = prior, class = c(0,6,7,8,9), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv,y.tes) # confusion matrix

#hist(qda$predv) # histogram of qda

# posterior dist'n
post <- table(qda$predv) / length(qda$predv) 
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y.tes)
wrong_idx

#####################################################-
##  Phi: SGCS                                    ####
#####################################################-
## rho function -------------------------------
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

## ka --------------------------------------
alpha <- 0.05
ka_vec <- rep(NA, nrow(qda$postp))
for(i in 1:length(ka_vec) ){
  ka_idx <- which.max(rho_list[[i]] <= 1-alpha)
  ka_vec[i] <- grid[ka_idx]
}

## phi --------------------------------------
phi_vec <- matrix(NA, nrow(qda$postp), ncol(qda$postp))
for(i in 1:nrow(qda$postp)){
  phi_vec[i,] <- phi(qda$postp[i,])
}



#####################################################-
##  Plot the prediction in 3d                    ####
#####################################################-
## test data prediction
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==0],color = "red")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==6],color = "green")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==7],color = "blue")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==8],color = "black")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==9],color = "orange")

dim(y.tes)
dim(qda$predv)
###############################################-
##  steering wheel plot in 3d              ####
###############################################-
## plot with test data prediction
n=length(y.tes)

colororder <-  c('orange','red','blue','black','green')
phi_vec_order <- phi_vec[,c(5,1,3,4,2)]

flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==0],
             phi_vec_order,r=0.1,'red',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==6],
             phi_vec_order,r=0.1,'green',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==7],
             phi_vec_order,r=0.1,'blue',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==8],
             phi_vec_order,r=0.1,'black',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==9],
             phi_vec_order,r=0.1,'orange',colororder)

###############################################-
##  selected prediction case               ####
###############################################-
preds <- cbind(pred, y.tes, qda$predv)
x0 <- cbind(pred[qda$predv==0,1], pred[qda$predv==0,2]) 
x6 <- cbind(pred[qda$predv==6,1], pred[qda$predv==6,2])
x7 <- cbind(pred[qda$predv==7,1], pred[qda$predv==7,2])
x8 <- cbind(pred[qda$predv==8,1], pred[qda$predv==8,2])
x9 <- cbind(pred[qda$predv==9,1], pred[qda$predv==9,2])

set.seed(434)
n0 <- sample(nrow(x0), 10)
n6 <- sample(nrow(x6), 10)
n7 <- sample(nrow(x7), 10)
n8 <- sample(nrow(x8), 10)
n9 <- sample(nrow(x9), 10)

x0 <- x0[n0, ]
x6 <- x6[n6, ]
x7 <- x7[n7, ]
x8 <- x8[n8, ]
x9 <- x9[n9,]

## test data
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==0][n0],color = "red")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==6][n6],color = "green")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==7][n7],color = "blue")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==8][n8],color = "black")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==9][n9],color = "orange")

## test data prediction
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==0][n0],color = "red")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==6][n6],color = "green")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==7][n7],color = "blue")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==8][n8],color = "black")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==9][n9],color = "orange")

## plot with test data prediction
colororder <-  c('orange','red','blue','black','green')
phi_vec_order <- phi_vec[,c(5,1,3,4,2)]

flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==0][n0],
             phi_vec_order,r=0.1,'red',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==6][n6],
             phi_vec_order,r=0.1,'green',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==7][n7],
             phi_vec_order,r=0.1,'blue',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==8][n8],
             phi_vec_order,r=0.1,'black',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==9][n8],
             phi_vec_order,r=0.1,'orange',colororder)

###############################################-
##  selected prediction case               ####
###############################################-
preds <- cbind(pred, y.tes, qda$predv)
x0 <- cbind(pred[qda$predv==0,1], pred[qda$predv==0,2]) 
x6 <- cbind(pred[qda$predv==6,1], pred[qda$predv==6,2])
x7 <- cbind(pred[qda$predv==7,1], pred[qda$predv==7,2])
x8 <- cbind(pred[qda$predv==8,1], pred[qda$predv==8,2])
x9 <- cbind(pred[qda$predv==9,1], pred[qda$predv==9,2])

set.seed(434)
n0 <- sample(nrow(x0), 10)
n6 <- sample(nrow(x6), 10)
n7 <- sample(nrow(x7), 10)
n8 <- sample(nrow(x8), 10)
n9 <- sample(nrow(x9), 10)

x0 <- x0[n0, ]
x6 <- x6[n6, ]
x7 <- x7[n7, ]
x8 <- x8[n8, ]
x9 <- x9[n9,]

## test data
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==0][n0],color = "red")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==6][n6],color = "green")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==7][n7],color = "blue")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==8][n8],color = "black")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==9][n9],color = "orange")

## test data prediction
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==0][n0],color = "red")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==6][n6],color = "green")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==7][n7],color = "blue")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==8][n8],color = "black")
data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==9][n9],color = "orange")

## plot with test data prediction
colororder <-  c('orange','red','blue','black','green')
phi_vec_order <- phi_vec[,c(5,1,3,4,2)]

flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==0][n0],
             phi_vec_order,r=0.1,'red',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==6][n6],
             phi_vec_order,r=0.1,'green',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==7][n7],
             phi_vec_order,r=0.1,'blue',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==8][n8],
             phi_vec_order,r=0.1,'black',colororder)
flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[qda$predv==9][n8],
             phi_vec_order,r=0.1,'orange',colororder)


###############################################-
##  selected prediction case               ####
###############################################-
#phi_vec_order <- phi_vec[,c(5,1,3,4,2)]
preds <- cbind(pred, y.tes, qda$predv, phi_vec)
x0 <- preds[y.tes==0,]
x6 <- preds[y.tes==6,]
x7 <- preds[y.tes==7,]
x8 <- preds[y.tes==8,]
x9 <- preds[y.tes==9,]

# x0 <- cbind(preds[y.tes==0,1], preds[y.tes==0,2]) 
# x6 <- cbind(preds[y.tes==6,1], preds[y.tes==6,2])
# x7 <- cbind(preds[y.tes==7,1], preds[y.tes==7,2])
# x8 <- cbind(preds[y.tes==8,1], preds[y.tes==8,2])
# x9 <- cbind(preds[y.tes==9,1], preds[y.tes==9,2])

set.seed(434)
x0 <- x0[sample(nrow(x0), 10), ]
x6 <- x6[sample(nrow(x6), 10), ]
x7 <- x7[sample(nrow(x7), 10), ]
x8 <- x8[sample(nrow(x8), 10), ]
x9 <- x9[sample(nrow(x8), 10),]
x <- rbind(x0,x6,x7,x8,x9)

nx <- nrow(x)

## test data
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,4]==0],color = "red")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,4]==6],color = "green")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,4]==7],color = "blue")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,4]==8],color = "black")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,4]==9],color = "orange")

## test data prediction
frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==0],color = "red")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==6],color = "green")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==7],color = "blue")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==8],color = "black")
data3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==9],color = "orange")

## plot with test data prediction
colororder <-  c('orange','red','blue','black','green')
x_phi_vec <- x[,c(6:10)]
x_phi_vec_order <- x_phi_vec[,c(5,1,3,4,2)]

flowerplot3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==0],
             x_phi_vec_order,r=0.1,'red',colororder)
flowerplot3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==6],
             x_phi_vec_order,r=0.1,'green',colororder)
flowerplot3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==7],
             x_phi_vec_order,r=0.1,'blue',colororder)
flowerplot3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==8],
             x_phi_vec_order,r=0.1,'black',colororder)
flowerplot3d(x,h1=1.3,h2=1.5,h3=1.2,th=20,(1:nx)[x[,5]==9],
             x_phi_vec_order,r=0.1,'orange',colororder)



## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   2D view                                     ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
## select 2 directions
pred12 <- x[,c(1,2)]
pred13 <- x[,c(1,3)]
pred23 <- x[,c(2,3)]
# 
# set.seed(434)
# n0 <- sample(nrow(x0), 10)
# n6 <- sample(nrow(x6), 10)
# n7 <- sample(nrow(x7), 10)
# n8 <- sample(nrow(x8), 10)
# n9 <- sample(nrow(x9), 10)

## 1 vs 2 view
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(x[x[,5]==0,1],x[x[,5]==0,2],col="red")
points(x[x[,5]==6,1],x[x[,5]==6,2],col="green")
points(x[x[,5]==7,1],x[x[,5]==7,2],col="blue")
points(x[x[,5]==8,1],x[x[,5]==8,2],col="black")
points(x[x[,5]==9,1],x[x[,5]==9,2],col="orange")

# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==0][n0],phi_vec,r=0.25,'red')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==6][n6],phi_vec,r=0.25,'green')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==7][n7],phi_vec,r=0.25,'blue')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==8][n8],phi_vec,r=0.25,'black')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==9][n9],phi_vec,r=0.25,'orange')

arrange_phi_vec <- x_phi_vec_order[,c(1,3,4,2,5)]
arrange_colorset <- c('red','blue','black','green','orange')
flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==0][n0],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==6][n6],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==7][n7],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==8][n8],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==9][n9],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)

## 1 vs 3 view
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="third DR direction",
     asp = 1)
points(pred[qda$predv==0,1][n0],pred[qda$predv==0,3][n0],col="red")
points(pred[qda$predv==6,1][n6],pred[qda$predv==6,3][n6],col="green")
points(pred[qda$predv==7,1][n7],pred[qda$predv==7,3][n7],col="blue")
points(pred[qda$predv==8,1][n8],pred[qda$predv==8,3][n8],col="black")
points(pred[qda$predv==9,1][n9],pred[qda$predv==9,3][n9],col="orange")


arrange_phi_vec <- phi_vec[,c(4,1,3,5,2)]
arrange_colorset <- c('black','red','blue','orange','green')
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==0][n0],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==6][n6],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==7][n7],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==8][n8],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==9][n9],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)

## 2 vs 3 view
plot(squaregrid,squaregrid,pch=" ",
     xlab="second DR direction", ylab="third DR direction",
     asp = 1)
points(pred[qda$predv==0,2][n0],pred[qda$predv==0,3][n0],col="red")
points(pred[qda$predv==6,2][n6],pred[qda$predv==6,3][n6],col="green")
points(pred[qda$predv==7,2][n7],pred[qda$predv==7,3][n7],col="blue")
points(pred[qda$predv==8,2][n8],pred[qda$predv==8,3][n8],col="black")
points(pred[qda$predv==9,2][n9],pred[qda$predv==9,3][n9],col="orange")

arrange_phi_vec <- phi_vec[,c(4,2,3,5,1)]
arrange_colorset <- c('black','green','blue','orange','red')
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==0][n0],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==6][n6],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==7][n7],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==8][n8],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==9][n9],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)



## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   2D view                                     ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
## select 2 directions
pred12 <- x[,c(1,2)]
pred13 <- x[,c(1,3)]
pred23 <- x[,c(2,3)]
# 
# set.seed(434)
# n0 <- sample(nrow(x0), 10)
# n6 <- sample(nrow(x6), 10)
# n7 <- sample(nrow(x7), 10)
# n8 <- sample(nrow(x8), 10)
# n9 <- sample(nrow(x9), 10)

## 1 vs 2 view
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(x[x[,5]==0,1],x[x[,5]==0,2],col="red")
points(x[x[,5]==6,1],x[x[,5]==6,2],col="green")
points(x[x[,5]==7,1],x[x[,5]==7,2],col="blue")
points(x[x[,5]==8,1],x[x[,5]==8,2],col="black")
points(x[x[,5]==9,1],x[x[,5]==9,2],col="orange")

# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==0][n0],phi_vec,r=0.25,'red')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==6][n6],phi_vec,r=0.25,'green')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==7][n7],phi_vec,r=0.25,'blue')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==8][n8],phi_vec,r=0.25,'black')
# flowerplot2d(cbind(pred12,y.tes),(1:n)[qda$predv==9][n9],phi_vec,r=0.25,'orange')

arrange_phi_vec <- x_phi_vec[,c(1,3,4,2,5)]
arrange_colorset <- c('red','blue','black','green','orange')
flowerplot2d(x,(1:nx)[x[,5]==0],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(x,(1:nx)[x[,5]==6],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(x,(1:nx)[x[,5]==7],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(x,(1:nx)[x[,5]==8],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(x,(1:nx)[x[,5]==9],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)

## 1 vs 3 view
plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="third DR direction",
     asp = 1)
points(pred[qda$predv==0,1][n0],pred[qda$predv==0,3][n0],col="red")
points(pred[qda$predv==6,1][n6],pred[qda$predv==6,3][n6],col="green")
points(pred[qda$predv==7,1][n7],pred[qda$predv==7,3][n7],col="blue")
points(pred[qda$predv==8,1][n8],pred[qda$predv==8,3][n8],col="black")
points(pred[qda$predv==9,1][n9],pred[qda$predv==9,3][n9],col="orange")


arrange_phi_vec <- phi_vec[,c(4,1,3,5,2)]
arrange_colorset <- c('black','red','blue','orange','green')
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==0][n0],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==6][n6],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==7][n7],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==8][n8],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(cbind(pred13,y.tes),(1:n)[qda$predv==9][n9],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)

## 2 vs 3 view
plot(squaregrid,squaregrid,pch=" ",
     xlab="second DR direction", ylab="third DR direction",
     asp = 1)
points(pred[qda$predv==0,2][n0],pred[qda$predv==0,3][n0],col="red")
points(pred[qda$predv==6,2][n6],pred[qda$predv==6,3][n6],col="green")
points(pred[qda$predv==7,2][n7],pred[qda$predv==7,3][n7],col="blue")
points(pred[qda$predv==8,2][n8],pred[qda$predv==8,3][n8],col="black")
points(pred[qda$predv==9,2][n9],pred[qda$predv==9,3][n9],col="orange")

arrange_phi_vec <- phi_vec[,c(4,2,3,5,1)]
arrange_colorset <- c('black','green','blue','orange','red')
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==0][n0],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==6][n6],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==7][n7],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==8][n8],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(cbind(pred23,y.tes),(1:n)[qda$predv==9][n9],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)



#####################################################-
##  runcode: 3d plot for DR for whole test set   ####
#####################################################-
#pred=center(x.tra)%*%method("dr");n=length(y.tes) # sufficient pred
# frame3d(h1=1.3,h2=1.5,h3=1.2,th=20,"dr")
# data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==0],color = "red")
# data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==6],color = "green")
# data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==7],color = "blue")
# data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==8],color = "black")
# data3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==9],color = "orange")
# 
# flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==0],
#              phi_vec,r=0.1,'red')
# flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==6],
#              phi_vec,r=0.1,'green')
# flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==7],
#              phi_vec,r=0.1,'blue')
# flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==8],
#              phi_vec,r=0.1,'black')
# flowerplot3d(cbind(pred,y.tes),h1=1.3,h2=1.5,h3=1.2,th=20,(1:n)[y.tes==9],
#              phi_vec,r=0.1,'orange')



## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   Belows are not used                         ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####

#####################################################-
##   runcode: 3d plot for DR for training set    ####
#####################################################-
# #pred=center(x.tra)%*%method("dr");n=length(y.tes) # sufficient pred
# frame3d(,,,,"dr")
# data3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==0],color = "red")
# data3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==6],color = "green")
# data3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==7],color = "blue")
# data3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==8],color = "black")
# data3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==9],color = "orange")

# flowerplot3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==0],phi_vec,r=0.1,'red')
# flowerplot3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==6],phi_vec,r=0.1,'green')
# flowerplot3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==7],phi_vec,r=0.1,'blue')
# flowerplot3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==8],phi_vec,r=0.1,'black')
# flowerplot3d(cbind(pred,y.tra),,,,,(1:n)[y.tra==9],phi_vec,r=0.1,'orange')

## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   Steering wheel plot                         ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
# plot(squaregrid,squaregrid,pch=" ", asp = 1)
# points(x.tes2[qda$predv==1,1], x.tes2[qda$predv==1,2], col="red")
# points(x.tes2[qda$predv==2,1], x.tes2[qda$predv==2,2], col="green")
# points(x.tes2[qda$predv==3,1], x.tes2[qda$predv==3,2], col="blue")
# points(x.tes2[qda$predv==4,1], x.tes2[qda$predv==4,2], col="black")
# points(x.tes2[qda$predv==5,1], x.tes2[qda$predv==5,2], col="orange")
# 
# n <- length(qda$predv)
# flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==1],phi_vec, 'red',r=0.3)
# flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==2],phi_vec, 'green',r=0.3)
# flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==3],phi_vec, 'blue',r=0.3)
# flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==4],phi_vec, 'black',r=0.3)
# flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==5],phi_vec, 'orange',r=0.3)
# 
# 
# flowerplot2d(cbind(x.tes2,y.tes2),(1:n)[qda$predv==1],phi_vec, 'red',r=0.3)
# flowerplot2d(cbind(c(5,5),1),(1:n)[qda$predv==1],t(c(1,1,0.5,1,0)), 'purple',r=0.3)

#save.image("credible-set.Rdata")
