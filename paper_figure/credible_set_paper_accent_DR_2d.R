#####################################################-
##  Packages                                     ####
#####################################################-
library(dplyr)
library(ggplot2)
library(mvtnorm)

## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   Real data: Speaker Accent                   ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####

#setwd('/Users/chaegeunsong/GitHub/credible-set')
## load the data
data <- read.csv('/Users/chaegeunsong/GitHub/credible-set/data/speakeraccent/accent-mfcc-data-1.csv',
                 header = TRUE)
str(data) 
table(data$language)
# 329 obs. of 13 varaibles. 
# 6 classes: ES, FR, GE, IT, UK, US.
#summary(data)

## rename the language variable as y
data <- rename(data, y=language)
table(data$y)

## convert y into numeric
data$y <- data$y %>% as.matrix() %>% as.factor() %>% as.numeric()
#data$y <- y %>% unlist(.) %>% as.numeric()

## data matrix x and y
y <- select(data, y) %>% as.matrix()
x <- select(data, -y) %>% as.matrix()

## split the data by accent class
k <- length(unique(data$y)) # number of classes

for(i in 1:k){
  assign(paste0("x", i), split(data, data$y)[[i]])
}
# x1,...,x6


###############################################-
##  DR on training to get 2d               ####
###############################################-
## beta
beta <- dr(x,y,3,4,"categorical")
dim(beta)

## sufficient predictor
pred <- center(x) %*% beta
dim(pred)

###############################################-
##  plot data in 2d                        ####
###############################################-
predtemp <- pred
squaregrid <- c(min(range(predtemp[,1])-0.5, range(predtemp[,2])-0.5), 
                max(range(predtemp[,1])+0.5, range(predtemp[,2])+0.5))

plot(squaregrid,squaregrid,pch=" ", asp = 1)
points(pred[y==1,1],pred[y==1,2],col="red")
points(pred[y==2,1],pred[y==2,2],col="green")
points(pred[y==3,1],pred[y==3,2],col="blue")
points(pred[y==4,1],pred[y==4,2],col="black")
points(pred[y==5,1],pred[y==5,2],col="orange")
points(pred[y==6,1],pred[y==6,2],col="purple")
#title(main='')

#####################################################-
##  QDA                                          ####
#####################################################-
## prior distribution -------------------------------
(prior <- table(y) / length(y))
prior <- prior %>% round(.,3)

## mu, sigma estimation------------------------------
## assume x follows multivariate normal(MVN) dist'n like QDA
## estimate MVN parameters: mu and sigma

## split the data conditioned on 1,2,3,4,5,6
for(i in 1:k){
  assign(paste0("pred", i), pred[y==i,])
}


## define dimension objects
p <- ncol(pred)
n <- nrow(pred)

n1 <- nrow(pred1)
n2 <- nrow(pred2)
n3 <- nrow(pred3)
n4 <- nrow(pred4)
n5 <- nrow(pred5)
n6 <- nrow(pred6)

## estimate mu
mu1 <- (1/n1) * ( pred1 %>% colSums() ) # apply(pred1, 2, mean)
mu2 <- (1/n2) * ( pred2 %>% colSums() ) # apply(pred2, 2, mean)
mu3 <- (1/n3) * ( pred3 %>% colSums() ) # apply(pred3, 2, mean)
mu4 <- (1/n4) * ( pred4 %>% colSums() ) # apply(pred4, 2, mean)
mu5 <- (1/n5) * ( pred5 %>% colSums() ) # apply(pred5, 2, mean)
mu6 <- (1/n6) * ( pred6 %>% colSums() ) # apply(pred6, 2, mean)
mu <- list(mu1, mu2, mu3, mu4, mu5, mu6)


## estimate sigma
sig1 <- matrix(1, p, p)
for(i in 1:n1){
  sig1 <- sig1+(pred1[i,]-mu1) %*% t((pred1[i,]-mu1))
}
sig1 <- (1/(n1-1)) * sig1

sig2 <- matrix(1, p, p)
for(i in 1:n2){
  sig2 <- sig2+(pred2[i,]-mu2) %*% t((pred2[i,]-mu2))
}
sig2 <- (1/(n2-1)) * sig2

sig3 <- matrix(1, p, p)
for(i in 1:n3){
  sig3 <- sig3+(pred3[i,]-mu3) %*% t((pred3[i,]-mu3))
}
sig3 <- (1/(n3-1)) * sig3

sig4 <- matrix(1, p, p)
for(i in 1:n4){
  sig4 <- sig4+(pred4[i,]-mu4) %*% t((pred4[i,]-mu4))
}
sig4 <- (1/(n4-1)) * sig4

sig5 <- matrix(1, p, p)
for(i in 1:n5){
  sig5 <- sig5+(pred5[i,]-mu5) %*% t((pred5[i,]-mu5))
}
sig5 <- (1/(n5-1)) * sig5

sig6 <- matrix(1, p, p)
for(i in 1:n6){
  sig6 <- sig6+(pred6[i,]-mu6) %*% t((pred6[i,]-mu6))
}
sig6 <- (1/(n6-1)) * sig6

sig <- list(sig1, sig2, sig3, sig4, sig5, sig6)

## prediction ------------------------------
qda <- br(x = pred, prior = prior, class = c(1,2,3,4,5,6), mu = mu, sig = sig)
head(qda$predv)

table(qda$predv)
table(qda$predv,y) # confusion matrix
diag(table(qda$predv,y)) / table(y) # accuracy rate

# posterior dist'n
post <- table(qda$predv) / length(qda$predv) 
post %>% round(.,3)

## posterior probability
qda$postp %>% head()
# Most of them cosists of 0 and 1 because it is easy to 
# do clustering. 

## wrong clustering
wrong_idx <- which(qda$predv != y)
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



## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   2D view                                     ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
# define objects needed
n <- nrow(pred)
predy <- cbind(pred,y)
v <- qda$predv
table(v)
datav <- cbind(data,v)
data$v <- datav$v %>% as.matrix() %>% as.factor() %>% as.numeric()

## 1 vs 2 view
squaregrid <- c(min(range(pred[,1])-0.5, range(pred[,2])-0.5), 
                max(range(pred[,1])+0.5, range(pred[,2])+0.5))

plot(squaregrid,squaregrid,pch=" ",
     xlab="first DR direction", ylab="second DR direction",
     asp = 1)
points(pred[y==1,1],pred[y==1,2],col="red")
points(pred[y==2,1],pred[y==2,2],col="green")
points(pred[y==3,1],pred[y==3,2],col="blue")
points(pred[y==4,1],pred[y==4,2],col="black")
points(pred[y==5,1],pred[y==5,2],col="orange")
points(pred[y==6,1],pred[y==6,2],col="purple")

arrange_colorset <- c('green','orange','black','blue','red','purple')
arrange_phi_vec <- phi_vec[,c(2,5,4,3,1,6)]
flowerplot2d(predy,(1:n)[v==1],arrange_phi_vec,r=0.25,
             'red',arrange_colorset)
flowerplot2d(predy,(1:n)[v==2],arrange_phi_vec,r=0.25,
             'green',arrange_colorset)
flowerplot2d(predy,(1:n)[v==3],arrange_phi_vec,r=0.25,
             'blue',arrange_colorset)
flowerplot2d(predy,(1:n)[v==4],arrange_phi_vec,r=0.25,
             'black',arrange_colorset)
flowerplot2d(predy,(1:n)[v==5],arrange_phi_vec,r=0.25,
             'orange',arrange_colorset)
flowerplot2d(predy,(1:n)[v==6],arrange_phi_vec,r=0.25,
             'purple',arrange_colorset)


## .-------------------------------------------  ####
## .-------------------------------------------  ####
#####################################################-
##   Total 20                                    ####
#####################################################-
## .-------------------------------------------  ####
## .-------------------------------------------  ####
## select 10 from each class except class=4=Italian
# for(i in 1:k){
#   assign(paste0("v", i), split(datav, datav$v)[[i]])
# }
table(v)

v1 <- split(datav, datav$v)[[1]]
v2 <- split(datav, datav$v)[[2]]
v3 <- split(datav, datav$v)[[3]]
v4 <- split(datav, datav$v)[[4]]
v4 <- split(datav, datav$v)[[5]]
v6 <- split(datav, datav$v)[[6]]

set.seed(434)
r <- 5
n1 <- sample(nrow(v1), r)
n2 <- sample(nrow(v2), r)
n3 <- sample(nrow(v3), r)
n4 <- sample(nrow(v4), r)
n5 <- sample(nrow(v5), r)
n6 <- sample(nrow(v6), r)

## 1 vs 2 view
squaregrid <- c(min(range(pred[,1]), range(pred[,2])), 
                max(range(pred[,1])+0.2, range(pred[,2])+0.2))

(max(squaregrid)-min(squaregrid))/0.3
cairo_ps("accent_DR.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 2) + 0.1)
plot(squaregrid,squaregrid,pch=" ",
     xlab="First direction", ylab="Second direction",
     asp = 1,
     cex=2, cex.lab=2, cex.axis=2)
# points(pred[v==1,1][n1],pred[v==1,2][n1],col="red", cex=1.25)
# points(pred[v==2,1][n2],pred[v==2,2][n2],col="green", cex=1.25)
# points(pred[v==3,1][n3],pred[v==3,2][n3],col="blue", cex=1.25)
# points(pred[v==4,1][n4],pred[v==4,2][n4],col="black", cex=1.25)
# points(pred[v==5,1][n5],pred[v==5,2][n5],col="orange", cex=1.25)
# points(pred[v==6,1][n6],pred[v==6,2][n6],col="purple", cex=1.25)
draw.circles(pred[v==1,1][n1],pred[v==1,2][n1],border="red", cex=1.25, ir=0.06)
draw.circles(pred[v==2,1][n2],pred[v==2,2][n2],border="green", cex=1.25, ir=0.06)
draw.circles(pred[v==3,1][n3],pred[v==3,2][n3],border="blue", cex=1.25, ir=0.06)
draw.circles(pred[v==4,1][n4],pred[v==4,2][n4],border="black", cex=1.25, ir=0.06)
draw.circles(pred[v==5,1][n5],pred[v==5,2][n5],border="orange", cex=1.25, ir=0.06)
draw.circles(pred[v==6,1][n6],pred[v==6,2][n6],border="purple", cex=1.25, ir=0.06)


arrange_colorset <- c('green','purple','black','blue','red','orange')
arrange_phi_vec <- phi_vec[,c(2,6,4,3,1,5)]
flowerplot2d(predy,(1:n)[v==1][n1],arrange_phi_vec,r=0.3,
             'red',arrange_colorset, cex=1.25, ir=0.06)
flowerplot2d(predy,(1:n)[v==2][n2],arrange_phi_vec,r=0.3,
             'green',arrange_colorset, cex=1.25, ir=0.06)
flowerplot2d(predy,(1:n)[v==3][n3],arrange_phi_vec,r=0.3,
             'blue',arrange_colorset, cex=1.25, ir=0.06)
flowerplot2d(predy,(1:n)[v==4][n4],arrange_phi_vec,r=0.3,
             'black',arrange_colorset, cex=1.25, ir=0.06)
flowerplot2d(predy,(1:n)[v==5][n5],arrange_phi_vec,r=0.3,
             'orange',arrange_colorset, cex=1.25, ir=0.06)
flowerplot2d(predy,(1:n)[v==6][n6],arrange_phi_vec,r=0.3,
             'purple',arrange_colorset, cex=1.25, ir=0.06)
dev.off()
dev.off()
# # Add a legend
# legend('topleft', 
#        legend=c('ES', 'FR', 'GE', 'IT', 'UK', 'US'),
#        col=c('red', 'green', 'blue', 'black', 'orange', 'purple'),
#        pch=16, cex=1.25)

