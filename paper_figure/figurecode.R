##############################################################-
## Figure 1 (a)                                           ####
##############################################################-
f1 <- function(x){x^2}
f2 <- function(x){1}
f3 <- function(x){-(x-3)^2+2}
f4 <- function(x){(x-6)^2}

ff <- function(x){
  ifelse(x<=1, f1(x), 
         ifelse(1<x & x<=2, f2(x),
                ifelse(2<x & x<=4, f3(x),
                       ifelse(3<x & x<=5, f2(x),
                              ifelse(5<x & x<=6, f4(x), 0)))))
}

kacol <- 'darkorchid1'

par(family = "Times New Roman", mar = c(5, 5, 4, 3) + 0.1)
curve(ff, from = 0, to = 6,  
      xlab = expression(theta), ylab = expression(pi(theta~"|"~x)), 
      ylim = c(0, 2), las=1,
      lwd=2, cex=1.5,cex.lab=1.5, cex.axis=1.5)

#axis(2, at = seq(0, 2), cex.axis = 1.25, las=1)

# title(ylab = expression(pi(theta~"|"~x)),
#       cex.lab=1.5, cex.axis=1.25, line = 1
#       )

# curve(ff, from = 0, to = 6, yaxt='n', 
#       xlab = expression(theta),
#       ylab = expression(pi(theta~"|"~x)),
#       lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.25)

# default margin size: mar = c(5.1, 4.1, 4.1, 2.1) 

#side=2: left, line=1: margin 1, at=1: location y=1, las=2: perdendicular to axis
lines(c(1,1), c(0,1), lty=3, lwd=2)
lines(c(2,2), c(0,1), lty=3, lwd=2)
lines(c(4,4), c(0,1), lty=3, lwd=2)
lines(c(5,5), c(0,1), lty=3, lwd=2)

x1 <- seq(1,5,0.01)
y1 <-  ff(x1)
x1 <- c(1,x1,5)
y1 <-  c(0,y1,0)
polygon(x1,y1,border=NA,  col="red", angle=45, density=10, lwd = 2)

x2 <- seq(2,4,0.01)
y2 <-  ff(x2)
x2 <- c(2,x2,4)
y2 <-  c(0,y2,0)
polygon(x2,y2,border=NA,  col="blue", angle=135, density=10, lty = 2, lwd = 2)

abline(h=1, lty=3, col=kacol, lwd=3)
mtext(expression(kappa[alpha]), side=4, line=0.6, at=1,
      cex=2, col=kacol, las=2)
##############################################################-
## Figure 1 (b)                                           ####
##############################################################-
n <- 5
p <- 0.5
x <- 0:n

set.seed(123)
data <- dbinom(x=x, size=n, prob=p)
names(data) <- 0:n

par(family = "Times New Roman", mar = c(5, 5, 4, 3) + 0.1)
plot(x, data, type = 'h', 
     xlab = expression(theta),
     ylab = expression(pi(theta~"|"~x)),
     xlim = c(-0.1, 5.1),
     ylim = c(0.01, 0.4),
     las = 1,
     lwd = 2, cex=1.5, cex.lab=1.5, cex.axis=1.5
     )
# title(ylab = expression(pi(theta~"|"~x)),
#       cex.lab=1.5, cex.axis=1.25, line = 2.2)
abline(h=0.03125, lty=3, col=kacol, lwd=3)
mtext(expression(kappa[alpha]), side=4, line=0.6, at=0.03125,  
      cex=2, col=kacol, las = 2)

library(dplyr)
data2 <- data %>% round(.,2) %>% as.character()
data2y <- data %>% round(.,2) + 0.02
text(x = x, y = data2y, labels = data2, cex=1.5)


# n <- 5
# p <- 0.5
# x <- 0:n
# 
# set.seed(123)
# data <- dbinom(x=x, size=n, prob=p)
# names(data) <- 0:n
# 
# barplot(data, ylim=c(0, 0.4), col = 'white', 
#         xlab = expression(theta),
#         ylab = expression(pi(theta~"|"~x))
#         )
# abline(h=0.03125, lty=2)
# mtext(expression(kappa[alpha]), side=2, line=1, at=0.03125, las=2)

##############################################################-
## Figure 2 (a)                                           ####
##############################################################-
## standard normal posterior
# curve(dnorm, from=-3.5, to=3.5,
#       xlab = expression(theta),
#       ylab = expression(pi(theta~"|"~x))
#       )
x <- seq(-3.5, 3.5, length=2000)
set.seed(123)
data <- dnorm(x=x)

par(family = "Times New Roman", mar = c(5, 5, 4, 3) + 0.1)
plot(x, data, type = 'l',
     xlab = expression(theta),
     ylab = expression(pi(theta~"|"~x)),
     las = 1,
     lwd=2, cex=1.5, cex.lab=1.5, cex.axis=1.5
)
# title(ylab = expression(pi(theta~"|"~x)),
#       cex.lab=1.5, cex.axis=1.25, line = 2.2)
## get ka
alpha <- 0.05
c <- qnorm(1-alpha/2)
ka <- dnorm(c)
ka
abline(h=ka, lty=3, col=kacol, lwd=3)
mtext(expression(kappa[alpha]), side=4, line=0.6, at=ka,
      cex=2, col=kacol, las = 1)

lines(c(-1.96,-1.96), c(0,ka), lty=3, lwd=2)
lines(c(1.96,1.96), c(0,ka), lty=3, lwd=2)



##############################################################-
## Figure 2 (b)                                           ####
##############################################################-
## phi
x <- seq(-c, c, length=2000)
y <- rep(1, length=2000)

par(family = "Times New Roman", mar = c(5, 5, 4, 3) + 0.1)
plot(x,y, type = 'l',
     xlab = expression(theta),
     ylab = expression(paste( phi^"*", (theta~"|"~x))),
     las = 1, 
     xlim = c(-3.5, 3.5),
     ylim = c(0, 1.1),
     lwd=2, cex=1.5, cex.lab=1.5, cex.axis=1.5
)
# title(ylab = expression(paste( phi^"*", (theta~"|"~x))),
#       cex.lab=1.5, cex.axis=1.25, line = 2.2)
lines(c(-3.5,-c), c(0,0),lwd=2)
lines(c(c,3.5), c(0,0),lwd=2)
lines(c(-c,-c), c(0,1),lwd=2)
lines(c(c,c), c(0,1),lwd=2)


##############################################################-
## Figure 3 (a)                                           ####
##############################################################-

# f1 <- function(x){x^2}
# f2 <- function(x){1}
# f3 <- function(x){-(x-3)^2+2}
# f4 <- function(x){(x-6)^2}
# 
# ff <- function(x){
#   ifelse(x<=1, f1(x), 
#          ifelse(1<x & x<=2, f2(x),
#                 ifelse(2<x & x<=4, f3(x),
#                        ifelse(3<x & x<=5, f2(x),
#                               ifelse(5<x & x<=6, f4(x), 0)))))
# }
# 
# kacol <- 'darkorchid1'
# 
# curve(ff, from = 0, to = 6, yaxt='n', 
#       xlab = expression(theta),
#       ylab = expression(pi(theta~"|"~x)),
#       lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.25
# )
# abline(h=1, lty=2, col=kacol, lwd=2)
# mtext(expression(kappa[alpha]), side=2, line=1, at=1, las=2, 
#       cex=2, col=kacol)
# #side=2: left, line=1: margin 1, at=1: location y=1, las=2: perdendicular to axis
# # lines(c(1,1), c(0,1), lty=2)
# # lines(c(2,2), c(0,1), lty=2)
# # lines(c(4,4), c(0,1), lty=2)
# # lines(c(5,5), c(0,1), lty=2)


## phi 
par(family = "Times New Roman", mar = c(5, 5, 4, 3) + 0.1)
plot(NULL, xlim=c(0,6), ylim=c(0,1.1),
     xlab = expression(theta),
     ylab = expression(paste( phi[f]^"*", (theta~"|"~x))),
     las = 1,
     lwd=2, cex=1.5, cex.lab=1.5, cex.axis=1.5
)
# title(ylab = expression(paste( phi[f]^"*", (theta~"|"~x))),
#       cex.lab=1.5, cex.axis=1.25, line = 2.2)

alpha <- 0.35
nc <- 3/20 # normalizing constant
gamma <- (1-alpha-nc*4)/(nc*2)
gamma

lines(c(0,1), c(0,0),lwd=2)
lines(c(1,2), c(gamma,gamma),lwd=2)
lines(c(2,4), c(1,1),lwd=2)
lines(c(4,5), c(gamma,gamma),lwd=2)
lines(c(5,6), c(0,0),lwd=2)

lines(c(1,1), c(0,gamma),lwd=2)
lines(c(2,2), c(gamma,1),lwd=2)
lines(c(4,4), c(gamma,1),lwd=2)
lines(c(5,5), c(0,gamma),lwd=2)

gammacol <- 'blue'
abline(h=gamma, lty=3, col=gammacol, lwd=3)
mtext(expression(gamma), side=4, line=0.8, at=gamma, 
      cex=2, col=gammacol, las = 1)
# expression(gamma(theta~"|"~x))

##############################################################-
## Figure 3 (b)                                           ####
##############################################################-


##############################################################-
## Figure 4 (a)                                           ####
##############################################################-
n <- 5
p <- 0.5
x <- 0:n
alpha <- 0.05

set.seed(123)
data <- dbinom(x=x, size=n, prob=p)
names(data) <- 0:n

# plot(x, data, type = 'h', lwd = 2,
#      xlab = expression(theta),
#      ylab = expression(pi(theta~"|"~x)),
#      xlim = c(-0.1, 5.1),
#      ylim = c(0.01, 0.4)
# )
# abline(h=0.03125, lty=2, col=kacol, lwd=2)
# mtext(expression(kappa[alpha]), side=2, line=1, at=0.03125, las=2,
#       cex=1.25, col=kacol)
# 
# library(dplyr)
# data2 <- data %>% round(.,2) %>% as.character()
# data2y <- data %>% round(.,2) + 0.02
# text(x = x, y = data2y, labels = data2)

# ##############################################################-
# ## Figure 4 (b)                                           ####
# ##############################################################-
ka <- 0.03125
gamma <- (1-alpha-sum(data[data > ka]))/(sum(data[data == ka]))

data[c(1,6)] <- gamma
data[c(2,3,4,5)] <- 1

plot(x, data, type = 'h',
     xlab = expression(theta),
     ylab = expression(paste( phi[f]^"*", (theta~"|"~x))),
     xlim = c(-0.1, 5.1),
     ylim = c(0, 1.1),
     las = 1, 
     lwd=2, cex=1.5, cex.lab=1.5, cex.axis=1.5
)
# title(ylab = expression(paste( phi[f]^"*", (theta~"|"~x))),
#       cex.lab=1.5, cex.axis=1.25, line = 2.2)
abline(h=gamma, lty=3, col=gammacol, lwd=3)
mtext(expression(gamma), side=4, line=0.8, at=gamma,
      cex=2, col=gammacol, las = 1)
