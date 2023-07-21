library(Cairo)
##############################################################-
## flatdistn                                           ####
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

## EPS
cairo_ps("flatdistn.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
curve(ff, from = 0, to = 6,  
      xlab = expression(italic('\u03B8')),
      ylab = expression(paste( italic('\u03c0'), (italic('\u03B8')~"|"~italic(x)))),
      ylim = c(0, 2), las=1,
      lwd=2.5, cex=2,cex.lab=2, cex.axis=2)
# default margin size: mar = c(5.1, 4.1, 4.1, 2.1) 
#side=2: left, line=1: margin 1, at=1: location y=1, las=2: perdendicular to axis

lines(c(1,1), c(0,1), lty=3, lwd=2.5)
lines(c(2,2), c(0,1), lty=3, lwd=2.5)
lines(c(4,4), c(0,1), lty=3, lwd=2.5)
lines(c(5,5), c(0,1), lty=3, lwd=2.5)

x1 <- seq(1,5,0.01)
y1 <-  ff(x1)
x1 <- c(1,x1,5)
y1 <-  c(0,y1,0)
polygon(x1,y1,border=NA,  col="red", angle=135, density=10, lty = 2, lwd = 2.5)

x2 <- seq(2,4,0.01)
y2 <-  ff(x2)
x2 <- c(2,x2,4)
y2 <-  c(0,y2,0)
polygon(x2,y2,border=NA,  col="blue", angle=45, density=10, lty = 1, lwd = 2.5)

abline(h=1, lty=3, col=kacol, lwd=3)
mtext(expression(italic('\u03BA')[italic('\u03B1')]), side=4, line=0.6, at=1,
      cex=2, col=kacol, las=2)
dev.off()
##############################################################-
## discdistn                                          ####
##############################################################-
n <- 5
p <- 0.5
x <- 0:n

set.seed(123)
data <- dbinom(x=x, size=n, prob=p)
names(data) <- 0:n

cairo_ps("discdistn.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
plot(x, data, type = 'h', 
     xlab = expression(italic('\u03B8')),
     ylab = expression(paste( italic('\u03c0'), (italic('\u03B8')~"|"~italic(x)))),
     xlim = c(-0.1, 5.1),
     ylim = c(0.01, 0.4),
     las = 1,
     lwd = 2.5, cex=2, cex.lab=2, cex.axis=2
     )
# title(ylab = expression(pi(theta~"|"~x)),
#       cex.lab=1.5, cex.axis=1.25, line = 2.2)
abline(h=0.03125, lty=3, col=kacol, lwd=3)
mtext(expression(italic('\u03BA')[italic('\u03B1')]), side=4, line=0.6, at=0.03125,  
      cex=2, col=kacol, las = 2)

library(dplyr)
data2 <- data %>% round(.,2) %>% as.character()
data2y <- data %>% round(.,2) + 0.02
text(x = x, y = data2y, labels = data2, cex=1.5)
dev.off()

##############################################################-
## normal                                           ####
##############################################################-
## standard normal posterior
# curve(dnorm, from=-3.5, to=3.5,
#       xlab = expression(theta),
#       ylab = expression(pi(theta~"|"~x))
#       )
x <- seq(-3.5, 3.5, length=2000)
set.seed(123)
data <- dnorm(x=x)

cairo_ps("normal.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
plot(x, data, type = 'l',
     xlab = expression(italic('\u03B8')),
     ylab = expression(paste( italic('\u03c0'), (italic('\u03B8')~"|"~italic(x)))),
     las = 1,
     lwd=2.5, cex=2, cex.lab=2, cex.axis=2
)
## get ka
alpha <- 0.05
c <- qnorm(1-alpha/2)
ka <- dnorm(c)
ka # 0.05844
abline(h=ka, lty=3, col=kacol, lwd=3)
mtext(expression(italic('\u03BA')[italic('\u03B1')]), side=4, line=0.6, at=ka,
      cex=2, col=kacol, las = 1)

lines(c(-1.96,-1.96), c(0,ka), lty=3, lwd=2.5)
lines(c(1.96,1.96), c(0,ka), lty=3, lwd=2.5)
dev.off()


##############################################################-
## normalGHPD                                           ####
##############################################################-
## phi
x <- seq(-c, c, length=2000)
y <- rep(1, length=2000)

cairo_ps("normalGHPD.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
plot(x,y, type = 'l',
     xlab = expression(italic('\u03B8')),
     ylab = expression(paste( italic('\u03D5')^"*", (italic('\u03B8')~"|"~italic(x)))),
     las = 1, 
     xlim = c(-3.5, 3.5),
     ylim = c(0, 1.1),
     lwd=2.5, cex=2, cex.lab=2, cex.axis=2
)
lines(c(-3.5,-c), c(0,0),lwd=2.5)
lines(c(c,3.5), c(0,0),lwd=2.5)
lines(c(-c,-c), c(0,1),lwd=2.5)
lines(c(c,c), c(0,1),lwd=2.5)
dev.off()

##############################################################-
## flatGHPD                                           ####
##############################################################-
## phi 
cairo_ps("flatGHPD.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
plot(NULL, xlim=c(0,6), ylim=c(0,1.1),
     xlab = expression(italic('\u03B8')),
     ylab = expression(paste( italic('\u03D5')[italic(f)]^"*", (italic('\u03B8')~"|"~italic(x)))),
     las = 1,
     lwd=2.5, cex=2, cex.lab=2, cex.axis=2
)

alpha <- 0.35
nc <- 3/20 # normalizing constant
gamma <- (1-alpha-nc*4)/(nc*2)
gamma # 0.167

lines(c(0,1), c(0,0),lwd=2.5)
lines(c(1,2), c(gamma,gamma),lwd=2.5)
lines(c(2,4), c(1,1),lwd=2.5)
lines(c(4,5), c(gamma,gamma),lwd=2.5)
lines(c(5,6), c(0,0),lwd=2.5)

lines(c(1,1), c(0,gamma),lwd=2.5)
lines(c(2,2), c(gamma,1),lwd=2.5)
lines(c(4,4), c(gamma,1),lwd=2.5)
lines(c(5,5), c(0,gamma),lwd=2.5)

gammacol <- 'blue'
abline(h=gamma, lty=3, col=gammacol, lwd=3)
mtext(expression(italic('\u03B3')), side=4, line=0.8, at=gamma, 
      cex=2, col=gammacol, las = 1)
dev.off()

# ##############################################################-
# ## discGHPD                                           ####
# ##############################################################-
n <- 5
p <- 0.5
x <- 0:n
alpha <- 0.05

set.seed(123)
data <- dbinom(x=x, size=n, prob=p)
names(data) <- 0:n

ka <- 0.03125
gamma <- (1-alpha-sum(data[data > ka]))/(sum(data[data == ka]))
gamma # 0.2

data[c(1,6)] <- gamma
data[c(2,3,4,5)] <- 1

cairo_ps("discGHPD.eps", family = "Times New Roman")
par(family = "Times New Roman", mar = c(5, 6, 4, 3) + 0.1)
plot(x, data, type = 'h',
     xlab = expression(italic('\u03B8')),
     ylab = expression(paste( italic('\u03D5')[italic(f)]^"*", (italic('\u03B8')~"|"~italic(x)))),
     xlim = c(-0.1, 5.1),
     ylim = c(0, 1.1),
     las = 1, 
     lwd=2.5, cex=2, cex.lab=2, cex.axis=2
)
abline(h=gamma, lty=3, col=gammacol, lwd=3)
mtext(expression(italic('\u03B3')), side=4, line=0.8, at=gamma,
      cex=2, col=gammacol, las = 1)
dev.off()
