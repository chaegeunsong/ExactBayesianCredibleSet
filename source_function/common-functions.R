##############################################################
#                   commonly used functions
############################################################## 
##############################################################
#                   symmtrize a matrix
############################################################## 
symmetry = function(a){
return((a + t(a))/2)}
##############################################################
#            function 4:  distance between subspaces
#                         used in Li, Zha, Chiaromonte (2005)           
#  dism if v1 v2 are matrices; disv if they are    
#                    vectors                      
##############################################################
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
##############################################################
#             function 5:                          
#          trace multiple correlation              
##############################################################
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
##############################################################
#     function 6: power of a matrix               
##############################################################
matpower = function(a,alpha){
a = (a + t(a))/2
tmp = eigen(a)
return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
t(tmp$vectors))}
##############################################################
#  function 8: Moore-Penrose type power            
#  Taking power ignoring 0 eigenvalues;            
#    ignoring criterion=ignore                     
##############################################################
mppower = function(matrix,power,ignore){
eig = eigen(matrix)
eval = eig$values
evec = eig$vectors
m = length(eval[abs(eval)>ignore])
tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
t(evec[,1:m])
return(tmp)
}
##############################################################
#      function 9: center X (n by p matrix)        
##############################################################
center = function(x){
return(t(t(x)-apply(x,2,mean)))}
#######################################################################
#                  function 2: standardize a matrix
#                              treating each row as a random vector
#                              in an iid sample
#######################################################################
standmat = function(x){
mu = apply(x,2,mean)
sig = var(x)
signrt = matpower(sig,-1/2)
return(t(t(x) - mu)%*%signrt)
}
#######################################################################
#                  function 4: standardize a vector
#######################################################################
standvec = function(x){
return((x - mean(x))/sd(x))}


