### DGP
#
# MCdata <- function(){
#
#
# func1 = function(x){
#   if(x==0) x=0 + rnorm(1, 1, 1)
#   else if(x==1) x=0 + rnorm(1, 2, 1)
#   else if(x==2) x=0 + rnorm(1, 3, 1)
#   else if(x==3) x=0 + rnorm(1, 4, 1)
#   else if(x==4) x=0 + rnorm(1, 5, 1)
# }
#
# func2 = function(x){
#   if(x==0) x=0 + rnorm(1, 6, 1)
#   else if(x==1) x=0 + rnorm(1, 7, 1)
#   else if(x==2) x=0 + rnorm(1, 8, 1)
#   else if(x==3) x=0 + rnorm(1, 9, 1)
#   else if(x==4) x=0 + rnorm(1, 1, 1)
# }
#
# j = floor(runif(10**3, 0,5))
#
# va1 = sapply(j, func1)
# va2 = sapply(j, func2)
# coefs = c(10, 15, 20, 35, 50)
# x = matrix(rnorm(10**3*5), nrow = 10**3)
#
# y1 = x %*% coefs
# y2 = x %*% coefs
#
# y = cbind(y1,y2)
#
#
# start = Sys.time()
#cpp = cppmvam(x,y,j)
# end = Sys.time()
# print(end-start)

# }
