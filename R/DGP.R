### DGP
#
# MCdata <- function(){


# func1 = function(x){
#   if(x==0) x=0 + rnorm(1, 100, 1)
#   else if(x==1) x=0 + rnorm(1, 200, 5)
#   else if(x==2) x=0 + rnorm(1, 300, 6)
#   else if(x==3) x=0 + rnorm(1, 400, 7)
#   else if(x==4) x=0 + rnorm(1, 500, 8)
#   else x = 1
# }
# 
# func2 = function(x){
#   if(x==0) x=0 + rnorm(1, 6, 1)
#   else if(x==1) x=0 + rnorm(1, 7, 1)
#   else if(x==2) x=0 + rnorm(1, 8, 1)
#   else if(x==3) x=0 + rnorm(1, 9, 1)
#   else if(x==4) x=0 + rnorm(1, 1, 1)
#   else x = 0
# }
# 
# nrow = 10000
# ncol = 5
# 
# set.seed(100)
# j = floor(runif(nrow, 0,30))
# 
# va1 = sapply(j, func1)
# va2 = sapply(j, func2)
# coefs = c(10, 15, 20, 35, 50)
# x = matrix(rnorm(nrow*ncol), nrow)
# 
# y1 = x %*% coefs + va1
# y2 = x %*% coefs + va2
# 
# y = cbind(y1,y2)
# 
# 
# start = Sys.time()
# cpp = cppmdvam(x,y,j)
# end = Sys.time()
# print(end-start)
# gc()
# }



# test3 = function(){
#   l = list();
#   mat = matrix(runif(1000**2),nrow=1000)
#   for (i in c(1:5**2)) {
#     l = append(l ,mat)
#    # print(i)
#   }
# }
# 
# 
# benchmark('1' = test1(),
#           '2' = test2(),
#           replications = 1)





