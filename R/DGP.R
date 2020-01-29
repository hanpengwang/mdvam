### DGP

# MCdata <- function(){


# func1 = function(x) {
#   if (x == 0)
#     x = 0 + rnorm(1, 0.1, 1)
#   else if (x == 1)
#     x = 0 + rnorm(1, 1, 5)
#   else if (x == 2)
#     x = 0 + rnorm(1, 10, 6)
#   else if (x == 3)
#     x = 0 + rnorm(1, 100, 7)
#   else if (x == 4)
#     x = 0 + rnorm(1, 1, 8)
#   else
#     x = 1
# }
# 
# func2 = function(x) {
#   if (x == 0)
#     x = 0 + rnorm(1, 6, 1)
#   else if (x == 1)
#     x = 0 + rnorm(1, 7, 1)
#   else if (x == 2)
#     x = 0 + rnorm(1, 8, 1)
#   else if (x == 3)
#     x = 0 + rnorm(1, 9, 1)
#   else if (x == 4)
#     x = 0 + rnorm(1, 1, 1)
#   else
#     x = 0
# }
# 
# nrow = 20000
# ncol = 10
# 
# set.seed(151)
# j = matrix(floor(runif(nrow, 0, 40)), ncol = 1)
# 
# va1 = sapply(j, func1)
# va2 = sapply(j, func2)
# coefs = runif(ncol, 1, 5)
# x = matrix(rnorm(nrow * ncol), nrow) * 1000
# 
# y1 = x %*% coefs * 100 + va1
# y2 = x %*% coefs * 100 + va2
# y3 = x %*% coefs + va1
# y4 = x %*% coefs + va2
# y = cbind(y1, y2)
# #
# #
# start = Sys.time()
# cpp = cppmdvam(x, y, j)
# end = Sys.time()
# print(end - start)
# gc()
# }






