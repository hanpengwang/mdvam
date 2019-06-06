# get_params <- function(x, y, col_category){
#
#                 idj <<- as.matrix(unlist(unique(schools))) ### schools' names
#
#                 M <<- length(y[1,])
#
#                 x <<- Z.var(x, M, col_category)
#                 y <<- Y.var(y, col_category)
#
#                 n <<- as.matrix(table(col_category))
#
#                 w <<-  Within(M, n)
#
#
#
#
#
#
#
#
# }
#
#
#
library(devtools);library(Matrix)
data <- read.csv('C:\\Users\\wangh\\OneDrive\\Desktop\\data_uni.csv')
schools <- data['idj']
X <- data[, c(1:6)]
Y <- data[, c(10:14)]

#test_result <- mvam(X,Y, schools)

#
#

#beta.w <- solve(crossprod(test[[1]])) %*% crossprod(test[[1]], test[[2]])

# idj <- as.matrix(unlist(unique(schools)))
# n <- as.matrix(table(schools))
#
#
#
# dim(svd(result))   # this can inverse singular matrix















































# Within2x.j <- function(M,j,data.x,namesvec)
# {
#   W.j <- Within.j(M,j,data.x)
#   X.j <- Xj.var(M,j,data.x,namesvec)
#   WXj <- W.j%*%X.j
#   return(WXj)
# }
#
#
# Within2x <- function(M,J,data.x,namesvec)
# {
#   wx <- 0
#   for (j in 1:J){
#     wx <- rBind(wx, Within2x.j(M,j,data.x,namesvec))  }
#   wx <- wx[-1,]
#   return(wx)
# }
#
#
# #--------------------------------------------------------------------------------------------------------------------
#
#
# Within2y.j <- function(M,j,data.y, namesvecspro)
# {
#   W.j <- Within.j(M,j,data.y)
#   Y.j <- Yj.var(M,j,data.y, namesvecspro)
#   WY.j <- W.j %*% Y.j
#   return(WY.j)
# }
#
# Within2y <- function(M,J,data.y, namesvecspro)
# {
#   wy <- 0
#   for (j in 1:J){  wy <- rBind(wy, Within2y.j(M,j,data.y, namesvecspro))  }
#   wy <- wy[-1,]
#   return(Matrix(wy))
# }










