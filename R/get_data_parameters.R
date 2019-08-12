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
#library(devtools);library(Matrix); library(MASS)
#data2 <- read.csv('data_uni.csv')
 #data <- read.csv('C:/Users/wangh/OneDrive/Desktop/simce.csv')
#  schools2 <- data2['idj']
#  X2 = data2[,c(2:7, 16:21)]
# #  Y2 = data2[names(data2) %in% c("mod_comunica_escrita_punt", "mod_ingles_punt", "mod_lectura_critica", "mod_razona_cuantitativo_punt", "mod_comp_ciudadanas_punt"),]
# # # 
# # test_result <- mvam(X2,Y2, schools2)
# # 
# # #
# #
# 
# # # data2 <- read.csv('data_uni.csv')
# data <- read.csv('C:/Users/wangh/OneDrive/Desktop/simce_2013-15_4b2b.csv')
#  schools <- data['rbd4b']
# #X <- new_data[, c(4:7,12:13)]
# # Y <- new_data[, c(9:11)]
# #schools <- data['rbd4b']
# X <- data[, c(2:7,12:32)]
# Y <- data[, c(9:11)]
# #beta.w <- solve(crossprod(test[[1]])) %*% crossprod(test[[1]], test[[2]])
# 
# # idj <- as.matrix(unlist(unique(schools)))
# # n <- as.matrix(table(schools))
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










