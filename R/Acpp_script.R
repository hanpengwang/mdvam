#' Estimate value added
#'
#'
#' @param x dependent variables
#' @param y indenpent variables, can be one or multiple
#' @param j this is a vector where each value corresponding each observation's id
#' @param comp_effcet default null. a vector contains variables names or index which you want to create composition effects
#' @param ... for future implementation
#'
#' @return function returns multivariate value added outcome for each specific category
#'
#' @examples
#' cppmdvam(x, y, j)
#' cppmdvam(x, y, j, c('X1','X2',..,'Xm'))
#' cppmdvam(x, y, j, c('1,2,..,n))
#'
#' @export


cppmdvam <-
  function(x, y, j, comp_effect=NULL, ...)
  {
    x <- as.matrix(x); y <- as.matrix(y); j <- as.matrix(j);
    N <- nrow(x)
    if (nrow(y)!=N || nrow(j)!=N) stop("x,y,j should have the same number of obeservations")
    if (ncol(j)!=1) stop("each obeservation should have only one category")

    if (is.null(colnames(x))) {Xnames <- paste('x', c(1:ncol(x)), sep = ""); colnames(x) = Xnames}
    else Xnames <- colnames(x)
    if (is.null(colnames(y))) Ynames <- paste('y', c(1:ncol(y)), sep = "")
    else Ynames <- colnames(y)
    if (is.null(colnames(j))) Jnames <- c("J")
    else Jnames <- colnames(j)

    M <- ncol(y)
    uniqueJ <- rownames(table(j))
    Nj <- as.vector(table(j))
    K <- ncol(x)
    J <- length(Nj)
    DF <- M*N - M * K - M * J

    data_jyx <- cbind(j, y, x); colnames(data_jyx) <- c(Jnames, Ynames, Xnames)

    if ( !is.null(comp_effect)) {

      comp_df <- cbind(j , x[,c(comp_effect)])
      comp_col <- comp_create(comp_df, Nj)
      data_jyx <- cbind(data_jyx, comp_col)
      Xnames <- c(Xnames, colnames(comp_col))


    }

    Xlist <- list()
    Ylist <- list()
    for (j in c(1:J))
    {
      jname <- uniqueJ[j]
      dfj <- as.matrix(data_jyx[data_jyx[,Jnames] == jname, ])
      Xlist[[j]] <- dfj[,Xnames]
      Ylist[[j]] <- dfj[,Ynames]

    }
    All <- list(Ylist, Xlist, M, N, Nj,
                K, J, DF)

    va = ValueAdded(All)
    
    return(va)


  }

