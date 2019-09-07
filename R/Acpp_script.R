#' Estimate value added
#'
#'
#' @param x dependent variables
#' @param y indenpent variables, can be one or multiple
#' @param j this is a vector where each value corresponding each observation's id
#' @param comp_effcet default null. a vector contains variables names or index which you want to create composition effects
#' @param ... for future implementation
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#' @examples
#' ValueAdded(cppmdvam(x, y, j))
#' ValueAdded(cppmdvam(x, y, j, c('X1','X2',..,'Xm')))
#' ValueAdded(cppmdvam(x, y, j, c('1,2,..,n)))
#'
#' @export

### variable user guide
# category : a vector corresponding each id for X and Y
# groups or n : a vector with number of observations for each category
# unique.category : set(category)

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


    return(ValueAdded(All))
    ##---1. split dataX, Y, M into list elements;
    ##---2. get params N, M etc.

    ##--- check data variations

  }

# library(usethis);library(devtools)
# data = read.csv('data_uni.csv')
# j = data['idj']
# x = data[,c('quimica_punt', 'fisica_punt', 'ingles_punt')]
# y = data[,c('lenguaje_punt', 'matematicas_punt', 'ciencias_sociales_punt')]



# listvar = cppmvam(x,y,j)
# test = ValueAdded(listvar)

# p1 = listvar[[1]]
# p2 = listvar[[2]]
# p3 = listvar[[3]]
# p4 = listvar[[4]]
# p5 = listvar[[5]]
# p6 = listvar[[6]]
# p7 = listvar[[7]]
# p8 = listvar[[8]]
# test = ValueAdded(p1,p2,p3,p4,p5,p6,p7,p8)
