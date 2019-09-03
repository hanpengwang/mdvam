comp_create <- function(comp_df, n) ### create composition effect regards to independent variables
{
  col_names <- colnames(comp_df[,-1])
  COMP <- aggregate(comp_df, list(idj = comp_df[,c(1)]),  FUN=function(x){return(mean(x,na.rm=T))} )
  COM_cols <- sapply(COMP, function(x){rep(x, n)}); COM_cols <- COM_cols[,c(-1,-2)]
  colnames(COM_cols) <- paste(col_names, 'Comp',sep = '_')
  return(COM_cols)

}
