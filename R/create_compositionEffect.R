comp_create <- function(comp_df, n) ### create composition effect regards to independent variables
{
  
  COMP <- aggregate(comp_df, list(idj = comp_df[,c(1)]),  FUN=function(x){return(mean(x,na.rm=T))} )
  COM_cols <- sapply(COMP, function(x){rep(x, n)})
  colnames(COM_cols) <- paste(colnames(COM_cols), 'Comp',sep = '_') 
  
  return(COM_cols)
  
}