
calcSpaceMat<- function (adjacent.mat,par.space=0.9){

  if( !is.matrix(adjacent.mat)) stop("space matrix must be a square matrix.")
  if(!is.numeric(par.space)) stop("par.space must be a number")
  if( par.space<=0 | par.space>1 )  stop("par.space must be a number between 0 and 1.")
  if( any(c(  !is.numeric(adjacent.mat) , length(table(adjacent.mat)) !=2, max(adjacent.mat)!=1, min(adjacent.mat)!=0  ))  ) stop("adjacent.mat must contain 0 or 1 ")
  if( is.null(rownames(adjacent.mat)  )  ) stop("row names of adjacent.mat is necessary!")
  if(   any(  is.na(
    as.numeric(rownames(adjacent.mat))
  )
  )
  )stop("row names of adjacent.mat must be matched with location(number)")


  diag(adjacent.mat) <- par.space
  adjacent.mat[which(adjacent.mat==1)] <- par.space*(1-par.space)
  return (adjacent.mat)
}

