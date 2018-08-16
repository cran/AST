\name{calcSpaceMat}
\alias{calcSpaceMat}
\title{
Space Matrix Calculation for AST
}
\description{
a function that calculates the space matrix for AST model by determining adjacent locations and give them equal weights 1 and 0 for non-neighbors.
}
\usage{
calcSpaceMat(adjacent.mat, par.space = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjacent.mat}{
a matrix for adjacent locations with weight 1 and 0 for non-neighbors. Rownames of this matrix must be match with location variable.
}
  \item{par.space}{
A parameter which specify the amount of correlation between adjacent locations.
}
}
\value{
a square matrix
}

\author{
Parinaz Mehdipour, Ali Ghanbari
}
\note{
rowname of space matrix should be numeric and the same as location variable in data.residual.

If a shapefile is used for making adjacent matrix, package "spdep" can be used to identify neighbors.
}


\examples{

library("AST")
data(data.adjacent.mat,package="AST")
calcSpaceMat(data.adjacent.mat,par.space=0.9)

}
