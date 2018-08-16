\name{AST}
\alias{AST}
\title{
Age-Spatio-Temporal Model
}
\description{
The age-spatial-temporal models utilizes additional prediction analysis to take into account how the dependent variable further varies across time, space, and age. This type of age-spatial-temporal model has been used in many applications. We do this by calculating the residual (predicted - observed dependent variable) for each data point and then run local regression in three dimensions on the residual. The process assumes that residuals contain valuable information that cannot be directly observed but nonetheless vary systematically across geographic region, time, and age group. This allows us to predict how much the observed dependent variable differs from the regression model's prediction and to account for these differences. This extra variation in prediction model defines by three matrices of age, time and spaces.
Once these weights have been calculated, weighting every observation in the data set relative to the one being predicted, it a simple matter of calculating a weighted average of the residuals from regression. This "predicted residual" is then added back onto the prediction, creating an estimate that more closely takes into account aspects of the data that cannot be captured by a simple covariate model.
This package also considers the data sources' weight in calculating. But this is optional to include.
}
\usage{
AST(data.residual, spaceMatrix, par.time=0.5, par.age=1,
    weight.coverage=0.9, agecat, minyear, maxyear)
}
\arguments{
  \item{data.residual}{
Data frame of residuals from any regression model. It should be included age, year,location and residual variables. The corresponding numeric variable names in data set are residuals of any prediction model, year, location, and age. coverage variable is also optional for data sources and data points which are more reliable and representative.
}
  \item{spaceMatrix}{
A square matrix of location which the number of rows and columns are equal to the number of spaces. Colnames and rownames must be the ordered number of locations and match with the name of location variable in data.residual.
}
  \item{par.time}{
A positive parameter which can be tuned to increase or decrease how much smoothing occurs across time. Usual range is 0.5-2.0, higher value results in less strength drawn across time.
}
  \item{par.age}{
A positive parameter which can be tuned to increase or decrease how much smoothing occurs across age groups. Usual range is 0.5-2.0, higher value results in less strength drawn across age.
}
  \item{weight.coverage}{
A parameter which can be specified to give greater weights to data sources that are more reliable. If data frame does not include coverage variable, this parameter will not be considered in AST function. This parameter assigns to data points that their corresponding coverage variable is 1. [1-weight.coverage] computes for 0 ones.
}
  \item{agecat}{
The number of age groups which were predicted. It sets to 1 when we do not have age variable in data frame.
}
  \item{minyear}{
The first year of data frame which is predicted.
}
  \item{maxyear}{
The last year of data frame which is predicted.
}
}
\details{
This package adjusts residuals of any prediction model like regression by time, age, and spaces. The weighted residuals then added back to predictions.

Space matrix specifies the amount of weights in location. It should be a square matrix or data frame with equal number of rows and columns with location variable in data.residual. This data set is define by users with desirable amount of weights with adjacent locations or even how smooth they can be in far from areas. If users do not have this information, calcSpaceMat is a function in this package which builds a matrix of adjacent locations from a shapefile and give them equal weights.

Space matrix can be handled by downloading a shapefile.  The shapefile format is a popular geospatial vector data format for geographic information system (GIS) software.The term "shapefile" is quite common, but is misleading since the format consists of a collection of files with a common filename prefix, stored in the same directory. The three mandatory files have filename extensions .shp, .shx, and .dbf. In example part, we used a package to read this file format and determine a matrix of adjacent locations.
}
\value{
a list including three parts as given by:

1. adj.rate    a data frame that includes 5 variables: ID, time, year, locations, and adjusted
residuals.

2. Age_weight    matrix

3. time_weight   matrix
}
\references{
Foreman et al.: Modeling causes of death: an integrated approach using CODEm. Population Health Metrics 2012 10:1.
}
\author{
Parinaz Mehdipour, Ali Ghanbari
}
\note{
The age, time, and space variables should be ordered.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library("AST")

data(data.residual.AST,package="AST")
data.residual <- data.residual.AST[data.residual.AST$sex == 0 ,]

data(data.adjacent.mat,package="AST")
spaceMat <- calcSpaceMat(data.adjacent.mat,par.space=0.9)

AST(data.residual=data.residual,spaceMatrix=spaceMat,par.time = 0.5,
     par.age=1, weight.coverage=0.9, agecat=3, minyear=2000, maxyear=2005)

}

