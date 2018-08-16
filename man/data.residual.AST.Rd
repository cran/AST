\name{data.residual.AST}
\alias{data.residual.AST}
\docType{data}
\title{
residual data set
}
\description{
a data frame of residuals excluded from a prediction model
}
\usage{data("data.residual.AST")}
\format{
  A data frame with 9512 observations on the following 6 variables.
  \describe{
    \item{\code{year}}{a numeric vector that specify the time dimension in our predictions}
    \item{\code{sex}}{a numeric vector}
    \item{\code{age}}{a numeric vector that specify the age dimension in our predictions. It is categorical and sorted. }
    \item{\code{residual}}{a numeric vector of residuals from Generalized Linear Mixed Model (GLMM) that is the selected prediction model for this data set. }
    \item{\code{coverage}}{a numeric vector that specify the amount of reliability for some data points, National data sources get coverage equal to 1 and other sub-national data sources and systematic reviews get 0. The coverage weight is applied to data sources with the amount of 1 in coverage variable. }
    \item{\code{location}}{a numeric vector that specify the space dimension in our predictions and in this data set is Iran's provinces. }
  }
}
\details{
This is a sub-sample of Body Mass Index (BMI) real data from National and Sub-national Burden of Diseases, Injuries, and Metabolic Risk Factors in Iran (NASBOD). There are nine different data sources for BMI that some of them are more reliable and representative and we give higher weights to them.
}
\source{
https://www.ncdrc.info/
}

\examples{
data(data.residual.AST)
}
