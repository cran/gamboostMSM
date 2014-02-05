\name{dof}
\alias{dof}
\title{Degrees of Freedom}
\description{This function calculates the degrees of freedom as part of
  the calculation of the Akaike Information Criterion (AIC).}
\usage{
dof(m)
}
\arguments{
  \item{m}{ a fitted multi state model using 
    \code{gamboost(..., family=multistate(...))}.}  
}
\details{This function calculates the degress of freedom as part of
  the calculation of the Akaike Information Criterion.}
\value{ 
A vector of length equal to the number of boosting iterations in the
plugged in model object. 
}
\author{Holger Reulen}
\examples{
\dontrun{dof(m)}
}
