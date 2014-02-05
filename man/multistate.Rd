\name{multistate}
\alias{multistate}
\title{Family for Multistate Models}
\description{This function implements a family for
  fitting multistate models using \code{mboost}.}
\usage{
multistate(trans)
}
\arguments{
\item{trans}{ a vector denoting respective transitions.}
} 
\details{This function implements a family for multistate models and
  will be used inside the \code{gamboost} or \code{glmboost} function.}
\value{ Functions to be used inside gamboost.
}
\author{Holger Reulen}
