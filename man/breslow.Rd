\name{breslow}
\alias{breslow}
\title{Breslow estimator for cumulative baseline hazard rate}
\description{This function calculates the Breslow estimator for the
  cumulative hazard rate, given a fitted multi state model.}
\usage{
breslow(m, mstop, entry, exit, trans, event)
}
\arguments{
\item{m}{ a fitted multi state model using \code{gamboost(...,
    family=multistateInR(...))} or \code{gamboost(..., family=multistate(...))}.}  
\item{mstop}{ an integer defing at which boosting iteration to evaluate.}
\item{entry}{ entry times.}
\item{exit}{ exit times.}
\item{trans}{ transition index.}
\item{event}{ observed event indicator.}
}
\details{This function calculates the Breslow estimator for the
  cumulative hazard rate, given a fitted multi state model and a
  given transition q.}
\value{ 
A list of length Q with each element including including elements
\item{times}{ a vector of observed event times,}
\item{bhr}{ a vector of calculated cumulative hazard rate values.}
}
\author{Holger Reulen}
\examples{
\dontrun{breslow(m, mstop, entry, exit, trans, event)}
}
