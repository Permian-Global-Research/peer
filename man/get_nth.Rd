% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgee-helpers.R
\name{get_nth}
\alias{get_nth}
\title{Get nth item from an ee collection}
\usage{
get_nth(collect, n)
}
\arguments{
\item{n}{}
}
\value{
The entire unclipped image at the nth position of the collection.
}
\description{
Get nth item from an ee collection
}
\examples{
s1.test <- s1_collect(kuamut, '2019-01-01', '2019-12-31' )
s1_nth(s1.test, 3)

}
