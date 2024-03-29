% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s1-utils.R
\name{s1_add_ratio}
\alias{s1_add_ratio}
\title{s1_add_ratio}
\usage{
s1_add_ratio(img, method = c("DB", "DN"), bands = c("VV", "VH"))
}
\arguments{
\item{img}{an ee image}

\item{method}{"DB" when image values have been converted to dB using, "DN"
otherwise.}

\item{bands}{default is c('VV', 'VH'). must be length 2. The bands to used to
calculate ratio. First value is either subtracted from or used to divide the
second dependng on the method arg.}
}
\value{

}
\description{
Add band ratio values to image.
}
