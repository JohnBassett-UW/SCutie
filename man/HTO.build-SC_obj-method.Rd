% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCobj.r
\name{HTO.build,SC_obj-method}
\alias{HTO.build,SC_obj-method}
\title{HTO.build
Build an HTO_obj to attach to SC_obj.
The HTO_obj contains the formatting and transformations necessary to perfrom
demultiplexing using MULTI-seq by Chris Mcginnis}
\usage{
\S4method{HTO.build}{SC_obj}(x)
}
\arguments{
\item{x}{SC_obj.}
}
\value{
an SC_obj containing the built HTO_obj
}
\description{
HTO.build
Build an HTO_obj to attach to SC_obj.
The HTO_obj contains the formatting and transformations necessary to perfrom
demultiplexing using MULTI-seq by Chris Mcginnis
}
\examples{
HTO.build(sc_obj)
}
