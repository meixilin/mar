% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fileparser.R
\name{lonlat_parser}
\alias{lonlat_parser}
\title{Parse longitude/latitude coordinates from text file}
\usage{
lonlat_parser(lonlat.fn)
}
\arguments{
\item{lonlat.fn}{Path to input file (txt/csv/tsv, can be gzipped) containing coordinates. No missing values allowed.}
}
\value{
A data frame containing sample IDs and their corresponding longitude/latitude coordinates.
Returns error if coordinates contain missing values or incorrect number of columns.
}
\description{
Reads a file containing sample IDs with their corresponding longitude and latitude coordinates.
The input file must have a header with columns: ID, LON/LONGITUDE, LAT/LATITUDE (in that order).
While coordinates do not need to be in "+proj=longlat +datum=WGS84" projection, the WGS84 projection was used in testing..
Sample IDs must be unique and in the same order as the Sample IDs provided in the genotype matrix.
}
\examples{
\dontrun{
# Read coordinates from file
coords <- lonlat_parser("sample_locations.txt")
}
}
