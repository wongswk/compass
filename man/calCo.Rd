% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calCo}
\alias{calCo}
\title{Calculate Cartesian coordinates for 1-4 bonded atoms}
\usage{
calCo(prev_atoms, length, bAngle, tAngle)
}
\arguments{
\item{prev_atoms}{a 3x3 matrix of coordinates for atoms 1-2-3 in dihedral, listed by row}

\item{length}{bond length between atoms 3-4 in dihedral}

\item{bAngle}{planar bond angle between atoms 2-3-4}

\item{tAngle}{dihedral angle formed by atoms 1-2-3-4}
}
\value{
Returns the vector of coordinates for the fourth atom in the dihedral
}
\description{
Given coordinates for the first three atoms forming a dihedral,
the dihedral angle, bond angle, and bond length, calculate the coordinates
of the fourth atom forming the dihedral.
}
