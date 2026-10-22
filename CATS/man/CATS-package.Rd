\name{CATS-package}
\alias{CATS-package}
\alias{CATS}
\docType{package}
\title{
Power calculation for a joint analysis of a two-stage case control design for SNP data
}
\description{
Power calculation for a joint analysis of a two-stage case control design for SNP data. 
}
\details{
\tabular{ll}{
Package: \tab CATS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-10-16\cr
License: \tab GPL
}
The function cats calculates the power
}
\author{
Anders Albrechtsen

Maintainer: Anders Albrechtsen <albrecht@binf.ku.dk>

}
\references{
 Skol AD, Scott LJ, Abecasis GR, Boehnke M: Joint analysis is more efficient than replication-based analysis for two-stage genome-wide association studies. Nat Genet 38: 209-213, 2006.
}
\keyword{ htest }
\seealso{
http://www.sph.umich.edu/csg/abecasis/CaTS/
}
\examples{
#calculate the power under a multiplicative model using a two stage design and 
#assuming a relative risk of 1.5
cats(freq=0.2,ncases=500,ncases2=500,ncontrols=1000,ncontrols2=1000,
risk=1.5,multiplicative=1)
}
