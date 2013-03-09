\name{H2vec}
\alias{H2vec}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
	data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
	A numeric list the same length as \code{corlist} of weights for each correlation matrix, sum should be between 0 and 1.
}
\usage{data(H2vec)}
\format{
  The format is:
 num [1:2] 0.5 0.5
}
\examples{
data(H2vec)
str(H2vec)
length(H2vec)
}
\keyword{datasets}
