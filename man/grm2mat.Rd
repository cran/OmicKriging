\name{grm2mat}
\alias{grm2mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
	read one set of grm files into a correlation matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	This function inputs the data from filename.grm.gz and filename.grm.id and converts it into an n x n correlation matrix appropriate for \code{\link{okriging}}.
}
\usage{
grm2mat(grmfullheader)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grmfullheader}{
%%     ~~Describe \code{grmfullheader} here~~
	path to the .grm.gz and grm.id files
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	matrix of size n x n, where n is the number of samples in filename.grm.id

}
\references{
%% ~put references to the literature/web site here ~
        Wheeler HE, et al. (2013) Poly-Omic Prediction of Complex Traits: OmicKriging. arXiv:1303.1788 
        http://arxiv.org/abs/1303.1788
}
\author{
%%  ~~who you are~~
	Hae Kyung Im

	Heather E. Wheeler
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
	\code{\link{readcorlist}} \code{\link{okriging}}
}
%\examples{
%}
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (grmfullheader) 
%{
%    grm = read.table(grmfullheader \%&\% ".grm.gz", header = F, 
%        as.is = T)
%    grmid = read.table(grmfullheader \%&\% ".grm.id", header = F, 
%        as.is = T)
%    names(grmid) = c("FID", "IID")
%    grmid$IID = as.character(grmid$IID)
%    nn = nrow(grmid)
%    grmat = matrix(0, nn, nn)
%    indexfun = function(i, j, nn) return(i + (j - 1) * nn)
%    grmat[indexfun(grm[, 1], grm[, 2], nn)] = grm[, 4]
%    grmat = t(lower.tri(grmat) * grmat) + grmat
%   rownames(grmat) = grmid$IID
%    colnames(grmat) = grmid$IID
%    return(grmat)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
