\name{readcorlist}
\alias{readcorlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
	read multiple sets of grm files into a list of correlation matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	This function takes a list of grm file paths and converts the data into a list of n x n correlation matrices appropirate for \code{\link{okriging}}.
}
\usage{
readcorlist(corfilelist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corfilelist}{
%%     ~~Describe \code{corfilelist} here~~
	list of paths to grm files (e.g. corfilelist <- c("/user/dir/genetic","/user/dir/geneexp"))	
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
	list of n x n correlation matrices, where n is the number of samples in filename.grm.id
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
	\code{\link{grm2mat}} \code{\link{okriging}}
}
%\examples{
%}
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (corfilelist) 
%{
%    nc = length(corfilelist)
%    corlist = list()
%    for (cc in 1:nc) {
%        corlist[[cc]] = grm2mat(corfilelist[cc])
%    }
%    return(corlist)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
