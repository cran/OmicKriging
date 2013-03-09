\name{okriging}
\alias{okriging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
	predict the phenotype of an individual using the remaining individuals	
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	Given "omics" correlation matrices, matrix weights, and covariates, \code{okriging} predicts the phenotype of an individual by using the phenotypes of the remaining individuals.
}
\usage{
okriging(idtest, idtrain = NULL, corlist, H2vec, 
pheno, phenoname, Xcova = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{idtest}{
%%     ~~Describe \code{idtest} here~~
	IID of individual(s) in the test set for phenotype prediction
}
  \item{idtrain}{
%%     ~~Describe \code{idtrain} here~~
	IID list of the remaining individuals in the training set for phenotype prediction
}
  \item{corlist}{
%%     ~~Describe \code{corlist} here~~
	list of correlation matrices (rownames and colnames should be the IID list and include \code{idtest} and \code{idtrain}
}
  \item{H2vec}{
%%     ~~Describe \code{H2vec} here~~
	list same length as \code{corlist} of weights for each correlation matrix, sum should equal 1
}
  \item{pheno}{
%%     ~~Describe \code{pheno} here~~
	data.frame of the phenotype (rownames should be the IID list and colnames should be: FID, IID, \code{phenoname})	
}
  \item{phenoname}{
%%     ~~Describe \code{phenoname} here~~
	name of the phenotype
}
  \item{Xcova}{
%%     ~~Describe \code{Xcova} here~~
	matrix of covariates  with the IID list as the rownames
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
%% ~put references to the literature/web site here ~
        Wheeler HE, et al. (2013) Poly-Omic Prediction of Complex Traits: OmicKriging. arXiv:1303.1788 
        http://arxiv.org/abs/1303.1788

	Cressie, N. (1993) Statistics for Spatial Data, Wiley, New York. p.154

	Stein, M. L. (1999) Interpolation of Spatial Data: Some Theory for Kriging, Springer, New York.
	
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
	\code{\link{readcorlist}} \code{\link{grm2mat}}
}
\examples{
data(iddata)
data(corlist)
data(H2vec)
data(pheno)
data(phenoname)
data(Xcovamat)

idfamilylist=unique(iddata$FID)
pred = data.frame()
for(idtrio in idfamilylist){
  idtest = iddata$IID[iddata$FID ==idtrio]
  idtrain = iddata$IID[!(iddata$IID \%in\% idtest)]
  res = okriging(idtest,idtrain,corlist,H2vec,pheno,phenoname,Xcova = Xcovamat)
  pred = rbind(pred,res)
  print(res)
}
summary(lm(Ypred~Ytest,data=pred))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
