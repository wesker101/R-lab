\name{minorAlleleFrequency}
\alias{minorAlleleFrequency}
\alias{filterSNPs}
\alias{mafHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation and Visualization of Minor Allele Frequency
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
minorAlleleFrequency(SNPmat)
filterSNPs(SNPmat, minMAF = 0.05)
mafHist(SNPmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPmat}{
%%     ~~Describe \code{SNPmat} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{minorAlleleFrequency }{returns a vector of allele frequencies}
%%  \item{filterSNPs }{returns a SNP matrix in the same format as SNPmat but with columns removed for SNPs with very rare alleles}
%%  \item{mafHis }{draws a plot and doesn’t return anything}
%% ...
}

\author{
Kevin Mitnick
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SNPmat) 
{
    SNPmat <- matrix(1:4, nrow = 4, ncol = 4, byrow = FALSE,
       dimnames = NULL)
       minorAlleleFrequency(SNPmat)
       filterSNPs(SNPmat)
       mafHist(SNPmat)
  }
}