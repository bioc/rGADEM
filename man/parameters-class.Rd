\name{parameters-class}
\Rdversion{1.1}
\docType{class}
\alias{parameters}
\alias{parameters-class}


\title{Class "parameters" }
\description{This object contains contains parameters of GADEM analysis}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("motif_gadem", ...)}.
}
\section{Slots}{
	 \describe{
    \item{numWordGroup}{:number of non-zero k-mer groups. }
	\item{numTop3mer}{:Number of top-ranked trimers for spaced dyads (default: 20). }
	\item{verbose}{:Print immediate results on screen [1-yes (default), 0-no].}
	\item{numTop4mer}{:Number of top-ranked tetramers for spaced dyads (default: 40).}
	\item{numTop5mer}{:Number of top-ranked pentamers for spaced dyads (default: 60).}
	\item{numGeneration}{:Number of genetic algorithm (GA) generations (default: 5).}
	\item{populationSize}{:GA population size (default: 100).}
	\item{pValue}{:P-value cutoff for declaring BINDING SITES (default: 0.0002). }
	\item{eValue}{:ln(E-value) cutoff for selecting MOTIFS (default: 0.0).}
	\item{extTrim}{:Base extension and trimming (1 -yes, 0 -no) (default: 1).}
	\item{minSpaceWidth}{:Minimal number of unspecified nucleotides in spaced dyads (default: 0).}
	\item{maxSpaceWidth}{:Maximal number of unspecified nucleotides in spaced dyads (default: 10).}
	\item{useChIPscore}{:Use top-scoring sequences for deriving PWMs.}
	\item{numEM}{:Number of EM steps (default: 40).}
	\item{fEM}{:Fraction of sequences used in EM to obtain PWMs in an unseeded analysis (default: 0.5).}
	\item{widthWt}{:For -posWt 1 or 3, width of central sequence region with large EM weights for PWM optimization (default: 50).}
	\item{fullScan}{:GADEM keeps two copies of the input sequences internally.}
	\item{userBackgModel}{:To run analysis in background (default : 0).}
	\item{slideWinPWM}{:sliding window for comparing pwm similarity (default : 6).}
	\item{stopCriterion}{}
	\item{MarkovOrder}{:Background Markov order,user-specified order highest order available in user-specified background indicator (default : 0).}
	\item{userMarkovOrder}{:Background Markov order,user-specified order highest order available in user-specified background indicator (default : 0).}
	\item{numBackgSets}{:Number of sets of background sequences (default: 10).}
	\item{weightType}{:Weight profile for positions on the sequence.}
	\item{pgf}{:By default, GADEM uses the Staden probability generating function (pgf) method to approximate the exact llr score null distribution.}
	\item{startPWMfound}{:Value for the PWM (default : 0).}
	\item{bOrder}{:The order of the background Markov model for computing llr scores: 0 - 0th 1 - 1st 2 - 2nd 8 - 8th}
	\item{bFileName}{:Reading user-specified background models.}
	\item{fpwm0}{:File name for the seed PWM, when a seeded approach is used.}
	
	
  }
}

\author{Arnaud Droit \email{arnaud.droit@ircm.qc.ca}}

\seealso{
	 \code{\linkS4class{gadem}} , \code{\linkS4class{align}}, \code{\linkS4class{motif}}
}
\examples{
showClass("parameters")
}
\keyword{classes}