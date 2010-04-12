\name{GADEM}
\alias{GADEM}

\title{Motif Analysis with rGADEM}
\description{It is an R implementation of GADEM, a powerful computational tools for de novo motif discovery.}

\usage{

gadem<-GADEM(Sequences,seed=1,genome=NULL,verbose=TRUE,numWordGroup=3,numTop3mer=20,numTop4mer=40,numTop5mer=60,numGeneration=5,
	populationSize=100,pValue=0.0002,eValue=0.0,extTrim=1,minSpaceWidth=0,maxSpaceWidth=10,useChIPscore=0,numEM=40,fEM=0.5,widthWt=80,fullScan=0,userBackgModel=0,slideWinPWM=6,stopCriterion="NUM_NO_MOTIF",MarkovOrder=0,userMarkovOrder=0,numBackgSets=10,weightType=0,pgf=1,startPWMfound=0,bOrder=-1,bFileName="NULL",Spwm="NULL")}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{Sequences}{Sequences from BED or FASTA file are converted into XString object view}
	\item{seed}{When a seed is specified, the run results are deterministic}
	\item{genome}{Specify the genome}
	\item{verbose}{Print immediate results on screen [TRUE-yes (default), FALSE-no]. These results include the motif consensus sequence, number of sites (in sequences subjected  to EM optimization, see -fEM, above), and ln(E-value).}
	\item{numWordGroup}{number of non-zero k-mer groups}
	\item{numTop3mer}{Number of top-ranked trimers for spaced dyads (default: 20).}
	\item{numTop4mer}{Number of top-ranked tetramers for spaced dyads (default: 40).}
	\item{numTop5mer}{Number of top-ranked pentamers for spaced dyads (default: 60).}
	\item{numGeneration}{Number of genetic algorithm (GA) generations (default: 5).}
	\item{populationSize}{GA population size (default: 100).  Both default settings should work well for most datasets (ChIP-chip and ChIP-seq).  The above two arguments are ignored in a seeded analysis, because spaced dyads and GA are no longer needed (-gen is set to 1 and -pop is set to 10 internally, corresponding to the 10 maxp choices).}
	\item{pValue}{P-value cutoff for declaring BINDING SITES (default: 0.0002). Depending on data size and the motif, you might want to assess more than one value. For ChIP-seq data (e.g., 10 thousand +/-200-bp max-center peak cores), p=0.0002 often seems appropriate. However, short motifs may require a less stringent setting.}
	\item{eValue}{ln(E-value) cutoff for selecting MOTIFS (default: 0.0). If a seeded analysis fails to identify the expected motif, run GADEM with  -verbose 1 to show motif ln(E-value)s on screen, then rerun with a larger ln(E-value) cutoff. This can help in identifying short and/or low abundance motifs, for which the default E-value threshold may be too low.}
	\item{extTrim}{Base extension and trimming (1 -yes, 0 -no) (default: 1).}
	\item{minSpaceWidth}{Minimal number of unspecified nucleotides in spaced dyads (default: 0).}
	\item{maxSpaceWidth}{Maximal number of unspecified nucleotides in spaced dyads (default: 10). -mingap and -maxgap control the lengths of spaced dyads, and, with -extrim, control motif lengths. Longer motifs can be discovered by setting -maxgap to larger values (e.g. 50).}
	\item{useChIPscore}{Use top-scoring sequences for deriving PWMs. Sequence (quality) scores are stored in sequence header (see documentation). 0 - no (default, randomly select sequences), 1 - yes.}
	\item{numEM}{Number of EM steps (default: 40). One might want to set it to a larger value (e.g. 80) in a seeded run, because such runs are fast.}
	\item{fEM}{Fraction of sequences used in EM to obtain PWMs in an unseeded analysis (default: 0.5). For unseeded motif discovery in a large dataset (e.g. >10 million nt), one might want to set -fEM to a smaller value (e.g., 0.3 or 0.4) to reduce run time.}
	 \item{widthWt}{For -posWt 1 or 3, width of central sequence region with large EM weights for PWM optimization (default: 50). This argument is ignored when -posWt is 0 (uniform prior) or 2 (Gaussian prior).}
	\item{fullScan}{GADEM keeps two copies of the input sequences internally: one (D) for discovering PWMs and one (S) for scanning for binding sites using the PWMs Once a motif is identified, its instances in set D are always masked by Ns. However, masking motif instances in set S is optional, and scanning unmasked equences allows sites of discovered motifs to overlap.}
	 \item{userBackgModel}{To run analysis in background (default : 0).}
	 \item{slideWinPWM}{sliding window for comparing pwm similarity (default : 6).}
 	 \item{stopCriterion}{Stop analysis.}
	\item{MarkovOrder}{Background Markov order,user-specified order highest order available in user-specified background indicator (default : 0).}
	\item{userMarkovOrder}{Background Markov order,user-specified order highest order available in user-specified background indicator (default : 0).}
	\item{numBackgSets}{Number of sets of background sequences (default: 10). The background sequences are simulated using the [a,c,g,t] frequencies in the input sequences, with length matched between the two sets. The background sequences are used as the random sequences for assessing motif enrichment in the input data. Another set (same default: 10) of background sequences is independently generated to approximate the empirical llr score distribution when -pgf is set to 0.}
	\item{weightType}{Weight profile for positions on the sequence. 0 - no weight (uniform spatial prior, default), 1 - small or zero weights for the ends and large weights for the center (e.g. the center 50 bp). If you expect strong central enrichment (as in ChIP-seq) and your sequences are long (e.g. >200 bp), choose type 1.}
	\item{pgf}{By default, GADEM uses the Staden probability generating function (pgf) method to approximate the exact llr score null distribution.}
	\item{startPWMfound}{Value for the PWM (default : 0).}
	\item{bOrder}{The order of the background Markov model for computing llr scores: 0 - 0th 1 - 1st 2 - 2nd 8 - 8th}
	\item{bFileName}{Reading user-specified background models.}
	\item{Spwm}{File name for the seed PWM, when a seeded approach is used. can be used as the starting PWM for the EM algorithm. This will help find an expected motif and is much faster than unseeded de novo discovery.Also, when a seed PWM is specified, the run results are deterministic, so only a single run is needed (repeat runs with the same settings will give identical results). In contrast, unseeded runs are stochastic, and we recommend comparingresults from several repeat runs.}
 }

\author{Arnaud Droit \email{arnaud.droit@ircm.qc.ca}}

\examples{

	library(BSgenome.Hsapiens.UCSC.hg18)
	pwd<-"" #INPUT FILES- BedFiles, FASTA, etc.
	path<- system.file("extdata/Test_100.bed",package="rGADEM")
	BedFile<-paste(pwd,path,sep="")
	BED<-read.table(BedFile,header=FALSE,sep="\t")
	BED<-data.frame(chr=as.factor(BED[,1]),start=as.numeric(BED[,2]),end=as.numeric(BED[,3]))
	#Create RD files 
	rgBED<-IRanges(start=BED[,2],end=BED[,3])
	Sequences<-RangedData(rgBED,space=BED[,1])
	
	
gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens)
}
\keyword{GADEM}
\keyword{MOTIFS}
 
