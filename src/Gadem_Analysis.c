#include "config.h"
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <Rmath.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <unistd.h>
#include "gadem.h"

#include "defines.h"
#include "evalue_meme.h"
#include <Rdefines.h>
#include <Rversion.h>

#if defined(__APPLE__) || defined(macintosh)
#ifdef __MAC_OS_X_VERSION_MIN_REQUIRED
#if __MAC_OS_X_VERSION_MIN_REQUIRED >= 1060
#include <dispatch/dispatch.h>
//#define HAVE_DISPATCH 1
#endif
#endif
#endif

#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif


  // last modification 9/07/2009
  // modifications:
  //   1) fixed a minor bug in scan_sites.c
  //   2) remove 6-mers, reduce search space
  //   3) added C function getpid (process id)
  //   4) added C function for cpu running time
  //   5) set equal mutation rate for maxp and spaced dyads.
  //      An optimal maxp may be important as it may affect how the EM converges
  //   6) some cosmetic changes in output 
  //   7) set default number of generations to 10
  //   8) allowed a user-specified "seed" PWM
  //   9) allowed a user-specified background model
  //  10) included enrichment analysis
  //  11) re-wrote pgf function (Staden's pgf method for llr null distribution)
  //  12) fixed a bug in computing marginal probabilities for subsequences containing non-[a,c,g,t]
  //  13) allow motif to overlap as an option


void populationCalculation(int maxSeqLen, int numEM, Fitness *fitness, int startPWMfound, int minminSites, double maxpFactor, int numSeq, int numSeqEM, char **seq, char **rseq, int *seqLen, char *Iseq, double *bfreq, double ** posWeight, int weightType, double pvalueCutoff, int* emSeqLen, double **pwm, int pwmLen, double **epwm, double **opwm,char *pwmConsensus, int *scoreCutoff,char *sdyad, int ii);

void populationCalculationSequential(int maxSeqLen, int numEM, Fitness *fitness, int startPWMfound, int minminSites, double *maxpFactor, int numSeq, int numSeqEM, char **seq, char **rseq, int *seqLen, char *Iseq, double *bfreq, double ** posWeight, int weightType, Pgfs *llrDist, int llrDim, double pvalueCutoff, int* emSeqLen, double ***pwm, int *pwmLen, double ***epwm, double ***opwm,char **pwmConsensus, int *scoreCutoff,char **sdyad, double **t1pwm, double **t2pwm, Sites *siteEM, double **logpwm, double** score, double** rscore, int **ipwm, int ii, int verbose, int jjj, int numCycle, int maxSpaceWidth);


SEXP GADEM_Analysis(SEXP sequence,SEXP sizeSeq, SEXP accession, SEXP Rverbose,SEXP RnumWordGroup,SEXP RnumTop3mer,SEXP RnumTop4mer,SEXP RnumTop5mer,SEXP RnumGeneration,SEXP RpopulationSize, SEXP RpValue,SEXP ReValue,SEXP RextTrim,SEXP RminSpaceWidth,SEXP RmaxSpaceWidth,SEXP RuseChIPscore,SEXP RnumEM,SEXP RfEM, SEXP RwidthWt,SEXP RfullScan,SEXP RuserBackgModel, SEXP RslideWinPWM,SEXP RstopCriterion,SEXP RMarkovOrder,SEXP RuserMarkovOrder,SEXP RnumBackgSets,SEXP RweightType,SEXP Rpgf,SEXP RstartPWMfound,SEXP RbOrder,SEXP RbFileName,SEXP RListPWM) 
{
  int jjj,ii,i,j,k;  
  char *bFileName;
  
  SEXP ResultsGadem;
  SEXP RSpwm;
  PROTECT(ResultsGadem=NEW_LIST(100));

#if defined(CSTACK_DEFNS) & defined(HAVE_DISPATCH)
  R_CStackLimit = (uintptr_t)-1;
#endif
  
#ifdef HAVE_DISPATCH
  dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
#endif
  
  
  int increment=0;
  
  GetRNGstate();
  double testrand;
  
    //Number of sequences
  int numSeq = INTEGER_VALUE(sizeSeq);
  // const
//  char *Fastaheader[size];
  int incr=0;
  
  int longueur=length(sequence);
  int IncrementTemp=0;
  
    // basic settings/info
  int maxSeqLen,*seqLen;          // sequence info
  double aveSeqLen;                      // sequence info
  char **seq,**rseq;
  char **oseq,**orseq;                   // copy of the original sequences
  char **pseq,**rpseq;                   // permuted seqs.
  double *bfreq;                         // base frequencies
  double *ChIPScore;                     // chip score
    
    // pwms
  double ***pwm;                         // initial population of PWMs from spaced dyads
  int *pwmLen;                           // initial pwm lengths 
  double **t1pwm,**t2pwm;                // two pwms before and after EM steps
  double **opwm2;                        // EM-derived PWM 
  double **logpwm;                       // log-transformed EM-derived EM
  int **ipwm;                            // integer pwm for computing llr score distribution
  double ***opwm;                        // observed PWMs from identified sites
  double ***epwm;                        // em-optimized PWMs
  double **logepwm;                      // log(em-optimized PWM)
  int *pwmnewLen;                        // final motif length after extending to both ends
  
    // llr score distr.
  Pgfs *llrDist;                         // llr distribution from pgf method
  int llrDim;                            // llr distribution dimension
  double *empDist;                       // empirical llr score distr.
  int empDim;                            // dimension of empirical llr score distribution 
  int pgf;                               // indicator for using pgf method or not
  
    // EM, motif, sites
  double pvalueCutoff;                   // user input, used to determine score cutoff based on ipwm
  int *scoreCutoff;                      // pwm score cutoff for the corresponding p-value cutoff
  double llrCutoff;
  double **score,**rscore;               // subsequence score, plus and minus strands
  double logev;                          // log of E-value of a motif;
  int useChIPscore;                      // indicator for using ChIP-seq score for seq. selection for EM
  int numEM;                             // number of EM steps
  double E_valueCutoff;                  // log E-value cutoff
  int minsitesEM;                        // minimal number of sites in a motif in EM sequences
  Sites *siteEM;                         // binding sites in EM sequences
  int *nsites;                           // number of binding sites in full data
  int minsites;                          // minimal number of sites in a motif in full data
  Sites **site;                          // binding sites in all sequences
  int motifCn;                           // number of motifs sought and found
  int extTrim;
  int noMotifFound;                      // none of the dyads in the population resulted in a motif
  char **pwmConsensus;                   // consensus sequences of motifs
  double pwmDistCutoff;                  // test statistic for motif pwm similarity
  char *uniqMotif;                       // motifs in a population unique or not
  int numUniq;                           // number of unique motifs in a population
  int slideWinPWM;                       // sliding window for comparing pwm similarity
  int widthWt;                           // window width in which nucleotides are given large weights for PWM optimization 
  int fullScan;                          // scan scan on the original sequences or masked sequences
  
    // background
  BACKGROUND_Model *back;                // background model either user-specified or computed from forward data
  double **bscore,**rbscore;             // likelihood scores using background model
  double *pscore;                        // score for permuted seqs.
  int numTopWmerInB,numWmerInB;
  int numBackgSets;
  int MarkovOrder,userMarkovOrder;       // background Markov order,user-specified order
  int maxUserMarkovOrder;                // highest order available in user-specified background 
  int userBackgModel;                    // indicator
  
    // weights 
  double **posWeight;                    // spatial weights
  int weightType;                        // four weight types 0, 1, 2, 3, or 4
  
    // words for spaced dyad
  Words *word;                           // top-ranked k-mers as the words for spaced dyads
  int numTop3mer,numTop4mer,numTop5mer;  // No. of top-ranked k-mers as words for dyads
  int maxWordSize;                       // max of the above three
  int numWordGroup;                      // number of non-zero k-mer groups
  int minSpaceWidth,maxSpaceWidth;       // min and max width of spacer of the spaced dyads
  Chrs **dyad;                           // initial population of "chromosomes"
  char **sdyad;                          // char of spaced dyads
  
    // GA
  int populationSize,numGeneration;      // GA parameters
  double maxpMutationRate;
  Fitness *fitness;                      // "chromosome" fitness
  Wheel *wheel;                          // roulette-wheel selection
  
    // to speed up only select a subset of sequences for EM algorithm
  double fEM;                            // percentage of sequences used in EM algorithm
  int numSeqEM;                          // number of sequences subject to EM
  char *Iseq;                            // Indicator if a sequence is used in EM or not
  int *emSeqLen;                         // length of sequences used in EM
  double *maxpFactor;
  
  int numCycle;                          // number of GADEM cycles
  int generationNoMotif;                 // maximal number of GA generations in a GADEM cycle resulted in no motifs
  
    // mis.
    //seed_t  seed;                          // random seed
  int motifCn2,id,numCycleNoMotif,verbose,bOrder,minminSites;
  int startPWMfound,stopCriterion;
  char *mFileName,*oFileName,*pwmFileName;
  time_t start;
  int cn[4],bcn[4],*seqCn,*bseqCn,avebnsites,avebnsiteSeq,totalSitesInput;
  
  
  mFileName=alloc_char(500);     mFileName[0]='\0';
  oFileName=alloc_char(500);     oFileName[0]='\0';
  pwmFileName=alloc_char(500);   pwmFileName[0]='\0';
  bFileName=alloc_char(500);     bFileName[0]='\0';
  seq=NULL; aveSeqLen=0; maxSeqLen=0; 
  minsites=-1; 
  
    
  maxSeqLen=0;
  for(incr=1;incr<longueur;incr=incr+2)
  {	
    if (length(STRING_ELT(sequence,(incr)))>maxSeqLen) maxSeqLen=length(STRING_ELT(sequence,(incr))); 
  }
//  printf("maxLength=%d",maxSeqLen);
//  exit(0);
  seq=alloc_char_char(numSeq,maxSeqLen+1);
  for(incr=1;incr<longueur;incr=incr+2)
  {	
    for (j=0; j<length(STRING_ELT(sequence,(incr))); j++)
    {
      seq[IncrementTemp][j]=CHAR(STRING_ELT(sequence,(incr)))[j];
    }
    IncrementTemp++;
  }
  
  
  verbose=LOGICAL_VALUE(Rverbose);
  numWordGroup=INTEGER_VALUE(RnumWordGroup);
  numTop3mer=INTEGER_VALUE(RnumTop3mer);
  numTop4mer=INTEGER_VALUE(RnumTop4mer);
  numTop5mer=INTEGER_VALUE(RnumTop5mer);
  numGeneration=INTEGER_VALUE(RnumGeneration);
  populationSize=INTEGER_VALUE(RpopulationSize);
  pvalueCutoff=NUMERIC_VALUE(RpValue);
  E_valueCutoff=NUMERIC_VALUE(ReValue);
  extTrim=INTEGER_VALUE(RextTrim);
  minSpaceWidth=INTEGER_VALUE(RminSpaceWidth);
  maxSpaceWidth=INTEGER_VALUE(RmaxSpaceWidth);
  useChIPscore=NUMERIC_VALUE(RuseChIPscore);
  numEM=INTEGER_VALUE(RnumEM);
  fEM=NUMERIC_VALUE(RfEM);
  widthWt=INTEGER_VALUE(RwidthWt);
  fullScan=INTEGER_VALUE(RfullScan);
  userBackgModel=INTEGER_VALUE(RuserBackgModel);
  slideWinPWM=INTEGER_VALUE(RslideWinPWM);
  numUniq=populationSize;
  stopCriterion=NUM_NO_MOTIF;  
  MarkovOrder=INTEGER_VALUE(RMarkovOrder);
  userMarkovOrder=INTEGER_VALUE(RuserMarkovOrder);
  numBackgSets=INTEGER_VALUE(RnumBackgSets);
  weightType=NUMERIC_VALUE(RweightType);
  pgf=INTEGER_VALUE(Rpgf);
  startPWMfound=INTEGER_VALUE(RstartPWMfound);
  bOrder=INTEGER_VALUE(RbOrder);
  const char *tempRbFileName[1];
  tempRbFileName[0]=CHAR(STRING_ELT(RbFileName,0));
    
  if(numSeq>MAX_NUM_SEQ)
  {
    printf("Error: maximal number of seqences reached! \n");
    printf("Please reset MAX_NUM_SEQ in gadem.h and rebuild (see installation)\n\n");
    exit(0);
  }
  
  strcpy(bFileName,tempRbFileName[0]);
  back=alloc_background();
  
  ChIPScore=alloc_double(MAX_NUM_SEQ);
  seqLen=alloc_int(MAX_NUM_SEQ); 
  
//  seq=sequences;
  
//  numSeq=size;
  int len; 
  
  for (i=0; i<numSeq; i++) {
    len=strlen(seq[i]); 
    seqLen[i]=len;
  }

  aveSeqLen=0; for (i=0; i<numSeq; i++) aveSeqLen +=seqLen[i]; aveSeqLen /=(double)numSeq;
  
  for (i=0; i<numSeq; i++) {
    if (seqLen[i]>maxSeqLen) maxSeqLen=seqLen[i]; 
  }
  
  rseq=alloc_char_char(numSeq,maxSeqLen+1);
  oseq=alloc_char_char(numSeq,maxSeqLen+1);
  orseq=alloc_char_char(numSeq,maxSeqLen+1);
  
  for (i=0; i<numSeq; i++)
  {
    if(seqLen[i]>maxSeqLen) maxSeqLen=seqLen[i]; 
  }
  
  reverse_seq(seq,rseq,numSeq,seqLen);
  
    // make a copy of the original sequences both strands
  for (i=0; i<numSeq; i++)
  {
    for (j=0; j<seqLen[i]; j++)
    {
      oseq[i][j]=seq[i][j];
      orseq[i][j]=rseq[i][j];
    }
    oseq[i][seqLen[i]]='\0'; orseq[i][seqLen[i]]='\0'; 
  }
    
  if (strcmp(bFileName,"NULL")!= 0)
  {
    maxUserMarkovOrder=read_userBackgModel(bFileName,back);
    if (userMarkovOrder>maxUserMarkovOrder)
    {
      userMarkovOrder=maxUserMarkovOrder;
    }
    userBackgModel=1;
  }
  else if (GET_LENGTH(RListPWM)!= 0)
  {
    startPWMfound=1; 
  }
  else { }
  
    // check for input parameters
  if(numGeneration<1)
  { 
    printf("\nError: numbe of generaton < 1.\n");
    exit(0);
  }
  if(populationSize<1)
  {
    printf("\nError: population size < 1.\n");
    exit(0);
  }
  if (minSpaceWidth<0)
  { 
    printf("\nError: minimal number of unspecified bases in spaced dyads <0.\n"); 
    printf("   check -mingap setting\n\n");
    exit(0);
  }
  if (maxSpaceWidth<0)
  { 
    printf("\nError: maximal number of unspecified bases in spaced dyads <0.\n"); 
    printf("   check -maxgap setting\n\n");
    exit(0);
  }
  if (minSpaceWidth>maxSpaceWidth)
  {
    printf("\nError: mingap setting must <= to maxgap setting.\n\n"); 
    exit(0); 
  }
  if (maxSpaceWidth+12>MAX_PWM_LENGTH)
  {
    printf("\nError: maxgap setting plus word lengths exceed <MAX_PWM_LENGTH>.\n");
    printf("   For very long motifs, please set <MAX_PWM_LENGTH> in 'defines.h' accordingly.\n\n");
    exit(0); 
  }
  if (numEM<0)
  {
    printf("\nError: number of EM steps is zero.\n");
    exit(0);
  }
  if (numEM==0)
  {
    printf("\nNote: number of EM steps = 0, no EM optimization is carried out.\n");
  }
  
  if (fullScan!=0 && fullScan!=1)
    fullScan=0;
  
  
  maxWordSize=0;
  if (numTop3mer>maxWordSize) maxWordSize=numTop3mer;
  if (numTop4mer>maxWordSize) maxWordSize=numTop4mer;
  if (numTop5mer>maxWordSize) maxWordSize=numTop5mer;
  
    // any one, two or three: tetramer, pentamer, hexamer
  if (numTop3mer==0 && numTop4mer==0 && numTop5mer==0)
  {
    printf("\nError: maxw3, maxw4, and maxw5 all zero - no words for spaced dyads.\n");
    exit(0);
  }
  
  if (startPWMfound && fEM!=0.5 && fEM!=1.0)
  {
    printf("\n***Note: -fEM argument is ignored in a seeded analysis***\n\n");
  }
  
  if (startPWMfound)
  {
    if(verbose)
    {
      if (populationSize!=10 && populationSize!=100) printf("\n***Note: -pop argument is ignored in a seeded analysis, -pop is set to 10.***\n\n");
      if (numGeneration!=1 && numGeneration!=5)      printf("\n***Note: -gen argument is ignored in a seeded analysis, -gen is set to 1.***\n\n");
    }
    fEM=1.0;
    populationSize=FIXED_POPULATION; numGeneration=1; 
  }
  
    // number of sequences for EM
  if (fEM>1.0 || fEM<=0.0)
  { 
    printf("\nError: the fraction of sequences subject to EM is %3.2f.\n",fEM);
    exit(0); 
  } 
  numSeqEM=(int)(fEM*numSeq);
  
  if (pgf!=0 && pgf!=1)   
  {
    printf("\nError: -pgf can only take 0 or 1.\n"); 
    exit(0); 
  }


    // memory callocations
  Iseq  =alloc_char(numSeq+1); 
  opwm2 =alloc_double_double(MAX_PWM_LENGTH,4);
  t1pwm =alloc_double_double(MAX_PWM_LENGTH,4);
  t2pwm =alloc_double_double(MAX_PWM_LENGTH,4);
  ipwm  =alloc_int_int(MAX_PWM_LENGTH,4);
  logpwm=alloc_double_double(MAX_PWM_LENGTH,4);
  score =alloc_double_double(numSeq,maxSeqLen);
  rscore=alloc_double_double(numSeq,maxSeqLen);
  logepwm=alloc_double_double(MAX_PWM_LENGTH,4);
  siteEM=alloc_site(MAX_SITES);
  emSeqLen=alloc_int(numSeqEM);
  scoreCutoff=alloc_int(1000);
  // scoreCutoff=alloc_int(populationSize);
  llrDist=alloc_distr(MAX_DIMENSION);
  posWeight=alloc_double_double(numSeq,maxSeqLen);
  bscore=alloc_double_double(numSeq,maxSeqLen);
  rbscore=alloc_double_double(numSeq,maxSeqLen);
  pscore=alloc_double(numSeq*maxSeqLen);
  empDist=alloc_double((int)(numSeq*2*maxSeqLen*numBackgSets/4));
  pseq=alloc_char_char(MAX_NUM_SEQ,maxSeqLen+1);
  rpseq=alloc_char_char(MAX_NUM_SEQ,maxSeqLen+1);
  bfreq=base_frequency(numSeq,seq,seqLen);
  

    // if minN not specified, set the defaults accordingly
  if (minsites==-1) 
  {
    minsites =max(2,(int)(numSeq/20)); 
  }
  minsitesEM=(int)(fEM*minsites);
  
  maxpMutationRate=MAXP_MUTATION_RATE;
  
    // determine the distribution and critical cut point
  pwmDistCutoff=vector_similarity();
  
  /*---------- select a subset of sequences for EM only --------------*/
  if (useChIPscore==1)
  {
    select_high_scoring_seq_for_EM (ChIPScore,numSeq,numSeqEM,Iseq,fEM);
  }
  else
  {
    sample_without_replacement(Iseq,numSeqEM,numSeq);
  }
  /*-------------------- end of selection --------------------------*/
  
  if (!userBackgModel & verbose)
  {
    if (!pgf && userMarkovOrder!=0 && aveSeqLen<=500)
    {
      printf("\n***Note: it is not recommended to use a non-zero-th Markov model estimated from\n");
      printf("   the input sequences as the background model, especially for short sequences as\n");
      printf("   the sequences generated by the resulting background model may be too similar to\n");
      printf("   the input sequences\n\n"); 
    }
    printf("\nestimating background Markov models using input sequences...\n");
    generate_background(numSeq,seq,rseq,seqLen,back,userMarkovOrder); 
    printf("done\n\n");
  }
  if (widthWt<20 & verbose)
  {
    printf("\n***Note: the window width of sequence centered on the nucleotides having large weights\n");
    printf("   in EM for PWM optimization is small: %d\n",widthWt);
    printf("   Motif longer than %d will not be discovered\n\n",widthWt); 
  }
  
  time(&start);
    
  if(verbose)
  {
    printf("==============================================================================================\n");
    printf("input sequence file:  %s\n",mFileName);
    printf("number of sequences and average length:\t\t\t\t%d %5.1f\n",numSeq,aveSeqLen);
    
    if (pgf) 
    {
      printf("\nUse pgf method to approximate llr null distribution\n");
      printf("background Markov order:\t\t\t\t\t0th\n");
      printf("parameters estimated from sequences in:  %s\n\n",mFileName);
    }
    else 
    {
      printf("\nUse an empirical approach to approximate llr null using background sequences\n");
      printf("Background sequences are simulated using [a,c,g,t] frequencies in input data\n");
      printf("Background Markov order (in llr calculation):\t\t\t\t");
      switch (userMarkovOrder) {
        case 0: printf("0th\n"); break; 
        case 1: printf("1st\n"); break; 
        case 2: printf("2nd\n"); break; 
        case 3: printf("3rd\n"); break; 
        case 4: printf("4th\n"); break; 
        case 5: printf("5th\n"); break; 
        case 6: printf("6th\n"); break; 
        case 7: printf("7th\n"); break; 
        case 8: printf("8th\n"); break; 
        default: break;
      }
      if (bFileName[0]!='\0') printf("parameters from:\t\t\t\t\t\t%s\n\n",bFileName);
      else                    printf("parameters estimated from sequences in:\t\t%s\n\n",mFileName);
    }
    if (weightType!=0) 
      printf("non-uniform weight applies to each sequence - type:\t\t%d\n",weightType);
    printf("number of GA generations & population size:\t\t\t%d %d\n\n",numGeneration,populationSize);
    printf("PWM score p-value cutoff for binding site declaration:\t\t%e\n",pvalueCutoff);
    printf("ln(E-value) cutoff for motif declaration:\t\t\t%f\n\n",E_valueCutoff);
//    printf("number (percentage) of sequences selected for EM:\t\t%d(%4.1f\%)\n",numSeqEM,100.0*(double)numSeqEM/(double)numSeq);
    printf("number of EM steps:\t\t\t\t\t\t%d\n",numEM);
    printf("minimal no. sites considered for a motif:\t\t\t%d\n\n",minsites);
    printf("[a,c,g,t] frequencies in input data:\t\t\t\t%f %f %f %f\n",bfreq[0],bfreq[1],bfreq[2],bfreq[3]);
    printf("==============================================================================================\n");
  }
  
  if (pgf) 
  {
    if (userMarkovOrder!=0 & verbose) 
    {
      printf("\n***The user-specified background Markov order (%d) is ignored when -pgf is set to 1***\n",userMarkovOrder);
    }
    if (bFileName[0]!='\0' & verbose)
    {
      printf("\n***The user-specified background models: %s are not used when -pgf is set to 1***\n\n",bFileName);
    }
  }
  if (startPWMfound && fEM!=1.0  & verbose)
  {
    printf("\n***Note: -fEM argument is ignored in a seeded analysis***\n\n");
  }
  
//  if(verbose)
//  {
//    printf("\nsSarting GADEM... this may take a few hours to complete\n");
//    printf("type: ctrl/z, then, bg, to run it in background\n\n");
//  }
  
    // determine seq length by counting only [a,c,g,t], seqLen is used in E-value calculation 
  effect_seq_length(seq,numSeq,seqLen,Iseq,emSeqLen);
    // determine the distribution and critical cut point
  pwmDistCutoff=vector_similarity();
  
  if      (weightType==0) assign_weight_uniform(seqLen,numSeq,posWeight);
  else if (weightType==1) assign_weight_rectangle(seqLen,numSeq,posWeight,widthWt);
  else if (weightType==2) assign_weight_triangular(seqLen,numSeq,posWeight);
  else if (weightType==3) assign_weight_triangular_uniform(seqLen,numSeq,posWeight,widthWt);
  else if (weightType==4) assign_weight_normal(seqLen,numSeq,posWeight);
  else
  {
    printf("Motif prior probability type not found - please choose: 0, 1, 2, 3, or 4\n");
    printf("Consider: -posWt 1 for strong central enrichment as in ChIP-seq\n");
    printf("          -posWt 0 for others\n\n");
    exit(0);
  }
  /*	if (startPWMfound) minminSites=minsites;
   else               minminSites=(int)(0.40*minsitesEM);*/
  
  motifCn=0; noMotifFound=0; numCycle=0; numCycleNoMotif=0; 
  int compt=0;
  int lengthList=GET_LENGTH(RListPWM);
  
  do
  {
    if(!startPWMfound)
    {
      
      if(verbose)
      {
        printf("\n|------------------------------------------------------------------|\n");
        printf("|                                                                  |\n");
        printf("|              *** Running an unseeded analysis ***                |\n");
        printf("|                                                                  |\n");
        printf("|------------------------------------------------------------------|\n\n");
      }
      populationSize=INTEGER_VALUE(RpopulationSize);
      numGeneration=INTEGER_VALUE(RnumGeneration);
      dyad  =alloc_chrs(populationSize,4);
      wheel =alloc_wheel(populationSize);
      fitness=alloc_fitness(populationSize);
      maxpFactor=alloc_double(populationSize);
      uniqMotif=alloc_char(populationSize+1);
      opwm  =alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      epwm=alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      pwmConsensus=alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
      pwm   =alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      pwmLen=alloc_int(populationSize);
      sdyad =alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
      word  =alloc_word(numWordGroup,maxWordSize);
      minminSites=(int)(0.40*minsitesEM);
      
        //printf("Population size %d \n",populationSize);
        //printf("####################################\n");
        //printf("####################################\n");
      
      
        // identify top-ranked k-mers (k=3,4,5) for spaced dyads
      if(verbose)
        printf("\nGADEM cycle %2d: enumerate and count k-mers...   ",numCycle+1);      
      numWordGroup=word_for_dyad(word,seq,rseq,numSeq,seqLen,bfreq,&numTop3mer,&numTop4mer,&numTop5mer);
      if(verbose)
        printf("done.\n");
      
        // generating a "population" of spaced dyads
      if(verbose)
        printf("\ninitializing GA...   ");

      initialisation(dyad,populationSize,numWordGroup,word,minSpaceWidth,maxSpaceWidth,maxpFactor);
      if(verbose)
        printf("done.\n\n");        
      
    }
    else
    {
      if(verbose)
      {
        printf("\n|------------------------------------------------------------------|\n");
        printf("|                                                                  |\n");
        printf("|               *** Running a seeded analysis ***                  |\n");
        printf("|                                                                  |\n");
        printf("|------------------------------------------------------------------|\n\n");
      }
      populationSize=FIXED_POPULATION; 
      dyad  =alloc_chrs(populationSize,4);
      pwm=alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      pwmLen=alloc_int(populationSize);
      maxpFactor=alloc_double(populationSize);
      uniqMotif=alloc_char(populationSize+1);
      opwm  =alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      epwm=alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      pwmConsensus=alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
      sdyad =alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
      word  =alloc_word(numWordGroup,maxWordSize);
      wheel =alloc_wheel(populationSize);
      fitness=alloc_fitness(populationSize);
      minminSites=minsites;
      int lengthMatrix;
      
      lengthMatrix=GET_LENGTH(VECTOR_ELT(RListPWM,compt));
      RSpwm=allocMatrix(REALSXP,4,(lengthMatrix/4));
      RSpwm=VECTOR_ELT(RListPWM,compt);
      
      
      pwmLen[0]=read_pwm0(RSpwm,pwm[0],lengthMatrix);
      
      for(i=1; i<populationSize; i++)
      {
        for (j=0; j<pwmLen[0]; j++)
        {
          for (k=0; k<4; k++)
          {
            pwm[i][j][k]=pwm[0][j][k];
          }
        }
        pwmLen[i]=pwmLen[0];
      }
      for (i=0; i<populationSize; i++)
      {
        maxpFactor[i]=FIXED_MAXPF*(i+1);
        standardize_pwm(pwm[i],pwmLen[i]);
        consensus_pwm(pwm[i],pwmLen[i],pwmConsensus[i]);
        strcpy(sdyad[i],pwmConsensus[i]);
      }
    }
    generationNoMotif=0;
    
    for (jjj=0; jjj<numGeneration; jjj++)
    {
        // convert spaced dyads to letter probability matrix
      if (!startPWMfound)
      {
        dyad_to_pwm(word,populationSize,dyad,pwm,pwmLen);
      }
      
#if defined(CSTACK_DEFNS) & defined(HAVE_DISPATCH)
      R_CheckUserInterrupt();
      dispatch_apply(populationSize, queue, ^(size_t iiii) {populationCalculation(maxSeqLen, numEM, fitness+iiii, startPWMfound, minminSites, maxpFactor[iiii], numSeq, numSeqEM, seq, rseq, seqLen, Iseq, bfreq, posWeight, weightType, pvalueCutoff, emSeqLen, pwm[iiii], pwmLen[iiii], epwm[iiii], opwm[iiii], pwmConsensus[iiii], scoreCutoff+iiii, sdyad[iiii], iiii);
      });
#else
      for (ii=0; ii<populationSize; ii++)
      {
        populationCalculationSequential(maxSeqLen, numEM, fitness, startPWMfound, minminSites, maxpFactor, numSeq, numSeqEM, seq, rseq, seqLen, Iseq, bfreq, posWeight, weightType, llrDist, llrDim, pvalueCutoff, emSeqLen, pwm, pwmLen, epwm, opwm,pwmConsensus,scoreCutoff, sdyad, t1pwm, t2pwm, siteEM, logpwm, score, rscore, ipwm, ii, verbose, jjj, numCycle, maxSpaceWidth);

      // if(llrDist)
      // {
      //   free(llrDist);
      //   llrDist=NULL;
      // }  
      // if(siteEM)
      // { 
      //   free(siteEM);
      //   siteEM=NULL;
      // }
      // if (logpwm[0])
      // { free(logpwm[0]);
      //   logpwm[0]=NULL;
      // }
      // if (logpwm)
      // { 
      //   free(logpwm);
      //   logpwm=NULL;
      // }      
      // if(ipwm[0])
      // { 
      //   free(ipwm[0]);
      //   ipwm[0]=NULL;
      // }
      // if (ipwm)
      // {
      //   free(ipwm);
      //   ipwm=NULL;
      // }
      // if (score[0])
      // {
      //   free(score[0]); 
      //   score[0]=NULL;
      // }
      // if (score)
      // { 
      //   free(score);
      //   score=NULL;
      // }
//      if (rscore[0])
//      { 
//        free(rscore[0]);
//        rscore[0]=NULL;
//      }
//      if (rscore)
//      { 
//        free(rscore);
//        rscore=NULL;
//      }

      // llrDist=alloc_distr(MAX_DIMENSION);
      // ipwm  =alloc_int_int(MAX_PWM_LENGTH,4);
      // logpwm=alloc_double_double(MAX_PWM_LENGTH,4);
      // score =alloc_double_double(numSeq,maxSeqLen);
//      rscore=alloc_double_double(numSeq,maxSeqLen);
      // siteEM=alloc_site(MAX_SITES);
      }

#endif
      if (populationSize>1)
      {
        sort_fitness(fitness,populationSize);
      }


      numUniq=check_pwm_uniqueness_dist(opwm,pwmLen,populationSize,fitness,pwmDistCutoff,E_valueCutoff,uniqMotif,slideWinPWM);
      if(verbose)
      {
        printf("\nGADEM cycle[%3d] generation[%3d] number of unique motif: %d\n",numCycle+1,jjj+1,numUniq);
        for (i=0; i<populationSize; i++)
        {
          if (uniqMotif[i]=='1')
          {
            printf("   spacedDyad: %s ",sdyad[fitness[i].index]);
            for (j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) printf(" ");
            printf("motifConsensus: %s ",pwmConsensus[fitness[i].index]);
            for (j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) printf(" ");
            printf(" %3.2f fitness: %7.2f\n",maxpFactor[fitness[i].index],fitness[i].value);
          }
        }
        printf("\n"); 
      }
      
      
      if (jjj<numGeneration-1)
      {
          // fitness based selection with replacement 
        roulett_wheel_fitness(fitness,populationSize,wheel);
          // mutation and crossover operations
        if (populationSize>1) 
        {
          testrand=runif(0,1);
          if (testrand>=0.5) 
          {
            mutation (dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif,
                      maxpFactor,maxpMutationRate); 
          }
          else 
          {
            crossover(dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif, maxpFactor,maxpMutationRate); 
          }
        }
        else
        {
          mutation (dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif, maxpFactor,maxpMutationRate);
        }
      }
    }
    
    if((numCycle+1)< lengthList)
    {
      compt++;
    }
    else
    {
      startPWMfound=0;
    }
    numCycle++;
    
    
    site=alloc_site_site(numUniq+1,MAX_SITES);
    nsites=alloc_int(numUniq+1);
    pwmnewLen=alloc_int(numUniq+1); // after base extension and trimming	
    seqCn=alloc_int(MAX_NUM_SEQ);
    bseqCn=alloc_int(MAX_NUM_SEQ);
    
      // final step user-specified background model is used
    motifCn2=0; // motifCn per GADEM cycle
    for (ii=0; ii<populationSize; ii++) 
    {
      
      id=fitness[ii].index;
      if (uniqMotif[ii]=='0')
      {
        continue;
      }
      
      MarkovOrder=min(pwmLen[id]-1,userMarkovOrder);
      
        // approximate the exact llr distribution using Staden's method
      if(pgf)
      {
        if(verbose)
        {
          printf("\nApproximate the exact pwm llr score distribution using the pgf method.\n");
        }
        log_ratio_to_int(epwm[id],ipwm,pwmLen[id],bfreq);
        
          // compute score distribution of the (int)PWM using Staden's method 
        llrDim=pwm_score_dist(ipwm,pwmLen[id],llrDist,bfreq);
        
          //printf("Avant ScoreCutoff %d \n",scoreCutoff[id]);
        scoreCutoff[id]=determine_cutoff(llrDist,llrDim,pvalueCutoff);
          //printf("Apres ScoreCutoff %d \n",scoreCutoff[id]);
        
        if(fullScan)
        {
          nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,oseq,orseq,seqLen,ipwm,pwmLen[id],scoreCutoff[id],bfreq);
        }
        else
        {
          nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,seq,rseq,seqLen,ipwm,pwmLen[id],scoreCutoff[id],bfreq);
        }
      }// determine the null llr distribution using background sequences
      else
      {
        log_pwm(epwm[id],logepwm,pwmLen[id]);
        /* -----------------compute the null distribtion------------------------------------*/
          // this generates N*(L-w+1)*numBackgSets w-mers compared to N*(L-w+1) in input data
        if(verbose)
        {
          printf("\nUse an empirical approach to approximate the llr score distribution\n");
        }
        numTopWmerInB=0; empDim=0;
        for (i=0; i<numBackgSets; i++)
        {
          simulate_background_seq(bfreq,numSeq,seqLen,pseq);
          ll_score_backg_model(numSeq, pseq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
          numWmerInB=llr_score(pscore,numSeq,pseq,seqLen,logepwm,pwmLen[id],bfreq,bscore) ;
          sort_double(pscore,numWmerInB);        // truncated null dist.
          
          for (j=0; j<numWmerInB/4; j++)
          {       // take only top25%
            empDist[numTopWmerInB]=pscore[j];   // pool all top25% from each permutation
            numTopWmerInB++;
          }
          empDim +=numWmerInB;
        }
        sort_double(empDist,numTopWmerInB);      // truncated null dist.
        /* -----------------end the null distribtion------------------------------------*/
        
        llrCutoff=empDist[(int)(pvalueCutoff*empDim)];
          // print_null(empDist,numTopWmerInB,empDim);
        if (fullScan)
        {
          ll_score_backg_model(numSeq, oseq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
          ll_score_backg_model(numSeq,orseq,rbscore,seqLen,pwmLen[id],back,MarkovOrder);
          nsites[motifCn2]=scan_llr_empirical(site[motifCn2],numSeq,oseq,orseq,seqLen,logepwm,pwmLen[id],bfreq,bscore,rbscore,llrCutoff,empDist,numTopWmerInB,empDim);
        }
        else
        {
          ll_score_backg_model(numSeq, seq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
          ll_score_backg_model(numSeq,rseq,rbscore,seqLen,pwmLen[id],back,MarkovOrder);
          nsites[motifCn2]=scan_llr_empirical(site[motifCn2],numSeq,seq,rseq,seqLen,logepwm,pwmLen[id],bfreq,bscore,rbscore,llrCutoff,empDist,numTopWmerInB,empDim);
        }
      }
      if (nsites[motifCn2]>=max(2,minsites))
      {
        for (j=0; j<numSeq; j++) seqCn[j]=0;
        for (j=0; j<nsites[motifCn2]; j++) seqCn[site[motifCn2][j].seq]++;
        
        for (j=0; j<4; j++) cn[j]=0;
        for (j=0; j<numSeq; j++)
        {
          if (seqCn[j]==0) cn[0]++;
          if (seqCn[j]==1) cn[1]++;
          if (seqCn[j]==2) cn[2]++;
          if (seqCn[j]>2)  cn[3]++;
        }
        totalSitesInput=nsites[motifCn2];
        if (extTrim)
        {
          if (fullScan)
          {
            extend_alignment(site[motifCn2],numSeq,oseq,orseq,seqLen,nsites[motifCn2],pwmLen[id],&(pwmnewLen[motifCn2]));
          }
          else
          {
            extend_alignment(site[motifCn2],numSeq,seq,rseq,seqLen,nsites[motifCn2],pwmLen[id],&(pwmnewLen[motifCn2]));
          }
        }
        else
        { 
          pwmnewLen[motifCn2]=pwmLen[id];
        } 
        
        if (fullScan)
        {
          align_sites_count(site[motifCn2],oseq,orseq,nsites[motifCn2],pwmnewLen[motifCn2],opwm2);
        }
        else
        {
          align_sites_count(site[motifCn2],seq,rseq,nsites[motifCn2],pwmnewLen[motifCn2],opwm2);
        }
        standardize_pwm(opwm2,pwmnewLen[motifCn2]);
        logev=E_value(opwm2,nsites[motifCn2],bfreq,pwmnewLen[motifCn2],numSeq,seqLen);
        
        if (logev<=E_valueCutoff)
        {
          consensus_pwm(opwm2,pwmnewLen[motifCn2],pwmConsensus[id]);
          if (fullScan)
          {
            SET_VECTOR_ELT(ResultsGadem,increment,print_result_2(site[motifCn2],nsites[motifCn2],numSeq,oseq,orseq,seqLen,logev,opwm2,pwmnewLen[motifCn2],motifCn+1,sdyad[id],pwmConsensus[id],numCycle,pvalueCutoff,maxpFactor[id]));
            increment++;           
            print_motif(site[motifCn2],nsites[motifCn2],oseq,orseq,seqLen,pwmnewLen[motifCn2],motifCn+1,opwm2);
          }
          else
          {
            SET_VECTOR_ELT(ResultsGadem,increment,print_result_2(site[motifCn2],nsites[motifCn2],numSeq,seq,rseq,seqLen,logev,opwm2,pwmnewLen[motifCn2],
                                                                 motifCn+1,sdyad[id],pwmConsensus[id],numCycle,pvalueCutoff,maxpFactor[id]));
            increment++;
            print_motif(site[motifCn2],nsites[motifCn2],seq,rseq,seqLen,pwmnewLen[motifCn2],motifCn+1,opwm2);
          }
          
          mask_sites(nsites[motifCn2],seq,rseq,seqLen,site[motifCn2],pwmnewLen[motifCn2]);
          
          /* ----------------------compute the average number of sites in background sequences ----------------------*/
          avebnsites=0; avebnsiteSeq=0;
          for (i=0; i<numBackgSets; i++)
          {
            simulate_background_seq(bfreq,numSeq,seqLen,pseq);
            reverse_seq(pseq,rpseq,numSeq,seqLen);
            
            if (pgf)
            {
              nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,pseq,rpseq,seqLen,ipwm,pwmLen[id],scoreCutoff[id],bfreq);
            }
            else
            {
              ll_score_backg_model(numSeq, pseq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
              ll_score_backg_model(numSeq,rpseq,rbscore,seqLen,pwmLen[id],back,MarkovOrder);
              nsites[motifCn2]=scan_llr_empirical(site[motifCn2],numSeq,pseq,rpseq,seqLen,logepwm,pwmLen[id],bfreq,bscore,rbscore,llrCutoff,empDist,numTopWmerInB,empDim);
            }
            
            for (j=0; j<numSeq; j++) bseqCn[j]=0;
            for (j=0; j<nsites[motifCn2]; j++) bseqCn[site[motifCn2][j].seq]++;
            
            for (j=0; j<4; j++) bcn[j]=0;
            for (j=0; j<numSeq; j++)
            {
              if (bseqCn[j]==0) bcn[0]++;
              if (bseqCn[j]==1) bcn[1]++;
              if (bseqCn[j]==2) bcn[2]++;
              if (bseqCn[j]>2)  bcn[3]++;
            }
            avebnsites+=nsites[motifCn2]; avebnsiteSeq+=(numSeq-bcn[0]);
          } 
          avebnsites/=numBackgSets; avebnsiteSeq/=numBackgSets;
          /* -----------------end compute the average number of sites in background sequences ----------------------*/
          motifCn++; motifCn2++; numCycleNoMotif=0;
        }
      }
    }
    
    for (i=0; i<motifCn2; i++)
    {
      mask_sites(nsites[i],seq,rseq,seqLen,site[i],pwmnewLen[i]); 
    }
    
    if (site[0])
    { 
      free(site[0]);
      site[0]=NULL;
    }
    if (site)      
    {
      free(site);      
      site=NULL;    
    }
    if (nsites)
    {
      free(nsites);
      nsites=NULL;
    }
    if (pwmnewLen) 
    {
      free(pwmnewLen);
      pwmnewLen=NULL;
    }
    
    if (motifCn2==0)
      numCycleNoMotif++;   
    if (numCycleNoMotif==stopCriterion)
      noMotifFound=1;
  }while (!noMotifFound);
  
   /*if (!startPWMfound) {  
   if (dyad[0])      { free(dyad[0]);         dyad[0]=NULL;    }
   if (dyad)         { free(dyad);            dyad=NULL;       }
   }*/
  if (seqLen)
  { 
    free(seqLen);
    seqLen=NULL;
  }
  if (pwm[0][0])       
  {
    free(pwm[0][0]);
    pwm[0][0]=NULL; 
  }
  if (pwm[0])
  { 
    free(pwm[0]);
    pwm[0]=NULL;     
  }
  if (pwm)             
  {
    free(pwm); 
    pwm=NULL;        
  }
  if (opwm2[0])  
  { 
    free(opwm2[0]); 
    opwm2[0]=NULL;
  }
  if (opwm2)     
  {
    free(opwm2); 
    opwm2=NULL;
  }
  if (opwm[0][0])      
  { 
    free(opwm[0][0]);
    opwm[0][0]=NULL;
  }
  if (opwm[0])    
  {
    free(opwm[0]);
    opwm[0]=NULL;
  }
  if (opwm)       
  {
    free(opwm);
    opwm=NULL;
  }
  if (t1pwm[0])        
  {
    free(t1pwm[0]);
    t1pwm[0]=NULL;
  }
  if (t1pwm)
  {
    free(t1pwm);
    t1pwm=NULL;
  }
  if (t2pwm[0])        
  { free(t2pwm[0]);       
    t2pwm[0]=NULL;  
  }
  if (t2pwm)  
  { free(t2pwm);    
    t2pwm=NULL;     
  }
  if (logpwm[0])
  { free(logpwm[0]);    
    logpwm[0]=NULL; 
  }
  if (logpwm) 
  { 
    free(logpwm); 
    logpwm=NULL;   
  }
  if(ipwm[0])
  { 
    free(ipwm[0]);     
    ipwm[0]=NULL;  
  }
  if (ipwm)
  {
    free(ipwm);   
    ipwm=NULL;
  }
  if (pwmLen)   
  { 
    free(pwmLen);    
    pwmLen=NULL; 
  }
  if (seq[0])          { free(seq[0]);          seq[0]=NULL;     }
  if (seq)             { free(seq);             seq=NULL;        }
    //  if (rseq[0])         { free(rseq[0]);         rseq[0]=NULL;    }
    // if (rseq)            { free(rseq);            rseq=NULL;       }
    // if (oseq[0])         { free(oseq[0]);         oseq[0]=NULL;    }
    // if (oseq)            { free(oseq);            oseq=NULL;       }
    // if (orseq[0])        { free(orseq[0]);        orseq[0]=NULL;   }
    // if (orseq)           { free(orseq);           orseq=NULL;      }
  if (score[0])
  {
    free(score[0]); 
    score[0]=NULL;
  }
  if (score)    
  { 
    free(score);   
    score=NULL; 
  }
  if (rscore[0])
  { 
    free(rscore[0]);   
    rscore[0]=NULL;
  }
  if (rscore)       
  { 
    free(rscore);    
    rscore=NULL;  
  }
  if (bfreq)    
  { 
    free(bfreq);    
    bfreq=NULL;  
  }
  if (wheel)    
  { 
    free(wheel);    
    wheel=NULL;    
  }
  if (fitness)    
  { 
    free(fitness); 
    fitness=NULL;
  }
  if (mFileName)  
  { 
    free(mFileName);    
    mFileName=NULL; 
  }
  if (oFileName)    
  { 
    free(oFileName);  
    oFileName=NULL;
  }
  if (pwmFileName)    
  {
    free(pwmFileName);
    pwmFileName=NULL;
  }
  if (sdyad[0]) 
  { 
    free(sdyad[0]); 
    sdyad[0]=NULL;
  }
  if (sdyad)    
  {
    free(sdyad);
    sdyad=NULL;
  }
  if (siteEM)   
  { 
    free(siteEM);
    siteEM=NULL;
  }
  if (pwmConsensus[0])
  { 
    free(pwmConsensus[0]);
    pwmConsensus[0]=NULL;
  }
  if (pwmConsensus)   
  {
    free(pwmConsensus);
    pwmConsensus=NULL;
  }
    //if (!startPWMfound && word) destroy_word(word,numWordGroup);
  
  PutRNGstate(); 
  UNPROTECT(1);
  return(ResultsGadem);
}

void print_ptable(Pgfs *llrDist,int llrDim) {
  
 
}

void print_empirical(double *empDist,int totalKmer) {
  
  
}

void select_high_scoring_seq_for_EM (double *ChIPScore,int numSeq,int numSeqEM,char *Iseq,double fEM) {
  
  register int i;
  int numSeqWithQualityScore,numSeqEMtmp1,numSeqEMtmp2;
  double *tmpScore;
  double ChIPscoreCutoff;
  
  tmpScore=alloc_double(numSeq);
  
  numSeqWithQualityScore=0;
  for (i=0; i<numSeq; i++)
  {
    if (ChIPScore[i]>0) numSeqWithQualityScore++;
  }
  
  tmpScore=alloc_double(numSeq);
  for (i=0; i<numSeq; i++) tmpScore[i]=ChIPScore[i];
  sort_double(tmpScore,numSeq);
  
  ChIPscoreCutoff=tmpScore[(int)(fEM*numSeq)];
  
  if (numSeqWithQualityScore<=(int)(fEM*numSeq))
  {
    for (i=0; i<numSeq; i++) Iseq[i]='0';
    numSeqEMtmp1=0;
    for (i=0; i<numSeq; i++)
    {
      if (ChIPScore[i]>0)
      {
        Iseq[i]='1'; numSeqEMtmp1++;
      }
    }
    numSeqEMtmp2=0;
    for (i=0; i<numSeq; i++)
    {
      if (ChIPScore[i]<=0)
      {
        Iseq[i]='1'; numSeqEMtmp2++;
        if (numSeqEMtmp1+numSeqEMtmp2==numSeqEM) break;
      }
    }
  }
  else
  {
    for (i=0; i<numSeq; i++) Iseq[i]='0';
    numSeqEMtmp1=0; numSeqEMtmp2=0;
    for (i=0; i<numSeq; i++) {
      if (ChIPScore[i]>=ChIPscoreCutoff)
      {
        Iseq[i]='1'; numSeqEMtmp1++;
        if (numSeqEMtmp1==numSeqEM) break;
      }
    }
  }
  if (tmpScore) 
  { 
    free(tmpScore);  
    tmpScore=NULL;  
  }
  if (ChIPScore) 
  { 
    free(ChIPScore);
    ChIPScore=NULL; 
  }
  
}

void print_null(double *empDist,int numTopWmerInB,int empDim)
{
 
}


void populationCalculation(int maxSeqLen, int numEM, Fitness *fitness, int startPWMfound, int minminSites, double maxpFactor, int numSeq, int numSeqEM, char **seq, char **rseq, int *seqLen, char *Iseq, double *bfreq, double ** posWeight, int weightType, double pvalueCutoff, int* emSeqLen, double **pwm, int pwmLen, double **epwm, double **opwm,char *pwmConsensus, int *scoreCutoff,char *sdyad, int ii)
{
  int jj=0,llrDim;
  double **t1pwm=alloc_double_double(MAX_PWM_LENGTH,4);
  double **t2pwm=alloc_double_double(MAX_PWM_LENGTH,4);

    //These need to be allocated when running in parallel
  double **logpwm, ** score, ** rscore;
  int **ipwm;
  Sites *siteEM;
  Pgfs *llrDist;

  llrDist=alloc_distr(MAX_DIMENSION);
  ipwm  =alloc_int_int(MAX_PWM_LENGTH,4);
  logpwm=alloc_double_double(MAX_PWM_LENGTH,4);
  score =alloc_double_double(numSeq,maxSeqLen);
  rscore=alloc_double_double(numSeq,maxSeqLen);
  siteEM=alloc_site(MAX_SITES);


  int maxp=0, nsitesEM=0;
  double pwmDiff=0;

  if (!startPWMfound)
  {
      // This function modifles sdyad[ii]
    pwm_profile(pwm, pwmLen, sdyad);
  }
    // Make a copy and then subject the copy to EM
  copy_pwm(pwm, t1pwm, pwmLen);

    // standarize pwm
  standardize_pwm(t1pwm, pwmLen);
    // EM on randomly selected sequences
  maxp=(int)(maxpFactor*numSeqEM);

    //Check if the user wants to interrupt

  for(jj=0; jj<numEM; jj++)
  {
      // Modifies logpwm
    log_pwm(t1pwm,logpwm, pwmLen);
      // Compute ll score of each w-mer | motif model
      // Modifies score and rscore
    ll_score_motif_model(numSeq,seq,rseq,seqLen,logpwm,pwmLen,score,rscore,Iseq,bfreq);

      // compute p(zij|y=1) probability of binding sites started at position j on seq i
      // modifies score and rscore
    normalize(score,rscore,seqLen,pwmLen,numSeq,Iseq,maxp,posWeight,weightType);

      // E-step
      // Modifies t2pwm
    construct_pwm(t2pwm,score,rscore,seq,rseq,seqLen,numSeq,pwmLen,Iseq);

      // M-step
      // Modifies t2pwm
    standardize_pwm(t2pwm,pwmLen);
      // Only compute pwmDiff
    pwmDiff=check_convergence(t1pwm,t2pwm,pwmLen);
      // copy t2pwm to t1pwm
    copy_pwm(t2pwm,t1pwm,pwmLen);
    if (pwmDiff<=PWM_CONVERGENCE)  break;
  }

    // Copy t1pwm to epwm[ii]
  copy_pwm(t1pwm,epwm,pwmLen);
    // Modifies ipwm
  log_ratio_to_int(epwm,ipwm,pwmLen,bfreq);
    // Compute score distribution of the (int)PWM using Staden's method 
    // Modifies llrDist
  llrDim=pwm_score_dist(ipwm,pwmLen,llrDist,bfreq);
    // Compute the score cutoff
    // Does not modify anything
  *scoreCutoff=determine_cutoff(llrDist,llrDim,pvalueCutoff);

    // test each w-mer to see if a motif site - test statistic: ll, distribution: Staden method, cutoff: user-specified
    // modifies siteEM    
  nsitesEM=scan_em_seq_ptable(llrDist,llrDim,siteEM,numSeq,seq,rseq,seqLen,ipwm,pwmLen,*scoreCutoff,bfreq,Iseq);
    // loose threshould at this step, as em only on a subset of sequences
  if (nsitesEM>=max(2,minminSites))
  { 
      // construct pwm from the identified sites
      // Modifies opwm
    align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen,opwm);
      // Modifies opwm
    standardize_pwm(opwm,pwmLen);
      // Modifies pwmConsensus
    consensus_pwm(opwm,pwmLen,pwmConsensus);
      // compute E-value of the relative entroy score of each motif, use it as fitness      
      // Only compute fitness.value
    (*fitness).value=E_value(opwm,nsitesEM,bfreq,pwmLen,numSeqEM,emSeqLen);
  }
  else 
  {
      // if too few sites in a motif
      // Modifies opwm
    align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen,opwm);
      // Modifies opwm
    standardize_pwm(opwm,pwmLen);
      // Modifies pwmConsensus
    consensus_pwm(opwm,pwmLen,pwmConsensus);
      // Only compute fitness.value
    (*fitness).value=DUMMY_FITNESS;
  }
  (*fitness).index=ii;


  if (llrDist)   
  { 
    free(llrDist);
    llrDist=NULL;
  }  
  if (siteEM)   
  { 
    free(siteEM);
    siteEM=NULL;
  }  
  if (t1pwm[0])        
  {
    free(t1pwm[0]);
    t1pwm[0]=NULL;
  }
  if (t1pwm)
  {
    free(t1pwm);
    t1pwm=NULL;
  }
  if (t2pwm[0])        
  { free(t2pwm[0]);       
    t2pwm[0]=NULL;  
  }
  if (t2pwm)  
  { free(t2pwm);    
    t2pwm=NULL;     
  }
  if (logpwm[0])
  { free(logpwm[0]);    
    logpwm[0]=NULL; 
  }
  if (logpwm) 
  { 
    free(logpwm); 
    logpwm=NULL;   
  }
  if(ipwm[0])
  { 
    free(ipwm[0]);     
    ipwm[0]=NULL;  
  }
  if (ipwm)
  {
    free(ipwm);
    ipwm=NULL;
  }
  if (score[0])
  {
    free(score[0]); 
    score[0]=NULL;
  }
  if (score)    
  { 
    free(score);   
    score=NULL; 
  }
  if (rscore[0])
  { 
    free(rscore[0]);   
    rscore[0]=NULL;
  }
  if (rscore)       
  { 
    free(rscore);    
    rscore=NULL;  
  }  
}

void populationCalculationSequential(int maxSeqLen, int numEM, Fitness *fitness, int startPWMfound, int minminSites, double *maxpFactor, int numSeq, int numSeqEM, char **seq, char **rseq, int *seqLen, char *Iseq, double *bfreq, double ** posWeight, int weightType, Pgfs *llrDist, int llrDim, double pvalueCutoff, int* emSeqLen, double ***pwm, int *pwmLen, double ***epwm, double ***opwm,char **pwmConsensus, int *scoreCutoff,char **sdyad, double **t1pwm, double **t2pwm, Sites *siteEM, double **logpwm, double** score, double** rscore, int **ipwm, int ii, int verbose, int jjj, int numCycle, int maxSpaceWidth)
{
  int maxp=0, nsitesEM=0;
  double pwmDiff=0;
  int jj=0,j=0;
    // to see from which spaced dyad a motif is derived
    // for (int kkk=0;kkk<10;kkk++)//{		
    //printf("Longeur %d \n",pwmLen[20]);//}

  if (!startPWMfound)
  {
    pwm_profile(pwm[ii],pwmLen[ii],sdyad[ii]);
  }
    // make a copy and then subject the copy to EM
  copy_pwm(pwm[ii],t1pwm,pwmLen[ii]);

    // standarize pwm
  standardize_pwm(t1pwm,pwmLen[ii]);
    // EM on randomly selected sequences
  maxp=(int)(maxpFactor[ii]*numSeqEM);

    //Check if the user wants to interrupt
  R_CheckUserInterrupt();

  for (jj=0; jj<numEM; jj++)
  {

    log_pwm(t1pwm,logpwm,pwmLen[ii]);

      // compute ll score of each w-mer | motif model
    ll_score_motif_model(numSeq,seq,rseq,seqLen,logpwm,pwmLen[ii],score,rscore,Iseq,bfreq);      
      // compute p(zij|y=1) probability of binding sites started at position j on seq i
    normalize(score,rscore,seqLen,pwmLen[ii],numSeq,Iseq,maxp,posWeight,weightType);

      // E-step
    construct_pwm(t2pwm,score,rscore,seq,rseq,seqLen,numSeq,pwmLen[ii],Iseq);

      // M-step
    standardize_pwm(t2pwm,pwmLen[ii]);

    pwmDiff=check_convergence(t1pwm,t2pwm,pwmLen[ii]);
      // copy t2pwm to t1pwm
    copy_pwm(t2pwm,t1pwm,pwmLen[ii]); 
    if (pwmDiff<=PWM_CONVERGENCE)  break;
  }

  copy_pwm(t1pwm,epwm[ii],pwmLen[ii]); // from to


  log_ratio_to_int(epwm[ii],ipwm,pwmLen[ii],bfreq);
    // compute score distribution of the (int)PWM using Staden's method 
  llrDim=pwm_score_dist(ipwm,pwmLen[ii],llrDist,bfreq);
  scoreCutoff[ii]=determine_cutoff(llrDist,llrDim,pvalueCutoff);

    // test each w-mer to see if a motif site - test statistic: ll, distribution: Staden method, cutoff: user-specified
  nsitesEM=scan_em_seq_ptable(llrDist,llrDim,siteEM,numSeq,seq,rseq,seqLen,ipwm,pwmLen[ii],scoreCutoff[ii],bfreq,Iseq);

    // loose threshould at this step, as em only on a subset of sequences
  if (nsitesEM>=max(2,minminSites))
  { 
      // construct pwm from the identified sites
    align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen[ii],opwm[ii]);
    standardize_pwm(opwm[ii],pwmLen[ii]);
    consensus_pwm(opwm[ii],pwmLen[ii],pwmConsensus[ii]);
      // compute E-value of the relative entroy score of each motif, use it as fitness

    fitness[ii].value=E_value(opwm[ii],nsitesEM,bfreq,pwmLen[ii],numSeqEM,emSeqLen);
  }
  else 
  {
      // if too few sites in a motif
    align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen[ii],opwm[ii]);
    standardize_pwm(opwm[ii],pwmLen[ii]);
    consensus_pwm(opwm[ii],pwmLen[ii],pwmConsensus[ii]);
    fitness[ii].value=DUMMY_FITNESS;
  }
  fitness[ii].index=ii;
//  if(verbose) 
//  { 
//    printf("cyc.[%3d] gen.[%3d] pop.[%3d] spacedDyad: %s ",numCycle+1,jjj+1,ii+1,sdyad[ii]);
//    for (j=strlen(sdyad[ii]); j<maxSpaceWidth+10; j++) printf(" ");
//    printf(" motifConsensus: %s",pwmConsensus[ii]);
//    for (j=strlen(sdyad[ii]); j<maxSpaceWidth+10; j++) printf(" ");
//    printf(" maxpf: %3.2f fitness: %7.2f nsitesEM: %d\n",maxpFactor[ii],fitness[ii].value,nsitesEM);
//  }
}
