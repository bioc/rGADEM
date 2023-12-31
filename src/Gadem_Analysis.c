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
#include <stdbool.h>
#include "defines.h"
#include "evalue_meme.h"
#include <Rdefines.h>
#include <Rversion.h>

#include "config.h"

/*---------------------------------------------------------------------------
// v1.3.1: last modifications 5/14/2011
// added user-input background model
// added masking simple repetitive elements and list them in .mask
// added a new argument -nmotifs: the maximal number of motifs sought
// removed high order background model (to be added back later)
// edited usage and info.txt printout
// simplified usage by elimating several options or setting them as default)
// ---------------------------------------------------------------------------*/

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


void populationCalculation(int maxSeqLen, int numEM,
                           Fitness *fitness,
                           int startPWMfound,
                           int minminSites,
                           double maxpFactor,
                           int numSeq, int numSeqEM,
                           char **seq, char **rseq,
                           int *seqLen, char *Iseq,
                           double *bfreq,
                           double ** posWeight, int weightType,
                           double pvalueCutoff, int* emSeqLen,
                           double **pwm, int pwmLen,
                           double **epwm, double **opwm,
                           char *pwmConsensus, int *scoreCutoff,
                           char *sdyad,
                           int ii);

char* convertRString2Char(SEXP);

SEXP GADEM_Analysis(SEXP sequence,SEXP sizeSeq, SEXP accession, SEXP Rverbose,SEXP RnumWordGroup,SEXP RnumTop3mer,SEXP RnumTop4mer,SEXP RnumTop5mer,SEXP RnumGeneration,SEXP RpopulationSize, SEXP RpValue,SEXP ReValue,SEXP RextTrim,SEXP RminSpaceWidth,SEXP RmaxSpaceWidth,SEXP RuseChIPscore,SEXP RnumEM,SEXP RfEM, SEXP RwidthWt,SEXP RfullScan, SEXP RslideWinPWM,SEXP RstopCriterion,SEXP RnumBackgSets,SEXP RweightType,SEXP RbFileName,SEXP RListPWM,SEXP RminSites,SEXP RmaskR,SEXP Rnmotifs) 
{
  char *bFileName;
  
  SEXP ResultsGadem;
  SEXP RSpwm;
  PROTECT(ResultsGadem=NEW_LIST(100));  
  
  int increment=0;
  
  double testrand;
  
  //Number of sequences
  int numSeq = INTEGER_VALUE(sizeSeq);
  // const
//  char *Fastaheader[size];
  int incr=0;
  
  int longueur=length(sequence);
  int IncrementTemp=0;
  
  // basic settings/info
  int maxSeqLen,*seqLen;       		 // sequence info	
  double aveSeqLen;                      // sequence info
  char **seq,**rseq;
  int *geneID;         			 // sequence info
  char **oseq,**orseq;                   // copy of the original sequences
  char **sseq,**rsseq;                   // simulated seqs.
  double *bfreq1, *bfreq0=NULL;                // base frequencies
  double *ChIPScore;                     // chip score
  int maskR;				 // mask simple repeats before running the algorithm
    
  // pwms
  double ***pwm;                         // initial population of PWMs from spaced dyads
  int *pwmLen;                           // initial pwm lengths
  double **opwm2;                        // EM-derived PWM
  double ***opwm;                        // observed PWMs from identified sites
  double ***epwm;                        // em-optimized PWMs
  double **logepwm;                      // log(em-optimized PWM)
  int *pwmnewLen;                        // final motif length after extending to both ends
  
  // llr score distr.
  Pgfs *llrDist;                         // llr distribution from pgf method
  int llrDim;                            // llr distribution dimension
  int **ipwm;                            // integer pwm for computing llr score distribution
  
  // EM, motif, sites
  double pvalueCutoff;                   // user input, used to determine score cutoff based on ipwm
  int *scoreCutoff;                      // pwm score cutoff for the corresponding p-value cutoff
  double logev;                          // log of E-value of a motif;
  int useChIPscore;                      // indicator for using ChIP-seq score for seq. selection for EM
  int numEM;                             // number of EM steps
  double E_valueCutoff;                  // log E-value cutoff
  //int nsitesEM;                          // number of binding sites in sequences subjected to EM
  int minsitesEM;                        // minimal number of sites in a motif in EM sequences
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
  int numBackgSets;
  
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
  int motifCn2,id,numCycleNoMotif,verbose,minminSites,nmotifs;
  int startPWMfound,stopCriterion;
  char *mFileName,*oFileName,*pwmFileName,*tempRbFileName;
  time_t start;
  int cn[4],bcn[4],*seqCn,*bseqCn,avebnsites,avebnsiteSeq,totalSitesInput;
  int i; 
  int ii=0;
  int jjj=0;
  
  GetRNGstate();
  
  mFileName=alloc_char(500);         mFileName[0]='\0';
  oFileName=alloc_char(500);         oFileName[0]='\0';
  pwmFileName=alloc_char(500);       pwmFileName[0]='\0';
  bFileName=alloc_char(500);         bFileName[0]='\0';
  //tempRbFileName=alloc_char(500);    tempRbFileName[0]='\0';
  seq=NULL; aveSeqLen=0; maxSeqLen=0; 
  //minsites=-1; 
  
  startPWMfound=0;    

  maxSeqLen=0;
  for(incr=1;incr<longueur;incr=incr+2)
  { 
    if (length(STRING_ELT(sequence,(incr)))>maxSeqLen) maxSeqLen=length(STRING_ELT(sequence,(incr))); 
  }
//  Rprintf("maxLength=%d",maxSeqLen);
//  exit(0);
  seq=alloc_char_char(numSeq,maxSeqLen+1);
  for(incr=1;incr<longueur;incr=incr+2)
  { 
    for (int j=0; j<length(STRING_ELT(sequence,(incr))); j++)
    {
      seq[IncrementTemp][j]=CHAR(STRING_ELT(sequence,(incr)))[j];
    }
    IncrementTemp++;
  }
  
  
  verbose=LOGICAL_VALUE(Rverbose);
  numWordGroup=INTEGER_VALUE(RnumWordGroup);
  minsites=INTEGER_VALUE(RminSites);
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
  slideWinPWM=INTEGER_VALUE(RslideWinPWM);
  numUniq=populationSize;
  stopCriterion=INTEGER_VALUE(RstopCriterion);  
  numBackgSets=INTEGER_VALUE(RnumBackgSets);
  weightType=NUMERIC_VALUE(RweightType);
  //const char *tempRbFileName[1];

 	tempRbFileName = convertRString2Char(RbFileName);	

  //tempRbFileName[0]=CHAR(STRING_ELT(RbFileName,0));
  nmotifs = INTEGER_VALUE(Rnmotifs);
  maskR = INTEGER_VALUE(RmaskR);

  

  if(numSeq>MAX_NUM_SEQ)
  {
    error("Error: maximal number of seqences reached!\nPlease reset MAX_NUM_SEQ in gadem.h and rebuild (see installation)\n");
  }
  
  strcpy(bFileName,tempRbFileName);

  ChIPScore=alloc_double(MAX_NUM_SEQ);
  seqLen=alloc_int(MAX_NUM_SEQ); 
  geneID=alloc_int(MAX_NUM_SEQ);

//  seq=sequences;
  
//  numSeq=size;
  int len; 
  
  for (i=0; i<numSeq; i++)
  {
    len=strlen(seq[i]); 
    seqLen[i]=len;
    geneID[i]=INTEGER(accession)[i];
  }

  aveSeqLen=0; 
  for (i=0; i<numSeq; i++) aveSeqLen +=seqLen[i]; aveSeqLen /=(double)numSeq;
  
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
    for (int j=0; j<seqLen[i]; j++)
    {
      oseq[i][j]=seq[i][j];
      orseq[i][j]=rseq[i][j];
    }
    oseq[i][seqLen[i]]='\0'; orseq[i][seqLen[i]]='\0'; 
  }
    
  if (strcmp(bFileName,"NULL")!= 0)
  {
    bfreq0=alloc_double(5);
    read_background(bFileName,bfreq0);
  }

  if (GET_LENGTH(RListPWM)!= 0)
  {
    startPWMfound=1; 
  }
  else { }
  
    // check for input parameters
  if(numGeneration<1)
  { 
    error("number of generaton < 1.\n");
  }
  if(populationSize<1)
  {
    error("population size < 1.\n");
  }
  if (minSpaceWidth<0)
  { 
    error("minimal number of unspecified bases in spaced dyads <0.\n"); 
  }
  if (maxSpaceWidth<0)
  { 
    error("maximal number of unspecified bases in spaced dyads <0.\n"); 
  }
  if (minSpaceWidth>maxSpaceWidth)
  {
    error("mingap setting must <= to maxgap setting.\n\n"); 
  }
  if (maxSpaceWidth+12>MAX_PWM_LENGTH)
  {
    error("maxgap setting plus word lengths exceed <MAX_PWM_LENGTH>.\n");
  }
  if (numEM<0)
  {
    error("number of EM steps is zero.\n");
  }
  if (numEM==0)
  {
    error("number of EM steps = 0, no EM optimization is carried out.\n");
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
    error("maxw3, maxw4, and maxw5 all zero - no words for spaced dyads.\n");
  }
  
  // if (startPWMfound && fEM!=0.5 && fEM!=1.0 & verbose)
  // {
  //   warning("fEM argument is ignored in a seeded analysis\n");
  // }
  
  if (startPWMfound)
  {
    // if(verbose)
    // {
    //   if (populationSize!=10 && populationSize!=100) warning("pop argument is ignored in a seeded analysis, -pop is set to 10.\n");
    //   if (numGeneration!=1 && numGeneration!=5)      warning("gen argument is ignored in a seeded analysis, -gen is set to 1.\n");
    // }
    fEM=1.0;
    populationSize=FIXED_POPULATION; numGeneration=1; 
  }
  
    // number of sequences for EM
  if (fEM>1.0 || fEM<=0.0)
  { 
    error("The fraction of sequences subject to EM is %3.2f.\n",fEM);
  } 
  numSeqEM=(int)(fEM*numSeq);
  


  // memory callocations
  Iseq  =alloc_char(numSeq+1); 
  opwm2 =alloc_double_double(MAX_PWM_LENGTH,4);
  ipwm  =alloc_int_int(MAX_PWM_LENGTH,4);
  logepwm=alloc_double_double(MAX_PWM_LENGTH,4);
  emSeqLen=alloc_int(numSeqEM);
  scoreCutoff=alloc_int(1000);
  // scoreCutoff=alloc_int(populationSize);
  llrDist=alloc_distr(MAX_DIMENSION);
  posWeight=alloc_double_double(numSeq,maxSeqLen);
  sseq=alloc_char_char(MAX_NUM_SEQ,maxSeqLen+1);
  rsseq=alloc_char_char(MAX_NUM_SEQ,maxSeqLen+1);
  bfreq1=base_frequency(numSeq,seq,seqLen);

  if (strcmp(bFileName,"NULL") == 0)
  {
    bfreq0=alloc_double(5);
    for (i=0; i<4; i++)
      {
	bfreq0[i]=bfreq1[i];
      }
  }
  

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
  
  if (maskR==1) mask_repetitive(geneID,seq,numSeq,seqLen,mFileName);

  if (widthWt<20)
  {
    warning("The window width of sequence centered on the nucleotides having large weights in EM for PWM optimization is small\n Motif longer than %d will not be discovered\n",widthWt);
  }
  
  time(&start);
  
    // if (weightType==1 || weightType==3) 
    //Rprintf("window width of sequence centered on the nucleotides having large weights for PWM optimization: %d\n",widthWt);
    //Rprintf("pwm score p-value cutoff for declaring binding site:\t%e\n",pvalueCutoff);
  
  if(verbose)
  {
    Rprintf("==============================================================================================\n");
    Rprintf("input sequence file:  %s\n",mFileName);
    Rprintf("number of sequences and average length:\t\t\t\t%d %5.1f\n",numSeq,aveSeqLen);
    
    Rprintf("Use pgf method to approximate llr null distribution\n");
    Rprintf("parameters estimated from sequences in:  %s\n\n",mFileName);

    if (weightType!=0) 
      Rprintf("non-uniform weight applies to each sequence - type:\t\t%d\n",weightType);
    Rprintf("number of GA generations & population size:\t\t\t%d %d\n\n",numGeneration,populationSize);
    Rprintf("PWM score p-value cutoff for binding site declaration:\t\t%e\n",pvalueCutoff);
    Rprintf("ln(E-value) cutoff for motif declaration:\t\t\t%f\n\n",E_valueCutoff);
//    Rprintf("number (percentage) of sequences selected for EM:\t\t%d(%4.1f\%)\n",numSeqEM,100.0*(double)numSeqEM/(double)numSeq);
    Rprintf("number of EM steps:\t\t\t\t\t\t%d\n",numEM);
    Rprintf("minimal no. sites considered for a motif:\t\t\t%d\n\n",minsites);
    Rprintf("[a,c,g,t] frequencies in input data:\t\t\t\t%f %f %f %f\n",bfreq1[0],bfreq1[1],bfreq1[2],bfreq1[3]);
    Rprintf("==============================================================================================\n");
  }
  
  // if (pgf) 
  // {
  //   if (userMarkovOrder!=0 & verbose) 
  //   {
  //     warning("The user-specified background Markov order (%d) is ignored when -pgf is set to 1\n",userMarkovOrder);
  //   }
  //   if (bFileName[0]!='\0' & verbose)
  //   {
  //     warning("The user-specified background models: %s are not used when -pgf is set to 1\n",bFileName);
  //   }
  // }
  // if (startPWMfound && fEM!=1.0  & verbose)
  // {
  //   warning("fEM argument is ignored in a seeded analysis\n");
  // }
  
    // determine seq length by counting only [a,c,g,t], seqLen is used in E-value calculation
  effect_seq_length(seq,numSeq,seqLen,Iseq,emSeqLen);
    // determine the distribution and critical cut point
  pwmDistCutoff=vector_similarity();
  
  if      (weightType==0) assign_weight_uniform(seqLen,numSeq,posWeight);
  else if (weightType==1) assign_weight_triangular(seqLen,numSeq,posWeight);
  else if (weightType==2) assign_weight_normal(seqLen,numSeq,posWeight);
  else
  {
    error("Motif prior probability type not found - please choose: 0, 1, or 2\n");
    // Rprintf("Consider: -posWt 1 for strong central enrichment as in ChIP-seq\n");
    // Rprintf("          -posWt 0 for others\n\n");
    // exit(0);
  }
  /*    if (startPWMfound) minminSites=minsites;
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
        Rprintf("*** Running an unseeded analysis ***\n");
        // Rprintf("\n|------------------------------------------------------------------|\n");
        // Rprintf("|                                                                  |\n");
        // Rprintf("|              *** Running an unseeded analysis ***                |\n");
        // Rprintf("|                                                                  |\n");
        // Rprintf("|------------------------------------------------------------------|\n\n");
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

        // identify top-ranked k-mers (k=3,4,5) for spaced dyads
      if(verbose)
        Rprintf("GADEM cycle %2d: enumerate and count k-mers... ",numCycle+1);
        
      numWordGroup=word_for_dyad(word,seq,rseq,numSeq,seqLen,bfreq1,&numTop3mer,&numTop4mer,&numTop5mer);
      
      if(verbose)
        Rprintf("Done.\n");
      
        // generating a "population" of spaced dyads
      if(verbose)
        Rprintf("Initializing GA... ");

      initialisation(dyad,populationSize,numWordGroup,word,minSpaceWidth,maxSpaceWidth,maxpFactor);
      if(verbose)
        Rprintf("Done.\n");
      
    }
    else
    {
      if(verbose)
      {
        Rprintf("*** Running an seeded analysis ***\n");
        // Rprintf("\n|------------------------------------------------------------------|\n");
        // Rprintf("|                                                                  |\n");
        // Rprintf("|               *** Running a seeded analysis ***                  |\n");
        // Rprintf("|                                                                  |\n");
        // Rprintf("|------------------------------------------------------------------|\n\n");
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
        for (int j=0; j<pwmLen[0]; j++)
        {
          for (int k=0; k<4; k++)
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

      DO_APPLY(populationCalculation(maxSeqLen, numEM, fitness+ii, 
                                     startPWMfound, minminSites, maxpFactor[ii], 
                                     numSeq, numSeqEM, seq, rseq, seqLen, Iseq, 
                                     bfreq0, posWeight, weightType, 
                                     pvalueCutoff, emSeqLen, 
                                     pwm[ii], pwmLen[ii], epwm[ii], opwm[ii], 
                                     pwmConsensus[ii], scoreCutoff+ii, sdyad[ii], ii),
               populationSize, ii);




      // for (i=0; i<5; i++)
      // {
      //   Rprintf("fitness.value=%lf\n",fitness[i].value);
      //   Rprintf("fitness.index=%d\n",fitness[i].index);
      //   Rprintf("maxpfactor=%lf\n",maxpFactor[i]);
      //   Rprintf("scoreCutoff=%d\n",scoreCutoff[i]);
      //   Rprintf("   spacedDyad: %s\n",sdyad[i]);
      //   
      //   for (l=0; l<pwmLen[i]; l++)
      //   {
      //     for (m=0; m<4; m++) 
      //     { 
      //       Rprintf("opwm[%d][%d][%d]=%lf ",i,l,m,opwm[i][l][m]);
      //       Rprintf("epwm[%d][%d][%d]=%lf ",i,l,m,epwm[i][l][m]);
      //       Rprintf("pwm[%d][%d][%d]=%lf ",i,l,m,pwm[i][l][m]);
      //     }
      //     Rprintf("\n");
      //   }
      //   Rprintf("\n");
      // }
      // 
      // testrand=runif(0,1);
      // Rprintf("testrand1=%lf\n",testrand);
      
      if (populationSize>1)
      {
        sort_fitness(fitness,populationSize);
      }



      // for (i=0; i<5; i++)
      // {
      //   Rprintf("fitness.value=%lf\n",fitness[i].value);
      //   Rprintf("fitness.index=%d\n",fitness[i].index);
      // }
      numUniq=check_pwm_uniqueness_dist(opwm, pwmLen,
                                        populationSize, fitness,
                                        pwmDistCutoff, E_valueCutoff,
                                        uniqMotif, slideWinPWM);


      // for (i=0; i<5; i++)
      // {
      //   Rprintf("fitness.value=%lf\n",fitness[i].value);
      //   Rprintf("fitness.index=%d\n",fitness[i].index);
      //   Rprintf("maxpfactor=%lf\n",maxpFactor[i]);
      //   Rprintf("scoreCutoff=%d\n",scoreCutoff[i]);
      //   Rprintf("   spacedDyad: %s\n",sdyad[i]);
      //   
      //   for (l=0; l<pwmLen[i]; l++)
      //   {
      //     for (m=0; m<4; m++) 
      //     { 
      //       Rprintf("opwm[%d][%d][%d]=%lf",i,l,m,opwm[i][l][m]); 
      //     }
      //     Rprintf("\n");
      //   }
      //   Rprintf("\n");
      // }
      
      if(verbose)
      {
        Rprintf("GADEM cycle[%3d] generation[%3d] number of unique motif: %d\n",numCycle+1,jjj+1,numUniq);
        for (i=0; i<populationSize; i++)
        {
          if (uniqMotif[i]=='1')
          {
            Rprintf("   spacedDyad: %s ",sdyad[fitness[i].index]);
            for (int j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) Rprintf(" ");
            Rprintf("motifConsensus: %s ",pwmConsensus[fitness[i].index]);
            for (int j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) Rprintf(" ");
            Rprintf(" %3.2f fitness: %7.2f\n",maxpFactor[fitness[i].index],fitness[i].value);
          }
        }
        Rprintf("\n");
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
            mutation(dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif,
                      maxpFactor,maxpMutationRate); 
          }
          else
          {
            crossover(dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif, maxpFactor,maxpMutationRate); 
          }
        }
        else
        {
          mutation(dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif, maxpFactor,maxpMutationRate);
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
      if(uniqMotif[ii]=='0')
      {
        continue;
      }


      // approximate the exact llr distribution using Staden's method
      // if(verbose)
      // {
      //   Rprintf("Approximate the exact pwm llr score distribution using the pgf method.\n");
      // }
      log_ratio_to_int(epwm[id],ipwm,pwmLen[id],bfreq0);

        // compute score distribution of the (int)PWM using Staden's method
      llrDim=pwm_score_dist(ipwm,pwmLen[id],llrDist,bfreq0);

        //Rprintf("Avant ScoreCutoff %d \n",scoreCutoff[id]);
      scoreCutoff[id]=determine_cutoff(llrDist,llrDim,pvalueCutoff);
        //Rprintf("Apres ScoreCutoff %d \n",scoreCutoff[id]);
        
      if(fullScan)
      {
        nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,oseq,orseq,seqLen,ipwm,pwmLen[id],scoreCutoff[id],bfreq0);
      }
      else
      {
        nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,seq,rseq,seqLen,ipwm,pwmLen[id],scoreCutoff[id],bfreq0);
      }
      if (nsites[motifCn2]>=max(2,minsites))
      {
      for (int j=0; j<numSeq; j++) seqCn[j]=0;
        for (int j=0; j<nsites[motifCn2]; j++) seqCn[site[motifCn2][j].seq]++;
        
        for (int j=0; j<4; j++) cn[j]=0;
        for (int j=0; j<numSeq; j++)
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
        logev=E_value(opwm2,nsites[motifCn2],bfreq0,pwmnewLen[motifCn2],numSeq,seqLen);

        if (logev<=E_valueCutoff)
        {
          consensus_pwm(opwm2,pwmnewLen[motifCn2],pwmConsensus[id]);
          if (fullScan)
          {
            SET_VECTOR_ELT(ResultsGadem,increment,print_result_R(site[motifCn2],nsites[motifCn2],numSeq,oseq,orseq,seqLen,logev,opwm2,pwmnewLen[motifCn2],motifCn+1,sdyad[id],pwmConsensus[id],numCycle,pvalueCutoff,maxpFactor[id],geneID));
            increment++;           
            print_motif(site[motifCn2],nsites[motifCn2],oseq,orseq,seqLen,pwmnewLen[motifCn2],motifCn+1,opwm2);
          }
          else
          {
            SET_VECTOR_ELT(ResultsGadem,increment,print_result_R(site[motifCn2],nsites[motifCn2],numSeq,seq,rseq,seqLen,logev,opwm2,pwmnewLen[motifCn2],
                                                                 motifCn+1,sdyad[id],pwmConsensus[id],numCycle,pvalueCutoff,maxpFactor[id],geneID));
            increment++;
            print_motif(site[motifCn2],nsites[motifCn2],seq,rseq,seqLen,pwmnewLen[motifCn2],motifCn+1,opwm2);
          }

          mask_sites(nsites[motifCn2],seq,rseq,seqLen,site[motifCn2],pwmnewLen[motifCn2]);

          /* ----------------------compute the average number of sites in background sequences ----------------------*/
          avebnsites=0; avebnsiteSeq=0;
          for (i=0; i<numBackgSets; i++)
          {
            simulate_background_seq(bfreq0,numSeq,seqLen,sseq);
            reverse_seq(sseq,rsseq,numSeq,seqLen);

            nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,sseq,rsseq,seqLen,ipwm,pwmLen[id],scoreCutoff[id],bfreq0);
            
            for (int j=0; j<numSeq; j++) bseqCn[j]=0;
            for (int j=0; j<nsites[motifCn2]; j++) bseqCn[site[motifCn2][j].seq]++;
            
            for (int j=0; j<4; j++) bcn[j]=0;
            for (int j=0; j<numSeq; j++)
            {
              if (bseqCn[j]==0) bcn[0]++;
              if (bseqCn[j]==1) bcn[1]++;
              if (bseqCn[j]==2) bcn[2]++;
              if (bseqCn[j]>2)  bcn[3]++;
            }
              //Rprintf("background set[%2d] Seqs with 0,1,2,>2 sites: %d %d %d %d\n",i+1,bcn[0],bcn[1],bcn[2],bcn[3]);
            avebnsites+=nsites[motifCn2]; avebnsiteSeq+=(numSeq-bcn[0]);
          } 
          avebnsites/=numBackgSets; avebnsiteSeq/=numBackgSets;
          /* -----------------end compute the average number of sites in background sequences ----------------------*/
          motifCn++; motifCn2++; 

			//if((numCycle+1) > lengthList & fixSeeded)
			//	{	
			//	  numCycleNoMotif=1;
			//		startPWMfound=1;
			//		} else {
					numCycleNoMotif=0;
			//	}

        }
      }
    }
    
    /* for (int i=0; i<motifCn2; i++)
    {
      mask_sites(nsites[i],seq,rseq,seqLen,site[i],pwmnewLen[i]); 
    } */
    
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
    if (motifCn==nmotifs)
      {
	Rprintf("Maximal number of motifs (%d) reached\n",nmotifs);
	break;
      }
    if (numCycleNoMotif==stopCriterion)
      noMotifFound=1;
  }while (!noMotifFound);
  
  
    // fclose(fp);
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
  //if (seq[0])          { free(seq[0]);          seq[0]=NULL;     }
  //if (seq)             { free(seq);             seq=NULL;        }
    //  if (rseq[0])         { free(rseq[0]);         rseq[0]=NULL;    }
    // if (rseq)            { free(rseq);            rseq=NULL;       }
    // if (oseq[0])         { free(oseq[0]);         oseq[0]=NULL;    }
    // if (oseq)            { free(oseq);            oseq=NULL;       }
    // if (orseq[0])        { free(orseq[0]);        orseq[0]=NULL;   }
    // if (orseq)           { free(orseq);           orseq=NULL;      }
  if (bfreq1)    
  { 
    free(bfreq1);    
    bfreq1=NULL;  
  }
  if (bfreq0)
  {
    free(bfreq0);
    bfreq0=NULL;
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



void populationCalculation(int maxSeqLen, int numEM, 
               Fitness *fitness, 
               int startPWMfound, 
               int minminSites, 
               double maxpFactor, 
               int numSeq, int numSeqEM, 
               char **seq, char **rseq, 
               int *seqLen, char *Iseq, 
               double *bfreq, 
               double **posWeight, int weightType, 
               double pvalueCutoff, int *emSeqLen, 
               double **pwm, int pwmLen, 
               double **epwm, double **opwm,
               char *pwmConsensus, int *scoreCutoff,
               char *sdyad, 
               int ii)
{


  int jj=0,llrDim;

  // two pwms before and after EM steps
  double **t1pwm=alloc_double_double(MAX_PWM_LENGTH,4);
  double **t2pwm=alloc_double_double(MAX_PWM_LENGTH,4);

  //These need to be allocated when running in parallel
  double **logpwm;              // log-transformed EM-derived EM
  double ** score, ** rscore; // subsequence score, plus and minus strands
  int **ipwm;      // integer pwm for computing llr score distribution
  Sites *siteEM;                // binding sites in EM sequences
  Pgfs *llrDist;                // llr distribution from pgf method

  llrDist=alloc_distr(MAX_DIMENSION);
  ipwm=alloc_int_int(MAX_PWM_LENGTH,4);
  logpwm=alloc_double_double(MAX_PWM_LENGTH,4);
  score=alloc_double_double(numSeq,maxSeqLen);
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

  // test each w-mer to see if a motif site
  // test statistic: ll 
  // distribution: Staden method 
  // cutoff: user-specified
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

void read_background(char *filename,double *bfreq) {

   FILE *fp;
   char *buffer,*tok,letter[2];
   int i,len,numTab;
   double sum;

   fp=fopen(filename,"r");
   if (!fp)
      { 
         error("Incorrect filename for background model\n"); 
      } 

   buffer=alloc_char(250);

   for (i=0; i<4; i++) bfreq[i]=-1;

   while (!feof(fp)) {
      if (fgets(buffer,250,fp)>0) {
         if (buffer[0]=='#') continue;
         len=strlen(buffer);
         buffer[len-1]='\0';
         numTab=0;
         for (i=0; i<len; i++) {
            if (buffer[i]=='\0') numTab++;
         }
         if (numTab>0) {
            tok=strtok(buffer,"\t");
            if (strlen(tok)>1) continue;
            letter[0]=tok[0];
            tok=strtok(0,"\t");
            if      (letter[0]=='A' || letter[0]=='a') { 
               if (bfreq[0]==-1)  bfreq[0]=atof(tok); 
            }
            else if (letter[0]=='C' || letter[0]=='c') { 
               if (bfreq[1]==-1)  bfreq[1]=atof(tok); 
            }
            else if (letter[0]=='G' || letter[0]=='g') {
               if (bfreq[2]==-1)  bfreq[2]=atof(tok); 
            }
            else if (letter[0]=='T' || letter[0]=='t') {
               if (bfreq[3]==-1)  bfreq[3]=atof(tok); 
            }
            else  { 
               error("Error reading %s: non-[A,C,G,T]\n",filename);
               /*Rprintf("Error reading %s: non-[A,C,G,T]\n",filename); exit(0);*/
            } 
         }
         else {
            tok=strtok(buffer," ");
            letter[0]=tok[0];
            if (strlen(tok)>1) continue;
            tok=strtok(0," ");
            if      (letter[0]=='A' || letter[0]=='a') { 
               if (bfreq[0]==-1)  bfreq[0]=atof(tok); 
            }
            else if (letter[0]=='C' || letter[0]=='c') { 
               if (bfreq[1]==-1)  bfreq[1]=atof(tok); 
            }
            else if (letter[0]=='G' || letter[0]=='g') {
               if (bfreq[2]==-1)  bfreq[2]=atof(tok); 
            }
            else if (letter[0]=='T' || letter[0]=='t') {
               if (bfreq[3]==-1)  bfreq[3]=atof(tok); 
            }
            else  {
               error("Error reading %s: non-[A,C,G,T]\n",filename);
               /*Rprintf("Error reading %s: non-[A,C,G,T]\n",filename); exit(0);*/
            }
         }
      }
   }
   fclose(fp);

   for (i=0; i<4; i++) {
      if (bfreq[i]==-1) {
         switch (i) {
            case 0: Rprintf("freq. for 'a' not found in %s\n",filename); break;
            case 1: Rprintf("freq. for 'c' not found in %s\n",filename); break;
            case 2: Rprintf("freq. for 'g' not found in %s\n",filename); break;
            case 3: Rprintf("freq. for 't' not found in %s\n",filename); break;
            default: break;
         }
         error("");
         /*exit(0);*/
      }
   }
   sum=0; for (i=0; i<4; i++) sum +=bfreq[i];
   if (fabs(sum-1.0)>0.001) {
      error("Warning: frequenices do not add to 1.0\nPlease check %s\n",filename);
      /*Rprintf("Warning: frequenices do not add to 1.0\n");*/
      /*Rprintf("Please check %s\n",filename);*/
      /*exit(0);*/
   }
   if (buffer) { free(buffer); buffer=NULL; }
}

char* convertRString2Char(SEXP rstring)
{
  char* charArray;
  int rstringlength=length(STRING_ELT(rstring, 0));
  //int i=0;
  int lengthcharArray=0;

  charArray=alloc_char(rstringlength+1);

  strcpy(charArray, CHAR(STRING_ELT(rstring, 0)));

  lengthcharArray = (int)(strlen(charArray)); 
  charArray[lengthcharArray]='\0';

  return charArray;
}

