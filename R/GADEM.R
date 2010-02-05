
GADEM<- function (Sequences,seed=1,genome=NULL,verbose=0,numWordGroup=3,numTop3mer=20,numTop4mer=40,numTop5mer=60,numGeneration=5,
	populationSize=100,pvalue=0.0002,EValue=0.0,extTrim=1,minSpaceWidth=0,maxSpaceWidth=10,useChIPscore=0,numEM=40,fEM=0.5,widthWt=80,
	fullScan=0,userBackgModel=0,slideWinPWM=6,stopCriterion="NUM_NO_MOTIF",MarkovOrder=0,userMarkovOrder=0,numBackgSets=10,weightType=0,pgf=1,startPWMfound=0,bOrder=-1,
	bFileName="NULL",Spwm=NULL) 
	{

		if(is(Sequences,"RangedData"))
		{
			spSeq<-space(Sequences)
			stSeq<-start(Sequences)
			edSeq<-end(Sequences)
			FastaSeq<-getSeq(genome,spSeq,start=stSeq,end=edSeq)
			FastaXstring <- XStringViews(FastaSeq,subjectClass="DNAString")
		}

		else if(is(Sequences,"XStringViews"))
		{
			FastaXstring<-Sequences
		}
		
		else if(is(Sequences,"DNAStringSet"))
		{
		FastaXstring<-Sequences	
		}

		else
		{
			stop("Object 'Sequences' should be of type 'XStringViews', 'DNAStringSet' or 'RangedData'")
		}

		FastaSequence<-DNAStringSet(FastaXstring)
		fastarecords<-XStringSetToFASTArecords(FastaSequence)
		sequenceFasta<-sapply(fastarecords,"tolower")
		accession<-as.character(1:length(FastaSequence))	

		Lengthfasta<-length(FastaSequence)

		cat("**Start Programm C **\n")

		if(!is.null(seed))
		{
			cat("*******Seed fixe--**********\n")
			set.seed(seed)	
		}
	
		# Calling C code with .Call
		obj<-.Call("GADEM_Analysis",sequenceFasta,Lengthfasta,accession,verbose,numWordGroup,numTop3mer,numTop4mer,numTop5mer,numGeneration,populationSize,
		pvalue,EValue,extTrim,minSpaceWidth,maxSpaceWidth,useChIPscore,numEM,fEM,widthWt,fullScan,userBackgModel,slideWinPWM,stopCriterion,
		MarkovOrder,userMarkovOrder,numBackgSets,weightType,pgf,startPWMfound,bOrder,bFileName,Spwm)

		i<-1
		j<-1
		
		parameter=list()
		parameter[[1]]<-new("parameters",numWordGroup=numWordGroup,numTop3mer=numTop3mer,verbose=verbose,numTop4mer=numTop4mer,numTop5mer=numTop5mer,numGeneration=numGeneration,
		populationSize=populationSize,pvalue=pvalue,EValue=EValue,extTrim=extTrim,minSpaceWidth=minSpaceWidth,maxSpaceWidth=maxSpaceWidth,useChIPscore=useChIPscore,
		numEM=numEM,fEM=fEM,widthWt=widthWt,fullScan=fullScan,userBackgModel=userBackgModel,slideWinPWM=slideWinPWM,stopCriterion=stopCriterion,MarkovOrder=MarkovOrder,
		userMarkovOrder=userMarkovOrder,numBackgSets=numBackgSets,weightType=weightType,pgf=pgf,startPWMfound=startPWMfound,bOrder=bOrder,bFileName=bFileName)
		
		list2=list()

		while(i<100&&(!is.null(obj[[i]])))
		{
			list=list()
			length(obj[[i]][[4]][[1]])

			for(j in 1:length(obj[[i]][[4]][[1]]))
			{

				if(is(Sequences,"RangedData"))
				{	
					ind<-as.numeric(obj[[i]][[4]][[5]][[j]])
					list[[j]]<-new("align",seq=obj[[i]][[4]][[1]][[j]],strand=obj[[i]][[4]][[2]][[j]],pos=obj[[i]][[4]][[3]][[j]],pval=obj[[i]][[4]][[4]][[j]],chr=spSeq[ind],start=stSeq[ind],end=edSeq[ind],seqID=obj[[i]][[4]][[6]][[j]],fastaHeader=obj[[i]][[4]][[5]][[j]])
				}

				else if(is(Sequences,"XStringViews"))
				{
					list[[j]]<-new("align",seq=obj[[i]][[4]][[1]][[j]],strand=obj[[i]][[4]][[2]][[j]],pos=obj[[i]][[4]][[3]][[j]],pval=obj[[i]][[4]][[4]][[j]],seqID=obj[[i]][[4]][[6]][[j]],fastaHeader=obj[[i]][[4]][[5]][[j]])			
				}

			}
			matrixPWM<-round(obj[[i]][[2]],4)
			rownames(matrixPWM)<-c("A","C","G","T")
			colnames(matrixPWM)<-seq(dim(matrixPWM)[2])
			list2[[i]]<-new("motif",alignList=list,consensus=obj[[i]][[1]],pwm=matrixPWM,name=obj[[i]][[5]])
			cat("i");
			i=i+1
		}
		gadem<-new("gadem",motifList=list2,parameters=parameter)
		return(gadem)
	}
 