TEnGExA<-function(x, ...)
{ UseMethod("TEnGExA")}

TEnGExA.default <- function(x, fpkm_flag=0, threshold=5, tissue_num= 5, min_fpkm= 1, ...) {
  
  finalfile="Output1.csv"
  logFile="logFile.txt"
  
  r1=nrow(x)
  c1=ncol(x)
  
  if (is.data.frame(x))
  {
  x1=data.matrix(x)
    
  if(!is.matrix(x1) & !is.numeric(x1))
  {
     stop("Data must be numeric ")
   #  write.csv("wrong input file provided ",finalfile)
     write("wrong input file provided... ", file=logFile, append=TRUE)
  }
  }
  else
  {
    stop("Data must be numeric ")
  #  write.csv("wrong input file provided ",finalfile)  
    write("wrong input file provided... ", file=logFile, append=TRUE)
  }

 
    
###################### Thershold ##########
#threshold=5
########  number of tissue need to be considered for group enrichedment class #
#tissue_num= 5

########  FPKM value minimum threshold defualt 1 #
#min_fpkm= 1

######## flag 0 for count file or 1 for FPKM values
#fpkm_flag=1


## minimum should tissue_num=3
##########################################

fpkmfile="FPKM.csv"


fpkm=x









# ****************** calculate **************
fpkm=x
tfactors=array()
maxfpkm=array()

u=2



if (fpkm_flag>1 || fpkm_flag<0)
{ stop("fpkm flag must be 0 or 1")
 # write.csv("wrong inputs,fpkm flag must be 0 or 1 ",finalfile)
  write("wrong inputs,fpkm flag must be 0 or 1... ", file=logFile, append=TRUE)
}

if (is.double(x1) && fpkm_flag==0)
{
 # stop("Data is double but you have choosen data is in counts")
 # write.csv("Data is double but you have choosen data is in counts",finalfile)
  write("Data is double but you have choosen data is in counts... ", file=logFile, append=TRUE)
}

if (is.integer(x1) && fpkm_flag==1)
{
 # stop("Data is integer (Count) but you have choosen data is in FPKM")
 # write.csv("Data is integer (Count) but you have choosen data is in FPKM",finalfile)
   write("Data is integer (Count) but you have choosen data is in FPKM",file=logFile, append=TRUE)
}

if (fpkm_flag==0)
{
  u=3
for (i in u:c1)
{tfactors[i]=sum(x[,i])/1000000}

for (i in 1:r1)
{
  for (j in u:c1)
  {
    fpkm[i,j]=(x[i,j]/tfactors[j])/(x[i,2]/1000)

    names(fpkm)[j]<-names(x[j])
  }
}
}
  
c11=c1-u+1
if (tissue_num>=c11 || tissue_num<=1)
{ stop("tissue_num must be more than 1 and less than columns(tissues)") 
 # write.csv("tissue_num must be more than 1 and less than columns(tissues)",finalfile)
  write.csv("Tissues consider must be more than 1 and less than columns(tissues)",finalfile)
}

tissue_num=tissue_num+1
t1=tissue_num 
c2=c1-u


for (i in 1:r1)
{
    
  maxfpkm=sort(fpkm[i,u:c1],decreasing = TRUE)
  sumt=sum(fpkm[i,u:c1])
  
  flageg=0
    
  egstr=""


  for (j in 1:(t1-1))
  {
    if (maxfpkm[j+1]!=0)
      fpkm[i,c1+(2*j-1)]=(maxfpkm[j])/(maxfpkm[j+1])
    else
      fpkm[i,c1+(2*j-1)]=(maxfpkm[j])/1;
    
   if (is.na(fpkm[i,c1+(2*j-1)])) fpkm[i,c1+(2*j-1)]=0
   
   
   
   if (fpkm[i,c1+(2*j-1)]>=threshold & maxfpkm[j]>=min_fpkm)
   {
     fpkm[i,c1+2*j]=colnames(maxfpkm)[j]
     if (j<=1)
     {flageg=1
      egstr<-paste("Tissue-Enriched:" ,colnames(maxfpkm)[j])
     }
     else
    {   if (flageg==0)
       egstr<-paste("Group-Enriched")
        flageg=1}
    }
   else
     fpkm[i,c1+2*j]="FALSE"
   
   names(fpkm)[c1+2*j]<-paste ("Fold",j,"/",(j+1), sep = "")
  # sumt=sumt+maxfpkm[j+1]
          
  }
  
  if ((sumt-maxfpkm[1])!=0)
   fpkm[i,c1+2*(t1-1)+1]=(maxfpkm[1])/((sumt-maxfpkm[1])/c2)
  else
    fpkm[i,c1+2*(t1-1)+1]=0
  
  fpkm[i,c1+2*(t1-1)+2]=colnames(maxfpkm)[1]
  
  if (fpkm[i,c1+2*(t1-1)+1]>=threshold  & flageg==0 & maxfpkm[1]>=min_fpkm )
   { egstr<-paste("Tissue-Enhanced:" ,colnames(maxfpkm)[1])
    flageg=1}


  	 if (flageg==0 & maxfpkm[1]>=min_fpkm  & maxfpkm[c2+1]<min_fpkm)
  	 { egstr<-paste("Mixed")
  	   flageg=1
    #	 cat ( "I=", i, "nchar=", nchar(egstr))
  	 }
 
  if (flageg==0 & maxfpkm[c2+1]>=min_fpkm)
  { egstr<-paste("Expressed in all")
    flageg=1
  }
  
  if (flageg==0)
    egstr<-paste("Not Expressed")
      
  fpkm[i,c1+2*(t1-1)+3]=egstr
  
  if (maxfpkm[2]==0 & maxfpkm[1]>=threshold)
    fpkm[i,c1+2*(t1-1)+4]="Tissue specific"
  else
    fpkm[i,c1+2*(t1-1)+4]="NA"
  names(fpkm)[c1+2*(t1-1)+4]= "Tissue Specific"
}


#write.csv(fpkm,paste("Details_",finalfile, sep = ""))
names(fpkm)[c1+2*(t1-1)+3]= "Tissue Enrichment"
names(fpkm)[c1+2*(t1-1)+4]= "Tissue Specific"


finalTable=data.frame(x = r1, y = (c1+2))
if (fpkm_flag==0)
{finalTable=fpkm[, c(1,3:c1,c1+2*(t1-1)+3,c1+2*(t1-1)+4)]}
else
{finalTable=fpkm[, c(1:c1,c1+2*(t1-1)+3,c1+2*(t1-1)+4)]} 

# write.csv(finalTable,finalfile)

# write.csv(finalTable,paste("Final_",finalfile, sep = ""))

fpkm_result<-model<-list(finalTable)
return(fpkm_result)

}
