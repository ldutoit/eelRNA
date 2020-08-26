  
#This short script find and output sequences to a fasta file. It can take a  transcript or a prefix to get all the isoforms in .e. TRINITY_DN701_c2_g1_i2 or TRINITY_DN701_c2). It can take "nucl" or "petide" to look in the MRNA or the peptide  records.

##Start by reading the files once ( a bit long so we might as well do it only once
library("seqinr") # can be installed using install.packages("seqinr")
nucl<-read.fasta(file="Trinity.fasta")
pep<-read.fasta(file="Trinity.fasta.transdecoder.pep")

#####Parameters
  geneofinterest <- "TRINITY_DN2157_c0"  ## Enter a ranscript or a prefix to get all the isoforms plotted (i.e. TRINITY_DN701_c2_g1_i2 or TRINITY_DN701_c2)
  type<-"DNA" # enter "DNA" or "peptide". It can print the matching peptide from the file  or the original DNA from the file ...

  
  
### Main
  if (type == "DNA"){
    sequences<-nucl[grep(geneofinterest,names(nucl))]
    names_of_sequences<-names(nucl)[grep(geneofinterest,names(nucl))]
     write.fasta(sequences = sequences , names = names_of_sequences, file.out=paste(geneofinterest,"_mRNA.fasta",sep=""))
     print(paste("Outputted",length(sequences) , "sequences to file",paste(geneofinterest,"_mRNA.fasta",sep="")))
  }
  if (type == "peptide"){
    sequences<-pep[grep(geneofinterest,names(pep))]
    names_of_sequences<-names(pep)[grep(geneofinterest,names(pep))]
    write.fasta(sequences = sequences , names = names_of_sequences, file.out=paste(geneofinterest,"_pep.fasta",sep=""))
    print(paste("Outputted",length(sequences) , "sequences to file",paste(geneofinterest,"_pep.fasta",sep="")))
    }
