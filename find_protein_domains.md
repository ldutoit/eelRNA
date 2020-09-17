# Find protein domains

This small script aims at finding protein domains of interest.


##Create a searchable peptide fasta file.

This aims to create Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep

First, I  only keep the genes that are considered expressed to some levels to avoid excessive number of matches from the peptide file.

I also transform the multiline seuqences in single line sequences so that I search for protein across the whole sequence regardless of line breaks.

```python
from Bio import SeqIO

###Expressed_genes

genes_to_keep = [ line.split()[1] for line in open("/Users/dutlu42p/repos/mahuika/eelRNA/denovoassembly_annotation_report_withDEPASSINGFILTER34606.xls")]

i=0
output = open("Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep","w")
with open("Trinity.fasta.transdecoder.pep", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
    	i+=1
  		if i%100 ==0:print(i)
    	transcript_name = record.id.split(".")[0]
    	print(transcript_name)
    	if transcript_name in genes_to_keep:
    		output.write(">"+record.id+"\n"+str(record.seq)+"\n")
output.close()

```
I check how many peptides we have:

```bash
grep -c  ">"  Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep
#26698 peptides out of 33359 unique transcript passing filter
``` 
## Finding patterns

```
grep -B 1 -o  -E PATTERN  Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep > exoressed_genes_match_[exp] # B show us the sequence name -o return only the matching pattern

```

The following four regexp patterns were prepared by Lucila.

```
>pattern1_EGF
[A-Z]{4}C[A-Z]{0,48}C[A-Z]{3,12}C[A-Z]{1,70}C[A-Z]{1,6}C[A-Z]{2}G[YWF][A-Z]{0,21}G[A-Z]{2}C[A-Z]


>pattern2_LDLa
C[A-Z]{1,12}C[A-Z]{1,8}C[A-Z]{1,8}C[A-Z]{1,10}C[A-Z]{1,21}C

>pattern3_LDLb
YWTD

>pattern4_internalisation
F[A-Z]NP[A-Z]Y
```

Search creating four files:


```
grep -B 1 -o -E  "[A-Z]{4}C[A-Z]{0,48}C[A-Z]{3,12}C[A-Z]{1,70}C[A-Z]{1,6}C[A-Z]{2}G[YWF][A-Z]{0,21}G[A-Z]{2}C[A-Z]" Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep > expressed_genes_matches_to_EGF.txt 

grep -B 1 -o -E "C[A-Z]{1,12}C[A-Z]{1,8}C[A-Z]{1,8}C[A-Z]{1,10}C[A-Z]{1,21}C"   Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep > expressed_genes_matches_to_LDLa.txt 

grep -B 1 -o -E  "YWTD"  Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep > expressed_genes_matches_to_LDLb.txt 

grep -B 1 -o -E  "F[A-Z]NP[A-Z]Y"  Trinity.fasta.transdecoder_withDEPASSINGFILTER.pep > expressed_genes_matches_to_internalisation.txt 

grep -c ">" *txt
#expressed_genes_matches_to_EGF.txt:150
#expressed_genes_matches_to_LDLa.txt:1467
#expressed_genes_matches_to_LDLb.txt:28
#expressed_genes_matches_to_internalisation.txt:44
```


sort file1 file2| uniq -d
