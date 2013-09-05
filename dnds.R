### Sébastien Renaut
### University of British Columbia, Vancouver 
### sebastien.renaut@gmail.com
### 2013

###please cite is "Renaut S, Nolte AW, Bernatchez L (2010) Mining transcriptome sequences towards identifying adaptive single nucleotide polymorphisms in lake whitefish species pairs (Coregonus spp Salmonidae). Molecular Ecology, 19, 115–131." if you are going to use this. This is were the paper for which I originally developped the code, although I have used it in several other papers since.

#This scripts does five main things. But use is at your own risks!
#It should work fine on Linux or Mac. Probably not on a PC.
#######1. It annotates all the reference sequences with "snp_table" (here the SNP table contains the two polymorphic alleles).
#######2. It finds all the best OpenReadingFrames (using getorf which you need to install yourself (http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html).
#######3. It places the annotated reference sequences in the right orientation and reading frame. It removes alternate stop codons. 
#######4. It prepares the PAML files for codeml (codeml is not supplied because it is not mine. You can easily get it yourself by installing PAML and copied from there in your working directory (Phylogenetic Analysis by Maximum Likelihood, Yang 2007) ).
#######5. It parses the output files of codeml and it creates a list of genes with their dn, ds and dn/ds values. (pnps.out) Given that these SNPs are not fixed among populations, these are actually pn (non-synonymous polymorphism per non-synonymous site), ps, pn/ps values  


###############################
###define some required parameters and file###
###############################
wd = "~/Documents/dnds_examples/dnds" # working directory where files are. 
ref = "reference.fa" #raw reference transcriptome
getorf = "getorf" #Were is getorf installed?
snp_table = "snp_table" #the SNP table. 
aa_code = "aa_code" #table with amino acid codes and codons correspondences.
setwd(wd)
		
#########################################
### UPLOAD AND ANNOTATE REFERENCE SEQUENCES#######
#########################################

###Replace names with unique numerical ID### 
reference_transcriptome = as.matrix(read.delim(ref, header = F,sep = "\t"))
reference_transcriptome_unique_ID = reference_transcriptome
reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1] = paste(">",c(1:length(reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1])), sep = "") 
write.table(reference_transcriptome_unique_ID,"unique_ID.fa", row.names = F, col.names = F, quote = F)

###create a objet reference_matrix which contains the annotated consensus sequences. 
reference_vector = c(1:nrow(reference_transcriptome))[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0)] #where do sequences start#
reference_vector = cbind(reference_transcriptome[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0),1], reference_vector,0)
reference_vector = rbind(reference_vector,c("null",nrow(reference_transcriptome),0)) #extra line

for(i in 1:(nrow(reference_vector)-1)) #add the sequence on one row
{reference_vector[i,3] = paste(reference_transcriptome[(as.numeric(reference_vector[i,2])+1): (as.numeric(reference_vector[(i+1),2])-1),1],collapse = "")} 

reference_matrix = cbind(gsub(">","",reference_vector[1:nrow(reference_vector),1]),reference_vector[,3],reference_vector[,3],reference_vector[,3]) #have two versions of the sequence for 2 pseudo haplotypes.
colnames(reference_matrix) = c("name","raw_sequence","ONE_haplo1","TWO_haplo1")

### annotate reference matrix with the consensus object (ie. a SNP table).
snp = read.delim(snp_table, header = T, stringsAsFactors = F)

for(i in 1:nrow(reference_matrix))
	{
	x = snp[,1] %in% reference_matrix[i,1]
	y = snp[x == T, ]

	if(nrow(y) == 1)  #only one SNP in contig
	{
	substring(reference_matrix[i,3], as.numeric(y[2]), as.numeric(y[2])) = substring(y[4],1,1);
	substring(reference_matrix[i,4], as.numeric(y[2]), as.numeric(y[2])) = substring(y[5],1,1)
	}

	if(nrow(y) > 1) for(j in 1:nrow(y)) #more than one SNP in contig
	{
	substring(reference_matrix[i,3], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,4],1,1)
	substring(reference_matrix[i,4], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,5],1,1)
	}

	if(i %% 500 == 0) print(paste(i, "of",nrow(reference_matrix),"sequences annotated, time is:",  Sys.time())) #progress report
}

#############################
### FIND BEST ORF using getorf#######
#############################
###Alternatively, if you already know were your ORF are, you can simply skip this section and supply a file which contains name of contig, start position and end position as in the object "orf_positions" created at the end of this section.

###RUN GETORF PACKAGE FROM THE EMBOSS PIPELINE (FOUND HERE: http://emboss.sourceforge.net/) ###
system(paste(getorf,"-sequence unique_ID.fa -minsize 300 -find 2 -outseq unique_ID_orf.fa"))

### PARSE THE OUTPUT OF GETORF###
out = as.matrix(read.delim("unique_ID_orf.fa", header = F, sep = " "))
out = out[,1:5]
out[,2] = gsub("^.","", out[,2]);out[,4] = gsub(".$","", out[,4]);out = gsub("SENSE)","", out, fixed = T);out = cbind(gsub("(REVERSE","", out, fixed = T), 0);out = rbind(out,">")

x = c(1:nrow(out))[(regexpr(">",out[,1],fixed = T) > 0)] #sequence starting lines
for(i in 1: (length(x)-1)) {out[x[i],5]  = paste(out[(x[i]+1): (x[(i+1)]-1),1], collapse = "")} #add full sequence to column 5

out_seq = out[(regexpr(">",out[,1],fixed = T) > 0),] #remove uneeded rows
out_seq[,3] =  abs(as.numeric(out_seq[,2]) - as.numeric(out_seq[,4]))  #sequence length
out_seq[,6] = apply(out_seq,2,substring,(regexpr(">", out_seq[,1],fixed = T)+1),(regexpr("_", out_seq[,1],fixed = T)-1))[,1] #unique ID 
unique_gene = unique(out_seq[,6])
unique_orf = NULL #These will be the longest ORF.

for(i in 1:(length(unique_gene)-1)) # this loop is to keep only the longest ORF in case there are several for one contig.
	{
	x = out_seq[,6] %in% unique_gene[i]
	temp = out_seq[x == T, ]
	if(length(temp) == 6) longest_temp = temp[as.numeric(temp[3]) == max(as.numeric(temp[3]))] else longest_temp = temp[as.numeric(temp[,3]) == max(as.numeric(temp[,3])),]
	if(length(longest_temp) == 6) unique_orf = rbind(unique_orf, longest_temp) else unique_orf = rbind(unique_orf,longest_temp[1,])
	}

### get the proper names back from the original file
	names = reference_transcriptome[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0),1] 
	for(i in 1:nrow(unique_orf)) {unique_orf[i,1] = names[as.numeric(unique_orf[i,6])]}
	unique_orf[,6] = gsub(">","",unique_orf[,1], fixed = T)
	colnames(unique_orf) = c("contig","start","length","stop","sequence","contigs")

#These are the position of the ORF in the reference contigs.
orf_positions = unique_orf[,c(6,2,4)]


##########################################
### create a new matrix which contains the ANNOTATED ORF ###
##########################################
symbols =      c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B")
replacements = c("t","g","c","a","k","y","w","s","r","m","b","d","h","v")

orf_consensus = cbind(reference_matrix[,1],0,0)#This object will contain the SNP annotated ORF
colnames(orf_consensus) = c("reference","ann","deb")

for(i in 1:nrow(reference_matrix))
{

x = unique_orf[,6] %in% reference_matrix[i,1]
x_positions = unique_orf[x == T, c(1,2,4)]

if((length(x_positions) == 3) & (as.numeric(x_positions[2]) < as.numeric(x_positions[3]))) orf_consensus[i,2:3] =  substring(reference_matrix[i,c(3,4)],as.numeric(x_positions[2]),as.numeric(x_positions[3]))
if((length(x_positions) == 3) & (as.numeric(x_positions[2]) > as.numeric(x_positions[3]))) # the orf for the reverse complement case.
	{
		orf_consensus[i,2:3] = substring(reference_matrix[i,c(3,4)],as.numeric(x_positions[3]),as.numeric(x_positions[2])) #substring the ORF
		orf_consensus[i,2] = paste(rev(strsplit(orf_consensus[i,2],"")[[1]]),collapse = "")  # ORF in the reverse order
		orf_consensus[i,3] = paste(rev(strsplit(orf_consensus[i,3],"")[[1]]),collapse = "")  # ORF in the reverse order
		for(s in 1:length(symbols)) orf_consensus[i,2:3] = gsub(symbols[s], replacements[s], orf_consensus[i,2:3]) #Complement sequence. 
	}
}
orf_consensus[,2:3] = toupper(orf_consensus[,2:3]) #gsub for small caps to capital letters. 

###This is to check for stop codons and cut the sequence until the STOP codon appears### ####STOP CODONS: TAA TRA TAG TGA#### 
premature_stop = NULL #list of premature STOP. If there are STOP, sequence is cut until the stop codon.

for(i in 1: nrow(orf_consensus))
{
	for(t in 1: (nchar(orf_consensus[i,3])/3))
		{
			for(k in 2:3)
			{
		if((substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TAA") | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TRA")  | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TAG") | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TGA")) {  premature_stop = c(premature_stop, orf_consensus[i,1]); orf_consensus[i,2:3] = substring(orf_consensus[i,2:3],1,((t-1)*3)) }
  			}
  		}	
  if(i %% 500 == 0) print(paste(i, "of",nrow(orf_consensus),"ORF checked, time is:",Sys.time())) #progress report
}

########################################
###FORMAT PAML INPUT FILE, RUN PAML, PARSE OUTPUT###
########################################
temp = c(rbind(colnames(orf_consensus)[2:3], orf_consensus[i,2:3]))

orf_length = nchar(orf_consensus[,3]) #number of character in each sequence
orf_paml = NULL # this will be the object which contains all the sequences
name_paml = NULL # this will contain all the names 

for(i in 1: nrow(orf_consensus))
	{
	if(orf_length[i] > 300) temp = c(paste(2,orf_length[i]),c(rbind(colnames(orf_consensus)[2:3], orf_consensus[i,2:3])))
	if(orf_length[i] > 300) orf_paml = c(orf_paml, temp)
	if(orf_length[i] > 300) name_paml = c(name_paml, orf_consensus[i,1])
	}

write.table(c(orf_paml,"//end", name_paml),"file_for_CODEML", row.names = F, col.names = F, quote = F)

###modify the control file for codeml
codeml.ctl = read.table("codeml.ctl", stringsAsFactors = F, header = F, sep = "\t")
codeml.ctl[1,1] = paste("ndata =",length(orf_paml[orf_paml == "ann"])) #how many sequences
codeml.ctl[2,1] = "seqfile = file_for_CODEML" #sequence file
codeml.ctl[3,1] = "treefile = ann_deb.tree" #tree file
codeml.ctl[4,1] = "outfile = ann_deb.out"#outfile
write.table(codeml.ctl,"codeml.ctl", row.names = F, col.names = F, quote = F)

###RUN PAML########
system(paste("./codeml codeml.ctl",sep = ""))

###Parse the output of PAML###
ds = as.matrix(read.delim("2NG.dS", header = F))
dn = as.matrix(read.delim("2NG.dN", header = F))
dnds = matrix(0, nrow = (length(ds)/3), ncol = 4)
colnames(dnds) = c("name","pn","ps","pnps")

for(i in 1: (length(ds)/3))
	{
	dnds[i,2] = tail(strsplit(dn[(i*3),1],split = " ")[[1]],1)
	dnds[i,3] =  tail(strsplit(ds[(i*3),1],split = " ")[[1]],1)
	}

	dnds[,4] = signif(as.numeric(dnds[,2]) / as.numeric(dnds[,3]),4)
	dnds[,1] = name_paml
	dnds[,4] = gsub("Inf","NaN", dnds[,4])

write.table(dnds,"pnps.out", row.names = F, col.names = T, quote = F) ###pnps results###

###tidy up
system("rm 2NG* 4fold.nuc rub rst1 rst lnf unique_ID_orf.fa unique_ID.fa file_for_CODEML")




