dnds
====
Sébastien Renaut
University of British Columbia, Vancouver 
sebastien.renaut@gmail.com
2013

please cite is "Renaut S, Nolte AW, Bernatchez L (2010) Mining transcriptome sequences towards identifying adaptive single nucleotide polymorphisms in lake whitefish species pairs (Coregonus spp Salmonidae). Molecular Ecology, 19, 115–131." if you are going to use this. This is were the paper for which I originally developped the code, although I have used it in several other papers since.

This scripts does five main things. But use is at your own risks!
It should work fine on Linux or Mac. Probably not on a PC.

1. It annotates all the reference sequences with "snp_table" (here the SNP table contains the two polymorphic alleles).

2. It finds all the best OpenReadingFrames (using getorf which you need to install yourself (http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html).

3. It places the annotated reference sequences in the right orientation and reading frame. It removes alternate stop codons. 

4. It prepares the PAML files for codeml (codeml is not supplied because it is not mine. You can easily get it yourself by installing PAML and copied from there in your working directory (Phylogenetic Analysis by Maximum Likelihood, Yang 2007) ).

5. It parses the output files of codeml and it creates a list of genes with their dn, ds and dn/ds values. (pnps.out) Given that these SNPs are not fixed among populations, these are actually pn (non-synonymous polymorphism per non-synonymous site), ps, pn/ps values  
