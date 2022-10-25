This project is about analysis Cas12f proteins.

a. Identify the Cas proteins.

a.1. Install the “prodigal.linux” , "bedtools" and “pilercr” in the default environment. 
prodigal.linux:https://github.com/yszhou2016/Cas13/tree/master/0.Cas-Finder
bedtools: https://sourceforge.net/projects/bedtools/ 
pilecr: http://www.drive5.com/pilercr/

a.2. run “perl 0.Cas-Finder $sample.fasta” to obtain the Cas proteins and generate "$sample.pep.cas.fasta" file.
a.3. run " perl 1.Cas12f-Finder.pl $sample.pep.cas.fasta" to otain the Cas12f proteins and generate "$sample.pep.cas.Cas12.fa" file.
a.4. Multiple alignment of Cas12f proteins with mafft: mafft --maxiterate 1000 --thread 12 --globalpair Cas12f.fa > Cas12f.mafft.fasta

b. TracrRNA prediction.
run "perl  2.Cas12f.tracrRNA.Finder.pl  Contig_name  Cas12f_start  Cas12f_end  Direction  CRISPR_array_start  CRISPR_array_end  Direct_repeat  output" to obtain the potential tracrRNA and generate "*tracrRNA.txt" file.

c. Indel analysis
Install the “bwa” in the default environment.
run "perl  3.Indel_Calculate.pl  Sample.R1.fastq   Target_sequence   Reference_sequence" and generate "*output.txt" file with Indel rate and Indel distrubution imformation.  
