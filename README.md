# RNA-seq_data_processing

## Prerequisites

```
fastqc
bowtie2
samtools
bedtools
```

1. perform a quality check on the paired-end raw sequence data (which are in gzip compressed fastq format; https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html) using the programme fastqc.
   
2. assess the numbers and quality of the raw sequence data based on the output of fastqc.

3. align the read pairs to the Trypanosoma congolense genome using installed programme bowtie2 OR hisat2 , converting the output to indexed "bam" format with samtools.

4. generate counts data: the number of reads that align to the regions of the genome that code for genes; this is to be done using the programme bedtools and the TriTrypDB- 46_TcongolenseIL3000_2019.bed "bedfile" that contains the
information about the gene locations in the genome that was assembled and annotated in 2019. The gene names are in the 4th column of the "bedfile".for the purposes of this analysis, you should assume, incorrectly, that all genes have no introns.

5. generate plain text tab-delimited output files that give the statistical mean (average) of the counts per gene (i.e. expression levels) for each group; as the gene names are pretty uninformative to a biologist, the gene descriptions (provided in the bed file) should also be included.

6. use the mean expression levels to generate "fold change" data for the "group-wise" comparisons; these data are indicative only, as you are not being asked to do statistical modelling, but should be output such that the fold-changes are in decreasing order; as the gene names are pretty uninformative to a biologist, the gene descriptions (provided in the bed file) should also be included.
