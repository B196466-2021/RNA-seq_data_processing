cd #Make sure we are in our home space 
mkdir f#Make a directory in my home directory called f
cd f#Change directory into f
cp -r /localdisk/data/BPSM/AY21/fastq/ ./f#Copy files 
cd ~/f/f #Change directory 
ls *gz |xargs fastqc -t 4#List all files with the suffix gz and use fastqc to control quality by 4 threads

cp /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz #Copy files 
gzip -d *.gz #Unzip file 
hisat2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta TriTrypDB-46_TcongolenseIL3000_2019_Genome #Build index

#mapping and translate format into bam
for i in *_1.fq
do
i=${i/_1.fq}
hisat2 -p 4 -x TriTrypDB-46_TcongolenseIL3000_2019_Genome -1 ${i}_1.fq -2${i}_2.fq | samtools sort -@ 4 -o ${i}_sort.bam
done

#Build index
ls *_sort.bam | while read i 
do
i=${i/_sort.bam}
samtools index -@ 4 ${i}_sort.bam 
done

cp /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed TriTrypDB-46_TcongolenseIL3000_2019.bed #Copy files
bedtools multicov -bams *.bam -bed TriTrypDB-46_TcongolenseIL3000_2019.bed > data.txt #Using bedtools to generate counts data and save to data.txt

cp /localdisk/data/BPSM/AY21/fastq/100k.fqfiles 100k.fqfiles #Copy files
cat 100k.fqfiles #View 100k.fqfiles
cat data.txt | less -S #View data.txt

awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6+$9+$12)/3"\t"($7+$10+$13)/3"\t"($8+$11+$14)/3"\t"($15+$17+$19)/3"\t"($16+$18+$20)/3"\t"($21+$24+$27)/3"\t"($22+$25+$28)/3"\t"($23+$26+$29)/3"\t"($30+$32+$34)/3"\t"($31+$33+$35)/3"\t"($36+$39+$42)/3"\t"($37+$40+$43)/3"\t"($38+$41+$44)/3"\t"($45+$47+$48)/3"\t"($46+$48+$50)/3}' data.txt >average.txt #Give the statistical mean (average) of the counts per gene 
cat average.txt | less -S #View average.txt

#Fold change and Single variable is time
awk -F "\t" '{if($8==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7/$8}' average.txt>C1_Induced_24_48.txt
sort -t $'\t' -k 6 -n -r C1_Induced_24_48.txt -o C1_Induced_24_48.txt
awk -F "\t" '{if($10==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9/$10}' average.txt>C1_Uninduced_24_48.txt
sort -t $'\t' -k 6 -n -r C1_Uninduced_24_48.txt -o C1_Uninduced_24_48.txt
awk -F "\t" '{if($9==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6/$9}' average.txt>C1_Uninduced_0_24.txt
sort -t $'\t' -k 6 -n -r C1_Uninduced_0_24.txt -o C1_Uninduced_0_24.txt
awk -F "\t" '{if($10==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6/$10}' average.txt>C1_Uninduced_0_48.txt
sort -t $'\t' -k 6 -n -r C1_Uninduced_0_48.txt -o C1_Uninduced_0_48.txt
awk -F "\t" '{if($13==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12/$13}' average.txt>C2_Induced_24_48.txt
sort -t $'\t' -k 6 -n -r C2_Induced_24_48.txt -o C2_Induced_24_48.txt
awk -F "\t" '{if($15==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$14/$15}' average.txt>C2_Uninduced_24_48.txt
sort -t $'\t' -k 6 -n -r C2_Uninduced_24_48.txt -o C2_Uninduced_24_48.txt
awk -F "\t" '{if($14==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11/$14}' average.txt>C2_Uninduced_0_24.txt
sort -t $'\t' -k 6 -n -r C2_Uninduced_0_24.txt -o C2_Uninduced_0_24.txt
awk -F "\t" '{if($15==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11/$15}' average.txt>C2_Uninduced_0_48.txt
sort -t $'\t' -k 6 -n -r C2_Uninduced_0_48.txt -o C2_Uninduced_0_48.txt
awk -F "\t" '{if($18==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$17/$18}' average.txt>WT_Induced_24_48.txt
sort -t $'\t' -k 6 -n -r WT_Induced_24_48.txt -o WT_Induced_24_48.txt
awk -F "\t" '{if($20==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$19/$20}' average.txt>WT_Uninduced_24_48.txt
sort -t $'\t' -k 6 -n -r WT_Uninduced_24_48.txt -o WT_Uninduced_24_48.txt
awk -F "\t" '{if($19==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$16/$19}' average.txt>WT_Uninduced_0_24.txt
sort -t $'\t' -k 6 -n -r WT_Uninduced_0_24.txt -o WT_Uninduced_0_24.txt
awk -F "\t" '{if($20==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$16/$20}' average.txt>WT_Uninduced_0_48.txt
sort -t $'\t' -k 6 -n -r WT_Uninduced_0_48.txt -o WT_Uninduced_0_48.txt

#Fold change and Single variable is Induction
awk -F "\t" '{if($9==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7/$9}' average.txt>C1_Induced_Uninduced_24.txt
sort -t $'\t' -k 6 -n -r C1_Induced_Uninduced_24.txt -o C1_Induced_Uninduced_24.txt
awk -F "\t" '{if($10==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8/$10}' average.txt>C1_Induced_Uninduced_48.txt
sort -t $'\t' -k 6 -n -r C1_Induced_Uninduced_48.txt -o C1_Induced_Uninduced_48.txt
awk -F "\t" '{if($14==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12/$14}' average.txt>C2_Induced_Uninduced_24.txt
sort -t $'\t' -k 6 -n -r C2_Induced_Uninduced_24.txt -o C2_Induced_Uninduced_24.txt
awk -F "\t" '{if($15==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$13/$15}' average.txt>C2_Induced_Uninduced_48.txt
sort -t $'\t' -k 6 -n -r C2_Induced_Uninduced_48.txt -o C2_Induced_Uninduced_48.txt
awk -F "\t" '{if($19==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$17/$19}' average.txt>WT_Induced_Uninduced_24.txt
sort -t $'\t' -k 6 -n -r WT_Induced_Uninduced_24.txt -o WT_Induced_Uninduced_24.txt
awk -F "\t" '{if($20==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$18/$20}' average.txt>WT_Induced_Uninduced_48.txt
sort -t $'\t' -k 6 -n -r WT_Induced_Uninduced_48.txt -o WT_Induced_Uninduced_48.txt

#Fold change and Single variable is clone
awk -F "\t" '{if($16==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6/$16}' average.txt>C1_WT_Uninduced_0.txt
sort -t $'\t' -k 6 -n -r C1_WT_Uninduced_0.txt -o C1_WT_Uninduced_0.txt
awk -F "\t" '{if($17==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7/$17}' average.txt>C1_WT_Induced_24.txt
sort -t $'\t' -k 6 -n -r C1_WT_Induced_24.txt -o C1_WT_Induced_24.txt
awk -F "\t" '{if($18==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8/$18}' average.txt>C1_WT_Induced_48.txt
sort -t $'\t' -k 6 -n -r C1_WT_Induced_48.txt -o C1_WT_Induced_48.txt
awk -F "\t" '{if($19==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9/$19}' average.txt>C1_WT_UNinduced_24.txt
sort -t $'\t' -k 6 -n -r C1_WT_UNinduced_24.txt -o C1_WT_UNinduced_24.txt
awk -F "\t" '{if($20==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10/$20}' average.txt>C1_WT_UNinduced_48.txt
sort -t $'\t' -k 6 -n -r C1_WT_UNinduced_48.txt -o C1_WT_UNinduced_48.txt
awk -F "\t" '{if($16==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11/$16}' average.txt>C2_WT_Uninduced_0.txt
sort -t $'\t' -k 6 -n -r C2_WT_Uninduced_0.txt -o C2_WT_Uninduced_0.txt
awk -F "\t" '{if($17==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12/$17}' average.txt>C2_WT_Induced_24.txt
sort -t $'\t' -k 6 -n -r C2_WT_Induced_24.txt -o C2_WT_Induced_24.txt
awk -F "\t" '{if($18==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$13/$18}' average.txt>C2_WT_Induced_48.txt
sort -t $'\t' -k 6 -n -r C2_WT_Induced_48.txt -o C2_WT_Induced_48.txt
awk -F "\t" '{if($19==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$14/$19}' average.txt>C2_WT_Uninduced_24.txt
sort -t $'\t' -k 6 -n -r C2_WT_Uninduced_24.txt -o C2_WT_Uninduced_24.txt
awk -F "\t" '{if($20==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$15/$20}' average.txt>C2_WT_Uninduced_48.txt
sort -t $'\t' -k 6 -n -r C2_WT_Uninduced_48.txt -o C2_WT_Uninduced_48.txt
awk -F "\t" '{if($11==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6/$11}' average.txt>C1_C2_Uninduced_0.txt
sort -t $'\t' -k 6 -n -r C1_C2_Uninduced_0.txt -o C1_C2_Uninduced_0.txt
awk -F "\t" '{if($12==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7/$12}' average.txt>C1_C2_Induced_24.txt
sort -t $'\t' -k 6 -n -r C1_C2_Induced_24.txt -o C1_C2_Induced_24.txt
awk -F "\t" '{if($13==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8/$13}' average.txt>C1_C2_Induced_48.txt
sort -t $'\t' -k 6 -n -r C1_C2_Induced_48.txt -o C1_C2_Induced_48.txt
awk -F "\t" '{if($14==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9/$14}' average.txt>C1_C2_Uninduced_24.txt
sort -t $'\t' -k 6 -n -r C1_C2_Uninduced_24.txt -o C1_C2_Uninduced_24.txt
awk -F "\t" '{if($15==0) print "NA"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10/$15}' average.txt>C1_C2_Uninduced_48.txt
sort -t $'\t' -k 6 -n -r C1_C2_Uninduced_48.txt -o C1_C2_Uninduced_48.txt
 
