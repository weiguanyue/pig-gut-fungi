# pig-gut-fungi

# De-redundancy of fungal genomes
dRep dereplicate drep99 -g *.fna -p 25  -d  -nc 0.3 -pa 0.9 -sa 0.99 --ignoreGenomeQuality --greedy_secondary_clustering  --run_tertiary_clustering   --skip_plots

# Quality Control
fastp -i ${fqs1} \
-o ${base2}_clean.1.fq.gz \
-I ${fqs2} \
-O ${base2}_clean.2.fq.gz \
-W 4 -M 20 -n 5 -c -l 50 -w 3 --cut_tail

# Remove host
bowtie2 -x sus_11.1 -p 5  -1 ${fqs1} -2 ${fqs2} | samtools view -@2 -bS -o ${base2}.bam
samtools view -@3  -bf 12 -F 256 ${base2}.bam |samtools sort  -@3 -n  -o ${base2}_bothEndsUnmapped_sorted.bam
wait
${bedtools} bamtofastq -i ${base2}_bothEndsUnmapped_sorted.bam -fq ${base2}_1.fastq -fq2 ${base2}_2.fastq

# Generate Reports
${fastp} -A -G -Q -L -i ${base2}_1.fastq  -I ${base2}_2.fastq  -o ${base2}_R1.fastq.gz   -O ${base2}_R2.fastq.gz   -h ${base2}_final.html  -j ${base2}_final.json

# Remove bacteria
bowtie2 -x /hl/weiguanyue/databases/pig_MAG/pig_bacteria  -p 5  -1 ${fqs1} -2 ${fqs2} |samtools view  -bS -o ${base2}.bam
# exect pairs nomapping 
samtools view   -bf 12 -F 256 ${base2}.bam >${base2}_nobac.bam
samtools sort -n ${base2}_nobac.bam -o ${base2}_nobacsort.bam
samtools fastq ${base2}_nobacsort.bam -1 ${base2}_1.fastq.gz -2 ${base2}_2.fastq.gz

# Fungal reads alignment
bowtie2 -x /hl/weiguanyue/genome/fungi/2_dereplicated_fungi/00_redundancy/00/00_databases/index/all_3156  -p 3  -1 ${fqs1} -2 ${fqs2} | samtools view  -bS -o ${base2}.bam
samtools view -b -h -F 4 -F  256  ${data}/${base2}.bam >${base2}_mapped.bam

# quilty MAPQ 30
samtools view -q 30 ${base2}_mapped.bam >${base2}_30mapped.sam
samtools view -bq 30 ${base2}_mapped.bam >${base2}_30mapped.bam
python coverage80filter.py -i  ${base2}_30mapped.sam -o ${base2}_list

# exact coverage 80%
samtools view  ${base2}_30mapped.bam |grep -f  ${base2}_list >${base2}_30mappedc80.sam
samtools view -bt  3156_rename.fa.fai  ${base2}_30mappedc80.sam -o ${base2}_30mappedc80.bam
samtools sort ${base2}_30mappedc80.bam -o  ${base2}_30mappedc80sort.bam
samtools index ${base2}_30mappedc80sort.bam

# exact genomes reads
samtools idxstats  ${base2}_30mappedc80sort.bam >${base2}_idx.txt

awk '$3 > 0'  ${base2}_idx.txt|awk '{print $1"_"$3}'|awk -F _ '{print$1"_"$2,$4}' >${base2}_hitall.txt

# sum same reads
awk '{a[$1]+=$2} END {for(i in a) print i" "a[i]b[i]}' ${base2}_hitall.txt >${base2}_hit.txt

# Taxonomy
buglist.py -i  ${base2}_hit.txt -t final_taxonomy.txt  -o taxonomy.txt

# Abundance calculation
coverm genome --bam-files   ${base2}_30mappedc80sort.bam  -m rpkm -o ${base2}_rpkm.tsv --min-covered-fraction 0 --genome-fasta-files   *.fa   -t 2
