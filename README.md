# Qiime2 pipeline
## Step 1. Data input
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path se-33-manifest \
  --output-path clean_data/single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33

## Step 2. Quality control, denoise, chimaera removal and feature clustering
qiime dada2 denoise-single \
  --i-demultiplexed-seqs clean_data/single-end-demux.qza \
  --p-trunc-len 291 \
  --p-trunc-q 2 \
  --o-representative-sequences clean_data/rep-seqs-dada2.qza \
  --o-table clean_data/table-dada2.qza \
  --output-dir clean_data/dada2_unspecified_re \
  --p-n-threads 30 

## Step 3. Assign taxonomy, selection of representative sequences
qiime feature-classifier classify-sklearn \
  --i-classifier clean_data/silva-138-99-F515R806-classifier.qza \
  --i-reads clean_data/rep-seqs-dada2.qza \
  --o-classification clean_data/taxonomy.qza \
  --p-n-jobs 30 

## Step 4. Barplot
qiime taxa barplot \
  --i-table clean_data/table-dada2.qza \
  --i-taxonomy clean_data/taxonomy.qza \
  --m-metadata-file sample.tsv \
  --o-visualization clean_data/taxa-bar-plots.qzv 

## Step 5. BLAST
makeblastdb -in reference.fasta -dbtype nucl
blastn -query feature_sequences.fasta -db reference.fasta -out blast_results.txt -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30

# Metagenomic analysis for single sample each time
## Step 1. Trimming sequences
java -jar ~/Trimmomatic-0.38/trimmomatic-0.38.jar PE raw_data1/2 output_forward/reverse_paired.fq.gz
output_forward/reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;

## Step 2. Assembly
python ~/SPAdes-3.12.0-Linux/bin/spades.py --meta -k 21,33,55,77 -t 10 -1 output_forward_paired.fq.gz -2 
output_reverse_paired.fq.gz -o spa_res;

## Step 3. Metadata geneprediction
prodigal -i scaffold.fasta -a prot.fasta -d nucl.fasta-p meta -q;

## Step 4. Binning
metawrap binning -o 04INITIAL_BINNING -t 60 -a 02ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct --universal --run-checkm 01READ_QC/final_pure_reads_1.fastq 01READ_QC/final_pure_reads_2.fastq

## Step 5. genome bins summary
checkm lineage_wf -x fa -t 40 checkm_res;

## Step 6. gtdbtk taxnomy
gtdbtk classify_wf --genome_dir gtdbtk_test/genomes --out_dir gtdbtk_test/output --cpus 10;

# Identification of genes coding for enzymes with alkane activation potentials
perl /scripts/SeqTools/limit2Length.pl -f assembled_scaffold.fasta -len 500 -o scaffold_500.fasta 
prodigal -f gff -i scaffold_500.fasta -o scaffold_500.fasta.gff -p meta -a scaffold_500.faa -d scaffold_500.fna
hmmsearch mcrA.hmm scaffold_500.faa > scaffold_500_mcrA.out
perl extract_fasta_from_list.pl -f scaffold_500.faa -l scaffold_500_mcrA.list -o mcrA.faa (used for phylogenetic analysis)
perl extract_fasta_from_list.pl -f scaffold_500.fna -l scaffold_500_mcrA.list -o mcrA.fna (used for relative abundance analysis)
/software/cdhit-4.6.2/cd-hit -i mcrA.fna -o ref_mcrA.fna -c 0.99 -d 20 -T 32

# Phylogenomic analyses
## Using concatenated ribosomal proteins (16 r-proteins)
checkm lineage_wf -x .fa /bins checkm_result -f checkm_result.txt -t 30 --pplacer_threads 30
perl /script/find_ribosome_faa_in_checkM2.0.pl -f checkm_result
cat Ribosomal*.txt > all_ribosomal.fas
blastp -query all_ribosomal.fas -db /new_tree_of_life/find_ribosomal/all_ribosomal.fas.clean -out all_ribosonmal_blast_new_tree.txt -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30
mkdir ribosomal_for_tree/
perl /script/pick_ribosomal_with_blast_v2.0.pl
Input the fasta file:all_ribosomal.fas
Input the blast result:all_ribosonmal_blast_new_tree.txt
Input the out file:ribosomal_for_tree
cd ribosomal_for_tree/
mkdir aligned_ribosomals/
muscle -in ribosomal_for_tree_Ribosomal_L*(S*).fas -out aligned_ribosomals/ribosomal_for_tree_Ribosomal_L*(S*).muscle.fas
cd aligned_ribosomals/
trimal -in ribosomal_for_tree_Ribosomal_L*(S*).muscle.fas -out ribosomal_for_tree_Ribosomal_L*(S*).muscle.trim.fas -gt 0.95
perl /script/combine_aligned_ribosomal_v2.0.pl
Input alignment folder:./
/tools/iqtree/bin/iqtree-omp -s ref_ribosomal_combine.phy -nt 30 -m WAG -bb 1000

## Phylogenetic trees for AcrA/McrA, AssA related, BcrB, and BcrC protein sequence
makeblastdb -in mcrA_reference_protein_seqs.faa -dbtype prot
blastp -query bin.fa.faa -db mcrA_reference_protein_seqs.faa -out bin.fa.txt -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30
muscle -in combined_mcrA.faa -out aligned_combined_mcrA.faa
/tools/iqtree/bin/iqtree-omp -s aligned_combined_mcrA.faa -nt 32 -m WAG -bb 1000

# Evaluation of the relative abundance and activity of ‘Ca. Methanoliparum’

## 16S rRNA gene sequences extraction and assignment from metagenomic data
kraken2-build --db /software/kraken2_db/silva_db --special silva
kraken2 --db /software/kraken2_db/silva_db/16S_SILVA138_k2db  --paired --classified-out cseqs.fq 1.clean.fq 2.clean.fq
kraken2 --db /software/kraken2_db/silva_db/16S_SILVA138_k2db  cseqs_1.fq cseqs_2.fq --output 16S_rRNA_copy.selection.fasta.out.txt --use-names --report sample_16S.txt --threads 32 

## Relative abundance and expression activity of dereplicated MAGs
dRep dereplicate outout_directory -g /bins/*.fa -p 32 -comp 50 -con 10 -sa 0.97 -genomeInfo ./checkM_gtdbtk_dRep.csv
/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build dRepMAGs.fa scaffold.fa
/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x scaffold.fa -1 1.clean.fq -2 2.clean.fq -S all.sam -p 128
/tools/samtools-1.3.1/samtools view -bS all.sam > all.bam -@ 128
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 128
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f dRepMAGs.fa -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat

## Expression activity of all annoted genes of four ‘Ca. Methanoliparum’ representative MAGs
prodigal -f gff -i four_rep_MAGs.fasta -o four_rep_MAGs.fasta.gff -p meta -a four_rep_MAGs.faa -d four_rep_MAGs.fna
bwa index four_rep_MAGs.fna
bwa mem -t 32 four_rep_MAGs.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f four_rep_MAGs.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat

## Expression activity of all annoted genes of assembled scaffold
perl /scripts/SeqTools/limit2Length.pl -f assembled_scaffold.fasta -len 500 -o scaffold_500.fasta 
prodigal -f gff -i scaffold_500.fasta -o scaffold_500.fasta.gff -p meta -a scaffold_500.faa -d scaffold_500.fna
bwa index scaffold_500.fna
bwa mem -t 32 scaffold_500.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f scaffold_500.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat

## Expression activity of AcrA/McrA, AssA-related and BcrB/BcrC 
bwa index ref_mcrA.fna
bwa mem -t 32 ref_mcrA.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f ref_mcrA.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat

# Metagenomic analysis for single sample each time

## Step 1. Dereplicate sequences

```
perl ~/dereplicate.pl -fq raw_data1/2.fq -out raw_data1/2.d.fastq;
```

## Step 2. Trimming sequences

```
sickle pe -f raw_data1.d.fastq -r raw_data2.d.fastq -t sanger -q 25 -o raw_data1.td.fastq -p raw_data2.td.fastq -s raw_data.s.fastq;
```

## Step 3. Assembly

```
perl ~/parseFastq.pl -i raw_data1/2.dt.fastq -f raw_data1/2.dt.fasta
perl ~/interleave.pl -fwd raw_data1.dt.fasta -rev raw_data2.dt.fasta -o raw_data_int.fasta
idba_ud -o k65-k145 --mink 65 --maxk 145 --step 10 -r raw_data_int.fasta;
```

## Step 4. Metadata gene prediction

```
perl ~/limit2Length.pl -f scaffold.fa -len 500 -o scaffold_500.fa
prodigal -i scaffold_500.fa -a prot.fasta -d nucl.fasta-p meta -q;
```

## Step 5. Binning

### Step 5.1. Mapping

```
~/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build scaffold.fa scaffold.fa
~/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x scaffold.fa -1 raw_data1.fq -2 raw_data2.fq -S scaffold.sam;
```

### Step 5.2. Sorting

```
~/samtools-1.3.1/samtools view -bS scaffold.sam >scaffold.bam
~/samtools-1.3.1/samtools sort scaffold.bam -o scaffold_sort.bam
~/samtools-1.3.1/samtools index scaffold_sort.bam;
```

### Step 5.3. Depth calculation

```
~/jgi_summarize_bam_contig_depths --outputDepth scaffold_depth.txt --pairedContigs paired.txt scaffold_sort.bam
```

### Step 5.4 binning

```
metabat -i scaffold.fa -a scaffold_depth.txt -o Binning/bin;
```

## Step 6. genome bins summary

```
checkm lineage_wf -x fa -t 40 checkm_res;
```

## Step 7. derep

```
dRep dereplicate dRepMAGs -g /*.fa -p 32 -comp 50 -con 10 -sa 0.97 -genomeInfo ./checkM_gtdbtk_dRep.csv;
```

## Step 8. gtdbtk taxnomy

```
gtdbtk classify_wf --genome_dir gtdbtk_test/genomes --out_dir gtdbtk_test/output --cpus 10;
```

# Identification of hgcA

```
hmmsearch hgcA.hmm scaffold_500.faa > scaffold_500_hgcA.out
perl extract_fasta_from_list.pl -f scaffold_500.faa -l scaffold_500_hgcA.list -o hgcA.faa (used for phylogenetic analysis)
perl extract_fasta_from_list.pl -f scaffold_500.fna -l scaffold_500_hgcA.list -o hgcA.fna (used for relative abundance analysis)
```

# Phylogenomic analyses

## Phylogenetic tree for hgcA carrying MAGs using concatenated ribosomal proteins (16 r-proteins)

```
checkm lineage_wf -x .fa /bins checkm_result -f checkm_result.txt -t 30 --pplacer_threads 30
perl /script/find_ribosome_faa_in_checkM2.0.pl -f checkm_result
cat Ribosomal*.txt > all_ribosomal.fas
blastp -query all_ribosomal.fas -db /new_tree_of_life/find_ribosomal/all_ribosomal.fas.clean -out all_ribosonmal_blast_new_tree.txt -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30
mkdir ribosomal_for_tree/
perl /script/pick_ribosomal_with_blast_v2.0.pl
Input the fasta file:all_ribosomal.fas
Input the blast result:all_ribosonmal_blast_new_tree.txt
Input the out file:ribosomal_for_tree
cd ribosomal_for_tree/
mkdir aligned_ribosomals/
muscle -in ribosomal_for_tree_Ribosomal_L*(S*).fas -out aligned_ribosomals/ribosomal_for_tree_Ribosomal_L*(S*).muscle.fas
cd aligned_ribosomals/
trimal -in ribosomal_for_tree_Ribosomal_L*(S*).muscle.fas -out ribosomal_for_tree_Ribosomal_L*(S*).muscle.trim.fas -gt 0.95
perl /script/combine_aligned_ribosomal_v2.0.pl
Input alignment folder:./
/tools/iqtree/bin/iqtree-omp -s ref_ribosomal_combine.phy -nt 30 -m WAG -bb 1000
```

## Phylogenetic trees for HgcA protein sequence

```
muscle -in combined_hgcA.faa -out combined_hgcA.muscle.faa
trimal -in combined_hgcA.muscle.faa -out combined_hgcA.muscle.trim.faa -gt 0.95
/tools/iqtree/bin/iqtree-omp -s combined_mcrA.muscle.trim.faa -nt 32 -m LG+F+G4 -bb 1000
```

# Evaluation of the relative abundance and activity of hgcA carrying MAGs

## Relative abundance and expression activity of hgcA carrying MAGs

```
/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build hgcAMAGs.fa scaffold.fa
/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x scaffold.fa -1 1.clean.fq -2 2.clean.fq -S all.sam -p 128
/tools/samtools-1.3.1/samtools view -bS all.sam > all.bam -@ 128
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 128
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f hgcAMAGs.fa -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
```

## Expression activity of hgcA genes 

```
bwa index hgcA.fna
bwa mem -t 32 hgcA.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f hgcA.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
```
