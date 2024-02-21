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
perl limit2Length.pl -f assembled_scaffold.fasta -len 500 -o scaffold_500.fasta 
prodigal -f gff -i scaffold_500.fasta -o scaffold_500.fasta.gff -p meta -a scaffold_500.faa -d scaffold_500.fna
hmmsearch mcrA.hmm scaffold_500.faa > scaffold_500_mcrA.out
perl extract_fasta_from_list.pl -f scaffold_500.faa -l scaffold_500_mcrA.list -o mcrA.faa (used for phylogenetic analysis)
perl extract_fasta_from_list.pl -f scaffold_500.fna -l scaffold_500_mcrA.list -o mcrA.fna (used for relative abundance analysis)
cd-hit -i mcrA.fna -o ref_mcrA.fna -c 0.99 -d 20 -T 32

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
iqtree-omp -s aligned_combined_mcrA.faa -nt 32 -m WAG -bb 1000

# Evaluation of the relative abundance and activity of ‘Ca. Methanoliparum’

## Relative abundance and expression activity of dereplicated MAGs
dRep dereplicate outout_directory -g /bins/*.fa -p 32 -comp 50 -con 10 -sa 0.97 -genomeInfo ./checkM_gtdbtk_dRep.csv
bowtie2-build dRepMAGs.fa scaffold.fa
bowtie2 -x scaffold.fa -1 1.clean.fq -2 2.clean.fq -S all.sam -p 128
samtools view -bS all.sam > all.bam -@ 128
samtools sort all.bam -o all_sort.bam -@ 128
samtools index all_sort.bam
perl length+GC.pl -f dRepMAGs.fa -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat

## Expression activity of AcrA/McrA, AssA-related genes
bwa index ref_mcrA.fna
bwa mem -t 32 ref_mcrA.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
samtools view -bS scaffold.sam > all.bam -@ 32
samtools sort all.bam -o all_sort.bam -@ 32
samtools index all_sort.bam
length+GC.pl -f ref_mcrA.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
