# nanoHiMe

Here we provide a tutorial about how to use nanoHiMe to detect 6mA and CpG methylations from nanopore sequencing data. More information could be seen at http://github.com/YinLab/nanoHiMe upon request.



Dependencies:

perl v5.26
Guppy v4.4.2
nanopolish v0.13.2
minimap2 v2.17
bedtools v2.26.0


Requirements:

human reference genome hg19 (UCSC)
samtools v1.9
Ecoli reference genome K12_MG1655


nanoHiMe modules:

nanoHiMe call-CpG-methylation: predict methylated and unmethylated CpG sites in the genome
nanoHiMe call-6mA-modification: call 6mA-containing regions from the genome
nanoHiMe train-model: learn parameters for the emission distributions for each k-mer



Examples of analysis workflow:

---------------PartI. Data preparation: from fast5 to nanopolish eventalign files ----------------------------

Sample information:

Name: example_hg19.fast5
Path: nanoHiMe/example/fast5
Source: HepG2 cell line
Experiment: nanoHiMe-seq using H3K27me3 antibody


In the first step, perform basecalling from FAST5 files using the following Guppy command:

cd nanoHiMe/example
GPU mode: guppy_basecaller --input_path fast5 --save_path output_fast5 -c dna_r9.4.1_450bps_hac.cfg --qscore_filtering -x cuda:0,1
CPU mode (16 thread): guppy_basecaller --input_path fast5 --save_path output_fast5 -c dna_r9.4.1_450bps_hac.cfg --qscore_filtering --cpu_threads_per_caller 16

After basecalling, combine all of the fastq files into a single file:

cat output_fast5/*/*.fastq > test.fastq



In the second step, use nanopolish to create an index file linking individual reads to their signal-level data stored in the FAST5 files:

nanopolish index -d fast5 test.fastq

The following new files were created: test.fastq.index, test.fastq.index.fai, test.fastq.index.gzi, test.fastq.index.readdb



In the third step, use minimap2 to align the basecalled reads to the reference genome (e.g. hg19):

wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz' -O hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
minimap2 -ax map-ont -t 8 hg19.fa test.fastq | /data/software/samtools-1.9/samtools sort -o test.sorted.bam -T test.tmp
samtools index test.sorted.bam
samtools view -b -q 20 -F 4 test.sorted.bam > test.sorted.q20.mapped.bam
samtools index test.sorted.q20.mapped.bam


In the fourth step, align the nanopore events to a reference:

nanopolish eventalign -n -t 16 --reads test.fastq --bam test.sorted.q20.mapped.bam --genome hg19.fa --scale-events > test.eventalign.txt


---------------PartII. NanoHiMe modification calling ----------------------------


Before modification calling, add htslib to LD_LIBRARY_PATH:

example:
export LD_LIBRARY_PATH=/user/nanoHiMe/bin/htslib
"/user/nanoHiMe" should be replaced by the program path on your own computer.


In the last step, use nanoHiMe to predict modified bases, 5-methylcytosine in CpG context and N6-methyladenine in all contexts:

6mA-containing region calling:


#ref.fa is the reference genome fasta file, such as hg19.fa and E.coli_K12_MG1655.fasta
perl nanoHiMe/perl_script/upper.pl ref.fa > ref_upper.fa  #convert lower case letters of bases to upper case
samtools faidx ref_upper.fa
nanoHiMe/nanoHiMe 6mA ref_upper.fa test.eventalign.txt(.gz) output.6mA 50 25 [input.peak.bed]  # test.eventalign.txt is created with nanopolish , could be in gzip format


If you set input.peak.bed , 6mA calling will be limited in these peak regions. Otherwise, 6mA calling will be althrough the whole genome and may be quite slow.

The following output file were created: output.6mA.methylation.txt

The information of predicted 6mA-containing regions was stored in test.win.6mA_region.txt file:

CHR     Chain   Start   End     Read_ID Log(LIKELIHOOD_RATIO)_mA        Ref_sequence
CP014225.1      +       361642  361691  2532edf7-19f5-4911-b4ab-1194e071e2ff    -2.21   TGATCCTGTTAGATCTGATGCTCCCTGGCACCGATGGCCTGACGCTGTG
CP014225.1      +       361667  361716  2532edf7-19f5-4911-b4ab-1194e071e2ff    -5.89   CTGGCACCGATGGCCTGACGCTGTGCCGGGAAATTCGTCGTTTTTCTGAC
CP014225.1      +       361692  361741  2532edf7-19f5-4911-b4ab-1194e071e2ff    -18.25  CCGGGAAATTCGTCGTTTTTCTGACATTCCGATCGTGATGGTGACGGCAA
CP014225.1      +       361717  361766  2532edf7-19f5-4911-b4ab-1194e071e2ff    -17.77  ATTCCGATCGTGATGGTGACGGCAAAAATCGAAGAGATCGATCGCCTGCT
CP014225.1      +       361742  361791  2532edf7-19f5-4911-b4ab-1194e071e2ff    10.11   AAATCGAAGAGATCGATCGCCTGCTGGGGCTGGAGATTGGCGCAGATGAT
CP014225.1      +       361767  361816  2532edf7-19f5-4911-b4ab-1194e071e2ff    -9.04   GGGGCTGGAGATTGGCGCAGATGATTATATCTGTAAGCCGTACAGCCCAC
CP014225.1      +       361792  361841  2532edf7-19f5-4911-b4ab-1194e071e2ff    13.63   TATATCTGTAAGCCGTACAGCCCACGGGAAGTGGTAGCGCGCGTCAAAAC
CP014225.1      +       361817  361866  2532edf7-19f5-4911-b4ab-1194e071e2ff    10.65   GGGAAGTGGTAGCGCGCGTCAAAACCATTTTGCGCCGTTGCAAACCGCAG
CP014225.1      +       361842  361891  2532edf7-19f5-4911-b4ab-1194e071e2ff    23.85   CATTTTGCGCCGTTGCAAACCGCAGCGCGAGTTGCAGCAACAGGATGCTG
...

A positive value in LOG(LIKELIHOOD_RATIO)_mA column shows the support for 6mA methylation.


mCG-calling mode:

./nanoHiMe mCG ref_upper.fa test.eventalign.txt(.gz) output.mCG

The output file, named output.mCG.eventalign.txt, contains CpG methylation information of CpG cotaining regions :

CHR     Chain   Start   End     Read_ID Log(LIKELIHOOD_RATIO)_mCG       Ref_sequence
CP014225.1      +       361664  361684  2532edf7-19f5-4911-b4ab-1194e071e2ff    10.90   TCCCTGGCACCGATGGCCTGA
CP014225.1      +       361675  361715  2532edf7-19f5-4911-b4ab-1194e071e2ff    44.25   GATGGCCTGACGCTGTGCCGGGAAATTCGTCGTTTTTCTGA
CP014225.1      +       361711  361735  2532edf7-19f5-4911-b4ab-1194e071e2ff    27.31   TCTGACATTCCGATCGTGATGGTGA
CP014225.1      +       361726  361769  2532edf7-19f5-4911-b4ab-1194e071e2ff    50.60   GTGATGGTGACGGCAAAAATCGAAGAGATCGATCGCCTGCTGGG
CP014225.1      +       361772  361792  2532edf7-19f5-4911-b4ab-1194e071e2ff    22.47   TGGAGATTGGCGCAGATGATT
CP014225.1      +       361795  361815  2532edf7-19f5-4911-b4ab-1194e071e2ff    6.10    ATCTGTAAGCCGTACAGCCCA
CP014225.1      +       361806  361826  2532edf7-19f5-4911-b4ab-1194e071e2ff    -0.53   GTACAGCCCACGGGAAGTGGT
CP014225.1      +       361819  361843  2532edf7-19f5-4911-b4ab-1194e071e2ff    28.19   GAAGTGGTAGCGCGCGTCAAAACCA
CP014225.1      +       361839  361879  2532edf7-19f5-4911-b4ab-1194e071e2ff    50.77   AACCATTTTGCGCCGTTGCAAACCGCAGCGCGAGTTGCAGC
CP014225.1      +       361888  361928  2532edf7-19f5-4911-b4ab-1194e071e2ff    33.43   GCTGAAAGCCCGTTGATTATCGACGAAGGTCGTTTTCAGGC
CP014225.1      +       361926  361948  2532edf7-19f5-4911-b4ab-1194e071e2ff    13.58   GGCTTCATGGCGCGGTAAAATG
CP014225.1      +       361948  361982  2532edf7-19f5-4911-b4ab-1194e071e2ff    49.51   CTTGACCTGACGCCTGCGGAATTTCGTCTGCTG
CP014225.1      +       361975  362005  2532edf7-19f5-4911-b4ab-1194e071e2ff    8.82    CTGCTGAAAACGCTCTCTCACGAACCAGGA
CP014225.1      +       362007  362029  2532edf7-19f5-4911-b4ab-1194e071e2ff    14.20   AGTGTTCTCCCGCGAGCAATT
...

A positive value in Log(LIKELIHOOD_RATIO)_mCG shows the support for mCG methylation, and a negative value indicates the support for unmethylation.


---------------PartIII. NanoHiMe parameters training ----------------------------

To train your own model, run the model training model :

Model_training \
        -m mcg \
        -c 5 \
        -r /.../ref_upper.fa \
        -e /.../training.eventalign.txt.gz \
        -o /.../output_dir \
        -g 3 \
        -p /.../EM_shift.txt \
        -i 500

----------------------------------------------------------------------
-m      <str>   Modification : ma, mcg or ma_mcg
-c      [int]   Cycles to re-alignment events :  [default: 5]
-r      <str>   Reference fasta with samtools index (abs_path)
-e      <str>   Input eventalign files from nanopolish (abs_path)
-o      [str]   Output directory (abs_path)
-g      [int]   Number of gaussian distribution in each 6-mer [default: 3]
-p      [str]   File with initial values of EM parameters (abs_path) [default: nanoHiMe/perl_script/EM_shift.txt]
-i      [int]   Max iterations of EM [default: 1000]

-r,-e,-o,-p : please use absolute path


The output file EM.out will contain the parameters of each gaussian distribution for each 6-mer:

6-mer   omega1  mu1     sigma1  omega2  mu2     sigma2  omega3  mu3     sigma3
AAAAAA  0.925549559849263       86.4533687505569        1.21282717235525        0.0372252200753692      85.7714165032309        1.53680334781428        0.0372252200753692      85.7714165032309        1.53680334781428
AAAAAC  0.606430807125636       83.8422378791128        1.20582910220601        0.0794899443967189      83.6039521077065        2.84518668032173        0.314079248477646       83.7509059178093        1.97676968685331
...
