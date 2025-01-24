# BiGEST

BiGEST (Biosynthetic Gene Eukaryotic Search Tool) is a BLAST-based tool that will search any genome fasta file and annotate them with Nonribosomal peptide synthetases (NRPS) and polyketide synthases (PKS). BiGEST uses a subset of the [NCBI CDD Database](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) to search for domains that are a part of NRPS/PKS Gene clusters, according to [Synthaser: a CD-Search enabled Python toolkit for analysing domain architecture of fungal secondary metabolite megasynth(et)ases](https://pubmed.ncbi.nlm.nih.gov/34763725/). 

BiGEST can also take an antiSMASH output and overlay those results with BiGEST annotated domains to get a fully comprehensive look at the genome given.

## Installation

- Clone the repo
- Install dependencies

The required dependencies are:

- python 3.6+
- blast+ 2.14.0+
- biopython 1.81+
- antiSMASH 7.1.0+ 

These can be downloaded individually or as a conda package with :
    BiGEST_environment.yml

Due to an issue	in the antiSMASH conda download, you need to downgrade the GFF parser	(bcbio-gff) in order for antismash to run correctly.

```
conda env create --name BiGEST --file=BiGEST_environment.yml
conda activate BiGEST
conda install bcbio-gff=0.7.0
```

## Usage

At its most basic, BiGEST requires an rpstblastn search, fasta file, and output directory. However, there are multiple ways BiGEST can be run. The blast input can be pre-run with the blast scripts, or at the same time as the entire BiGEST run. 

The BiGEST annotation workflow is typically: 

- **BLAST** analysis on the given fasta file against the full CDD database, or the subset of Sythaser-identified-domains included in this repository can be used instead (seen at db/BiGEST). 

    This can be done beforehand but it is recommended to run the rpstblastn search simultaenously with the BiGEST run to ensure the correct formatting of the blast output  see section "Reminders" before running BiGEST.

- **antiSMASH** analysis on the given fasta file

    This can be done beforehand and the output GenBank file can be used as BiGEST input. 

- **Filtering and Compiling** the BLAST and antiSMASH results as necessary to only include hits relevant to NRPS and PKS domains. 

- **Visualizing** with a PDF and GenBank output. 

Fasta input can be gzipped or not.

### To run BiGEST with blast and antiSMASH simultaneously

```
./BiGEST.py\
    -d /home/user/BiGEST/db/BiGEST_DB
    -f /home/user/fasta_file/my_fa.fasta
    -o /home/user/BiGEST/output/
    -a True
```

### To run BiGEST with blast and antiSMASH already done

```
./BiGEST.py\
    -i /home/user/BiGEST/blast_output/my_fasta_rpstblastn.txt
    -f /home/user/fasta_file/my_fasta.fa
    -o /home/user/BiGEST/output/
    -a /home/user/BiGEST/antiSMASH_output/my_fasta/my_fasta.gbk
```

### Reminders before running BiGEST:

 - Most importantly, if BLAST has already been run and will be an input to BiGEST, there are several positions in the output required for BiGEST to run. The BLAST command should use output format 6 with options in this order: 
    
        "-outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle gaps qseq sseq sacc slen"

### Expected Output

Depending on the size of your input fasta, BiGEST can take anywhere from a few seconds to many minutes to run. If blast or antiSMASH needs to be run, it will take longer. Using the included subsetted blast database (db/BIGEST_DB) in place of the full CDD blast database will decrease your run time incredibly.

When the script is done, there should be mutliple files, all named with the same basename as your fasta file. These files include: 

 - a GenBank file with the BGC matches found by BiGEST
 - a text file with the rpstblastn hits that were qualified as BiGEST BGC hits
 - a BED file with the locations of each blast (and antiSMASH if given) domain related to BGCs
 - a GFF3 file with the locations of each blast domain related to BGCs
 - a folder, 'figures', with pdf images displaying each BGC found
 
 If antiSMASH was included in the run:

 - a "combined" GenBank file with the BGC matches found by both BiGEST and antiSMASH to allow for easy viewing with overlaid results
 

Sample scripts of these files can be found in Example/output of this github repository. In the Example directory, there are instructions to test the installation with an example algal genome, [Chlorella sorokiniana UTEX 1230](https://phycocosm.jgi.doe.gov/Chloso1230_1_1/Chloso1230_1_1.home.html). 