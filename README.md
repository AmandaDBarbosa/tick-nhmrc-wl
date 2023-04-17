# tick-nhmrc-wl-(In Qiime2)


#Analysis of 16S rRNA (hypervariable region v1-2) metabarcoding.

#Raw Illumina MiSeq `.fastq.gz` reads analysed using QIIME2-2020.11 pipeline using dada2 denoising to create ASVs.

**Background**

This workflow is written for analysing amplicon data in QIIME2. Input data is Illumina MiSeq paired-end data prepared using Nextera XT indexes (i.e. no additional demultiplexing steps are needed in this case however should your data require demultiplexing it can easily be added in).

## 0. Install & activate QIIME2 environment (commandline) 

This workflow utilsing commandline interface with QIIME2.

Requires miniconda/conda, see [here](https://docs.qiime2.org/2020.11/install/native/#installing-miniconda)

Latest version = [QIIME2-2022.11](https://https://docs.qiime2.org/2022.8/), see QIIME2 documentation for install based on your platform.

**Activate qiime2 environment**
```bash
conda activate qiime2-2022.11
```

## 1. Input data

> Assumes paired-end data that does not require demultiplexing

Place raw data files in zipped format (i.e. `.fastq.gz` in a directory called `raw-seqs/WL-ticks/`).

### File naming conventions

In [Casava 1.8 demultiplexed (paired-end)](https://docs.qiime2.org/2020.11/tutorials/importing/#casava-1-8-paired-end-demultiplexed-fastq) format, there are two `.fastq.gz` files for each sample in the study, each containing the forward or reverse reads for that sample. The file name includes the sample identifier. The forward and reverse read file names for a single sample might look like `XXXX_L001_R1_001.fastq.gz` and `XXXX_L001_R2_001.fastq.gz`, respectively. The underscore-separated fields in this file name are:

- the sample identifier,
- the barcode sequence or a barcode identifier,
- the lane number,
- the direction of the read (i.e. R1 or R2), and
- the set number.

Depending on sequencing facility you may need to add the `_001` prefix to sample files.

Note however that you do **not** need to unzip fastq data to analyse.

Navigate into the directory with raw data files:
```bash
for file in ~/raw-seqs/WL-ticks/*.fastq.gz;
do
newname=$(echo "$file" | sed 's/0_BPDNR//' | sed 's/.fastq/_001.fastq/')
mv $file $newname
done
```

### Import as QIIME2 artefact

Import `.fastq.gz` data into QIIME2 format using [Casava 1.8 demultiplexed (paired-end)](https://docs.qiime2.org/2020.11/tutorials/importing/#casava-1-8-paired-end-demultiplexed-fastq) option. Remember assumes raw data is in directory labelled `raw_seqs/LoD/` and file naming format as above.

```bash
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ~/raw-seqs/WL-skin-seqs \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path 16S_demux_seqs.qza
```

In this case we are using Nextera Indexes which mean they are demultiplexed automatically by basespace and therefore we can skip over any reference to demultiplexing steps.

**Inspect reads for quality**
To inspect raw reads

```bash
qiime demux summarize \
  --i-data 16S_demux_seqs.qza \
  --o-visualization 16S_demux_seqs.qzv
```

View this output by importing into [QIIME2 view](https://view.qiime2.org/). Use this output to choose your parameters for QC such as trimming low quality sequences and truncating sequence length.

### Sample metadata

This holds you associated metadata related to your samples (e.g. host information, sampling data, etc). [Tutorial here](https://docs.qiime2.org/2020.11/tutorials/moving-pictures/#sample-metadata)

The metadata needs to be in `.tsv` format, the best way to do this is to access the QIIME2 googlesheet example. Save a copy and edit/add in your sample details. Then select `File > Download as > Tab-separated values`. Alternatively, the  command `wget "https://data.qiime2.org/2020.11/tutorials/moving-pictures/sample_metadata.tsv"` will download the sample metadata as tab-separated text and save it in the file `sample-metadata.tsv`. It is import you don't change the header for the first column `sample-id`.


## 2. Sequence quality control and feature table construction

> Denoise using dada2

Based on quality plot in the above output `16S_demux_seqs.qza` adjust trim length to where quality falls.

Then you can also trim primers. In this case working with 16S v1-2 data with the following primers

Example data - amplicon NGS data targeting bacteria using 16S rRNA hypervariable region 1-2 with the following primers:

- 27F-Y (20 nt): AGAGTTTGATCCTGGCTYAG #16S v1-2 primer, ref Gofton et al. Parasites & Vectors (2015) 8:345
- 338R (19 nt): TGCTGCCTCCCGTAGGAGT #16S v1-2 primer, ref Turner et al. J Eukaryot Microbiol (1999) 46(4):32

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S_demux_seqs.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 19 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --o-table 16S_denoise_table.qza \
  --o-representative-sequences 16S_denoise_rep-seqs.qza \
  --o-denoising-stats 16S_denoise-stats.qza
```

At this stage, you will have artifacts containing the feature table, corresponding feature sequences, and DADA2 denoising stats. You can generate summaries of these as follows.

```bash
qiime feature-table summarize \
  --i-table 16S_denoise_table.qza \
  --o-visualization 16S_denoise_table.qzv \
  --m-sample-metadata-file metadata-tick.csv # Can skip this bit if needed.

qiime feature-table tabulate-seqs \
  --i-data 16S_denoise_rep-seqs.qza \
  --o-visualization 16S_denoise_rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file 16S_denoise-stats.qza \
  --o-visualization 16S_denoise-stats.qzv
```

### Merging denoised artefacts (SKIP - AB-XB 160123)

To merge denoised data sets and generate one `FeatureTable[Frequency]` and `FeatureData[Sequence]`  artifacts

```bash
qiime feature-table merge \
  --i-tables table-1.qza \
  --i-tables table-2.qza \
  --o-merged-table table.qza
qiime feature-table merge-seqs \
  --i-data rep-seqs-1.qza \
  --i-data rep-seqs-2.qza \
  --o-merged-data rep-seqs.qza
```

### Export ASV table

To produce an ASV table with number of each ASV reads per sample that you can open in excel. Use tutorial [here](https://rstudio-pubs-static.s3.amazonaws.com/489645_5fff8a6a02d84084a55e3b5b6ff960a4.html#)

Need to make biom file first

```bash
qiime tools export \
--input-path 16S_denoise_table.qza \
--output-path feature-table

biom convert \
-i feature-table/feature-table.biom \
-o feature-table/feature-table.tsv \
--to-tsv
```

### Phylogenetic tree

Several downstream diversity metrics, available within QIIME 2, require that a phylogenetic tree be constructed using the Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs) being investigated. Documentation [here](https://docs.qiime2.org/2020.11/tutorials/phylogeny/)

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 16S_denoise_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

**Export**

Covert unrooted tree output to newick formatted file

```bash
qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path exported-tree
```

## 3.Taxonomy

Assign taxonomy to denoised sequences using a pre-tarined naive bayes classifier and the [q2-feature-classifier](https://docs.qiime2.org/2020.11/plugins/available/feature-classifier/) plugin. Details on how to create a classifier are available [here](link).

Note that  taxonomic classifiers perform best when they are trained based on your specific sample preparation and sequencing parameters, including the primers that were used for amplification and the length of your sequence reads.

```bash
qiime feature-classifier classify-sklearn \
--i-classifier ../training-feature-classifiers/gg_13_8_otus/classifier.qza \
--i-reads 16S_denoise_rep-seqs.qza \
--o-classification taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv
```


In order to be able to download the sample OTU table need to do the taxonomy assignment and then make the taxa barplot. Then can download csv file with sequence number, samples and taxonomy.
[see here](https://docs.qiime2.org/2020.11/tutorials/moving-pictures/#taxonomic-analysis)

```
qiime taxa barplot \
  --i-table 16S_denoise_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata-tick.tsv \
  --o-visualization taxa-bar-plots.qzv
```

Details on sample metadata available [here](https://docs.qiime2.org/2020.11/tutorials/moving-pictures/#sample-metadata)


Extra bit of code to generate a taxonomy table table to tsv from the commandline

```bash
qiime tools export \
--input-path taxonomy.qza \
--output-path exports
```

#Alpha Rarefaction Plot

qiime diversity alpha-rarefaction \
--i-table 16S_denoise_table.qza \
--i-phylogeny rooted-tree.qza \
--m-metadata-file metadata-tick.tsv \
--p-max-depth 40000 \
--o-visualization alpha-rarefaction-plot.qzv

#----------------------------------------------------------------------------

###See if I have all outputs as in this template, including alpha rarefaction plot (https://view.qiime2.org/)
### Now do Contaminants and Phyloseq (R)

