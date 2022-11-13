#!/bin/bash
#SBATCH --partition=standard --time=05:00:00 --mem=30G --output=qiime2.log

basepath=/scratch/fmorandi/fula_test/Fula-microbiome/02_pipeline_output
cd $basepath

# Import fastqs
qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path ../01_fastqs/manifest_new.tsv \
	--output-path reads.qza \
	--input-format PairedEndFastqManifestPhred33V2
# Summarize fastqs (Read count, quality information)
qiime demux summarize \
	--i-data reads.qza \
	--o-visualization reads_summary.qzv
  
# Make ASVs, representative sequences and denoising stats
qiime dada2 denoise-paired \
	--i-demultiplexed-seqs reads.qza \
	--p-trim-left-f 10 \
	--p-trim-left-r 10 \
	--p-trunc-len-f 0 \
	--p-trunc-len-r 0 \
	--o-representative-sequences rep-seqs-dada2.qza \
	--o-table table-dada2.qza \
	--o-denoising-stats stats-dada2.qza
# Summarize ASV table
qiime feature-table summarize \
	--i-table table-dada2.qza \
	--o-visualization table-dada2.qzv
# Tabulate representative seqs
qiime feature-table tabulate-seqs \
	--i-data rep-seqs-dada2.qza \
	--o-visualization rep-seqs-dada2.qzv
# Summarize denoising stats
qiime metadata tabulate \
	--m-input-file stats-dada2.qza \
	--o-visualization stats-dada2.qzv

# Phylogeny tree
qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences rep-seqs-dada2.qza \
	--o-alignment aligned-rep-seqs.qza \
	--o-masked-alignment aligned-rep-seqs-masked.qza \
	--o-tree tree-unrooted.qza \
	--o-rooted-tree tree-rooted.qza
	
# Download pre-trained classifier
wget -q https://data.qiime2.org/2022.8/common/silva-138-99-nb-classifier.qza
# Assign taxonomical annotations
qiime feature-classifier classify-sklearn \
	--i-classifier silva-138-99-nb-classifier.qza \
	--i-reads rep-seqs-dada2.qza \
	--o-classification taxonomy.qza \
	--p-reads-per-batch 1000
qiime metadata tabulate \
	--m-input-file taxonomy.qza \
	--o-visualization taxonomy.qzv

# Table collapsed at species level
qiime taxa collapse \
	--i-table table-dada2.qza \
	--i-taxonomy taxonomy.qza \
	--o-collapsed-table table-species.qza \
	--p-level 7
qiime feature-table summarize \
	--i-table table-species.qza \
	--o-visualization table-species.qzv
	
# Table collapsed at genus level
qiime taxa collapse \
	--i-table table-dada2.qza \
	--i-taxonomy taxonomy.qza \
	--o-collapsed-table table-genus.qza \
	--p-level 6
qiime feature-table summarize \
	--i-table table-genus.qza \
	--o-visualization table-genus.qzv

# Table collapsed at phylum level
qiime taxa collapse \
	--i-table table-dada2.qza \
	--i-taxonomy taxonomy.qza \
	--o-collapsed-table table-phylum.qza \
	--p-level 2
qiime feature-table summarize \
	--i-table table-phylum.qza \
	--o-visualization table-phylum.qzv
	
# Export abundance tables in biom format
mkdir -p exports
qiime tools export \
	--input-path table-dada2.qza \
	--output-path exports/asv
qiime tools export \
	--input-path table-species.qza \
	--output-path exports/species
qiime tools export \
	--input-path table-genus.qza \
	--output-path exports/genus
qiime tools export \
	--input-path table-phylum.qza \
	--output-path exports/phylum
# Export rooted tree
qiime tools export \
	--input-path tree-rooted.qza \
	--output-path exports
# Export dada2 stats
qiime tools export \
	--input-path stats-dada2.qza \
	--output-path exports

# Move somes files and tidy up
mv exports/asv/feature-table.biom exports/table-asv.biom
mv exports/species/feature-table.biom exports/table-species.biom
mv exports/genus/feature-table.biom exports/table-genus.biom
mv exports/phylum/feature-table.biom exports/table-phylum.biom
rm -r exports/asv
rm -r exports/species
rm -r exports/genus
rm -r exports/phylum

# Export taxonomic assignment of ASVs
qiime tools export \
	--input-path taxonomy.qza \
	--output-path exports

# Add header row to taxonomy
sed -i "1s/.*/#OTUID\ttaxonomy\tconfidence/" exports/taxonomy.tsv  
# Add taxonomy assignment in ASV level abundance table
biom add-metadata \
	-i exports/table-asv.biom \
	-o exports/table-asv.biom \
	--observation-metadata-fp exports/taxonomy.tsv \
	--sc-separated taxonomy