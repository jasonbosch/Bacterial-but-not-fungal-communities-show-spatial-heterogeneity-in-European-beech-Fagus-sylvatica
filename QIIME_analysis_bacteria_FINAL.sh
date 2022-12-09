#Bacterial analysis
conda activate qiime2-2022.2

#Bacterial analysis
conda activate qiime2-2022.2
mkdir ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED
cd ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED

#Determine the Phred encoding
#https://wiki.bits.vib.be/index.php/Identify_the_Phred_scale_of_quality_scores_used_in_fastQ
#https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py
awk 'NR % 4 == 0' ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/01_Raw_Data/Bacteria/EMA_R1.fq | ~/PostDoc/04_Software/guess-encoding.py -n 10000
#They are Illumina 1.8 (Phred+33)

#Import data into QIIME2 object with manifest
qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/01_Raw_Data/Bacteria/EMA_split_by_sample/bacteria_manifest \
	--input-format 'PairedEndFastqManifestPhred33V2' \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/demux-paired-end.qza

#Created a summary of the multiplexed files to check for quality
qiime demux summarize --i-data ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/demux-paired-end.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/demux.qzv

#The high trimming showed the highest number of usable reads in sample-wise processing and will be used here too the rest of the analysis
#Try to keep quality at the lower end of the box >=20
qiime dada2 denoise-paired --i-demultiplexed-seqs ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/demux-paired-end.qza \
	--p-n-threads 7 \
	--p-trim-left-f 20 --p-trim-left-r 20 \
	--p-trunc-len-f 235 --p-trunc-len-r 0 \
	--p-pooling-method "pseudo" \
	--o-table ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/table_high_trim.qza \
	--o-representative-sequences ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rep-seqs_high_trim.qza \
	--o-denoising-stats ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/denoising-stats_high_trim.qza

#Made a visualisation of the stats
qiime metadata tabulate --m-input-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/denoising-stats_high_trim.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/denoising-stats_high_trim.qzv

#Made a Feature Table summary
qiime feature-table summarize --i-table ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/table_high_trim.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/table.qzv \
	--m-sample-metadata-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/06_Miscellaneous/Metadata.tsv

#And feature data summary
qiime feature-table tabulate-seqs --i-data ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rep-seqs_high_trim.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rep-seqs.qzv

#Created a tree for phylogenetic analysis
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rep-seqs_high_trim.qza \
 	--o-alignment ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/aligned-rep-seqs.qza \
 	--o-masked-alignment ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/masked-aligned-rep-seqs.qza \
 	--o-tree ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/unrooted-tree.qza \
 	--o-rooted-tree ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rooted-tree.qza

#Created rarefaction curves to show how well we had got all the diversity
qiime diversity alpha-rarefaction --i-table ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/table_high_trim.qza \
	--i-phylogeny ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rooted-tree.qza \
	--p-max-depth 50000 \
	--m-metadata-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/06_Miscellaneous/Metadata.tsv \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/alpha-rarefaction.qzv

# #Extract the Dada2 ASVs to annotate with Kraken
# mkdir dada2ASVs
# qiime tools extract --input-path rep-seqs_high_trim.qza --output-path dada2ASVs

#Had a pre-trained classifier which was used to classify samples with the Naive Bayes classifier.
qiime feature-classifier classify-sklearn --i-classifier ~/PostDoc/04_Software/QIIME2-classifier/silva-138.1-ssu-nr99-515F_806R-classifier.qza \
	--i-reads ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rep-seqs_high_trim.qza \
	--o-classification ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/taxonomy_high_trim_silva_138.1.qza
	
#To see the classification scores
qiime metadata tabulate --m-input-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/taxonomy_high_trim_silva_138.1.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/taxonomy_high_trim_silva_138.1.qzv

#Put all the stuff we want to export for use in a separate folder and export to there
mkdir export
qiime tools export --input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/taxonomy_high_trim_silva_138.1.qza \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export
qiime tools export --input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/rooted-tree.qza \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export
qiime tools export --input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/table_high_trim.qza \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export
biom convert -i ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export/feature-table.biom \
	-o ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export/OTU.tsv --to-tsv
sed -i s/"#OTU"/"OTU"/ ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export/OTU.tsv
