#Fungal analysis
conda activate qiime2-2022.2
mkdir ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED
cd ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED

#Determine the Phred encoding
#https://wiki.bits.vib.be/index.php/Identify_the_Phred_scale_of_quality_scores_used_in_fastQ
#https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py
awk 'NR % 4 == 0' ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/01_Raw_Data/Fungi/NEM_R1_fungi.fq | ~/PostDoc/CZECH_POSTDOC/04_Software/guess-encoding.py -n 10000
#They are Illumina 1.8 (Phred+33)

#Import data into QIIME2 object with manifest
qiime tools import \
	--type 'SampleData[SequencesWithQuality]' \
	--input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/01_Raw_Data/Fungi/NEM_split_by_sample/fungi_manifest \
	--input-format 'SingleEndFastqManifestPhred33V2' \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/demux.qza

#Created a summary of the multiplexed files to check for quality
qiime demux summarize --i-data ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/demux.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/demux.qzv
	
#Some Dada2 variables adjusted according to recommendations in: 
#Rolling, T., Zhai, B., Frame, J., Hohl, T. M., and Taur, Y. (2022). Customization of a DADA2-based pipeline for fungal internal transcribed spacer 1 (ITS1) amplicon data sets. JCI Insight 7, e151663. doi:10.1172/jci.insight.151663.

#The high trimming showed the highest number of usable reads in sample-wise processing and will be used here too the rest of the analysis
#Try to keep quality at the lower end of the box >=20
qiime dada2 denoise-single --i-demultiplexed-seqs ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/demux.qza \
	--p-n-threads 7 \
	--p-max-ee 8 \
	--p-trunc-q 8 \
	--p-trunc-len 0 \
	--p-trim-left 20\
	--p-pooling-method "pseudo" \
	--o-table ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/table_high_trim.qza \
	--o-representative-sequences ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rep-seqs_high_trim.qza \
	--o-denoising-stats ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/denoising-stats_high_trim.qza

#Made a visualisation of the stats
qiime metadata tabulate --m-input-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/denoising-stats_high_trim.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/denoising-stats_high_trim.qzv

#Made a Feature Table summary
qiime feature-table summarize --i-table ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/table_high_trim.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/table.qzv \
	--m-sample-metadata-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/06_Miscellaneous/Metadata.tsv

#And feature data summary
qiime feature-table tabulate-seqs --i-data ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rep-seqs_high_trim.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rep-seqs.qzv

# #Created a tree for phylogenetic analysis
# qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rep-seqs_high_trim.qza \
# 	--o-alignment ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/aligned-rep-seqs.qza \
# 	--o-masked-alignment ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/masked-aligned-rep-seqs.qza \
# 	--o-tree ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/unrooted-tree.qza \
# 	--o-rooted-tree ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rooted-tree.qza

#Created rarefaction curves to show how well we had got all the diversity
	#--i-phylogeny ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rooted-tree.qza \
qiime diversity alpha-rarefaction --i-table ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/table_high_trim.qza \
	--p-max-depth 50000 \
	--m-metadata-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/06_Miscellaneous/Metadata.tsv \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/alpha-rarefaction.qzv

# #Extract the Dada2 ASVs to annotate with Kraken
# mkdir dada2ASVs
# qiime tools extract --input-path rep-seqs_high_trim.qza --output-path dada2ASVs

#Used my pre-trained ITS classifier to classify samples with the Naive Bayes classifier.
qiime feature-classifier classify-sklearn --i-classifier ~/PostDoc/04_Software/QIIME2-classifier/UNITE_8.3-dynamic-classifier.qza \
	--i-reads ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rep-seqs_high_trim.qza \
	--o-classification ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/taxonomy_high_trim_UNITE_8.3.qza

#To see the classification scores
qiime metadata tabulate --m-input-file ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/taxonomy_high_trim_UNITE_8.3.qza \
	--o-visualization ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/taxonomy_high_trim_UNITE_8.3.qzv

#Put all the stuff we want to export for use in a separate folder and export to there
mkdir export
qiime tools export --input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/taxonomy_high_trim_UNITE_8.3.qza \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/export
# qiime tools export --input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/rooted-tree.qza \
# 	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/export
qiime tools export --input-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/table_high_trim.qza \
	--output-path ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/export
biom convert -i ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/export/feature-table.biom \
	-o ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/export/OTU.tsv --to-tsv
sed -i s/"#OTU"/"OTU"/ ~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED/export/OTU.tsv
# cp dada2ASVs/*/data/dna-sequences.fasta export/


