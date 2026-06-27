#Load the files to qiime
qiime tools import \
      --type 'SampleData[SequencesWithQuality]' \
      --input-path manifest.tsv \
      --output-path joined_sequences_clean.qza \
      --input-format SingleEndFastqManifestPhred33V2

#Visualization
qiime demux summarize \
  --i-data joined_sequences_clean.qza \
  --o-visualization joined_sequences_clean_lote*.qzv

#Denoised
qiime dada2 denoise-single \
--i-demultiplexed-seqs joined_sequences_clean.qza \
--p-trim-left 0 \
--p-trunc-len 440 \
--p-max-ee 5.0 \
--p-n-threads 16 \
--output-dir resultados_dada2_joined_440 \
--verbose

#Visualization
qiime metadata tabulate \
  --m-input-file denoising_stats.qza \
  --o-visualization stats_440_lote5.qza
#Repeat this steps for all the lotes that we want to analyse

#Join the count tables
qiime feature-table merge \
  --i-tables lote1/resultados_dada2_joined_440/table.qza \
  --i-tables lote2/resultados_dada2_joined_440/table.qza \
  --i-tables lote3/resultados_dada2_joined_440/table.qza \
  --i-tables lote4/resultados_dada2_joined_440/table.qza \
  --i-tables lote5/resultados_dada2_joined_440/table.qza \
  --o-merged-table merge_data/table-merged.qza

#Join the representative features
qiime feature-table merge-seqs \
      --i-data lote1/resultados_dada2_joined_440/representative_sequences.qza \
      --i-data lote2/resultados_dada2_joined_440/representative_sequences.qza \
      --i-data lote3/resultados_dada2_joined_440/representative_sequences.qza \
      --i-data lote4/resultados_dada2_joined_440/representative_sequences.qza \
      --i-data lote5/resultados_dada2_joined_440/representative_sequences.qza \
      --o-merged-data merge_data/representative_sequences_merge.qza

#Creation of phylogenetic three
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merge_data/representative_sequences_merge.qza \
  --o-alignment merge_data/phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment merge_data/phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree merge_data/phylogeny/unrooted-tree.qza \
  --o-rooted-tree merge_data/phylogeny/rooted-tree.qza
  --p-n-threads 4

#Classification taxonomy, we use the Silva 138-99 naive bayes classifier
qiime feature-classifier classify-sklearn \
  --i-classifier ../classifier/silva-138-99-nb-classifier.qza \
  --i-reads merge_data/representative_sequences_merge.qza \
  --o-classification merge_data/taxonomy.qza \
  --p-n-jobs 16

#Visualization
qiime taxa barplot \
  --i-table merge_data/table-merged.qza \
  --i-taxonomy merge_data/taxonomy.qza \
  --m-metadata-file ../metadata/sample-metadata.tsv \
  --o-visualization merge_data/taxa-bar-plots.qzv

#Data cleaness, without  mitochondrias and chloroplast
qiime taxa filter-table \
  --i-table merge_data/table-merged.qza \
  --i-taxonomy merge_data/taxonomy.qza\
  --p-exclude mitochondria,chloroplast \
  --p-mode contains \
  --o-filtered-table merge_data/taxonomy_clean.qza

#Discard ASVs with low abundancy
#feauture table
qiime feature-table filter-features \
  --i-table merge_data/taxonomy_clean.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table merge_data/final_table_clean.qza

#Representative sequences
qiime feature-table filter-seqs \
  --i-data representative_sequences_merge.qza \
  --i-table final_table_clean.qza\
  --o-filtered-data representative_sequences_filtered.qza
