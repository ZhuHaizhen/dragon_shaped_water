#!/bin/bash

cd /media/dell/data2/zhuhz/Results/dragon_shape_water/20200323/

ln -s /media/dell/data2/zhuhz/Results/dragon_shape_water/20200316_ehbio/seq ./

mkdir -p temp

mkdir -p result

cp /media/dell/data2/zhuhz/Results/dragon_shape_water/20200316_ehbio/result/otus.fa result/

ln -s /media/dell/data2/zhuhz/Results/dragon_shape_water/20200316_ehbio/temp/filtered.fa temp/

## remove non-bacterial sequences

source activate qiime1

time align_seqs.py -i result/otus.fa -o temp/aligned/  # 15m15.673s

grep -c '>' temp/aligned/otus_failures.fasta  #559

grep '>' temp/aligned/otus_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > temp/aligned/otus_failures.id

filter_fasta.py -f result/otus.fa -o result/otus_filtered.fa -s temp/aligned/otus_failures.id -n

grep '>' -c result/otus_filtered.fa  #24515

## make OTU table

usearch11 -usearch_global temp/filtered.fa -db result/otus_filtered.fa -otutabout result/otu_table.txt -biomout result/otu_table.biom -strand plus -id 0.97 -threads 24

## add taxonomy to OTU table

assign_taxonomy.py -i result/otus_filtered.fa -r /media/dell/data3/database/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /media/dell/data3/database/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt -o result/

biom add-metadata -i result/otu_table.biom --observation-metadata-fp result/otus_filtered_tax_assignments.txt -o result/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy

biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy

biom summarize-table -i result/otu_table_tax.biom -o result/otu_table_tax.sum

## OTU table filtration

filter_otus_from_otu_table.py --min_count_fraction 0.00001 -i result/otu_table_tax.biom -o result/otu_table_tax_filtered.biom

biom summarize-table -i result/otu_table_tax_filtered.biom -o result/otu_table_tax_filtered.sum

biom convert -i result/otu_table_tax_filtered.biom -o result/otu_table_tax_filtered.txt --to-tsv --header-key taxonomy

sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_tax_filtered.txt

## after removing mitochondria and Chloroplast with customized R script

biom convert -i result/otu_table_filtered_clean.txt -o result/otu_table_filtered_clean.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata=taxonomy

biom summarize-table -i result/otu_table_filtered_clean.biom -o result/otu_table_filtered_clean.sum

#filter_fasta.py -f temp/filtered.fa -b result/otu_table_filtered_clean.biom -o result/otu_final.fa

## make phylogenetic tree

align_seqs.py -i result/otus_filtered.fa -o temp/pynast_aligned_seqs

filter_alignment.py -i temp/pynast_aligned_seqs/otus_filtered_aligned.fasta -o temp/filtered_alignment

make_phylogeny.py -i temp/filtered_alignment/otus_filtered_aligned_pfiltered.fasta -o result/otus.tree

## rarefaction

single_rarefaction.py -i result/otu_table_filtered_clean.biom -o result/otu_table_rare.biom -d 42367

biom convert -i result/otu_table_rare.biom -o result/otu_table_rare.txt --to-tsv --header-key taxonomy

sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_rare.txt

## sumarize OTU table

summarize_taxa.py -i result/otu_table_tax_filtered.biom -o result/sum_taxa

sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/sum_taxa/*.txt

## alpha diversity

alpha_diversity.py -i result/otu_table_rare.biom -o result/alpha.txt -t result/otus.tree -m shannon,simpson,goods_coverage,chao1,observed_otus,PD_whole_tree

## beta diversity

normalize_table.py -i result/otu_table_rare.biom -o result/otu_table_css.biom -a CSS

biom convert -i result/otu_table_css.biom -o result/otu_table_css.txt --to-tsv --table-type="OTU table" --header-key=taxonomy

sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt

beta_diversity.py -i result/otu_table_css.biom -o result/beta -t result/otus.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac

sed -i 's/^\t//g' result/beta/*
