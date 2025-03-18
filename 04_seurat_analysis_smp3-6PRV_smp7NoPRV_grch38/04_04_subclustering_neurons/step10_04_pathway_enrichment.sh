#!/usr/bin/env bash

# Wrapper to run pseudobulk DE analysis on all cells

data_dir="$HOME/Dropbox (Gladstone)/GB-LZ-1373/results/10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/02_neuron_subclustering"
script_dir="$HOME/Dropbox (Gladstone)/github/GB-LZ-1373/10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF"
cd "$data_dir"

for subdir in ./*; do 
	echo "Working on: $(basename "$subdir")" 
	Rscript "$script_dir"/10_04_pathway_enrichment.R \
	--input "$data_dir"/$subdir \
	--output "$data_dir"/$subdir \
	--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
	--species 'human' \
	--p_val_cutoff 0.01 \
	--fold_change_cutoff 1
done


# END -----------------------------------------------------