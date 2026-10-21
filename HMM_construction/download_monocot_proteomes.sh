#!/bin/bash -e

datasets download genome \
    taxon 4447 \
    --dehydrated \
    --reference \
    --include protein \
    --filename monocot_proteomes.zip

unzip monocot_proteomes.zip -d monocot_proteomes

datasets rehydrate --directory monocot_proteomes

cat monocot_proteomes/ncbi_dataset/data/*/*.faa > monocot_proteome.faa