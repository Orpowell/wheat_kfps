data="OUTDIR"

datasets download genome taxon 4447 \
	--dehydrated \
	--reference \
	--include genome \
	--filename $data/monocot_4447_reference_genomes.zip

wait -n

unzip $data/monocot_4447_reference_genomes.zip -d $data/monocot_4447_reference_genomes

wait -n

datasets rehydrate --directory $data/monocot_4447_reference_genomes