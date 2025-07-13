data="/ibex/project/c2210/powellor/kfp_paper/monocot_presence_absence/data"

conda activate ncbi

datasets download genome accession GCA_020804685.1 \
	--include genome \
	--filename $data/s_angustfolia.zip

wait -n

unzip $data/s_angustfolia.zip -d $data/s_angustfolia