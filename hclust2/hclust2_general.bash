#!/bin/bash

### Vars
make_species_table=false

input_file=all_results_merged_flo.tsv
output_file=flo_abundance_heatmap_species.png
flabel=10
slabel=25
dpi=600
image_size=20
ftop=30

### MAIN
if [ $make_species_table ]
then
    grep -E "s__|SampleID" $input_file | sed 's/^.*s__//g'\ > 'merged'${input_file}
fi

hclust2.py -i 'merged'${input_file} \
    -o $output_file \
    --out_table ${output_file}.txt \
    --legend_file ${output_file}_legend.png \
    --f_dist_f braycurtis \
    --s_dist_f braycurtis \
    --cell_aspect_ratio 0.5 \
    -l \
    --flabel_size $flabel \
    --slabel_size $slabel \
    --max_flabel_len 100 \
    --max_slabel_len 100 \
    --minv 0.1 \
    --dpi $dpi \
    --image_size $image_size \
    --ftop $ftop
