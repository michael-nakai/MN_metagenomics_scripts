#!/bin/bash

# For use on my desktop

### Vars here
medium='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/2-trimmomatic/outputs/'

flopaired='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/2-trimmomatic/outputs/paired/flo/'
hamdipaired='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/2-trimmomatic/outputs/paired/hamdi/'

floCatFiles='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/HUMAnN3/outputs/catfiles/flo/'
hamdiCatFiles='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/HUMAnN3/outputs/hamdi/'

flooutFolder='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/HUMAnN3/outputs/flo/'
hamdioutFolder='/home/michael/Data_Storage/Bioinformatics/Metagenome_Seq/2021-05-27/HUMAnN3/outputs/hamdi/'

uniref_database='/home/michael/Data_Storage/Bioinformatics_Tools/shared_databases/uniref'
chocophlan_database='/home/michael/Data_Storage/Bioinformatics_Tools/shared_databases/chocophlan'


### MAIN
mkdir -p $flooutFolder $hamdioutFolder $floCatFiles $hamdiCatFiles

flo=$(ls "${flopaired}" | grep "1_paired.fq.gz")
hamdi=$(ls "${hamdipaired}" | grep "1_paired.fq.gz")

# Flo's stuff
for x in ${flo[@]}
do
        y=${x%_*}
        z=${y%_*}
        fw=${flopaired}${z}_1_paired.fq.gz
        rv=${flopaired}${z}_2_paired.fq.gz
        catname=${z}"_cat.fq.gz"

        # Cat the paired end files together (https://github.com/biobakery/humann#humann-30-and-paired-end-sequencing-data)
        echo 'Catting' ${z}
        cat $fw $rv > ${floCatFiles}${catname}

        #Run HumanN3 on the catfile
        humann --input ${floCatFiles}${catname} --output ${flooutFolder}${z} --protein-database $uniref_database --nucleotide-database $chocophlan_database
done

# Hamdi's stuff
for x in ${hamdi[@]}
do
        y=${x%_*}
        z=${y%_*}
        fw=${hamdipaired}${z}_1_paired.fq.gz
        rv=${hamdipaired}${z}_2_paired.fq.gz
        catname=${z}"_cat.fq.gz"

        # Cat the paired end files together (https://github.com/biobakery/humann#humann-30-and-paired-end-sequencing-data)
        echo 'Catting' ${z}
        cat $fw $rv > ${hamdiCatFiles}${catname}

        #Run HumanN3 on the catfile
        humann --input ${hamdiCatFiles}${catname} --output ${hamdioutFolder}${z} --protein-database $uniref_database --nucleotide-database $chocophlan_database
done
