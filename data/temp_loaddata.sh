#!/bin/bash 

dtpath="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/mutationarearelationship/mar/tmpobjects/"
source_species=('arabidopsis' 'alyrata' 'amaranthus' 'eucalyptus' 'joshua')
sink_species=('Arabidopsis_thaliana' 'Arabidopsis_lyrata' 'Amaranthus_tuberculatus' 'Eucalyptus_melliodora' 'Yucca_brevifolia')

for ((i = 0; i < ${#source_species[@]}; i++)); do
    echo "Source: ${source_species[i]}, Sink: ${sink_species[i]}"
    rsync -ahv ${dtpath}/genome-${source_species[i]}.rda ./geno-${sink_species[i]}.rda
    rsync -ahv ${dtpath}/coords-${source_species[i]}.rda ./coords-${sink_species[i]}.rda
done
