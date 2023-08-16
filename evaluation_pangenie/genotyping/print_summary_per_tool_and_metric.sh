#!/bin/bash

metric="concordance"
#choose among concordance, precision-recall-typable and fscore

mode="graph"
#choose among external or graph, i.e. leave-one-out or external evaluation

region="nonrep-simple"
#choose variants from external, repeats-simple, repeats-complex, nonrep-simple, nonrep-complex

variant="indel"
#choose variants from indel, sv, large, large-deletion, large-insertion

#Note that not every combination exist, i.e. for region=external, only sv, indel, large-deletion and large-insertion are available

tools=("pangenie" "pangenie" "bayestyper" "graphtyper")
qualities=("qual_0" "qual_200" "qual_0" "qual_0")

# Iterate through each index of the tools array
for index in "${!tools[@]}"
do
    tool="${tools[index]}"
    quality="${qualities[index]}"

    echo "---------------------------------------------------------------"
    echo "metric=${metric} mode=${mode} tool=${tool} region=${region} variant=${variant} quality=${quality} "
    echo "---------------------------------------------------------------"

    echo ""


    # Construct the folder name based on the tool and quality
    folder_name="results/NA24385/evaluation/${metric}/${mode}/${tool}-genotyping-biallelic/coverage-full_${region}_${variant}/${quality}"

    # Check if the folder exists
    if [ -d "$folder_name" ]; then
        # Print the contents of the summary.txt file
        cat "$folder_name/summary.txt"
    else
        echo "Folder not found: $folder_name"
    fi

    echo ""
done
