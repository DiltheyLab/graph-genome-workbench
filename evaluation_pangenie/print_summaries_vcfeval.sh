#!/bin/bash

variants=("large" "midsize" "small")
tools=("graphtyper" "bayestyper" "pangenie")

for variant in "${variants[@]}"
do
    for tool in "${tools[@]}"
    do
        folder="evaluation/NA24385/${variant}/${tool}/vcfeval"
        echo "Folder: $folder"
        file_path="$folder/summary.txt"
        genotyping_file_path="genotyping/NA24385/${tool}/genotyping_${variant}.vcf.gz" 

        if [ -f "$file_path" ]; then
            if [ -f "$genotyping_file_path" ]; then
                line_count=$(bcftools view -H "$genotyping_file_path" | wc -l)
                echo "#lines in $genotyping_file_path: $line_count"
            fi
            cat "$file_path"
        else
            echo "summary.txt not found in $folder"
        fi

        echo 
        echo "-----------------------------"
        echo ""
        
    done
done
