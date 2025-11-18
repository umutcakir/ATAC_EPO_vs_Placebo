#!/bin/bash
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 2-00:00:00
#SBATCH --mem=120G

export PATH=/scratch1/users/umut.cakir/scATAC_analysis/software/cellranger-atac-2.1.0/bin:$PATH
raw_path="/scratch1/users/umut.cakir/scATAC_analysis/raw_data"

for folder in $raw_path/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        echo "$folder"
        cellranger-atac count --id="$folder_name" \
                        --reference=/scratch1/users/umut.cakir/scATAC_analysis/software/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                        --fastqs=$raw_path/$folder_name \
                        --sample="$folder_name" \
                        --localcores=24 \
                        --localmem=96
    fi
done

