#!/bin/bash

OUTPUT="/scratch/dbihnam/lsc585"

while IFS= read -r sra_id; do
echo "Downloading and converting $sra_id..."

fasterq-dump $sra_id --outdir $OUTPUT

done < sra_ids1.txt

