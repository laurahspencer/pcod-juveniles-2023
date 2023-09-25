#!/bin/bash

IN=/share/afsc/pcod-lcwgs-2023/raw-data
OUT=/share/afsc/pcod-lcwgs-2023/concat

SAMPLES=$(ls -1d ${IN}/G* | awk -F "/" '{ print $NF }')

for sample in ${SAMPLES}
do
  echo ${sample}
  zcat ${sample}/${sample}_*_L*_1.fq.gz \
      >> ${OUT}/${sample}_R1.fastq
  zcat ${sample}/${sample}_*_L*_2.fq.gz \
      >> ${OUT}/${sample}_R2.fastq
done
