BarDecodeRGv26.pl
A Perl script that detects and extracts 49–52 bp barcode structures formatted as:
bc1 – offset – constant1 – bc2 – constant2
The script identifies reads containing this structure and writes the parsed barcode information to an output file.

CollapseBarcodeTags.py
A Python script that processes the extracted 49–52 bp barcode sequences, corrects barcode segments against a predefined whitelist, and collapses them into a final 18 bp barcode consisting of the combined bc1–bc2 sequence.

```
# Barcode decode
time perl ${baseDir}/BarDecodeRGv26.pl --input-1 ${workDir}/${prefix}_step1_R1_trimmed.fastq.gz --input-2 ${workDir}/${prefix}_step1_R2_trimmed.fastq.gz --output-1 ${workDir}/${prefix}_starInput_1.fastq --output-2 ${workDir}/${prefix}_starInput_2.fastq --discards ${tempDir}/${prefix}_bardecode_discards.sam

# Collapse barcodes
python ${baseDir}/CollapseBarcodeTags.py --input ${tempDir}/${prefix}Aligned.sortedByCoord.out.bam --barcodes1 ${baseDir}/bc1_v26.txt --barcodes2 ${baseDir}/bc2_v26.txt --output ${workDir}/${prefix}_collapsed.bam --perl-diff
```
