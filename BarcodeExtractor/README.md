# Barcode Extractor

### About
A set of tools for decoding data obtained from the Mission Bio Tapestri instrument. It designed to allow barcode extraction and correction from single-cell sequencing fastq/bam files in Linux. 

### Installation

Barcode Extractor is available for installation using [installer](https://dl.missionbio.io/tools/Barcode-v2.0.1_10_g3614e85-Linux-x86_64.sh)

- Download above installer, then install in Linux
```
bash Barcode-v2.0.1_10_g3614e85-Linux-x86_64.sh
```
- Activate the conda environment
```
conda activate barcode
```

### Usage

- Get help information
```
tapestri barcode extract -h
```

- Extract barcodes from fastq files
```
tapestri barcode extract --r1 /path-to-R1/JSN-59IIRNA_S2_L003_R1_001.fastq.gz --r2 /path-to-R2/JSN-59IIRNA_S2_L003_R2_001.fastq.gz --output-r1 debarcode.output.R1.fq --output-r2 debarcode.output.R2.fq --barcode-version dna.v2 --output-discards discard.fq --output-barcodes detected_barcodes 
```


Please contact Mission Bio support for any questions relating to provided files
