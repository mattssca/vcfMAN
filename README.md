# vcfMAN
Pipeline consists of scripts (R and bash) for processing Variant Call Format (VCF) files. The purpose of this pipeline is to standardize output and perform downstream plotting for overview of relevant variant metrics and variant distribution. Pipeline takes VCF files related to both Structural Variants (SV) and small variants (variants ≤ 50 bp). The main outputs are tables, BED formatted variant annotations and figures.

## Instructions on how to execute
1. Downlaod repo to local machine with: `wget https://github.com/mattssca/vcfMAN/archive/refs/heads/main.zip`
2. Unpack with: `unzip -a main.zip`
3. Set directory as current directory: `cd vcfMAN-main/`
4. Install dependencies: `sh install_dep.sh`
5. Migrate VCF files (SVs/mall variants) to corresponding directory (pipeline also takes compressed VCF files in .gzip format
6. Execute master script with: `sh 01_vcf_man.sh`
7. Output files (figures, tables, summaries and reports) are generated and saved to corresponding output folder (out/SVs or out/small_variants)

## Flowchart
Overview of associated processes and workflow described in vcfMAN.

![flowchart](https://github.com/mattssca/vcfMAN/blob/main/dep/flowchart.png)
