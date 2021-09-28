# vcfMAN
Pipeline consists of scripts (R and bash) for processing Variant Call Format (VCF) files. The purpose of this pipeline is to standardize output and perform downstream plotting for overview of relevant variant metrics and variant distribution. Pipeline takes VCF files related to both Structural Variants (SV) and small variants (variants ≤ 50 bp). The main outputs are tables, BED formatted variant annotations and figures.

## Instructions on how to execute
1. Download repo to local machine with: `wget https://github.com/mattssca/vcfMAN/archive/refs/heads/main.zip`
2. Unpack with: `unzip -a main.zip`
3. Set directory as current directory: `cd vcfMAN-main/`
4. Install dependencies: `sh install_dep.sh`
5. Migrate VCF files (SVs/mall variants) to corresponding directory (pipeline also takes compressed VCF files in .gzip format)
6. Execute master script with: `sh 01_vcf_man.sh`
7. Output files (figures, tables, summaries and reports) are generated and saved to corresponding output folder (out/SVs or out/small_variants)

## Flowchart
Overview of associated processes and workflow described in vcfMAN. Program prompts user to specify if input VCF are SVs, small variants or if both types of VCFs are present. Valid commands are; `SVs`, `small_variants` and `both`.

![flowchart](https://github.com/mattssca/vcfMAN/blob/main/example_figures/flowchart.png)

## Example Output
Here goes example figures with explanations and comments on expected output, with references for interpretations etc.

### Structural Variants
#### SV Size Distribution
![sv_size_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_sv_size_violin.png)

#### SV Distribution per Chromosome (stacked)
![sv_chrdist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_sv_chrdist.png)

#### Binned SV Sizes
![sv_binned](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_sv_binned.png)

#### Circos Plot
![sv_circos](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_circos.png)

#### Tables and Summaries

### Small Variants
#### Small variants Size Distribution
![small_variants_size_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_size_violin.png)

#### SNV Distribution per Chromosome
![small_variants_chrdist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_snv_chr_dist.png)

#### SNV Distance
![small_variants_snv_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_snv_distance.png)

#### Loss of Heterozygosity
![]()

#### Tables and Summaries

## Dependencies
Pipeline is designed to work on MacOSX systems. Disclaimer, pipeline has not been tested on either Windows or Linux systems.

| Package | Enviroment | Version |
| ------- | ---------- | ------- |
| Brew | MacOSX | 3.2.0 |
| wget | MacOSX | 1.21.2 |
| imagemagick | C | 7.1.0 |
| PhantomJS | C | 2.1.1 |
| Webshot | R | 0.5.1 |
| stringr | R | 1.4.0 |
| table1 | R  | 1.4.2 |
| dplyr | R | 2.1.1 |
| knitr | R | 1.3.4 |
| devtools | R | 2.4.2 |
| gridExtra | R | 2.3 |
| ggthemr | R | 1.1.0 |
| BiocManager | R | 1.30.16 |
| karyoploteR | R | 1.18.0 |
| openxlsx | R | 4.2.4 |
| RCircos | R | 1.2.1 |