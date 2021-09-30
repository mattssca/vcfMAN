# vcfMAN
Pipeline consists of scripts (R and bash) for processing Variant Call Format (VCF) files. The purpose of this pipeline is to standardize output and perform downstream plotting for overview of relevant variant metrics and variant distribution. Pipeline takes VCF files related to both Structural Variants (SV) and small variants (variants ≤ 50 bp). The main outputs are tables, BED formatted variant annotations and figures.

## Instructions on how to execute
1. Download repo to local machine with: `wget https://github.com/mattssca/vcfMAN/archive/refs/heads/main.zip`
2. Unpack with: `unzip -a main.zip`
3. Set directory as current directory: `cd vcfMAN-main/`
4. Install dependencies: `sh install_dep.sh`
5. Migrate VCF files (SVs/mall variants) to corresponding directory (pipeline also takes compressed VCF files in .gzip format)
6. Execute master script with: `sh 01_vcf_man.sh`
7. Software asks user to specify VCF input (SVs and small variants). Valid commands are; `SVs`, `small_variants` and `both`.
8. Output files (figures, tables, summaries and reports) are generated and saved to corresponding output folder (out/SVs or out/small_variants)

## Flowchart
Overview of associated processes and workflow described in vcfMAN.

1. vcf_man.sh acts as a master script and calls appropiate scripts based on the user input (SVs, small variants or both)
2. 01-gunzip.sh (unpacksmallvariants.sh and unpackSVs.sh) checks if input VCFs are compressed with Gzip, if so, VCF files are extracted from compressed format.
3. 02-read_vcf.R (read_vcf_smallvariants.R and read_vcf_structuralvaiants.R) initially performs data wrangling associated tasks in order to extract relevant information from input files. This is done a few different steps. The steps are;

   * List VCf files located in /in folder. Strips the path and saves the sample name as a variable.
   * Read listed VCF files into the R environment, skipping header (header is being extracted and saved in /out folder).
   * Subset input VCF on relevant variables.
   * Format genotype field.
   * Format INFO field, to get variant-type (e.g deletion, duplication, etc.).
   * Rename variable-names (e.g chromosome, start, end, genotype, etc.),
   * Creates two new variables based on n nucleotides present in the alt and ref column.
   * Compute SV length (alt nucleotide count - ref nucleotide count). Deletions are defined as negative sv_lenght values.
   * Negative values are transformed into absolute values.
   * Compute variant end-coordinates (start coordinate + sv_lenght).
   * Genotype information is converted into character-string.
   * End-coordinates for deleterious variants are transformed into start coordinates + 1.
   * SV sub-types are defined (i.g SIMPLEDEL = del, SUBSDEL = del, CONTRAC = del, DUP = dup, SIMPLEINS = ins).
   * Variants not belonging to specified input VCF are subset and removed from downstream analysis (i.e for small variants, variants > 50 bp are removed, and for SV calls variants ≤ 50 bp) Non hard coded genotypes are also subset and removed in a similar approach. Metrics related to removed variants are printed to console and saved as individual .txt files allowing for more in-depth interrogation of such variants.
   * Data frame is sorted and exported as .tsv (main input for plotting).
   * Summary metrics are generated and exported.
   * Summary tables is printed to console and exported as .png.
   
04. 03-plot.R (plot_smallvariants.R and plot_structuralvariants.R) main plotting scripts utilizing R base functions as well as thirdparty R packages (listed in dependencies) to generate all associated plots.
05. 04-img.sh (img_man_smallvariants.sh, img_man_structuralvariants.sh and img_man_combine.sh) are called to format variant reports in pdf. Scripts are combining tables and figures to a standardized report that can be used for interrogating call-set quality and varaint distributions.
  
![flowchart](https://github.com/mattssca/vcfMAN/blob/main/example_figures/flowchart.png)

## Example Output
Brief overview and comments on output figures and tables.

### Structural Variants
#### SV Size Distribution
Violin plot visualizing size distributions of SVs (deletions, duplications and insertions). Variant sizes in log10 scale on y-axis and sub types of SVs (deletions, duplications and insertions) on x-axis. Black dot annotates mean variant size. 
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
![small_variants_snv_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_snv_distance.png)¨

#### SNV Ideogram
![snv_ideogram](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_ideogram.png)

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
| psych | R | 2.1.6 |
| data.table | R | 1.14.0 |
