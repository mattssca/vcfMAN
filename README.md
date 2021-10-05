# vcfMAN
Pipeline consists of scripts (R and bash) for processing Variant Call Format (VCF) files. The purpose of this pipeline is to standardize output and perform downstream plotting for overview of relevant variant metrics and variant distribution. Pipeline takes VCF files related to both Structural Variants (SV) and small variants (variants ≤ 50 bp). The main outputs are tables, BED formatted variant annotations and figures.

## Instructions on how to execute
1. Download repo to local machine with: `wget https://github.com/mattssca/vcfMAN/archive/refs/heads/main.zip`
2. Unpack with: `unzip -a main.zip`
3. Set directory as current directory: `cd vcfMAN-main/`
4. Install dependencies: `sh install_dep.sh`
5. Migrate VCF files (SVs/mall variants) to corresponding directory (pipeline also takes compressed VCF files in .gzip format)
6. Execute master script with: `sh vcf_man.sh` follwed by command line argument for input VCF file. Valid commands are; `SVs`, `small_variants` and `both`
7. Example for running pipeline on small variant VCF `sh vcf_man.sh small_variants`
8. Output files (figures, tables, summaries and reports) are generated and saved to corresponding output folder (out/SVs or out/small_variants)

## Flowchart
Overview of associated processes and workflow described in vcfMAN.

1. vcfMAN.sh ([vcf_man.sh](https://github.com/mattssca/vcfMAN/blob/main/vcf_man.sh)) acts as a master script and calls appropiate scripts based on the user input (SVs, small variants or both)
2. 01-gunzip.sh ([unpacksmallvariants.sh](https://github.com/mattssca/vcfMAN/blob/main/scripts/unpacksmallvariants.sh) and [unpackSVs.sh](https://github.com/mattssca/vcfMAN/blob/main/scripts/unpackSVs.sh)) checks if input VCFs are compressed with Gzip, if so, VCF files are extracted from compressed format.
3. 02-read_vcf.R ([read_vcf_smallvariants.R](https://github.com/mattssca/vcfMAN/blob/main/scripts/read_vcf_smallvariants.R) and [read_vcf_structuralvaiants.R](https://github.com/mattssca/vcfMAN/blob/main/scripts/read_vcf_structuralvariants.R)) initially performs data wrangling associated tasks in order to extract relevant information from input files. This is done a few different steps. The steps are;

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
   
04. 03-plot.R ([plot_smallvariants.R](https://github.com/mattssca/vcfMAN/blob/main/scripts/plot_smallvariants.R) and [plot_structuralvariants.R](https://github.com/mattssca/vcfMAN/blob/main/scripts/plot_structuralvariants.R)) main plotting scripts utilizing R base functions as well as thirdparty R packages (listed in dependencies) to generate all associated plots.
05. 04-img.sh ([img_man_smallvariants.sh](https://github.com/mattssca/vcfMAN/blob/main/scripts/img_man_smallvariants.sh), [img_man_structuralvariants.sh](https://github.com/mattssca/vcfMAN/blob/main/scripts/img_man_structuralvariants.sh) and [img_man_combine.sh](https://github.com/mattssca/vcfMAN/blob/main/scripts/img_man_combine.sh)) are called to format variant reports in pdf. Scripts are combining tables and figures to a standardized report that can be used for interrogating call-set quality and varaint distributions. Individual plots are also avaialble in out/fig folders.
  
![flowchart](https://github.com/mattssca/vcfMAN/blob/main/example_figures/flowchart.png)

## Example Output
Brief overview and comments on output figures and tables.

### Structural Variants
#### SV Size Distribution
Violin plot visualizing size distributions of SVs (deletions, duplications and insertions). Variant sizes in log10 scale on y-axis and sub types of SVs (deletions, duplications and insertions) on x-axis. Black dot annotates mean variant size. 
![sv_size_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_sv_size_violin.png)

#### SV Distribution per Chromosome (stacked)
Stacked histogram depicting variant distributions in chromosome-dependent manner. Y-axis shows number of variants (n) and chromosomes are arranged on the x-axis. Note, certain chromosomes (e.g chr2, chr16 and chr19) typically shows higher fractions of SVs compared to other chromosomes. This is typically related to complex and difficult to map regions (segmental duplications, homypolymeres and long repetitive sequences) that are shown to be enriched for SVs (Chaisson, Mark J P et al. “Multi-platform discovery of haplotype-resolved structural variation in human genomes.” Nature communications vol. 10,1 1784. 16 Apr. 2019, doi:10.1038/s41467-018-08148-z)
![sv_chrdist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_sv_chrdist.png)

#### Binned SV Sizes
Histogram with set bin sizes showing size distributions of SVs. Bin sizes can easily be configured to match specific desires (lines 84 and 87 of [plot_structuralvariants.R](https://github.com/mattssca/vcfMAN/blob/main/scripts/plot_structuralvariants.R)). y-axis depicts the actual number (n) of variants residing in each bin and bins are shown along the x-axis.
![sv_binned](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_sv_binned.png)

#### Circos Plot
Ideogram with cytogentic bands visualized as a circos plot. Each variant sub type is shown as its own individual track (colors follow the same pattern as for violin and chromosome distribution plots). User can also supply additional BED-tracks to visualize any genomic feature of interest (e.g specific genomic regions, medically relevant genes, etc). In order to use this feature, the user needs to specify additional BED-tracks on line 179 and 215 of [plot_structuralvariants.R](https://github.com/mattssca/vcfMAN/blob/main/scripts/plot_structuralvariants.R). For additional functions and features please see [Rciorcos documentation](https://cran.r-project.org/web/packages/RCircos/vignettes/Using_RCircos.pdf)
![sv_circos](https://github.com/mattssca/vcfMAN/blob/main/example_figures/SV/example_SV_figure_circos.png)

#### Tables and Summaries
Summaries are exported as png files and located in out/SVs/fig. In addition to summary figures, tables are also generated. These include: 

   * VCF header (previously exported).
   * non-structural variants - text file annotating all variants < 50 bp.
   * Spreadsheet (.xlsx) with different pages for each metric and subset
   * BED formatted txt file with **chr | start | end | lenght | sv type | sv family | genotype**

### Small Variants
#### Small variants Size Distribution
Violin plot visualizing size distributions of small variants (i.g deletions and duplications) with SNVs excluded. Variant sizes in log10 scale on y-axis and small variant sub type on x-axis. Black dot annotates mean variant size. 
![small_variants_size_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_size_violin.png)

#### SNV Distribution per Chromosome
Histogram showing the number of SNVs sorted on chromosome with chromosomes on the y-axis and number of SNVs on the x-axis.
![small_variants_chrdist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_snv_chr_dist.png)

#### SNV Distance
Plot showing the distribution of distances between neighboring SNVs with SNV distances above 3rd quantile excluded, to compensate for variants with exceptionally long distance between each other (e.g variants on opposite sides of the same centromere). Mean SNV distance for each chromosome is shown with black line inside each box. Metric can be used to understand the breadth (coverage) of called SNVs compared to an expected output (SNV occurs on average every 1000 nucleotide). Plot also shows if any chromosome shows an increase or decrease of SNVs compared to other chromosmes (i.e even distribution of variants).
![small_variants_snv_dist](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_snv_distance.png)¨

#### SNV Ideogram
Horizontally alligned ideogram highlighting SNVs in a genomic context. Plot makes use of chromosome lenghts located in dep/ folder. Currently all coordinates are in respect to grch38 and regions excluded (centromeres) are also in refernce to the same build. Tables can be customized to accomodate for other versions of the reference genome, as well as blacklsited regions can also be added to further exclude specific genomic regions.
![snv_ideogram](https://github.com/mattssca/vcfMAN/blob/main/example_figures/small_variants/example_small_variants_figure_ideogram.png)

#### Tables and Summaries
Summaries are exported as png files and located in out/amll_variants/fig. In addition to summary figures, tables are also generated. These include: 

   * VCF header (previously exported).
   * non-small variants variants - text file annotating all variants > 50 bp.
   * non-hardcoded genotypes (i.e 1|2, 2|1) - text file
   * Spreadsheet (.xlsx) with different pages for each metric and subset
   * BED formatted txt file with **chr | start | end | lenght | sv type | genotype**

## Dependencies
Pipeline is designed to work on MacOSX systems. Disclaimer, pipeline has not been tested on either Windows or Linux systems.
In order to install all dependencies, execute [install.dep.sh](https://github.com/mattssca/vcfMAN/blob/main/install_dep.sh)

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
