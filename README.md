# Disentangling the genetic basis of rhizosphere microbiome assembly in tomato

This is a repository of code used by Oyserman, et al. "Disentangling the genetic basis of rhizosphere microbiome assembly in tomato." Nature communications (2022).

The code has been organized below based on the sections outline in the manuscript. Numbers have been added for convenience.

## Organization of codes
* 1.1 Mapping files, MAG Bins, & Contigs greater than 10 KB
* 2.1 Baseline analyses of the tomato Recombinant Inbred Line population
* 2.2 Taxonomic microbiome features as quantitative traits
    * 1. DADA2 rRNA amplicon sequence processing
    * 2. ASV statistics
    * 3. 16S QTL analysis
* 2.3 Functional microbiome features as quantitative traits 
    * 1. Metagenomics analysis
         * Assembly
         * Mapping
    * 2. Binning of metagenomic contigs
    * 3 Making phenotype files based on contig depth
    * 4 Feature selection
    * 5 Taxonomic and functional annotation of the metagenome
         * Functional annotation
         * Taxonomic annotation
    * 6 Single Nucleotide Variant analysis
         * Bin QTL
         * Contig QTL
    * 7 Quantitative Trait Locus Analysis
* 2.4 Amplicon-based bulk segregant analysis of Streptomyces and Cellvibrio abundance 
    * 1. Independent validation of QTLs through bulk segregant analysis
    * 2. Bulk Segregation Analysis
* 2.5 Illuminating metagenomic traits in Cellvibrio and Streptomyces
* 2.6 Genomic structure in Cellvibrio and Streptomyces provides insights into adaptations for differential recruitment
    * 1. Single Nucleotide Variant analysis
        * InStrain
        * SNV extraction
        * SNV selection and annotation
    * 2. Quantitative Trait Loci Analysis
