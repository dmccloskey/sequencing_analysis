# sequencing_analysis
methods for analyzing sequencing data

##Douglas McCloskey

To install, please follow the [instructions](INSTALL.md).

##Library overview:

###genome_annotations.py
methods to annotate a genome based on the orgnism and strand position

methods to extract genomic information from database genomic records

methods to convert gene identifiers for a given orgnanism

###genome_diff.py
methods to parse a breseq .gd file into tables for metadata, mutations, evidence, and validation

methods for annotating and filtering breseq data prior to further analysis

###genome_diff_mutations.py
methods to handle multiple mutation tables for analysis of multiple .gd samples

###mutations_lineage.py
methods to track and visualize genomic changes that occur accross time in a given sample

###mutation_endpoints
methods to compare and contrast multiple samples to identify unique and conserved features between samples

###mutations_heatmap
methods for clustering and visualizing genomic differences between multiple samples

###general_feature_format.py
methods to extract read information on the plus and minus strands from .gff files

methods to perform basic visualization and analyses on the plus and minus strands

###gff_coverage.py
methods to perform basic statistical analyses and visualize the read coverage of a sequenced sample based on the read count of the plus and minus strands

methods to identify, annotate, and visualize genomic amplifications
