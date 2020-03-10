# E052_exome_pipeline
A whole exome somatic variant calling pipeline.

This project was a collaboration between the labs of Trevor Dale (Cardiff University) and Wales Gene Park (c.2018-2020). In this project, there were 10 organoids derived from primary colorectal cancers, matched with 10 germline blood samples. A somatic variant calling pipeline was used to generate call sets for the 10 somatic (organoid) samples. A custom perl script was used to filter the resulting call sets for ~1605 genes of interest.

The pipeline employed commonly used software, including BWA, SAMtools, GATK and VEP. The perl script 'exome_pipeline_raven.pl' is the wrapper script used to parse samples to these software. This is a generic script used for multiple projects and so some commands may have been commented out, but can be uncommented to recreate the workflow for this project.

The perl script 'gene_extractor.pl' was used to extract variants in genes of interest. The gene_list.CSV file contains the ~1600 genes of interest.
