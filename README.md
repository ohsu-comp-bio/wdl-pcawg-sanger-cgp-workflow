# ICGC-TCGA-PanCancer BWA Workflow

## Overview

This is the workflow derived from Sanger's Cancer Genome Project core somatic 
calling pipeline for the TCGA/ICGC PanCancer Analysis of Whole Genomes (PCAWG) project.

For more information about the project overall see the
[PanCancer wiki space](https://wiki.oicr.on.ca/display/PANCANCER/PANCANCER+Home).

More detailed documentation about the production use of this workflow can be
found in the [PanCancer-Info](https://github.com/ICGC-TCGA-PanCancer/pancancer-info)
project where we maintain our production documentation and SOPs.

## Building the Worklfow Docker Image

You can also build a Docker image that has the workflow ready to run in it.

    docker build -t pcawg-cgp-somatic-workflow


## Running the Workflow with Cromwell

    java -jar cromwell-19.3.jar run pcawg-cgp-somatic-workflow.wdl pcawg-cgp-somatic-workflow.json


## Sample Data

Some sample data.

* https://s3-eu-west-1.amazonaws.com/wtsi-pancancer/testdata/HCC1143_ds.tar

## Reference Data

* https://s3-eu-west-1.amazonaws.com/wtsi-pancancer/reference/GRCh37d5_CGP_refBundle.tar.gz
* https://s3-eu-west-1.amazonaws.com/wtsi-pancancer/reference/GRCh37d5_battenberg.tar.gz
