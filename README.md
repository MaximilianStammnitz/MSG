MSV
=========

### Wrapper script for structural variant calling and genotyping

__Max Stammnitz 
Transmissible Cancer Group, University of Cambridge__

MSV is an open-source wrapper that makes use of the powerful structural variant caller Manta ([![DOI](10.1093/bioinformatics/btv710)](10.1093/bioinformatics/btv710)) for identifying germline and somatic SVs in sequencing data coming from a set of unpaired samples.

It has been designed to offer great sensitivity and specificity, even for low-frequency somatic mutations found in cancer genomes. Although this pipeline can be applied to any cancer sequence data, it is particularly useful for the processing of data obtained from many tumour samples that are closely related to each other and lack matched normal samples, such as transmissible cancer or metastasis data.

Somatypus has been tested on Ubuntu (14.04.4) systems, and it should work well on any Linux distribution. It has not been tested on Mac systems, although it may work after some minor code modifications.
