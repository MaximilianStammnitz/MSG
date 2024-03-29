MSG - Wrapper for SV calling and genotyping
===========================================

![MSG](./MSG.png)

MSG is an open-source wrapper that makes use of the powerful structural variant caller [**M**anta](https://github.com/Illumina/manta) ([Chen et al., 2016](https://academic.oup.com/bioinformatics/article/32/8/1220/1743909)) for identifying germline and somatic SVs in sequencing data of unpaired samples, as well as the [**S**vimmer](https://github.com/DecodeGenetics/svimmer) and [**G**raphTyper](https://github.com/DecodeGenetics/graphtyper) frameworks for SV merging and genotyping, by use of pangenome graphs ([Eggertson et al., 2017](https://www.nature.com/articles/ng.3964); [Eggertson et al., 2019](https://www.nature.com/articles/s41467-019-13341-9)).

This workflow:
* facilitates high-resolution SV discovery, reassembly, genotyping and variant allele frequency inference from Illumina PE reads (> 30X)
* is suitable for the genomic studies of asexually evolving genomes, such as the somatic tracing of cancer and normal tissue clones
* has been benchmarked for >512 Mb chromosomes, via CSI indexing support
* follows the recommendations by [Eggertson et al., 2019](https://www.nature.com/articles/s41467-019-13341-9)
* has been recently employed to study structural variation across hundreds of humans ([Almarri et al., 2020](https://doi.org/10.1016/j.cell.2020.05.024))

In the Transmissible Cancer Group, we rely on MSG because of our observation of comparably low false-negative genotyping rates with respect to input SV sets across multiple clonally related samples. In addition, we find that GraphTyper's SV variant allele frequency estimates are exceptionally reliable, in the sense that they match well with tumour purity/cellularity estimates derived from orthogonal SNV analyses.

MSG is a wrapper script written in GNU Bash, and has been tested extensively on Ubuntu (18.04) systems. It should work well on other Linux distributions, and has a checkpoint flag implementation for resumptions after steps (1) and (2). Moreover, thanks to both Manta's and GraphTyper's parallelised implementation, runtimes can be substantially reduced by the specification of multiple computing nodes (if available).


---

## Dependencies

The following tools/scripts need to be installed and placed in your $PATH environment, in addition to Python v2.6+, tabix, vcf-sort, vcf-concat, bcftools and bgzip:
* [Manta v1.6.0 (dev branch, contains convertIntervion.py)](https://github.com/MaximilianStammnitz/MSG/blob/master/manta-1.6.0.tar.bz2)
* [svimmer v0.1](https://github.com/DecodeGenetics/svimmer/releases/tag/v0.1)
* [Graphtyper v2.6.1 (dev branch)](https://github.com/MaximilianStammnitz/MSG/blob/master/graphtyper2.6.1-dev)


---

## Running MSG

A run of MSG_v1.3.sh requires the following inputs:
* list with absolute paths to indexed BAMs for SV discovery (-d)
* list with absolute paths to indexed BAMs for SV genotyping (-g)
* reference genome FASTA (-r)
* path to the global output folder (-o)
* region window file (-w)
* number of CPUs (-c)


---

## Contact

This script was written by Max Stammnitz, Transmissible Cancer Group (2020), with UNIX steals from [Adrian Baez-Ortega](https://github.com/baezortega), Wellcome Sanger Institute and svimmer/GraphTyper optimisation support from [Hannes Eggertson](https://github.com/hannespetur), deCODE Genetics.

Please get in touch if there are any issues: maxrupsta@gmail.com


---

## Citation

If you can make good use of these functions in your work, I would be grateful for your citation of our associated preprint: **The evolution of two transmissible cancers in Tasmanian devils ([Stammnitz _et al._ 2023, Science 380:6642](https://doi.org/10.1126/science.abq6453))**
