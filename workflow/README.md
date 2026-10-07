# Read processing

Raw reads: SRA BioProject PRJNA1134421 (24 libraries, paired-end Illumina NovaSeq, Novogene).
`samples.tsv` maps each sample to its group (Ki = Kisumu, Ng = Ngousso, MC/MS = coluzzii control/survivor,
SC/SS = gambiae control/survivor).

Reads were processed with the RNA-Seq-Pop Snakemake workflow (Nagi et al. 2023,
https://github.com/sanjaynagi/rna-seq-pop): fastp, Kallisto (AgamP4.12), DESeq2, HISAT2 and FreeBayes
variant calling, AIM and karyotype analyses, gene-set enrichment. Its outputs are the inputs to `analysis/rnaseq/`.

**Missing from this repository:** the `config.yaml` and workflow version used for the run (the run directory,
`hp_scratch/rna-seq-bouake` on the LSTM cluster, is not currently readable). Add them here before release.
