# Readme

The example includes three individuals, each with multiple chromosome arms. Each individual can be separated into reconstructions for chromosome 2L, 2R, 3L, and 3R (we'll exclude X for these male individuals, because they're haploid).

The files `#.counts` include the expected file format for input. As-is, it's currently counts across all chromosomes, even though you'll only be working with one at a time. This means you'll either need to separate the file using bash commands, or to parse one chromosome at a time in Julia.

The file `founders.txt` lists the total set of potential parents any individual could have., while the files `#.mlp` include the expected file format for reducing search space by ranking most likely parents (mlp). N represents the 'score' where higher is better. Though it ranks all founders, we may only take the top N per chromosome (like top 12 or top 16).

The file `population.haps` includes the true path used to generate the simulated individuals. All chromosomes for all three individuals are included in this single file.

# Retrieving haplotype data:

`wget -L https://virginia.box.com/shared/static/4r7uqcrk00dlinsc6134yte22vnjvy53.gz -O haplotypes.polarized.vcf.gz`

# Note to self on generating the data:

* `cut -f 1-9 haplotypes.polarized.vcf | sed 's@$@\t0/0@g' > dm3.intervals.vcf`
* Ensure `dm3.intervals.vcf` includes header
* extract `.bam` files from `.zip` archive
* index `.bam` files with `samtools index`
* `module load gatk
* `gatk CollectAllelicCounts -I 1.bam -R ./IndoorData/data_external/all_dmel.fasta -O counts.tmp -L ../dm3.intervals.vcf`


