# \(More\) Assembly QC

Joe Fass - jnfass@ucdavis.edu  
14 December 2016

## N50 and Related Metrics

test N50; list sizes, add cumulative length column, print at desired cumulative length cutoff

```bash
cat pet_fungus.draft.fa | \
    ./fa.sizes.pl | sort -rn | \
    perl -ne 'chomp; $c+=$_; print "$_\t$c\n"' | \
    perl -ane 'if ($F[1] >= 500000000) {print; exit}'
```

---

## BUSCO

First download a relevant lineage tarball from [Evgeny Zdobnov's lab](http://busco.ezlab.org/).

```bash
tar xzvf fungi_odb9.tar.gz
```

We don't have a module for the most recent BUSCO yet, but the old module sets up the dependencies. So, load the module, but use the newer BUSCO script. Also we need to have a local Augustus config directory, so BUSCO has permission to muck around in it and make changes.

```bash
module load busco/1.1b1
# download augustus-3.0.3 locally, to point at legitimate config directory
wget http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.0.3.tar.gz
tar xzvf augustus-3.0.3.tar.gz
export AUGUSTUS_CONFIG_PATH=/share/biocore/jfass/Bioinformatics-Genome-Assembly-and-Analysis-Workshop-Pac-Bio-and-10x-Genomics-/augustus-3.0.3/config/
# download BUSCO 2.0 locally to use updated script
./busco-master-623f5a65e8467daea4dd1eedf13653a52de25897/BUSCO.py \
    --in pet_fungus.draft.fa --out pet_fungus_BUSCO --lineage \
    /share/biocore/jfass/Bioinformatics-Genome-Assembly-and-Analysis-Workshop-Pac-Bio-and-10x-Genomics-/fungi_odb9 \
    --mode genome -f
```

## Variant viz (?)

Now how to address variants found by falcon\_unzip? Try mauve-contig-reorder on all #F\_### contigs collected together (haplotigs) versus collected (or collected and reordered) primary contigs. Then progressive-mauve. Alternately separate haplotigs by each primary 'tig they align to (use LAST, or Zev's modified BLASR), to align each separately.





