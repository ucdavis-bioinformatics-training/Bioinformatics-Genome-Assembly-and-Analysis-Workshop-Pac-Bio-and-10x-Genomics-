# \(More\) Assembly QC

Joe Fass - jnfass@ucdavis.edu  
14 December 2016

## Setup

Copy my data and scripts into your home somewheres:

```bash
cp -av /share/biocore/workshops/Genome-Assembly-Workshop/examples/assemblyQC.joe/setup \
    ~/workshopStuff/assemblyQCfiles
cd ~/workshopStuff/assemblyQCfiles/
ls -ltrha  # or alias ll='ls -ltrha'; ll
```

---

## N50 and Related (Better!) Metrics

Test *NG50*; list sizes, add cumulative length column, print at desired cumulative length cutoff  

```bash
cat p_ctg.fa | \
    ./fa.sizes.pl | sort -rn | \
    perl -ne 'chomp; $c+=$_; print "$_\t$c\n"' \
    > seqs.lengths  # cumulative len vs len data points
cat p_ctg.fa | \
    ./fa.sizes.pl | sort -rn | \
    perl -ne 'chomp; $c+=$_; print "$_\t$c\n"' | \
    perl -ane 'if ($F[1] >= 30000000) {print; exit}'
gnuplot  # enters gnuplot interactive mode; ctrl-d to quit!
set term dumb  # not an insult, just plots to the terminal instead of an image file
plot "seqs.lengths" with dots
# hit ctrl-d when done looking
```

Pop-quiz: how would you compare the lengths of the alternative contigs to the lengths of the primary contigs? Hint: "replot" is a gnuplot way to plot more information without erasing previoius data points.   

What else do we want to know about our contigs? Base content?

```bash
cat p_ctg.fa | grep -v ">" | grep -o . | \
    perl -ne 'chomp; $h{$_}++ }{ print "$_\t$h{$_}\n" for (keys %h)'
```

Pop-quiz: what about *di-base* content?

---

## BUSCO

Ian Korf's lab published a **C**ore set of (highly conserved) **E**ukaryotic **G**enes, and a tool to look for these in draft genomes ([CEGMA](http://korflab.ucdavis.edu/datasets/cegma/)). Though CEGMA [isn't maintained anymore](http://www.acgt.me/blog/2015/5/18/goodbye-cegma-hello-busco), there's a new tool ... BUSCO. First download a relevant lineage tarball from [Evgeny Zdobnov's lab](http://busco.ezlab.org/) ... actually, nevermind; I've put it in your working directory. You'll also need an appropriate database - in this case a fungal database - which I've also put in your working directory.  

We don't have a module for the most recent BUSCO yet, but the old module sets up the dependencies. So, load the module, but use the newer BUSCO script. Also we need to have a local Augustus (a gene prediction tool) config directory, so BUSCO has permission to muck around in it and make changes. I've put the Augustus tarball in your working directory; see the commented 'wget' command for where I got it.  

```bash
module load busco/1.1b1
tar xzvf fungi_odb9.tar.gz
# download augustus-3.0.3 locally, to point at local config directory
# wget http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.0.3.tar.gz
tar xzvf augustus-3.0.3.tar.gz
export AUGUSTUS_CONFIG_PATH= {~~~~~full/path/to/your/working/directory~~~~~} /augustus-3.0.3/config/
# download BUSCO 2.0 locally to use updated script
tar xzvf busco-master-623f5a65e8467daea4dd1eedf13653a52de25897.tar.gz
# ... and run BUSCO on our FALCON primary contigs!:
./busco-master-623f5a65e8467daea4dd1eedf13653a52de25897/BUSCO.py \
    --in p_ctg.fa --out pet_fungus_BUSCO --lineage \
    /share/biocore/workshops/Genome-Assembly-Workshop/examples/assemblyQC.joe/test/fungi_odb9 \
    --mode genome -f  # beware; this will run for about a half hour
```

---

## Variant Visualization

How can we look at the variants found by falcon\_unzip (found in 4-quiver/cns\_output/)? Mauve, a whole genome aligner, is one possibility. With draft genomes, we'll often reorder contigs to more closely match a trusted reference (see command below). But in this case, with discrete variants, this can create false LCBs that cross contig boundaries and (visually) imply missing sequence between them. Just run progressiveMauve on the consensus haplotigs versus the consensus primotigs.

```bash
# wget http://darlinglab.org/mauve/snapshots/2015/2015-02-13/linux-x64/mauve_linux_snapshot_2015-02-13.tar.gz
tar xzvf mauve_linux_snapshot_2015-02-13.tar.gz
cd mauve_snapshot_2015-02-13/
# java -Xmx5g -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output ../mauve.reorder -ref ../published.fasta -draft ../p_ctg.fa  # example for draft genome versus a published relative
mkdir ../mauve.progMauve
# progressiveMauve runs on this data in about 5 minutes, but takes maybe 15 minutes to fully write and close the output!
./linux-x64/progressiveMauve --output=../mauve.progMauve/aln.xmfa ../cns_p_ctg.fasta ../cns_h_ctg.fasta
cd ..
```

We'll need to pull the alignment files down to a local machine to view them:  

```bash
# for ~~example~~:
cd ~/Documents/
scp -r class8@cabernet.genomecenter.ucdavis.edu:/home/class8/mauve.progMauve .
# and we need the sequence files so Mauve can show nucleotides and sequence boundaries:
cd ..
scp -r class8@cabernet.genomecenter.ucdavis.edu:/home/class8/cns*fasta .
```

Inspect the alignment by launching Mauve, then File --> Open alignment.  

Back on the server, pull out the cumulative lengths, as above, and plot at will. We can pull out the sizes of the alignments ... in this case, they represent the fraction of the genome that contains structural variation, according to FALCON\_unzip. Look in the backbone file. Pull out the sizes of alignment blocks, and cumulative lengths, this way:   

```bash
tail -n +2 mauve.progMauve/aln.xmfa.backbone | \
    perl -ane 'if ($F[2]!=0 and $F[3]!=0) {$b=abs($F[1]-$F[0]); print "$b\n"}' | \
    sort -rn | perl -ne '$b=$_; chomp $b; $c+=$b; \
    print "$b\t$c\n"' > blocks.lengths
gnuplot
set term dumb
plot "blocks.lengths" with points
# to plot contigs and alignment blocks together:
plot "seqs.lengths" with lines
replot "blocks.lengths" with dots
# can't see 'em? down left, near the origin!
plot [0:1e5] [0:4e6] "blocks.lengths" with dots
# ctrl-d to quit
```





