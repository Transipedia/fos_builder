
# fos_builder

**companon tools for Transipedia website and Reindeer.**

Until to version 1.0.2 of Reindeer (current version), the fos_builder script aims to add some info to Reindeer output.

the aim is to add an additional file to the index containing the names of the unitigs files used by Reindeer when creating an index. This will enable column names to be displayed in count tables.

In addition, if quality control has been carried out with fastqc and multiqc, fos_buider will add the total number of kmer for each of the files, enabling counts to be normalized at a later date.

## usage

```
fos_builder.py [-h] [-k KMER_LEN] [-o OUTPUT] [-v] fof_unitigs.txt [multiqc_dir]

positional arguments:
  fof_unitigs.txt       fof_unitigs.txt file.
  multiqc_dir           Directory of multiqc results.

options:
  -h, --help            show this help message and exit
  -k KMER_LEN, --kmer-len KMER_LEN
                        kmer length (default: 31)
  -o OUTPUT, --output OUTPUT
                        output file (default: stdin)
  -v, --version         show program's version number and exit
```

  Example:

```
  fos_builder.py  unitigs/fof_unitigs.txt  quality/multiqc/ -o reindeer-index/fos.txt
```

  The fof_unitigs file contain the list of paths to unitig files, as example:

```
output/bcalm/unitigs/SRR00000001.unitigs.fa
output/bcalm/unitigs/SRR00000002.unitigs.fa
output/bcalm/unitigs/SRR00000003.unitigs.fa
output/bcalm/unitigs/SRR00000004.unitigs.fa
output/bcalm/unitigs/SRR00000005.unitigs.fa
output/bcalm/unitigs/SRR00000006.unitigs.fa
output/bcalm/unitigs/SRR00000007.unitigs.fa
```

## Important

The next version of Reindeer will automatically add column headers and integrate a kmer counting system for normalization purposes. fos_builder has been created on a temporary basis.