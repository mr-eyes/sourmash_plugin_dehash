# sourmash_plugin_dehash

## Motivation

In some cases, we need to reverse-hash sourmash signatures to retreive the ACGT k-mers for further investigations or downstream analaysis. This plugin do it effeciently, but might make use of some refactoring and enhancements in the future.


## Description

The plugin is backed with C++ for performance. To use it, you can run the following command:

```
sourmash scripts dehash --sig-paths sample1.sig sample2.sig --fasta-paths raw_reads_1.fa raw_reads_2.fa -k 31 -c 1 -o dehashed_kmers.tsv 
```

## Usage

```
usage:  dehash [-h] [-q] [-d] --sig-paths SIG_PATHS [SIG_PATHS ...] --fasta-paths FASTA_PATHS [FASTA_PATHS ...] -k KSIZE [-c CHUNK_SIZE] -o OUT

Dehash sourmash signatures to ACGT k-mers

options:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -d, --debug           provide debugging output
  --sig-paths SIG_PATHS [SIG_PATHS ...]
                        signature files to process
  --fasta-paths FASTA_PATHS [FASTA_PATHS ...]
                        FASTA files to process
  -k KSIZE, --ksize KSIZE
                        k-mer size
  -c CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        Number of reads per chunk
  -o OUT, --out OUT     output signature path
```


## Installation

This plugin has a C++ extension, so you will need a C++ compiler to install it. On Linux, you can install the `build-essential` package to get a C++ compiler. On macOS, you can install the Xcode command line tools with `xcode-select --install`.

```
pip install git+https://github.com/mr-eyes/sourmash_plugin_dehash
```


## Technical details

We use `simdjson` for parsing the json file, and we use the `parallel_hashmap` to count the kmers in parallel.