# LoHaMMer: Local Li-Stephens Hidden Markov Modeler for Phased Genotype Imputation

This repository contains the source code and documentation for LoHaMMer -- local hidden Markov model-based genotype imputation evaluation.

## Build ##
Run following commands to build LoHaMMer:
```
make clean
make
```

The executable is built at 'bin/LoHaMMer'.

## Example Run ##
We created a bash script that can be used to parse and use Illumina Duo v3 typed variants using 1000 Genomes chromosome 20 dataset. This script is located under EXAMPLE/example.sh

This script extracts the untyped variants whose MAF is greater than 1%.

The script can be run by running each code blocks in order.