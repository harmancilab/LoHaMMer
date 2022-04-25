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

It is necessary to have *tabix*, *bcftools*, and *bgzip* installed on the system.

This script extracts the untyped variants whose MAF is greater than 5%, performs imputation using LoHaMMer with forward-backward algorithm and calculates R2 statistics.

The script can be run using following commands:
```
./example.sh -get_typed_vars
./example.sh -get_parse_1kG
./example.sh -parse_tag_target
./example.sh -build_query_reference_panels
./example.sh -eagle_phase_typed_variants
./example.sh -run_vicinity_HMM
```

After this, run the imputation jobs using following:
```
./q_foreback_submission_script.csh
```
By default, 40 parallel jobs are run.

Now, the genotype R2 can be calculated by following:
```
./example.sh -get_R2_stats
``` 