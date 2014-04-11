#!/usr/bin/env bash
rm -r test-output/

Rscript test.R -o test-output -l old_exome.interval_list test.maf
