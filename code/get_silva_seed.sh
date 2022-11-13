#!/usr/bin/env bash
# Author Richard Stack
#Inputs none
#outputs places SILVA SEED reference alignment into data/references/dilva_seed
#
#Download this version of the SILVA reference to help with aligning our sequence data. This is version 138 which was released in 2020. Because the tgz contains a README file we extracted to a directory within data/references

wget -nc -P data/references/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138.tgz
mkdir data/references/silva_seed
tar xvzf data/references/silva.seed_v138.tgz -C data/references/silva_seed
