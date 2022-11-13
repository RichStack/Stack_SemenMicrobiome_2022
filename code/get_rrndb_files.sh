#!/usr/bin/env bash
#Author: Richard Stack
# Inputs: Name of the file extracted from the archive (without the path)
#Outputs: The appropriate file into data/raw

archive=$1

wget -nc -P data/raw https://rrndb.umms.med.umich.edu/static/download/"$archive".zip
unzip -n -d data/raw data/raw/"$archive".zip
