Added SILVA reference files to align our sequences against.

https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz

We used wget, mkdir, and tar to download and extract silva seed files to data/references/silva_seed

wget -nc -P data/references/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz
mkdir data/references/silva_seed
tar xvzf data/references/silca_seed.v138.tgz -C data/references/silva_seed/
