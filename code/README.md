Downloaded version 1.48.0 of mothur from mothur wiki on github

https://github.com/mothur/mothur/releases/tag/v1.48.0

We used `wget` `mkdir` and `tar` to download and extract files to `data/references/silva_seed`

```
wget -nc -P data/references/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138.tgz
mkdir data/references/silva_seed
tar xvzf data/references/silva.seed_v138.tgz -C data/references/silva_seed
```
