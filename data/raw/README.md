Obtained files from the rrnDB located at their downloads page

https://rrndb.umms.med.umich.edu/static/download/ 

These are files from version 5.8 released in 2022
<<<<<<< HEAD
https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.8_pantaxa_stats_RDP.tsv.zip
https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.8_pantaxa_stats_NCBI.tsv.zip
https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.8_16S_rRNA.fasta.zip
https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.8.tsv.zip

We automated downloading and extracting the tsv file with wget wget -P data/raw -nc filename
unzip -n -d data/raw data/raw/rrnDN-5.8.tsv.zip
=======

We automated downloading the tsv file with wget and unzip

wget --no-clobber --directory-prefix=data/raw https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.8_pantaxa_stats_RDP.tsv.zip

unzip -n -d data/raw/rrnDB-5.8.tsv.zip
>>>>>>> master
