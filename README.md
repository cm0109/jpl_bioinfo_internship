# All Pipelines Run in Loop

[![hackmd-github-sync-badge](https://hackmd.io/fpv-E78CSwWopDzmIiLSfw/badge)](https://hackmd.io/fpv-E78CSwWopDzmIiLSfw)


<br>

## System Description
Running a **m1.xxlarge** (CPU: 44, Mem: 120 GB, Disk: 60 GB) cloud instance, with additional 320 Gb storage volume (/vol_b). This storage volume will be used to store the downloaded databases for all 3 pipelines.

Connection: ssh cmicro@<IP>

<br>

## Downloading the metagenomic dataset

#### From the home directory (/home/cmicro/)

```
mkdir data
cd data

curl -L -o metagenomic-read-files.tar.gz \
https://ndownloader.figshare.com/files/24079451

tar -xzvf metagenomic-read-files.tar.gz
```

# Make a copy of the files and then unzip R1 to count how many reads

find . -name "*_R1*"  -exec gunzip {} \;
find . -type f -name "*_R1*" | sort -n | xargs grep -c "@" > fastq_read_counts.txt


<br>

## Kraken2-Bracken

### Created a script (kraken2.sh) on local computer

```
#!/bin/bash

# This script runs Kraken2 on paired FASTQ files & Bracken on K2 report file

# Required: Kraken2 & Bracken installed, trimmed FASTQs/links in CWD

# arg1: number of threads
# arg2: location of db

# Make script executable: chmod +x kraken2_bracken.sh

# To run: 
# <path>/kraken2.sh <number_of_threads> <location_of_db>
# Example: /kraken2.sh 42 /vol_b/kraken2-db/

for f in *_R1_trimmed.fastq.gz # for each sample F

do
    n=${f%%_R1_trimmed.fastq.gz} # strip part of file name

	kraken2 --threads $1 --db $2 \
	--paired ${n}_R1_trimmed.fastq.gz ${n}_R2_trimmed.fastq.gz \
	--output ${n}_kraken2_out.txt --report ${n}_kraken2_report.txt

done

for k in *_kraken2_report.txt # for each sample k2 report file

do
    m=${k%%_kraken2_report.txt} # strip part of file name

	bracken -r 150 -d $2 -i ${m}_kraken2_report.txt -o ${m}_bracken_out.txt

done
```


#### Transferred BASH script from local computer to cloud instance

```
scp -C /Users/cm/JPL_Google_Drive/scripts/metagenomics/kraken2.sh \
cmicro@129.114.17.69:/home/cmicro/scripts
```

<br>

### On the cloud instance


#### Make bash script executable
```
chmod +x /home/cmicro/scripts/kraken2.sh
```


#### Create new directory for running Kraken2 + Bracken
```
mkdir k2b_all
cd k2b_all
```

#### Make Symbolic Links for all files in the sample directory

```
ln -s /home/cmicro/data/*_trimmed.fastq.gz .
```


#### Create new conda environment & install Kraken2 & Bracken
```
conda create -y -n kraken2 -c conda-forge -c bioconda -c defaults kraken2=2.0.9beta bracken=2.6.0
conda activate kraken2
```


#### Calling script to execute Kraken2 command over multiple FASTQ files
```
/home/cmicro/scripts/kraken2.sh 42 /vol_b/kraken2-db/
```

#### Moving output report files
```
mkdir bracken_out
find . -type f -iname "*bracken_out*" -exec mv '{}' bracken_out/ \;
```

#### Transfer report files back to local computer for analysis
```
scp -r cmicro@129.114.17.69:/home/cmicro/k2b_all/bracken_out \
/Users/cm/JPL_Google_Drive/hbcu_data_analysis
```

<br>
<br>


## Centrifuge

#### Created BASH script for running centrifuge on multiple fastq files

```
#!/bin/bash

# This script runs Centrifuge on paired FASTQ files & generates Kraken style report

# Required: Centrifuge installed, trimmed FASTQs/links in CWD, Centrifuge db built

# arg1: number of threads
# arg2: db prefix

# Make script executable: chmod +x kraken2_bracken.sh

# To run: 
# <path>/centrifuge.sh <number_of_threads> <location_of_db>
# Example: /centrifuge.sh 42 #/vol_b/centrifuge-db/centrifuge-complete-genomes-arc-bac-human-viral-fungi


for f in *_R1_trimmed.fastq.gz # for each sample F

do
    n=${f%%_R1_trimmed.fastq.gz} # strip part of file name
	
	centrifuge -p $1 -x $2 -k 1 -1 ${n}_R1_trimmed.fastq.gz -2 ${n}_R2_trimmed.fastq.gz -S ${n}_centrifuge_out.txt --report-file \
	${n}_centrifuge_report.txt

done

for r in *_centrifuge_out.txt # for each sample output file

do
    m=${r%%_centrifuge_out.txt} # strip part of file name
	
	centrifuge-kreport -x $2 ${m}_centrifuge_out.txt > ${m}_centrifuge_reformatted_out.txt

done
```

<br>

#### Transferred BASH script from local computer to cloud instance

```
scp -C /Users/cm/JPL_Google_Drive/scripts/metagenomics/centrifuge.sh \
cmicro@129.114.17.69:/home/cmicro/scripts
```

<br>

### On the cloud instance

#### Make script executable
```
chmod +x /home/cmicro/scripts/centrifuge.sh
```

#### Download Centrifuge database
```
scp -r astrobio@149.165.170.83:/vol_b/centrifuge-db/ /vol_b/
```

#### Create & move to working directory
```
mkdir centrifuge_all
cd centrifuge_all
```


#### Create links to sample FASTQ files

```
ln -s /home/cmicro/data/*_trimmed.fastq.gz .
```


#### Installing Centrifuge


```
conda create -y -n centrifuge -c conda-forge -c bioconda -c defaults centrifuge=1.0.4_beta
conda activate centrifuge
```



#### Calling script to execute Centrifuge command over multiple FASTQ files

```
/home/cmicro/scripts/centrifuge.sh 42 /vol_b/centrifuge-db/centrifuge-complete-genomes-arc-bac-human-viral-fungi
```

#### Moving report files for transfer
```
mkdir cent_out
find . -type f -iname "*reformatted*" -exec mv '{}' cent_out/ \;
```

#### Transfer report files to local computer for downstream analysis
```
scp -r cmicro@129.114.17.69:/home/cmicro/centrifuge_all/cent_out \
/Users/cm/JPL_Google_Drive/hbcu_data_analysis
```

<br>
<br>



## Ganon

#### Created a BASH script to run Ganon

```
#!/bin/bash

# This script runs Ganon on paired FASTQ files & generates report file

# Required: Ganon installed, trimmed FASTQs/links in CWD, Ganon db built

# arg1: number of threads
# arg2: db prefix

# Make script executable: chmod +x kraken2_bracken.sh

# To run: 
# <path>/ganon.sh <number_of_threads> <db prefix>
# Example: ganon.sh 42 /vol_b/ganon-db/ganon-complete-genomes-arc-bac-human-viral-fungi

for f in *_R1_trimmed.fastq.gz # for each sample F

do
    n=${f%%_R1_trimmed.fastq.gz} # strip part of file name

	echo "ganon classify --threads $1 --db-prefix $2 \
	--paired-reads ${n}_R1_trimmed.fastq.gz ${n}_R2_trimmed.fastq.gz \
	-o ${n}_ganon"

done

for r in *_ganon.rep # for each sample ganon report file

do
    m=${r%%_ganon.rep} # strip part of file name
	
	echo "ganon report --db-prefix $2 --rep-file ${m}_ganon.rep --ranks species \
	--output-report ${m}_ganon_species.txt"

done
```

#### Transferred BASH script from local computer to cloud instance

```
scp -C /Users/cm/JPL_Google_Drive/scripts/metagenomics/ganon.sh \
cmicro@129.114.17.69:/home/cmicro/scripts
```


<br>


### On the cloud instance:

#### Make script executable
```
chmod +x /home/cmicro/scripts/ganon.sh
```


#### Download the pre-built Ganon database (~100 GB)

```
scp -r astrobio@149.165.170.83:/vol_b/ganon-db/ /vol_b/
```



#### Create & move to working directory
```
mkdir ganon_all
cd ganon_all
```


#### Create links to sample FASTQ files

```
ln -s /home/cmicro/data/*_trimmed.fastq.gz .
```


#### Install & Activate Conda env. for Ganon

```
conda create -y -n ganon -c conda-forge -c bioconda -c defaults ganon=0.2.3
conda activate ganon
```


#### Calling script to execute Ganon command over multiple FASTQ files

```
/home/cmicro/scripts/ganon.sh 42 /vol_b/ganon-db/ganon-complete-genomes-arc-bac-human-viral-fungi
```

#### Moving report files for transfer
```
mkdir ganon_out
find . -type f -iname "*species*" -exec mv '{}' ganon_out/ \;
```

#### Transfer report files to local computer for downstream analysis
```
scp -r cmicro@129.114.17.69:/home/cmicro/ganon_all/ganon_out \
/Users/cm/JPL_Google_Drive/hbcu_data_analysis
```

<br>
<br>