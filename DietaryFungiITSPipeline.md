### Setup of working directory - change this path to match your file system
```
WD_path=/anvil/scratch/x-jnash12/DietaryFungi
mkdir ${WD_path}
cd ${WD_path}
```

### Generates filelist to loop through, including only Mengyi's samples which were also on this run
```
ls /anvil/scratch/x-jnash12/Nash_8872_24011001 | grep "R1_001.fastq.gz" | grep "Mengyi" | sed 's/_R1_001.fastq.gz//g' > ${WD_path}/filelist
```

### Uses PEAR to merge paired reads
```
mkdir ${WD_path}/PEARReads/ && cd ${WD_path}/PEARReads/
for i in $(cat ${WD_path}/filelist)
do
sbatch -o slurm-%j-PEAR.out --partition=shared --account=BIO230020 --export=ALL -t 2:00:00 -c 1 --wrap="pear \
	-f /anvil/scratch/x-jnash12/Nash_8872_24011001/${i}_R1_001.fastq.gz \
	-r /anvil/scratch/x-jnash12/Nash_8872_24011001/${i}_R2_001.fastq.gz \
	-o ${WD_path}/PEARReads/${i}"
done
```

### Uses ITSxpress to extract ITS2 region from merged reads
```
conda activate ITSxpress
mkdir ${WD_path}/ITSxpressReads/ && cd ${WD_path}/ITSxpressReads/
for i in $(cat ${WD_path}/filelist)
do
sbatch -o slurm-%j-ITSxpress.out --partition=shared --account=BIO230020 --export=ALL -t 24:00:00 -c 128 --wrap="itsxpress --fastq ${WD_path}/PEARReads/${i}.assembled.fastq --single_end --region ITS2 --taxa Fungi \
--log ${WD_path}/ITSxpressReads/${i}.logfile.txt --outfile ${WD_path}/ITSxpressReads/${i}.ITS.fastq.gz --threads 128"
done
```

### Prepares QIIME2 manifest
```
printf "%s\t%s\n" "sample-id" "absolute-filepath" > ${WD_path}/QIIMEManifest.tsv
for i in $(cat ${WD_path}/filelist)
do
SampleID=$(echo "$i" | grep -oP 'Mengyi.')
printf "%s\t%s\n" "${SampleID}" "${WD_path}/ITSxpressReads/${i}.ITS.fastq.gz" >> ${WD_path}/QIIMEManifest.tsv
done
```

### Imports data into QIIME2. Submit this as a batch script titled DataImport
```
#!/bin/bash
#SBATCH -o slurm-%j-QIIME_Import.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 4:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/DietaryFungi

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ${WD_path}/QIIMEManifest.tsv \
  --output-path ${WD_path}/ITS2_demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data ${WD_path}/ITS2_demux.qza \
  --o-visualization ${WD_path}/ITS2_demux.qzv

 ```

### Uses Dada2 to denoise, quality filter, chimera filter, and ASV call. Submit this as a batch script titled Dada2SE
```
#!/bin/bash
#SBATCH -o slurm-%j-dada2.out
#SBATCH -c 128
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 48:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/DietaryFungi

qiime dada2 denoise-single \
	--i-demultiplexed-seqs ${WD_path}/ITS2_demux.qza \
	--p-trunc-len 0 \
	--output-dir ${WD_path}/dada2out \
	--p-n-threads 0 \
	--o-table ${WD_path}/ITS2_dada2table.qza \
	--o-representative-sequences ${WD_path}/ITS2_dada2seqs.qza \
	--o-denoising-stats ${WD_path}/ITS2_dada2denoising.qza

qiime metadata tabulate \
	--m-input-file ${WD_path}/ITS2_dada2denoising.qza \
	--o-visualization ${WD_path}/ITS2_dada2denoising.qzv

qiime feature-table summarize \
	--i-table ${WD_path}/ITS2_dada2table.qza \
	--o-visualization ${WD_path}/ITS2_dada2table.qzv

qiime feature-table tabulate-seqs \
	--i-data ${WD_path}/ITS2_dada2seqs.qza \
	--o-visualization ${WD_path}/ITS2_dada2seqs.qzv
```


### OPTIONAL: The Dada2 ASVs can be used as-is, or you can run Vsearch clustering to generate 97% OTUs. Submit this as a batch script titled ClusterDada2
```
#!/bin/bash
#SBATCH -o slurm-%j-vsearch_cluster.out
#SBATCH -c 128
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 6:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/DietaryFungi

qiime vsearch cluster-features-de-novo \
	--i-sequences ${WD_path}/ITS2_dada2seqs.qza\
	--i-table ${WD_path}/ITS2_dada2table.qza \
	--p-perc-identity 0.97 \
	--p-threads 0 \
	--o-clustered-table ${WD_path}/ITS2_Dada2_table97.qza \
	--o-clustered-sequences ${WD_path}/ITS2_Dada2_repseqs97.qza

qiime feature-table summarize \
	--i-table ${WD_path}/ITS2_Dada2_table97.qza \
	--o-visualization ${WD_path}/ITS2_Dada2_table97.qzv

qiime feature-table tabulate-seqs \
	--i-data ${WD_path}/ITS2_Dada2_repseqs97.qza \
	--o-visualization ${WD_path}/ITS2_Dada2_repseqs97.qzv
```

### Taxonomic Classification - Below uses a pre-trained taxonomic classifier for the UNITE database of full length ITS sequences. I am working on training my own classifier for the UNITE database that has been trimmed to only include the ITS2 region. The trimmed ITS2 classifier should perform better, but the pretrained version below 
```
#!/bin/bash
#SBATCH -o slurm-%j-sklearn_classify.out
#SBATCH -c 128
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 6:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/DietaryFungi
cd ${WD_path}

wget https://github.com/colinbrislawn/unite-train/releases/download/v9.0-v25.07.2023-qiime2-2023.9/unite_ver9_dynamic_all_25.07.2023-Q2-2023.9.qza -O unite_ver9_dynamic_all_25.07.2023-Q2-2023.9.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ${WD_path}/unite_ver9_dynamic_all_25.07.2023-Q2-2023.9.qza \
  --i-reads ${WD_path}/ITS2_Dada2_repseqs97.qza \
  --o-classification ${WD_path}/ITS2_Dada2_repseqs97_taxonomy.qza \
  --p-n-jobs -1
```

### ExportData
```
#!/bin/bash
#SBATCH -o slurm-%j-ExportQIIME.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 1:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/DietaryFungi

qiime tools extract \
  --input-path ${WD_path}/ITS2_Dada2_repseqs97_taxonomy.qza \
  --output-path ${WD_path}/ITS2_Dada2_repseqs97_taxonomy

qiime tools export \
  --input-path ${WD_path}/ITS2_Dada2_table97.qza \
  --output-path ${WD_path}/QIIME_exported_files

biom convert -i ${WD_path}/QIIME_exported_files/feature-table.biom -o ${WD_path}/QIIME_exported_files/ITS2_OTUTable_97.tsv --to-tsv

cp ${WD_path}/ITS2_Dada2_repseqs97_taxonomy/*/data/taxonomy.tsv ${WD_path}/QIIME_exported_files/ITS2_Dada2_repseqs97_taxonomy.tsv

qiime tools export \
  --input-path ${WD_path}/ITS2_Dada2_repseqs97.qza \
  --output-path ${WD_path}/QIIME_exported_files

mv ${WD_path}/QIIME_exported_files/dna-sequences.fasta ${WD_path}/QIIME_exported_files/ITS2_Dada2_repseqs97.fasta

```