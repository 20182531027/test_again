1. [fastp quality control](#Fastp)
2. [bwa alignment](#BWA)
	2.1 samtools view the result file
3. genome assemble
    3.1 [Spades](https://github.com/ablab/spades)
	3.2 [Novoplasty](https://github.com/ndierckx/NOVOPlasty)
4. Blast to identify the type of HPV 
# 1. [Fastp](https://github.com/OpenGene/fastp/blob/master/README.md)
## 1.1 Install

```bash
conda install -c bioconda fastp
```
Already installed, path: /home/chenqi5/miniconda3/bin/fastp
## 1.2 Usage

``look  the  useage of fastp :``
```bash
fastp -h 
```

`common options:`
```
fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
-i, --in1  read1 input file name (string)
-o, --out1  read1 output file name (string [=])
-I, --in2   read2 input file name (string [=])
-O, --out2  read2 output file name (string [=])
# reporting options
-j, --json the json format report file name (string [=fastp.json])
-h, --html the html format report file name (string [=fastp.html])
```

`example:`
```bash
# the directory info
pwd
```

```bash
#output
/home/chenqi5/result
```

```bash
fastp -i BGIHPVTEST1_L01_1-9.hpv.1.fq.gz \
-o BGIHPVTEST1_L01_1-9.hpv.qc.1.fq.gz \
-I BGIHPVTEST1_L01_1-9.hpv.2.fq.gz \
-O BGIHPVTEST1_L01_1-9.hpv.qc.2.fq.gz \
-h BGIHPVTEST1_L01_1-9.hpv.qc.html \
-j BGIHPVTEST1_L01_1-9.hpv.qc.json
```

`notes`:
-i, -o info are paired, so are the -I and -O;
output fastq after fastp were sent for next alignment; while the .html or .json file was the quality control report. 
## 1.3 look result 

use [**Seqkit**](https://bioinf.shenwei.me/seqkit/) to look the changes:
```bash 
seqkit stats BGIHPVTEST1_L01_1-9.hpv.1.fq.gz
```
![[Pasted image 20241019105944.png]]

```bash 
seqkit stats BGIHPVTEST1_L01_1-9.hpv.qc.1.fq.gz
```
![[Pasted image 20241019110134.png]]

use **WinSCP** to see the .html file:
	*insatll the AccessClient.exe*
	*install the WinSCP-6.3.5-Setup.exe*

# 2. [BWA](https://bio-bwa.sourceforge.net/bwa.shtml#3)
## 2.1 Install

```bash
conda install bwa
```
## 2.2 Usage
### 2.2.1 bulid index
```bash
bwa index pave_hpv193ref.fasta
```
### 2.2.2 mem alignmnet
```bash
/home/chenqi5/miniconda3/bin/bwa/ mem -t 4 /home/chenqi5/file/test_again/pave_hpv193ref.fasta \ # the reference
            BGIHPVTEST1_L01_1-9.hpv.qc.1.fq.gz \ # reads1 path
            BGIHPVTEST1_L01_1-9.hpv.qc.2.fq.gz \ # reads2 path
            > BGIHPVTEST1_L01_1-9.sam \ # output sam file
            2> bwa.log # the record of process
```

**transform  *.sam to *.bam**
```bash
samtools view -bS BGIHPVTEST1_L01_1-9.sam > BGIHPVTEST1_L01_1-9.bam
```

![[Pasted image 20241019110803.png]]
The [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map) format and the BAM (Binary Alignment/Map) format, use `less` to see SAM directly, while `samtools view` to see BAM. BAM could save more storage space.

**extract aligned reads**
```bash
samtools view -bS -F 12 BGIHPVTEST1_L01_1-9.bam| samtools sort -n | samtools fastq -1 BGIHPVTEST1_L01_1-9.hpv.map.1.fq -2 BGIHPVTEST1_L01_1-9.hpv.map.2.fq
```


# 3. [Spades](https://ablab.github.io/spades/getting-started.html) for genome Assemble
```bash
python /home/chenqi5/miniconda3/bin/spades.py -o ./ --pe1-1 BGIHPVTEST1_L01_1-9.hpv.map.1.fq --pe1-2 BGIHPVTEST1_L01_1-9.hpv.map.2.fq -t 4
```


# 4. [Blast](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a) to identify the type of HPV


![[Pasted image 20241022104258.png]]

### [Online blast ](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
![[Pasted image 20241022104506.png]]



### Local blast
contrust the database, already constructed
```bash
makeblastdb \
 -dbtype nucl \
 -in pave_hpv193ref.fasta \
 -input_type fasta \
 -parse_seqids \
 -out hpv193ref
```
Parameter Description:
- `-dbtype` is a mandatory parameter, choose either "nucl" or "prot" for nucleic acid sequence or protein sequence database.
- `-in` is used for the FASTA file to build the search database.
- `-input_type` defaults to FASTA files, other supported file formats include asn1_bin, asn1_txt, and blastdb.
- `-parse_seqids` with this parameter, the result file will generate three additional files; without it, the default is three result files.
- `-out` the prefix for the output file.
After executing the above code, six files with the suffixes nin, nhr, nsq, nsi, nsd, nog will be generated, indicating that the custom search database files have been successfully created.

```bash
blastn \
 -query NODE_1_length_8015_cov_315.231156.fasta \
 -db ~/file/test_again/hpv193ref \
 -out blastn_results.xls \
 -task blastn
```


```bash
blastn \
 -query scaffolds.fasta \
 -db ~/file/test_again/hpv193ref \
 -out scaffolds_blastn_results.xls \
 -task blastn
```

```bash
blastn -query scaffolds.fasta -db ~/file/test_again/hpv193ref -outfmt 7 -out scaffolds_blastn_results.txt
```

