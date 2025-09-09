# ShigaPass 

ShigaPass is a new *in silico* tool used to predict *Shigella* serotypes and to differentiate between *Shigella*, EIEC (Enteroinvasive *E. coli*), and non *Shigella*/EIEC using assembled whole genomes.

## Dependencies
 ShigaPass is a command-line tool written in **Bash (v4.4.20)** and requires **awk** and **BLAST+ (v2.12.0)**.   
 - BLAST+ can be installed using the following command line (for Ubuntu users) : sudo apt-get install ncbi-blast+ 
 - For more information on how to install Blast+, please refer to https://www.ncbi.nlm.nih.gov/books/NBK279690/

## Installation
**1.** Clone this repository with the following command line:
```
git clone https://github.com/imanyass/ShigaPass.git
```
**2.** Give the execute permission to the file ShigaPass.sh:
```
chmod +x ShigaPass.sh
```
**3.** Execute ShigaPass  with the following command line model:
```
./ShigaPass.sh  [options]
```
## Usage 
Run ShigaPass without option to read the following documentation:
````
###### This tool is used to predict Shigella serotypes  #####
        Usage : ShigaPass.sh [options]
   
        options :
        -l	List of input file(s) (FASTA) with their path(s) (mandatory)
        -o	Output directory (mandatory)
        -p	Path to databases directory (mandatory)
        -t	Number of threads (optional, default: 2)
        -u	Call the makeblastdb utility for databases initialisation (optional, but required when running the script for the first time)
        -k	Do not remove subdirectories (optional)
       	-v	Display the version and exit
        -h	Display this help and exit
        Example: ShigaPass.sh -l list_of_fasta.txt -o ShigaPass_Results -p ShigaPass/ShigaPass_DataBases -t 4 -u -k
        Please note that the -u option should be used when running the script for the first time and after databases updates
````

> [!NOTE]  
>When you use ShigaPass for the first time, you need to **initialise the internal databases**.  
>This is done with the `-u` (update) flag, which **must be run together with a full command**.  
>For example:  
>```
>ShigaPass.sh -l examples/example_list.txt -o results/ -p ShigaPass/ShigaPass_DataBases -u
>```
>This will:
 >- Build/update the required databases.
>- Run ShigaPass on the example dataset.
>
>⚠️ Running **ShigaPass.sh -u** *alone* is not sufficient — the update step is only triggered when combined with a full run. Once the databases are initialised, you can run ShigaPass normally without -u.
>You will need to run with -u again if you update ShigaPass or any of its databases.

## Quick Start Example
Example FASTA sequence files are provided in `Example/Input`.  

   * Please unzip the sequences (using gunzip) before running ShigaPass
### Input File
Create a list file containing the paths to the FASTA files. You can use `readlink` or `ls`.  For example: 
   ```bash
   readlink -e *.fasta > ShigaPass_test.txt
  ```
**Running ShigaPass**

```
 ShigaPass.sh -l ShigaPass_test.txt -o ShigaPass_Results -p ShigaPass_DataBases -u -k
```
- All output files are available in the directory `Example/ShigaPass_Results`
  
Here's an example of ShigaPass summary file
|Name |rfb|rfb_hits,(%)|MLST|fliC|CRISPR|ipaH|Predicted_Serotype|Predicted_FlexSerotype|Comments|
| :------------- |:-------------|:---------------- |:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------| 
|ERR5888634|	C2|	79,(48.2%)|	ST145|	ShH57(ShH3cplx)|	A-var2|	ipaH+|	SB2||
|ERR5952732|	B1-5|	139,(93.3%)|	ST245|	ShH2(ShH2cplx)|	A-var3,x,16|	ipaH+|	SF1-5|	1b||
|ERR5976293|	D|	202,(70.6%)|	ST152|	ShH25(ShH1cplx)|	A-var0,27|	ipaH+|	SS||
|ERR5982186|	A2|	100,(61.7%)|	ST147|	none|	A-var1,12,3,5,11-var1|	ipaH+|	SD2||

"none" means that no allele/profile is detected (in the ERR5982186 example no *fliC* allele was detected)

SB: *S. boydii*; SD: *S. dysenteriae*; SF: *S. flexneri*; SS: *S. sonnei*

### Output Files
* In the output directory, two files will be written:
  1. ShigaPass_summary.csv: semicolon-delimited file with one row per genome inclinding the sample name; type of *rfb*; number of *rfb* hits, (% of *rfb* coverage); MLST profile; type of *fliC*; CRISPR spacers; the presence of *ipaH*; the predicted serotype and *S. flexneri* subserotype; comments to show the number of *rfb* when more than one is detected 
  2. ShigaPass_Flex_summary.csv: semicolon-delimited file detailing the phage and plasmid-encoded O-antigen modification (POAC) genes detected for the predicted *S. flexneri* genomes
 * In case -k option is used, a directory will be created for every assembled genome and will contain the following files:
 
| Extension       | Description |
| :------------- |:-------------| 
|blastout.txt  |Blast results in tabular format|
|allrecords.txt     |Blast hits that passed the selected thresholds| 
|records.txt     |The best blast hit that passed the selected thresholds| 
|hits.txt     |Name and number of hits that passed the selected thresholds (only for k-mers databases: *rfb*, *ipaH* and POAC genes)| 
|hitscoverage.txt|This file displays in addition to the name and the number of hits detected present in hits.txt, the total hits number for the identified gene (3rd column) and the percentage of the hits detected (number of hits detected/total number of hits) (4th column) |

### Notes
The Fasta sequences were assembled using SPAdes version 3.15 (Bankevich et al. Journal of Computational Biology, 2012) with the following options: -k 21,33,55,77  --only-assembler --careful --cov-cutoff auto 

You can download the short reads using the following command lines:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR588/004/ERR5888634/ERR5888634_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR588/004/ERR5888634/ERR5888634_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR595/002/ERR5952732/ERR5952732_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR595/002/ERR5952732/ERR5952732_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR597/003/ERR5976293/ERR5976293_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR597/003/ERR5976293/ERR5976293_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR598/006/ERR5982186/ERR5982186_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR598/006/ERR5982186/ERR5982186_2.fastq.gz
```
All reads were filtered with FqCleanER version 3.0 (https://gitlab.pasteur.fr/GIPhy/fqCleanER) with options -q 15 -l 50 

## Citation

If you use **ShigaPass** in your research or analysis pipeline, please cite:
|**Publication**|**DOI**|
|-------------------------------|------|
|[ShigaPass ](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000961)|![doi](https://img.shields.io/badge/DOI-https://doi.org/10.1099/mgen.0.001280-blue)|  
|[BLAST+ ](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421)|![doi](https://img.shields.io/badge/DOI-https://doi.org/10.1099/mgen.0.001280-blue)|

