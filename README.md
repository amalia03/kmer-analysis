# kmer-analysis
Includes Perl and R scripts that perform various kmer analyses from FASTA files. 

### Programs required

Programs described here are running using the following versions for: 
- Perl: (v5.16.3)
- R:  3.6.0 

For the rapid kmer assesment, use the following C code from lmj: https://github.com/lmjakt/R_c_plugins

### Input data format

Use files in FASTA format. The Identifiers are not used, only the sequence content. Please make sure that the file does not contain any Ns or lowercase letters as the C code cannot operate with multisemantic letters.

**
Please refer to the Wiki for more information on the R functions.

