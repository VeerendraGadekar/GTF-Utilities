### NAME

<pre><code>dissect_gtf.R </code></pre>

### SYNOPSIS

A utility to readily extract gene and transcript level bed files, statistics on transcripts per gene and exons per transcript, and gene and transcript counts per annotated biotype from a gtf file. The script also generates a SQLite database file (which could be accessed from
programs other than R if desired).

#### Usage example:

<pre><code>Rscript dissect_gtf.R --gtf path_to_standard.gtf --outdir path_to_desired_output_directory
Rscript dissect_gtf.R --gtf path_to_standard.gtf --outdir path_to_desired_output_directory --biotype gene_biotype,transcript_biotype
</code></pre>



#### Options:

<pre><code>Rscript dissect_gtf.R --help

Options:
	--gtf=CHARACTER
		     path_to_gtf_file

	--TxDB=CHARACTER
		     path_to_txdb_file

	--biotype=CHARACTER
		     Use this flag to get gene and transcript counts per biotype. 
                     (Only when input GTF contains a biotype meta data). 
                     This requires you to enter the right keys, which could be 
	             either "gene_type,transcript_type" or "gene_biotype,transcript_biotype"
		     depending upon your gtf file

	--outdir=CHARACTER
		     output folder name [default= NULL]

	-h, --help
		     Show this help message and exit
</code></pre>

If a SQLite annotation database already exists, it can be used as input instead of a GTF file


### INPUT
A standard GTF file or a SQLite database file with gene annotation

### OUTPUT

The script generates following files -
* SQLite database file (if .gtf file is used as input)
* gene level .bed file
* transcript level .bed file
* per gene transcript count .tsv file
* per transcript exon count .tsv file
* gene and transcript counts per annotated biotype (only when --biotype flag is used)
* log file with any warnings generated while creating the TxDB object from .gtf file

### DEPENDENCIES

R packages - [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html), [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html)


### AUTHOR
Veerendra Gadekar, IBCH PAS, contact [gpveerendra09@gmail.com](mailto:gpveerendra09@gmail.com)
