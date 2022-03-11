### NAME

<pre><code>dissect_gtf.R </code></pre>

## SYNOPSIS

A utility to quickly extract gene and transcript bed files, and get the statistics on transcripts per gene and exons per transcript. The script also generates a SQLite database file (which could be accessed from
programs other than R if desired).

#### Usage example:

<pre><code>Rscript dissect_gtf.R --gtf path_to_standard.gtf --outdir path_to_desired_output_directory
</code></pre>

#### Options:

<pre><code>Rscript dissect_gtf.R --help

Options:
    --gtf=CHARACTER
         path_to_gtf_file

    --TxDB=CHARACTER
         path_to_txdb_file

    --outdir=CHARACTER
        output folder name [default= NULL]

    -h, --help
        Show this help message and exit
</code></pre>


If a SQLite annotation database already exists, it can be used as input instead of a GTF file

## INPUT
A standard GTF file or a SQLite database file with gene annotation

## OUTPUT

The script generates following files -
* SQLite database file (if .gtf file is used as input)
* gene level .bed file
* transcript level .bed file
* per gene transcript count .tsv file
* per transcript exon count .tsv file
* log file with any warnings generated while creating the TxDB object from .gtf file

## DEPENDENCIES

R packages - [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html)


## AUTHOR
Veerendra Gadekar, IBCH PAS, contact [gpveerendra09@gmail.com](mailto:gpveerendra09@gmail.com)


```python

```
