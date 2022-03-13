#!/usr/bin/env Rscript

# AUTHOR: Veerendra Gadekar
# CONTACT: gpveerendra09@gmail.com

args = commandArgs(TRUE)

library("optparse")

option_list = list(

  make_option(
     c("--gtf"),
     type="character",
     default=NULL,
     help=" path_to_gtf_file",
     metavar="character"),

  make_option(
     c("--TxDB"),
     type="character",
     default=NULL,
     help=" path_to_txdb_file",
     metavar="character"),

  make_option(
     c("--biotype"),
     type="character",
     default=NULL,
     help=" Use this flag to get gene and transcript counts per biotype. 
                 (Only when input GTF contains a biotype meta data). 
                 This requires you to enter the right keys, which could be 
	         either \"gene_type,transcript_type\" or \"gene_biotype,transcript_biotype\"
		 depending upon your gtf file",
     metavar="character"),

  make_option(
     c("--outdir"),
     type="character",
     default=NULL,
     help="output folder name [default= %default]",
     metavar="character")
 )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test if there is at least one argument: if not, return an error
if(length(args)==0) {
  print_help(opt_parser)
  stop("At least one argument must be supplied.", call.=FALSE)
 }


if(length(opt$biotype != 0)){

       cat("\nLoading R packages ... ")

       suppressMessages(library(GenomicFeatures, warn.conflicts = F, quietly = T))
       suppressMessages(library(data.table, warn.conflicts = F, quietly = T))

       file <- opt$gtf
       fname <-  gsub("\\..*","", basename(file))
       wd <- opt$outdir

       warning_file = file(paste(wd,"/",fname,"_dissect_gtf_biotypes.log", sep=""), open = "wt")
       sink(warning_file, type = "message")

       biotype_col <- unlist(strsplit(opt$biotype, ","))
       gene_biotype_col <- biotype_col[1]
       tx_biotype_col <- biotype_col[2]

       cat("\nloading gtf file ...")
       gtffile <- rtracklayer::import(file)
 
       cat("\ngenerating biotype statistics ... ")

       gene_type <- setDT(unique((as.data.frame(gtffile)[,c("gene_id", get("gene_biotype_col"))])))
       setkeyv(gene_type,names(gene_type)[2])
       gene_type_counts <- gene_type[, .(gene_counts = uniqueN(gene_id)), key(gene_type)]
       gene_type_counts <- gene_type_counts[order(gene_counts, decreasing = TRUE)]

       write.table(
         gene_type_counts, 
         file=paste(wd,"/",fname,"_gene_biotype_counts.tsv", sep=""), 
         quote=F, sep="\t", row.names=F)

       tx_type <- setDT(unique((as.data.frame(gtffile)[,c("transcript_id", get("tx_biotype_col"))])))
       setkeyv(tx_type,names(tx_type)[2])    
       tx_type_counts <- tx_type[, .(transcript_counts = uniqueN(transcript_id)), key(tx_type)]
       tx_type_counts <- tx_type_counts[order(transcript_counts, decreasing = TRUE)]
    
       write.table(
         tx_type_counts, 
         file=paste(wd,"/",fname,"_transcript_biotype_counts.tsv", sep=""), 
         quote=F, sep="\t", row.names=F)

       gtffile <- NULL

     cat("\nAll done!\n")
     cat("\nCheck the log file for warnings, if any!\n")

 }else{

    if(length(opt$TxDB) != 0){

        cat("\nLoading R packages ... ")
        suppressMessages(library(GenomicFeatures, warn.conflicts = F, quietly = T))
        suppressMessages(library(data.table, warn.conflicts = F, quietly = T))

        file <- opt$TxDB
        fname <-  gsub("\\..*","", basename(file))
        wd <- opt$outdir

        warning_file = file(paste(wd,"/",fname,"_dissect_gtf.log", sep=""), open = "wt")
        sink(warning_file, type = "message")

        cat("\nLoading TxDB object ... ")
        TxDb  <- loadDb(file)
      }

    if(length(opt$gtf) != 0){

        cat("\nLoading R packages ... ")
        suppressMessages(library(GenomicFeatures, warn.conflicts = F, quietly = T))
        suppressMessages(library(data.table, warn.conflicts = F, quietly = T))

        file <- opt$gtf
        fname <-  gsub("\\..*","", basename(file))
        wd <- opt$outdir

        warning_file = file(paste(wd,"/",fname,"_dissect_gtf.log", sep=""), open = "wt")
        sink(warning_file, type = "message")

        cat("\nCreating TxDB object ... ")
        TxDb  <- makeTxDbFromGFF(file)
        saveDb(TxDb, file=paste(wd,"/",fname,".sqlite", sep=""))
      }

    cat("\nGenerating results ... ")

    genes <- genes(TxDb)
    transcripts <- transcripts(TxDb)

    convert_GRanges_to_bed <- function(gr, names){
        df <- data.frame(
                seqnames = seqnames(gr),
                starts = start(gr)-1,
                ends = end(gr),
                names = mcols(gr)[,names],
                scores = c(rep(".", length(gr))),
                strands = strand(gr))
        fname2 <- deparse(substitute(gr))   
        write.table(df, file=paste(wd,"/",fname,"_",fname2,".bed", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
        df
     }

    gb <- convert_GRanges_to_bed(genes, 'gene_id')
    tb <- convert_GRanges_to_bed(transcripts, 'tx_name')

    getx <- transcriptsBy(TxDb)
    getxdt <- setDT(as.data.frame(as.data.frame(getx)))
    setnames(getxdt, old = c("group_name","tx_name"), new = c("gene_id","transcript_id"))
    getxdt_counts <- getxdt[, .(tx_count=.N),by=gene_id]
    write.table(getxdt_counts, paste(wd,"/",fname,"_gene_transcripts_counts.tsv", sep=""), quote=F, sep="\t", row.names=F, col.names=F)

    txex <- exonsBy(TxDb)
    txexdt <- as.data.frame(as.data.frame(txex))
    txexdt <- setDT(merge(txexdt, subset(as.data.frame(transcripts), select = c(tx_id, tx_name)), by.x = "group_name", by.y = "tx_id"))
    txexdt_counts <- txexdt[, .(ex_count=.N),by=tx_name]
    write.table(txexdt_counts, paste(wd,"/",fname,"_transcripts_exon_counts.tsv", sep=""), quote=F, sep="\t", row.names=F, col.names=F)


    cat("\nAll done!\n")
    cat("\nCheck the log file for warnings, if any!\n")

}
