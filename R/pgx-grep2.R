## https://cran.r-project.org/web/packages/GREP2/vignettes/vignette.html
if(0) {
    install.packages(c("devtools", "XML", "parallel", "utils", "rentrez", "RCurl"))
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("GEOquery", "Biobase", "tximport", "AnnotationDbi",
               "EnsDb.Hsapiens.v86",  "EnsDb.Mmusculus.v79", "EnsDb.Rnorvegicus.v79",
               "org.Rn.eg.db", "org.Hs.eg.db", "org.Mm.eg.db"))
    devtools::install_github("uc-bd2k/GREP2")    

    ## also need install FastQC and salmon:
    ## >> sudo apt install fastqc salmon sra-toolkit
}

options(warn=-1)
library(GREP2)
## Sys.setenv(PATH="/bin:/usr/bin:/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin:/opt/sratoolkit/bin")

if(0) {
    ##------------------------------------------------
    ## How to use
    ##------------------------------------------------

    FASTQ  = "/home/kwee/Projects/Data/fastq-examples/"  ## folder where FASTQ file are
    OUTPUT = "/tmp/fastq-out/"  ## output folder 
    INDEX = "/data/Projects/Data/salmon_index"     ## existing or where to save Salmon index

    fastq_dir=FASTQ;destdir=OUTPUT;indexdir=INDEX;nthread=4;do.qc=FALSE
    instrument="HiSeq";library_layout="SINGLE"
    nthread=24
    pgx.fastq2counts(fastq_dir=FASTQ, destdir=OUTPUT, indexdir=INDEX, nthread=4, do.qc=FALSE)

    dir(OUTPUT)
    dir(paste0(OUTPUT,"/files_csv"))  ## counts are here
    
}


pgx.fastq2counts <- function(fastq_dir, destdir, indexdir, nthread=4, do.qc=FALSE,
                             instrument="HiSeq", library_layout="SINGLE")
{
    ## 
    ##
    ##
    require(GREP2)

    run_fastqc2 <- function (destdir, fastq_dir, n_thread) 
    {
        cat(paste("Running FastQC... ", Sys.time(), "\n", sep = ""))
        if (!dir.exists(paste0(destdir,"/fastqc"))) {
            system(paste0("mkdir -p ", destdir, "/fastqc"))
        }
        fastq_files = list.files(fastq_dir, pattern = "[.]fastq$", full.names = TRUE)
        length(fastq_files)
        cmd <- paste0("fastqc -o ", destdir, "/fastqc/ --threads ", 
                      n_thread, " ", paste(fastq_files,collapse=" "))
        system(cmd)
    }
    
    ## ----------- Unzip FASTQ files (if necessary)
    dir(fastq_dir)
    dir(fastq_dir,"[.]gz$",include.dirs=FALSE)
    if(length(dir(fastq_dir,"[.]gz$",include.dirs=FALSE))) {
        ##cmd <- paste0("gunzip ",fastq_dir,"/*gz")
        cmd <- paste0("(cd ",fastq_dir," && gunzip *gz)")
        cmd
        system(cmd)
    }
    
    ## ----------- Run FastQC on each fastq file to generate quality control (QC) reports.
    if(do.qc) {
        if(dir.exists(paste0(fastq_dir,"/fastqc"))) system(paste0("rm -fr ", fastq_dir,"/fastqc"))
        run_fastqc2(destdir=destdir, fastq_dir=fastq_dir, n_thread=nthread)
    }
    
    ## ----------- Before running Salmon, you will have to build index first.
    dir.exists(indexdir)
    if(!dir.exists(indexdir)) {
        system(paste("mkdir -p",indexdir))
        build_index(species="human", kmer=31, ens_release=92, destdir=indexdir)
    }
    
    file_id <- sub(".fastq$","",dir(fastq_dir, pattern="001.fastq$"))
    file_id
    
    ## ----------- Run Trimmomatic
    system(paste0("rm -f ",fastq_dir,"/*_trimmed.fastq"))
    i=1
    for(i in 1:length(file_id)){
        trim_fastq(srr_id=file_id[i], fastq_dir=fastq_dir,
                   instrument=instrument, library_layout=library_layout,
                   destdir=fastq_dir, n_thread=nthread)
    }
    
    ## ----------- Run Salmon
    index2 = file.path(indexdir,"human_transcripts_release92_index")
    i=1
    for(i in 1:length(file_id)) {
        run_salmon(srr_id=file_id[i], library_layout=library_layout,
                   index_dir=index2, destdir=destdir,
                   fastq_dir=fastq_dir, use_trimmed_fastq=TRUE,
                   other_opts="--validateMappings", n_thread=nthread)
    }
    
    ## ----------- Run MultiQC
    if(do.qc) {
        run_multiqc(fastqc_dir=fastq_dir,salmon_dir=destdir, destdir=destdir)
    }
    
    ## ----------- Run tximport
    txi <- run_tximport(srr_id=file_id, species="human",
                        salmon_dir=paste0(destdir,"/salmon"),
                        countsFromAbundance="lengthScaledTPM")
    names(txi)
    
    ## ----------- Extract counts
    genes  <- txi$gene_counts[,2:3]
    counts <- as.matrix(txi$gene_counts[,4:ncol(txi$gene_counts),drop=FALSE])
    rownames(genes) <- rownames(counts) <- txi$gene_counts$ENSEMBL

    ## ----------- Collapse multiple probes to single gene by summing up counts
    counts <- apply(counts, 2, function(x) tapply(x,genes$SYMBOL,sum))
    genes <- genes[match(rownames(counts),genes$SYMBOL),]
    rownames(genes) <- genes$SYMBOL
    counts <- round(counts, digits=3)
    
    ## ----------- Save counts as text files for Omics Playground
    colnames(genes) <- c("gene_name","gene_title")
    system(paste0("mkdir -p ",destdir,"/files_csv"))
    cat("writing files to",paste0(destdir,"/files_csv/ \n"))
    write.csv(counts, file=paste0(destdir,"/files_csv/counts.csv"))
    write.csv(genes,  file=paste0(destdir,"/files_csv/genes.csv"))
    samples <- data.frame(sample=colnames(counts), group=NA, phenotype1=NA)  ## empty template
    write.csv(samples,  file=paste0(destdir,"/files_csv/samples.csv"), row.names=FALSE)
}


