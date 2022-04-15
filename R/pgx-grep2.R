##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


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
                             species = c("human","mouse","rat"),
                             quant.method = "salmon",
                             trimming=TRUE, trimmethod="trimmomatic",
                             instrument="HiSeq", library_layout="SINGLE")
{
    
    species <- species[1]
    
    run_fastqc1 <- function (destdir, fastq_dir, nthread) 
    {
        cat(paste("Running FastQC... ", Sys.time(), "\n", sep = ""))
        dir.exists(paste0(destdir,"/fastqc"))        
        if (!dir.exists(paste0(destdir,"/fastqc"))) {
            system(paste0("mkdir -p ", destdir,"/fastqc"))
        } else {
            system(paste0("rm -fr ", destdir,"/fastqc"))
            system(paste0("mkdir -p ", destdir,"/fastqc"))            
        }
        fastq_files = list.files(fastq_dir, pattern = "[.]fastq$", full.names = TRUE)
        length(fastq_files)
        cmd <- paste0("fastqc -o ", destdir, "/fastqc/ --threads ", 
                      nthread, " ", paste(fastq_files,collapse=" "))
        system(cmd)
    }

    run_fastqc2 <- function (destdir, fastq_dir, nthread) 
    {
        cat(paste("Running FastQC... ", Sys.time(), "\n", sep = ""))
        dir.exists(paste0(destdir,"/fastqc_trimmed"))        
        if (!dir.exists(paste0(destdir,"/fastqc_trimmed"))) {
            system(paste0("mkdir -p ", destdir,"/fastqc_trimmed"))
        } else {
            system(paste0("rm -fr ", destdir,"/fastqc_trimmed"))
            system(paste0("mkdir -p ", destdir,"/fastqc_trimmed"))            
        }
        fastq_files = list.files(fastq_dir, pattern = "_trimmed.fastq$", full.names = TRUE)
        length(fastq_files)
        cmd <- paste0("fastqc -o ", destdir, "/fastqc_trimmed/ --threads ", 
                      nthread, " ", paste(fastq_files,collapse=" "))
        system(cmd)
    }
    
    ## ----------- Unzip FASTQ files (if necessary)
    dir(fastq_dir)
    dir(fastq_dir,"[.]gz$",include.dirs=FALSE)
    if(length(dir(fastq_dir,"[.]gz$",include.dirs=FALSE))) {
        cat(">>> Unzipping FASTQ files... \n")
        ##cmd <- paste0("gunzip ",fastq_dir,"/*gz")
        cmd <- paste0("(cd ",fastq_dir," && gunzip *gz)")
        cmd
        system(cmd)
    }

    ## remove any previously trimmed sequences
    system(paste0("rm -f ",fastq_dir,"/*_trimmed.fastq"))
    system(paste0("rm -f ",fastq_dir,"/*trimmed*"))
    system(paste0("rm -f ",fastq_dir,"/*trimming*"))    
    
    ## ----------- Run FastQC on each fastq file to generate quality control (QC) reports.
    if(do.qc) {
        cat(">>> Running FastQC ... \n")
        run_fastqc1(destdir=destdir, fastq_dir=fastq_dir, nthread=nthread)
    }
            
    ## ----------- Run Trimmomatic (some say it is not necessary for salmon...)
    file_id <- sub(".fastq$","",dir(fastq_dir, pattern=".fastq$"))
    file_id <- grep("trimmed$",file_id,value=TRUE,invert=TRUE)
    file_id
    ## trimmethod="trimmomatic"
    if(trimming && trimmethod=="trimmomatic") {
        cat(">>> Running Trimmomatic... \n")
        i=1
        for(i in 1:length(file_id)){
            GREP2::trim_fastq(srr_id=file_id[i], fastq_dir=fastq_dir,
                              instrument=instrument, library_layout=library_layout,
                              destdir=fastq_dir, n_thread=nthread)
        }
    }
    if(trimming && trimmethod=="trimgalore") {
        cat(">>> Running TrimGalore... \n")
        fastq1=fastq2=NULL
        if(library_layout=="SINGLE") {
            fastq1 <- file.path(fastq_dir,paste0(file_id,".fastq"))
            fastq2 = NULL
        } else {
            fastq1 <- file.path(fastq_dir,paste0(file_id,"_1.fastq"))
            fastq2 <- file.path(fastq_dir,paste0(file_id,"_2.fastq"))            
        }

        trimgalore_fastq(fastq1, fastq2 = fastq2,
                         adapter1 = NULL, adapter2 = NULL,
                         ## illumina = FALSE, nextera = FALSE, small_rna = FALSE,
                         ## minlength = 20, minqual = 20, trimN = TRUE,
                         ## retainUnpaired = TRUE, retain1length = 35,
                         ## retain2length = 35, clipR1 = NULL, clipR2  = NULL,
                         ## clip3primeR1 = NULL, clip3primeR2 = NULL,
                         ## robust_check = FALSE,
                         dest.dir = fastq_dir,
                         threads = nthread,
                         ## do.fastqc = TRUE,
                         trimgalore = "trim_galore") 

        ##system(paste("",))
        fq <- dir(fastq_dir,pattern="_trimmed.fq$",full.names=TRUE)
        fq
        for(f in fq) file.rename(from = f, to = sub("fq$","fastq",f))
        
    }
    
    if(0 && do.qc && trimming) {
        cat(">>> Running FastQC (trimmed) ... \n")
        run_fastqc2(destdir=destdir, fastq_dir=fastq_dir, nthread=nthread)
    }

    if(quant.method=="salmon") {

        ## ----------- Before running Salmon, you will have to build index first.
        dir.exists(indexdir)
        index2 = file.path(indexdir,paste0(species,"_transcripts_release99_index"))
        index2
        dir.exists(index2)
        if(!dir.exists(index2)) {
            cat(">>> Building index for species",species,"... \n")
            system(paste("mkdir -p",indexdir))
            ##build_index(species="human", kmer=31, ens_release=92, destdir=indexdir)
            build_index(species=species, kmer=31, ens_release=99, destdir=indexdir)
        } else {
            cat(">>> Found index folder at ",index2,"\n")
        }
        dir(indexdir)

        ## ----------- Run Salmon
        cat(">>> Running Salmon... \n")
        i=1
        for(i in 1:length(file_id)) {
            pgx.run_salmon(srr_id=file_id[i], library_layout=library_layout,
                           index_dir=index2, destdir=destdir,
                           fastq_dir=fastq_dir, use_trimmed_fastq=trimming,
                           other_opts="--validateMappings", nthread=nthread)
        }
    }
    if(quant.method=="kallisto") {
        ## ----------- Before running Kallisto, you will have to build index first.
        indexdir
        dir.exists(indexdir)
        kallisto_index = file.path(indexdir,paste0(species,"_transcripts.idx"))
        if(!file.exists(kallisto_index)) {
            cat(">>> Building index for species",species,"... \n")
            system(paste("mkdir -p",indexdir))
            ##build_index(species="human", kmer=31, ens_release=92, destdir=indexdir)
            ##build_index(species=species, kmer=31, ens_release=99, destdir=indexdir)
            if(species=="human") ref.genome = file.path(indexdir,"Homo_sapiens.GRCh38.cdna.all.fa.gz")
            if(species=="mouse") ref.genome = file.path(indexdir,"Mus_musculus.GRCm38.cdna.all.fa.gz")
            kallisto_index
            ref.genome
            system( paste("kallisto index -i",kallisto_index, ref.genome) )
        } else {
            cat(">>> Found index at ",kallisto_index,"\n")
        }

        ## ----------- Run Kallisto
        cat(">>> Running Kallisto... \n")
        i=1
        fq <- dir(fastq_dir,pattern=".fastq$",full.names=TRUE)
        fq <- grep("trimmed",fq,invert=TRUE,value=TRUE)
        if(trimming) fq <- dir(fastq_dir,pattern="_trimmed.fastq$",full.names=TRUE)
        kallisto_outdir <- file.path(destdir,"kallisto")
        system(paste("mkdir -p",kallisto_outdir))

        cmd = "kallisto quant"
        if(library_layout=="SINGLE") cmd = "kallisto quant --single -l 250 -s 50" ## NEED RETHINK!!!       
        cmd <- paste(cmd, "-i", kallisto_index, "-b 100 -t 32")
        cmd
        i=1
        for(i in 1:length(fq)) {
            f1 <- gsub(".*/|.fastq$|_trimmed.fastq$","",fq[i])
            out1 <- file.path(kallisto_outdir,f1)
            cmd1 <- paste(cmd,"-o",out1,fq[i] )
            cmd1
            system(cmd1)
        }
    }
    
    ## ----------- Run MultiQC
    if(do.qc) {
        cat(">>> Running MultiQC... \n")        
        ##run_multiqc(fastqc_dir=fastq_dir, salmon_dir=destdir, destdir=destdir)
        run_multiqc(fastqc_dir = file.path(destdir,"fastqc"),
                    salmon_dir = file.path(destdir,quant.method),
                    destdir = destdir)        
    }
    
    ## ----------- Run tximport
    if(quant.method=="salmon") {
        cat(">>> Running TxImport on Salmon files... \n")        
        txi <- run_tximport(srr_id = file_id,
                            species = species,
                            salmon_dir = paste0(destdir,"/salmon"),
                            countsFromAbundance = "lengthScaledTPM")
    }

    if(quant.method=="kallisto") {
        cat(">>> Running TxImport on Kallisto files... \n")               
        txi <- run_tximport_kallisto (
            srr_id = file_id,
            species = species,
            kallisto_dir = paste0(destdir,"/kallisto")
        )
    }

    ## ----------- Extract counts
    names(txi)    
    genes  <- txi$gene_counts[,2:3]
    counts <- as.matrix(txi$gene_counts[,4:ncol(txi$gene_counts),drop=FALSE])
    rownames(genes) <- rownames(counts) <- txi$gene_counts$ENSEMBL
    
    ## ----------- Collapse multiple probes to single gene by summing up counts
    counts <- apply(counts, 2, function(x) tapply(x,genes$SYMBOL,sum))
    genes  <- genes[match(rownames(counts),genes$SYMBOL),]
    rownames(genes) <- genes$SYMBOL
    counts <- round(counts, digits=3)
    
    ## ----------- Save counts as text files for Omics Playground
    colnames(genes) <- c("gene_name","gene_title")
    system(paste0("mkdir -p ",destdir,"/files_csv"))
    cat(">>> Writing CSV files to",paste0(destdir,"/files_csv/ \n"))
    write.csv(counts, file=paste0(destdir,"/files_csv/counts.csv"))
    write.csv(genes,  file=paste0(destdir,"/files_csv/genes.csv"))
    samples <- data.frame(sample=colnames(counts), group=NA, phenotype1=NA)  ## empty template
    write.csv(samples,  file=paste0(destdir,"/files_csv/samples.csv"), row.names=FALSE)

    cat(">>> done!\n")            
}


pgx.run_salmon <- function (srr_id, library_layout = c("SINGLE", "PAIRED"), index_dir, 
    destdir, fastq_dir, use_trimmed_fastq = FALSE, other_opts = NULL, nthread) 
{
    if (!dir.exists(paste0(destdir, "/salmon"))) {
        system(paste0("mkdir ", destdir, "/salmon"))
    }
    library_layout <- match.arg(library_layout, c("SINGLE", "PAIRED"))
    if (library_layout == "SINGLE") {
        if (use_trimmed_fastq) {
            system(paste0("salmon quant -i ", index_dir, " -p ", 
                nthread, " ", other_opts, " -l A -r ", fastq_dir, 
                "/", srr_id, "_trimmed.fastq -o ", destdir, "/salmon/", 
                srr_id, "_transcripts_quant"))
        }
        else {
            system(paste0("salmon quant -i ", index_dir, " -p ", 
                nthread, " ", other_opts, " -l A -r ", fastq_dir, 
                ##"/", srr_id, "_pass.fastq -o ", destdir, "/salmon/",
                "/", srr_id, ".fastq -o ", destdir, "/salmon/",                 
                srr_id, "_transcripts_quant"))
        }
    }
    else {
        if (use_trimmed_fastq) {
            system(paste0("salmon quant -i ", index_dir, " -p ", 
                nthread, " ", other_opts, " -l A -1 ", fastq_dir, 
                "/", srr_id, "_trimmed_1.fastq ", "-2 ", fastq_dir, 
                "/", srr_id, "_trimmed_2.fastq -o ", destdir, 
                "/salmon/", srr_id, "_transcripts_quant"))
        }
        else {
            system(paste0("salmon quant -i ", index_dir, " -p ", 
                nthread, " ", other_opts, " -l A -1 ", fastq_dir, 
                ##"/", srr_id, "_pass_1.fastq ", "-2 ", fastq_dir, 
                ##"/", srr_id, "_pass_2.fastq -o ", destdir, "/salmon/", 
                "/", srr_id, "_1.fastq ", "-2 ", fastq_dir, 
                "/", srr_id, "_2.fastq -o ", destdir, "/salmon/", 
                srr_id, "_transcripts_quant"))
        }
    }
    if (file.exists(paste0(destdir, "/salmon/", srr_id, "_transcripts_quant/quant.sf"))) {
        system(paste("cat ", destdir, "/salmon/", srr_id, "_transcripts_quant/quant.sf", 
            "| sed -E 's/\\.[0-9]+//' > ", destdir, "/salmon/", 
            srr_id, "_transcripts_quant", "/", srr_id, "_quant_new.sf", 
            sep = ""))
    }
    else {
        cat("quant.sf doesn't exist. Processing next sample.")
    }
}

run_tximport_kallisto <- function (srr_id, species = c("human", "mouse", "rat"), kallisto_dir)
{
    species <- match.arg(species, c("human", "mouse", "rat"))
    edb = function(species) {
        if (species == "human") {
            GenomicFeatures::transcripts(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, 
                columns = c("tx_id", "gene_id", "gene_name"), 
                return.type = "DataFrame")
        }
        else if (species == "mouse") {
            GenomicFeatures::transcripts(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79, 
                columns = c("tx_id", "gene_id", "gene_name"), 
                return.type = "DataFrame")
        }
        else if (species == "rat") {
            GenomicFeatures::transcripts(EnsDb.Rnorvegicus.v79::EnsDb.Rnorvegicus.v79, 
                columns = c("tx_id", "gene_id", "gene_name"), 
                return.type = "DataFrame")
        }
        else {
            return(NULL)
        }
    }
    gene_ensembl = function(species) {
        if (species == "human") {
            require(org.Hs.eg.db)
            return(org.Hs.eg.db::org.Hs.eg.db)
        }
        else if (species == "mousee") {
            require(org.Mm.eg.db)
            return(org.Mm.eg.db::org.Mm.eg.db)
        }
#        else if (species == "rat") {
#            require(org.Rn.eg.db)            
#            return(org.Rn.eg.db::org.Rn.eg.db)
#        }
        else {
            return(NULL)
        }
    }
    assign("Tx.ensemble", edb(species))
    Tx.ensemble <- get("Tx.ensemble")
    tx2gene <- Tx.ensemble[, c(1, 2)]
    files <- file.path(paste0(kallisto_dir, "/", srr_id, "/abundance.tsv"))
    names(files) <- srr_id
    file.exists(files)
    
    cat("generating counts table\n")
    txi.t <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, 
                               txOut = TRUE,
                               ##importer = utils::read.delim,
                               dropInfReps = TRUE)
    txi.g <- tximport::summarizeToGene(txi.t, tx2gene,
                                       ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
    gene_counts <- txi.g$counts
    gene_counts[is.na(gene_counts)] <- 0
    colnames(gene_counts) <- srr_id
    transcript_counts <- txi.t$counts
    transcript_counts[is.na(transcript_counts)] <- 0
    colnames(transcript_counts) <- srr_id
    annot_genes <- AnnotationDbi::select(gene_ensembl(species), 
        keys = rownames(gene_counts), columns = c("SYMBOL", "SYMBOL", 
            "GENENAME"), keytype = "ENSEMBL")
    annot_genes2 <- annot_genes[match(rownames(gene_counts), 
        annot_genes[, 1]), , drop = FALSE]
    gene_counts <- cbind(annot_genes2, gene_counts)
    counts <- list(gene_counts = gene_counts, transcript_counts = transcript_counts, 
        tximport_gene_data = txi.g, tximport_transcript_data = txi.t)
    return(counts)
}



##-------------------------------------------------------------------------
## https://rdrr.io/github/anilchalisey/rseqR/src/R/trim_fastq.R
##-------------------------------------------------------------------------

#' An R-based wrapper for Trim Galore!
#'
#' @description Run the Trim Galore! tool
#'
#' @details This script runs the Trim Galore! tool and requires installation
#' of both Cutadapt and Trim Galore!  It is essential that Cutadapt is in the
#' executable path otherwise this tool will not work.
#'
#' @param fastq1 a character vector indicating the read files to be trimmed.
#' @param fastq2 (optional) a character vector indicating read files to be
#' trimmmed.  If specified, it is assumed the reads are paired, and this vector
#' MUST be in the same order as those listed in \code{fastq1}.  If \code{NULL}
#' then it is assumed the reads are single-end.
#' @param adapter1 a character string specifying the adapter sequence to be
#' trimmed. If not specified explicitly, Trim Galore will try to auto-detect
#' whether the Illumina universal, Nextera transposase or Illumina small RNA
#' adapter sequence was used. Also see \code{illumina}, \code{nextera} and
#' \code{small_rna} options. If no adapter can be detected within the first 1
#' million sequences of the first file specified Trim Galore defaults to
#' \code{illumina}.
#' @param adapter2 a character string specifying an optional adapter sequence to
#' be trimmed off read 2 of paired-end files. This option requires paired-end
#' reads.
#' @param illumina a logical specifying that the adapter sequence to be trimmed
#' is the first 13bp of the Illumina universal adapter AGATCGGAAGAGC instead of
#' the default auto-detection of adapter sequence.  Default: \code{FALSE}
#' @param nextera adapter sequence to be trimmed is the first 12bp of the
#' Nextera adapter CTGTCTCTTATA instead of the default auto-detection of adapter
#' sequence.
#' @param small_rna a logical specifying that the adapter sequence to be trimmed
#' is the first 12bp of the Illumina Small RNA 3' Adapter TGGAATTCTCGG instead
#' of the default auto-detection of adapter sequence.  Selecting to trim
#' smallRNA adapters will also lower the \code{length} value to 18bp. If the
#' smallRNA libraries are paired-end then \code{adapter2} will be set to the
#' Illumina small RNA 5' adapter automatically (GATCGTCGGACT) unless
#' \code{adapter2} had been defined explicitly.
#' @param minlength an integer value; reads that become shorter than this length
#' as a result of either quality or adapter trimming are discarded. A value of 0
#' effectively disables this behaviour.  Default: 20 bp.  For paired-end files,
#' both reads of a read-pair need to be longer than bp to be printed out to
#' validated paired-end files. If only one read became too short there is the
#' possibility of keeping such unpaired single-end reads (see
#' \code{retain_unpaired}). Default pair-cutoff: 20 bp.
#' @param minqual an integer value specifying the quality threshold below which
#' to trim low-quality ends from reads in addition to adapter removal. Default
#' Phred score: 20.
#' @param trimN a logical specifying whether to remove Ns from the end of reads.
#' @param retainUnpaired a logical.  If only one of the two paired-end reads
#' become too short, the longer read will be written to either .unpaired_1.fq or
#' .unpaired_2.fq output files. The length cutoff for unpaired single-end reads
#' is governed by the parameters \code{retain1length} and \code{retain2length}.
#' Default: ON.
#' @param retain1length an integer.  Unpaired single-end read length cutoff
#' needed for read 1 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp.
#' @param retain2length an integer.  Unpaired single-end read length cutoff
#' needed for read 2 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp
#' @param clipR1 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 1 (or single-end reads). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clipR2 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 2 (paired-end reads only). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clip3primeR1 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param clip3primeR2 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param robust_check a logical indicating whether to check that the paired
#' files specified are matching and have equal numbers of reads.  Default:
#' \code{FALSE}
#' @param trimgalore a character string specifying the path to the trimgalore executable.
#' On Unix systems, if the executable is in \code{$PATH}, then it may be left as
#' the default. If it is not in \code{$PATH}, then the absolute path should be given.
#` If using the WSL on Windows 10, then the path must be the absolute path in WSL,
#' unless the system has been set up as described in the vignette.
#' @param dest.dir a character string specifying the output directory.  If NULL
#' a directory named "TRIMMED_FASTQC" is created in the current working directory
#' [DEFAULT = NULL].
#' @param threads an integer value indicating the number of parallel threads to
#' be used by FastQC. [DEFAULT = maximum number of available threads - 1].
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel

trimgalore_fastq <- function(fastq1, fastq2 = NULL, adapter1 = NULL, adapter2 = NULL,
                       illumina = FALSE, nextera = FALSE, small_rna = FALSE,
                       minlength = 20, minqual = 20, trimN = TRUE,
                       retainUnpaired = TRUE, retain1length = 35,
                       retain2length = 35, clipR1 = NULL, clipR2  = NULL,
                       clip3primeR1 = NULL, clip3primeR2 = NULL,
                       robust_check = FALSE, dest.dir = NULL,
                       threads = NULL, do.fastqc=FALSE,
                       trimgalore = "trim_galore")
{
    
    doParallel::registerDoParallel(cores=2)
    ##foreach(i=1:3) %dopar% sqrt(i)

    cmd <- trimgalore

    if (is.null(fastq2)) {
        paired = FALSE
    } else {
        if(length(fastq1) != length(fastq2)) {stop("The number of forward and reverse reads do not match")}
        if (robust_check) {
            fq1lengths <- lapply(fastq1, function (x) {sprintf("gzip -cd %s | wc -l", x)})
            fq2lengths <- lapply(fastq1, function (x) {sprintf("gzip -cd %s | wc -l", x)})
            fq1lengths <- unlist(lapply(fq1lengths, run_cmd, intern = TRUE))
            fq2lengths <- unlist(lapply(fq2lengths, run_cmd, intern = TRUE))
            if (!identical(fq1lengths, fq2lengths)) {
                stop("One or more of the forward and reverse reads pairs have differing number of reads.\n",
                     "Are you sure the two lists are in the correct paired order?")
            }
        }
        paired = TRUE
    }

    if (is.null(dest.dir)) dest.dir <- "TRIMMED_FASTQC"
    dir.create(dest.dir, showWarnings = FALSE)

    cmd <- paste(cmd,
                 "-q", minqual, "--length", minlength, "-o", dest.dir)
    if(do.fastqc) cmd <- paste(cmd, "--fastqc" )
    if(!is.null(adapter1)) cmd <- paste(cmd, "--adapter", adapter1)
    if(!is.null(adapter2)) cmd <- paste(cmd, "--adapter2", adapter2)
    if(illumina) cmd <- paste(cmd, "--illumina")
    if(nextera) cmd <- paste(cmd, "--nextera")
    if(small_rna) cmd <- paste(cmd, "--small_rna")
    if(trimN) cmd <- paste(cmd, "--trim-n")
    if(!is.null(clipR1)) cmd <- paste(cmd, "--clip_R1", clipR1)
    if(!is.null(clipR2)) cmd <- paste(cmd, "--clip_R2", clipR2)
    if(!is.null(clip3primeR1)) cmd <- paste(cmd, "--three_prime_clip_R1", clip3primeR1)
    if(!is.null(clip3primeR2)) cmd <- paste(cmd, "--three_prime_clip_R2", clip3primeR2)
    if(paired) {
        cmd <- paste(cmd, "--paired")
        if(retainUnpaired) {
            cmd <- paste(cmd, "--retain_unpaired",
                         "-r1", retain1length,
                         "-r2", retain2length)
        }
    }

    ##give_note("\nRemoving adapters and performing quality trimming...\n\n")
    cat("\nRemoving adapters and performing quality trimming...\n\n")  

    if (is.null(threads)) {threads <- parallel::detectCores() - 1}

    if (threads > 1) {
        cl <- parallel::makeCluster(threads)
        doParallel::registerDoParallel(cl)

        if (paired) {
            foreach (i = seq_along(fastq1)) %dopar% {
                tgcmd <- sprintf("%s %s %s", cmd, fastq1[[i]], fastq2[[i]])
                run_cmd <- function(cmd, intern = FALSE) {
                    if (.Platform$OS.type != "windows") {
                        system(command = cmd, intern = intern)
                    } else {
                        shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
                    }
                }
                run_cmd(tgcmd)
            }
        } else {
            foreach (i = seq_along(fastq1)) %dopar% {
                tgcmd <- sprintf("%s %s", cmd, fastq1[[i]])
                run_cmd <- function(cmd, intern = FALSE) {
                    if (.Platform$OS.type != "windows") {
                        system(command = cmd, intern = intern)
                    } else {
                        shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
                    }
                }
                run_cmd(tgcmd)
            }

        }
        parallel::stopCluster(cl)
    } else {
        if (paired) {
            for (i in seq_along(fastq1)) {
                tgcmd <- sprintf("%s %s %s", cmd, fastq1[[i]], fastq2[[i]])
                run_cmd(tgcmd)
            }
        } else {
            for (i in seq_along(fastq1)) {
                tgcmd <- sprintf("%s %s", cmd, fastq1[[i]])
                run_cmd(tgcmd)
            }
        }
    }

    if (paired) {
        trimmed.files <- list.files(path = dest.dir, pattern = "*val", full.names = TRUE)
        lapply(trimmed.files, function(x) {
            file.to = sub(".fq",".fastq",gsub("_val_[0-9]", "_trimmed", x))
            file.rename(from = x, to = file.to)
        })
        unpaired.files <- list.files(path = dest.dir, pattern = "*unpaired", full.names = TRUE)
        lapply(unpaired.files, function(x) {
            file.to = sub(".fq",".fastq",gsub("_unpaired_[0-9]", "_unpaired", x))
            file.rename(from = x, to = file.to)
        })
    } else {
        trimmed.files <- list.files(path = dest.dir, pattern = "*val", full.names = TRUE)
        lapply(trimmed.files, function(x) {
            file.to = sub(".fq",".fastq",gsub("_val_[0-9]", "_trimmed", x))
            file.rename(from = x, to = file.to)
        })
    }
}
