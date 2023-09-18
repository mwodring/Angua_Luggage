import_packages <- function(){
    suppressPackageStartupMessages({
        library(ORFik)
        library(GenomicRanges)
        #library(GenomicFeatures)
        library(Biostrings)
        library(jsonlite)
        library(dplyr)
        library(Gviz)
        library(plyranges)
        library(stringr)
	    library(gridExtra)
	    library(data.table)
    })
    }

import_packages()

setup_files <- function(dir, filetype) {
    file_pattern <- paste0( "\\" , ".", filetype, "$")
    all_files <- list.files(dir, pattern = file_pattern, ignore.case=TRUE)
    setwd(dir)
    return(all_files)
    }

ORF_from_fasta <- function(contig_dir, aa_dir, nt_dir, ORF_min_len) {
    import_packages()
    all_fastas <- setup_files(contig_dir, "fasta")
    if (file.exists("ORFs.rdata") && file.exists("grl.rdata")) {
    log <- "ORFs already found, skipping ORF step."
    return(list(log))
    } else {
    log_list <- list()
    all_grls <- GRangesList()
    for (fasta in all_fastas){
        wd <- getwd()
        filename <- tools::file_path_sans_ext(fasta)
        fa_filepath <- file.path(wd, tools::file_path_sans_ext(fasta))
        new_filename <- paste(filename, "ORF", ".fasta", sep = "_")
        aa_filename <- paste(filename, "ORF", "aa", ".fasta", sep = "_")
        fa <- FaFile(fasta)
        seq <- readDNAStringSet(fasta, format="fasta", use.names=TRUE) 
        ORFs <- findORFsFasta(seq, startCodon = "ATG", stopCodon = stopDefinition(1), longestORF = TRUE, minimumLength = ORF_min_len)
        if (length(ORFs) >= 1){
            gr <- GRanges(ORFs)
            extracted_ORFs <- getSeq(fa, gr)
            names(gr) <- paste0("ORF_", seq.int(length(gr)), "_", seqnames(gr))
            names(extracted_ORFs) <- names(gr)
            setwd(nt_dir)
            writeXStringSet(extracted_ORFs, new_filename, append=FALSE,
                            compress=FALSE, format="fasta")
            grl_ORFs <- GRangesList(gr)
	        suppressWarnings({all_grls <- append(all_grls, grl_ORFs, after = length(all_grls))})
            export.bed12(grl_ORFs, paste0(filename, ".bed12"))
            ORFs_aa <- Biostrings::translate(extracted_ORFs)
            setwd(aa_dir)
            writeXStringSet(ORFs_aa, aa_filename, append=FALSE, compress=FALSE, format="fasta")
            setwd(contig_dir)
            } else {
            log_list <- append(log_list, paste0("No ORFs of sufficient length found for ", filename, "."))
            all_fastas <- setdiff(all_fastas, fasta)
            next
            }
        }
    all_ORFs <- unlistGrl(all_grls)
    save(all_ORFs, file = "ORFs.rdata")
    save(all_grls, file = "grl.rdata")
    return(list("log"))
    }
    }

parse_pfam_json <- function(dir, ORFs_file) {
    all_jsons <- setup_files(dir, filetype = "json")
    if(file.exists("pfam_grl.rdata") && file.exists("pfam_df.rdata")) {
    return(list("Skipping pfam-parse step: already complete."))
    } else {
    ORFs <- load(file = ORFs_file)
    pfam_grl <- GRangesList()
    for(filename in all_jsons){
        pfam_df <- fromJSON(filename, simplifyDataFrame = TRUE)
        seq_names <- pfam_df$seq$name
        if(!(is.null(seq_names))){
            tsv_df <- data.frame(orf = seq_names, 
                                protein = pfam_df$name, 
                                accession = pfam_df$acc)
            write.csv(tsv_df, paste0(tools::file_path_sans_ext(filename), ".csv"), row.names = FALSE)
            warn <-options(warn=-1)
                seq_to <- as.numeric(unlist(pfam_df$seq$to))
                seq_from <- as.numeric(unlist(pfam_df$seq$from))
                pfam_gr <- GRanges(seqnames = Rle(seq_names),
                           ranges = IRanges(seq_from, end = seq_to, names = pfam_df$name))
                new_pfam_grl <- GRangesList(pfam_gr)
                pfam_grl <- append(pfam_grl, new_pfam_grl, after = length(pfam_grl))
            options(warn) } else {
            next }
            }
    save(pfam_grl, file = "pfam_grl.rdata")
    save(pfam_df, file = "pfam_df.rdata")
    }
    }

#Coverage with aid of https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/#convert-sequencing-depth-to-bedgraph-format
generate_orf_plots <- function(grl_file, fasta_dir, out_dir, pfam_file, pfam_df_file, bedgraph_dir) { 
    load(file = grl_file)
    load(file = pfam_file)
    load(file = pfam_df_file)
    pfam <- unlist(pfam_grl)
    listed_contigs <- as.list(all_grls)
    file_end <- "_plot.jpg"
    file_names <- setup_files(fasta_dir, "fasta")
    log <- list()
    setwd(out_dir)
    
    i <- 1
    for(filename in file_names){
        sample_name_vec <- filename %>%
                           tools::file_path_sans_ext() %>%
                           str_split("_") %>%
                           unlist()
        sample_name <- paste(sample_name_vec[1:2], collapse = "_") 
        all_con_bios <- fasta_dir %>%
                        paste(filename, sep = "/") %>%
                        FaFile() %>%
                        scanFa()
        bg_file <- paste0(bedgraph_dir, paste0("/", sample_name, ".bedGraph.gz"))
        bedgraph_dt <- fread(bg_file, col.names = c('chromosome', 'start', 'end', 'value'))     
        
        for(grange in listed_contigs){
            current_contig_dir <- paste(tools::file_path_sans_ext(filename), "plots", sep = "_")
            dir.create(current_contig_dir, showWarnings = FALSE)
            setwd(current_contig_dir)
            
            seq_names <- as.list(seqnames(grange))
            for(seq_name in seq_names) {
                bedgraph_dt_contig <- filter(bedgraph_dt, chromosome == seq_name)
                current_contig <- filter(grange, seqnames == seq_name)
                get_prots <- filter(pfam, seqnames %in% names(current_contig))
                orig_contig <- all_con_bios[grepl(seq_name, names(all_con_bios))]
                
                if(!(is.null(orig_contig)) & (length(get_prots) > 0)) {
                    orig_gr <- GRanges(seqnames = Rle(seq_name, 1),
                                       ranges = IRanges(start = 1, width = width(orig_contig), names = c("orig")))
                    seq_shorter <- seq_name %>%
                                   str_split("_", n= Inf, simplify = FALSE) %>%
                                   unlist() %>%
                                   setdiff(c("TRINITY"))
                    seq_title <- paste(seq_shorter[2:length(seq_shorter)], collapse = "")
                    ORF_names <- character()
                    for(name in names(current_contig)) {
                        current_name <- name %>%
                                        str_split("_", n= Inf, simplify = FALSE) %>%
                                        unlist()
                        current_name_str <- paste(current_name[1:2], collapse = "_")
                        ORF_names <- append(ORF_names, current_name_str, after = length(ORF_names))
                    }
                    
                    options(ucscChromosomeNames=FALSE)
                    details <- function(identifier, ...) {
                        proteins <- get_prots[grepl(identifier, seqnames(get_prots))]
                        if(length(proteins) <= 0) {
                            d <- data.frame(protein = c("NA"))
                            } else {
                            d <- data.frame(protein = names(proteins))
                            }
                        grid.text(paste(d$protein, collapse = "\n"), draw = TRUE)
                    }
                    
                    dtrack <- DetailsAnnotationTrack(range = current_contig, 
                                                     name = seq_title, 
                                                     id = ORF_names, 
                                                     fun = details)
                                                     
                    displayPars(dtrack) <- list(fontcolor.item = "black", 
                                            col = "darkblue", fill = "lightblue", detailsBorder.col = "blue",
                                            showFeatureId = TRUE, background.title = "darkgray")
        
                    gtrack <- GenomeAxisTrack(orig_gr, littleTicks = TRUE, cex = 1)
                
                    datrack <- DataTrack(range = bedgraph_dt_contig, genome = orig_gr,
                                         chromosome = seq_name,
                                         name = "Coverage")

                    datrack2 <- DataTrack(range = bedgraph_dt_contig, genome = orig_gr,
                                         chromosome = seq_name,
                                         name = "Line") 
                                         
                    displayPars(datrack) <- list(type = "gradient", 
                                                         gradient = 
                                                         c("mintcream", "lightskyblue1", "paleturquoise3", "lightsalmon", "orange", "orangered1"),
                                                         background.title = "darkgray", cex.axis = 1)
                    displayPars(datrack2) <- list(type = "a", alpha.title = 0, col= "black")                 
                    otrack <- OverlayTrack(trackList=list(datrack, datrack2), 
                                           name="Coverage", background.title = "darkgray")
                    
                    jpeg_name <- paste0(seq_name, file_end)
                    jpeg(jpeg_name, width = 700, height = 500)
                    
                    tryCatch( {
                        plotTracks(list(dtrack, gtrack, otrack), add53= TRUE, 
                                   stacking = "squish", stackHeight = 0.9, add = TRUE)
                    },
                        error=function(e) {
                        message(paste0('One of the plots broke: ', filename, ": ", seq_name))
                        print(e)
                    },
                        warning = function(w) {
                        print("Warning")
                        print(w)
                    }
                    )
                    dev.off()
                } else {
                log <- append(log, paste0("Unable to plot ", seq_name, " suggest manual review."), after=length(log)) 
                log <- append(log, toString(get_prots), after=length(log))
                log <- append(log, toString(bedgraph_dt_contig), after=length(log))
                }
            }
        setwd("..")
        }
    }
    return(list(log))
    }