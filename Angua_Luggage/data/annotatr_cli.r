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
    setwd(aa_dir)
    if (file.exists("ORFs.rdata") && file.exists("grl.rdata")) {
    print("ORFs already found, skipping ORF step.")
    } else {
    all_fastas <- setup_files(contig_dir, "fasta")
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
            print(paste0("No ORFs of sufficient length found for ", filename, "."))
            all_fastas <- setdiff(all_fastas, fasta)
            next
            }
        }
    all_ORFs <- unlistGrl(all_grls)
    setwd(aa_dir)
    save(all_ORFs, file = "ORFs.rdata")
    save(all_grls, file = "grl.rdata")
    success_file <- file("hits.txt")
    writeLines(all_fastas, success_file)
    close(success_file)
    }
    }

parse_pfam_json <- function(dir, ORFs_file) {
    all_jsons <- setup_files(dir, filetype = "json")
    if(file.exists("pfam_grl.rdata") && file.exists("pfam_dfs.rdata")) {
    return(list("Skipping pfam-parse step: already complete."))
    } else {
    ORFs <- load(file = ORFs_file)
    pfam_grl <- list()
    pfam_dfs <- list()
    for(filename in all_jsons){
        contigs_name_vec <- filename %>%
                           tools::file_path_sans_ext() %>%
                           str_split("_") %>%
                           unlist()
        contigs_file_name <- paste(contigs_name_vec, collapse = "_")
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
                pfam_grl[[contigs_file_name]] <- pfam_gr
                pfam_dfs[[contigs_file_name]] <- pfam_df
            options(warn) } else {
            next }
            }
    save(pfam_grl, file = "pfam_grl.rdata")
    save(pfam_dfs, file = "pfam_dfs.rdata")
    }
    }

#Coverage with aid of https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/#convert-sequencing-depth-to-bedgraph-format
generate_orf_plots <- function(grl_file, fasta_dir, out_dir, pfam_file, pfam_df_file, bedgraph_dir) { 
    print(grl_file)
    print(pfam_file)
    print(pfam_df_file)
    load(file = grl_file)
    orfs <- unlist(all_grls)
    load(file = pfam_file)
    load(file = pfam_df_file)
    file_end <- "_plot.jpg"
    file_names <- setup_files(fasta_dir, "fasta")
    log <- list()
    setwd(out_dir)
    
    i <- 1
    for(contig_file_name in names(pfam_dfs)){ 
        grange <- all_grls[[i]]
        get_prots <- pfam_dfs[[contig_file_name]]
        contig_full_file_name <- paste0(contig_file_name, ".fasta")
        current_contig_dir <- paste(contig_file_name, "plots", sep = "_")
        dir.create(current_contig_dir, showWarnings = FALSE)
        setwd(current_contig_dir)
        all_con_bios <- fasta_dir %>%
                        paste(contig_full_file_name, sep = "/") %>%
                        FaFile() %>%
                        scanFa()
        bg_file <- paste0(bedgraph_dir, paste0("/", contig_file_name, ".bedGraph.gz"))
        bedgraph_dt <- fread(bg_file, col.names = c('chromosome', 'start', 'end', 'value'))        
        current_orfs <- pfam_grl[[contig_file_name]]

        for(orig_contig_name in unique(bedgraph_dt[["chromosome"]])){
                current_contig <- grange[grepl(orig_contig_name, seqnames(grange))]  
                if(length(current_contig) <= 0) {
                    print(paste("No pfam hits for:", contig_file_name, orig_contig_name, sep=" "))
                    next
                    }
                orig_contig <- all_con_bios[grepl(orig_contig_name, names(all_con_bios))]
                bedgraph_dt_contig <- filter(bedgraph_dt, chromosome == orig_contig_name)
                if((length(orig_contig) > 0) & (length(get_prots) > 0)) {
                    orig_gr <- GRanges(seqnames = Rle(orig_contig_name, 1),
                                       ranges = IRanges(start = 1, width = width(orig_contig), names = c("orig")))
                    seq_shorter <- orig_contig_name %>%
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
                        proteins <- filter(get_prots, grepl(identifier, seq$name))
                        if(length(proteins) <= 0) {
                            d <- data.frame(protein = c("NA"))
                            } else {
                            d <- data.frame(protein = proteins$desc)
                            }
                        grid.text(paste(d$protein, collapse = "\n"), draw = TRUE)
                    }
                    
                    contig_chr <- GRanges(seqnames = Rle(orig_contig_name, length(seqnames(current_contig))),
                                  ranges(current_contig, use.mcols=TRUE), strand=strand(current_contig))
                    
                    dtrack <- DetailsAnnotationTrack(range = contig_chr, 
                                                     name = seq_title, 
                                                     id = ORF_names,
                                                     fun = details)
                    
                    displayPars(dtrack) <- list(fontcolor.item = "black", 
                                                col = "darkblue", fill = "lightblue", detailsBorder.col = "blue",
                                                showFeatureId = TRUE, background.title = "darkgray")
        
                    gtrack <- GenomeAxisTrack(orig_gr, littleTicks = TRUE, cex = 1)
                    
                    datrack <- DataTrack(range = bedgraph_dt_contig, genome = orig_gr,
                                         chromosome = orig_contig_name,
                                         name = "Coverage")
                    
                    datrack2 <- DataTrack(range = bedgraph_dt_contig, genome = orig_gr,
                                         chromosome = orig_contig_name,
                                         name = "Line") 
                    
                    displayPars(datrack) <- list(type = "gradient", 
                                                         gradient = 
                                                         c("mintcream", "lightskyblue1", "paleturquoise3", "lightsalmon", 
                                                           "orange", "orangered1"),
                                                         background.title = "darkgray", cex.axis = 1)
                    displayPars(datrack2) <- list(type = "a", alpha.title = 0, col= "black")                 
                    
                    otrack <- OverlayTrack(trackList=list(datrack, datrack2), 
                                           name="Coverage", background.title = "darkgray")
                    
                    jpeg_name <- paste0(orig_contig_name, file_end)
                    jpeg(jpeg_name, width = 700, height = 500)
                    
                    tryCatch( {
                        plotTracks(list(dtrack,
                                        gtrack,
                                        otrack), 
                                        add53= TRUE, 
                                   stacking = "squish", stackHeight = 0.9, add = TRUE)
                    },
                        error=function(e) {
                        message(paste0('One of the plots broke: ', contig_file_name, ": ", orig_contig_name))
                        print(e)
                    },
                        warning = function(w) {
                        print("Warning")
                        print(w)
                    }
                    )
                    dev.off()
                } else {
                print(paste0("Unable to plot ", orig_contig_name, " suggest manual review."))
                }
        }
        i <- i + 1
        setwd("..")
    }
    }