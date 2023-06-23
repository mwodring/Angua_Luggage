import_packages <- function(){
    suppressPackageStartupMessages({
        library(ORFik)
        library(GenomicRanges)
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

setup_files <- function(dir, new_dir, filetype) {
    file_pattern <- paste0( "\\" , ".", filetype, "$")
    all_files <- list.files(dir, pattern = file_pattern, ignore.case=TRUE)
    setwd(dir)
    if(new_dir != ""){
    dir.create(new_dir, showWarnings = FALSE)
        }
    return(all_files)
    }

loop_fasta_files <- function(dir, new_dir, ORF_min_len) {
    import_packages()
    all_fastas <- setup_files(dir, new_dir, "fasta")
    log_list <- list()
    all_grls <- GRangesList()
    for (fasta in all_fastas)
        {
        filename <- tools::file_path_sans_ext(fasta)
        new_filepath <- paste(filename, "ORF", ".fasta", sep = "_")
        aa_filepath <- paste(filename, "ORF", "aa", ".fasta", sep = "_")
        fa <- FaFile(fasta)
        ORFs <- findORFsFasta(fasta, startCodon = "ATG", stopCodon = stopDefinition(1), longestORF = TRUE, minimumLength = ORF_min_len, is.circular = FALSE)
        #file.remove(list.files(pattern = "*.fai"))
        if (length(ORFs) >= 1){
            names(ORFs) <- paste("ORF", seq.int(length(ORFs)), seqnames(ORFs), sep = "_")
            extracted_ORFs <- getSeq(fa, ORFs)
            names(extracted_ORFs) <- names(ORFs)
            setwd(new_dir)
            contig_dir <- paste(filename, "ORFs", sep = "_")
            dir.create(contig_dir, showWarnings = FALSE)
            setwd(contig_dir)
            writeXStringSet(extracted_ORFs, new_filepath, append=FALSE,
                            compress=FALSE, format="fasta")
            grl_ORFs <- GRangesList(ORFs)
	        suppressWarnings({all_grls <- append(all_grls, grl_ORFs, after = length(all_grls))})
            export.bed12(grl_ORFs, paste0(filename, ".bed12"))
            ORFs_aa <- Biostrings::translate(extracted_ORFs)
            writeXStringSet(ORFs_aa, aa_filepath, append=FALSE, compress=FALSE, format="fasta")
            setwd("..")
            setwd("..")
            } else {
            log_list <- append(log_list, paste0("No ORFs of sufficient length found for ", filename, "."))
            all_fastas <- setdiff(all_fastas, fasta)
            next
            }
        }
    all_ORFs <- unlistGrl(all_grls)
    return(list("log" = log_list, "ORFs" = all_ORFs, "grl" = all_grls, "files" = all_fastas))
    }

loop_json_files <- function(dir, ORFs) {
    all_jsons <- setup_files(dir, "", filetype = "json")
    pfam_grl <- GRangesList()
    for(filename in all_jsons)
        {
        pfam_df <- fromJSON(filename, simplifyDataFrame = TRUE)
        seq_names <- pfam_df$seq$name
	if(!(is.null(seq_names))){
        tsv_df <- data.frame(orf = seq_names, 
                             protein = pfam_df$name, 
                             accession = pfam_df$acc)
	    write.csv(tsv_df, paste0(tools::file_path_sans_ext(filename), ".txt"), row.names = FALSE)
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
    return(list("pfam" = pfam_grl, "df" = pfam_df))
    }

#Coverage with aid of https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/#convert-sequencing-depth-to-bedgraph-format
generate_orf_plots <- function(grl, fasta_dir, file_names, out_dir, pfam_full, pfam_df, bedgraph_dir) { 
    pfam <- unlist(pfam_full)
    listed_contigs <- as.list(grl)
    file_end <- "_plot.jpg"
    dir.create(out_dir, showWarnings = FALSE)
    setwd(out_dir)
    
    i <- 1
    for(grange in listed_contigs)
        {
        filename <- file_names[[i]]
        sample_name_vec <- filename %>%
                       str_split("_") %>%
                       unlist()
        sample_name <- paste(sample_name_vec[1], collapse = "")
                       
        bg_file <- filename %>%
                   tools::file_path_sans_ext() %>%
                   paste0(bedgraph_dir, "/", ., "/", sample_name, "/", paste0(sample_name, "_sort", ".bedGraph.gz"))
        bedgraph_dt <- fread(bg_file, col.names = c('chromosome', 'start', 'end', 'value'))     
        
        all_con_bios <- fasta_dir %>%
                        paste(filename, sep = "/") %>%
                        FaFile() %>%
                        scanFa()
        
        seq_names <- as.list(seqnames(grange))
        i <- i + 1
        
        current_contig_dir <- paste(tools::file_path_sans_ext(filename), "plots", sep = "_")
        dir.create(current_contig_dir, showWarnings = FALSE)
        setwd(current_contig_dir)
        
        for(seq_name in seq_names) 
            {
            bedgraph_dt_contig <- filter(bedgraph_dt, chromosome == seq_name)
            current_contig <- filter(grange, seqnames == seq_name)
            get_prots <- filter(pfam, seqnames %in% names(current_contig))
            orig_contig <- all_con_bios[grepl(seq_name, names(all_con_bios))]
            if(!(is.null(orig_contig)) & length(get_prots) > 0) {
                orig_gr <- GRanges(seqnames = Rle(seq_name, 1),
                                   ranges = IRanges(start = 1, width = width(orig_contig), names = c("orig")))
                seq_shorter <- seq_name %>%
                               str_split("_", n= Inf, simplify = FALSE) %>%
                               unlist() %>%
                               setdiff(c("TRINITY"))
                seq_title <- paste(seq_shorter[2:length(seq_shorter)], collapse = "")
                ORF_names <- character()
                for(name in names(current_contig))
                    {current_name <- name %>%
                                     str_split("_", n= Inf, simplify = FALSE) %>%
                                     unlist()
                     current_name_str <- paste(current_name[1:2], collapse = "_")
                     ORF_names <- append(ORF_names, current_name_str, after = length(ORF_names))
                    }

                details <- function(identifier, ...) {
                proteins <- get_prots[grepl(identifier, seqnames(get_prots))]
                if(length(proteins) <= 0) {
                       d <- data.frame(protein = c("NA"))} else {
		               d <- data.frame(protein = names(proteins))}
                grid.text(paste(d$protein, collapse = "\n"), draw = TRUE)
		}
                
                options(ucscChromosomeNames=FALSE)
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
                                                     gradient = c("mintcream", "lightskyblue1", "paleturquoise3", "lightsalmon", "orange", "orangered1"),
                                                     background.title = "darkgray", cex.axis = 1)
                displayPars(datrack2) <- list(type = "a", alpha.title = 0, col= "black")                 
                otrack <- OverlayTrack(trackList=list(datrack, datrack2), 
                                       name="Coverage", background.title = "darkgray")
                
                jpeg_name <- paste0(seq_name, file_end)
                jpeg(jpeg_name, width = 700, height = 500)
                plotTracks(list(dtrack, gtrack, otrack), add53= TRUE, 
                           stacking = "squish", stackHeight = 0.9, add = TRUE)
                dev.off()
                }  else {
		print(paste0("No proteins found for ", seq_name)) }
             }
        setwd("..")
        }
    }
