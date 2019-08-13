# Make a regional plot
####

# Notes:
# Function for fixing chromosome name convention is needed (chr1 vs 1)
# Check that the number of titles matches the number of bed files
# Assign default background colors to bed tracks
# Move from biomart query to tabix/htslib for gene regions

# Required packages
library(GenomicRanges)
library(Gviz)
library(biomaRt)
# library(RColorBrewer)
library(data.table)

makeMart <- function(genome_build = "hg19"){
  genome_build <- tolower(genome_build)
  if (genome_build == "hg19" | genome_build == "grch37"){
      mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    } else if (genome_build == "hg38" | genome_build == "grch38") {
      mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    } else {
      warning("Unrecognized genome build, defaulting to hg19.")
      mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    }
  return(mart)
}

makeGeneTrack <- function(chr, start, end, mart, 
                          genome_build = "hg19", 
                          highlight = NULL, 
                          title = "Genomic Context",
                          background_color = "#31a354",
                          background_frame = T){
  # query biomart for the gene positions and symbols
  gene.data <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),filters = c("chromosome_name"), values=list(chr), mart=mart)

  # subset to genes in the region
  gene.data <- gene.data[(gene.data$end_position >= start) & (gene.data$start_position <= end),]
  # gene.data[gene.data$hgnc_symbol == "", "hgnc_symbol"] <- gene.data[gene.data$hgnc_symbol == "", "ensembl_gene_id"]
  gene.data <- gene.data[gene.data$hgnc_symbol != "",]

  # define a color for non highlighted genes (tan)
  gene.data$color <- "#ffd470"

  # change the color of highlighted genes
  if (!is.null(highlight)){
    highlight <- c(highlight)
    gene.data[gene.data$hgnc_symbol %in% highlight, "color"] <- "#de2d26"
  }
  
  gene.track <- AnnotationTrack(
    name = title, 
    start = gene.data$start_position,
    end = gene.data$end_position,
    id = gene.data$hgnc_symbol,
    feature = gene.data$color,
    chromosome = gene.data$chromosome_name[1], 
    genome = genome_build,
    col="transparent",
    background.title = background_color,
    col.frame = background_color,
    frame = background_frame,
    groupAnnotation="id",
    "#ffd470" = "#ffd470",
    "#0fa3ff" = "#0fa3ff",
    "#de2d26" = "#de2d26"
  )
  # displayPars(gene.track) <- list(fontsize.group = 10)
  return(gene.track)
}

makeBedTrack <- function(bed_data, 
                         title = "", 
                         genome_build = "hg19", 
                         background_color = "#485467",
                         background_frame = T){
  # Assume that colors are stored in column 5
  # Need to make a better way of adding color to these annotation tracks
  
  # Make everything the same color if given a 4 column bed file
  if (ncol(bed_data) == 4){
    bed_data[,5] <- "#878787"
  } 

  bed.track <- AnnotationTrack(
      name = title, 
      start = bed_data[,2],
      end = bed_data[,3],
      id = bed_data[,4],
      feature = bed_data[,5],
      chromosome = bed_data[1,1], 
      genome = genome_build,
      stacking = "dense",
      col = "transparent",
      background.title = background_color,
      col.frame = background_color,
      frame = background_frame,
      groupAnnotation = "feature",
      "255,255,255" = "#ffffff",
      "255,195,77" = "#ffc34d",
      "255,0,0" = "#ff0000",
      "255,69,0" = "#ff4500",
      "0,100,0" = "#006400",
      "255,255,0" = "#ffff00",
      "192,192,192" = "#c0c0c0",
      "205,92,92" = "#cd5c5c",
      "0,128,0" = "#008000",
      "128,128,128" = "#808080",
      "255,195,77" = "#ffc34d",
      "255,69,0" = "#ff4500",
      "194,225,5" = "#c2e105",
      "#ffffff" = "#ffffff",
      "#ffc34d" = "#ffc34d",
      "#ff0000" = "#ff0000",
      "#ff4500" = "#ff4500",
      "#006400" = "#006400",
      "#ffff00" = "#ffff00",
      "#c0c0c0" = "#c0c0c0",
      "#cd5c5c" = "#cd5c5c",
      "#008000" = "#008000",
      "#808080" = "#808080",
      "#ffc34d" = "#ffc34d",
      "#ff4500" = "#ff4500",
      "#c2e105" = "#c2e105"
    )
  return(bed.track)
}

makeVariantTrack <- function(variant_data, chr_column, pos_column, y_column,
                              marker_column = NULL,
                              ld_data = NULL,
                              ld_ref = NULL,
                              genome_build = "hg19", 
                              should_log = T,
                              horizontal_value = 5e-8,
                              horizontal_color = "red",
                              horizontal_style = "dashed",
                              background_color = "#e6550d",
                              background_frame = T,
                              point_color = "#252525",
                              title = ""){
  
  # log the data if we need to
  if (should_log){
    variant_data[,y_column] <- -log10(variant_data[,y_column])
    horizontal_value <- -log10(horizontal_value)
  }

  # merge with LD information if we need to
  if (!is.null(marker_column) & !is.null(ld_data) & !is.null(ld_ref)){
    # get to right format
    ld_data <- t(ld_data[,2:ncol(ld_data)])
    ld_data <- data.frame(
      marker = row.names(ld_data),
      r = ld_data[,1],
      stringsAsFactors = F)

    # first make sure the markername column is correct in the ld data file
    names(ld_data)[1] <- marker_column

    # make sure that we have the same sep character in the marker names
    ld_data[,marker_column] <- gsub("[[:punct:]]", "_", ld_data[,marker_column])
    ld_ref <- gsub("[[:punct:]]", "_", ld_ref)
    variant_data[,marker_column] <- gsub("[[:punct:]]", "_", variant_data[,marker_column])

    # then merge
    variant_data <- merge(variant_data,
                          ld_data,
                          by = marker_column,
                          all.x = T)

    # column number with ld information
    ld_col <- ncol(variant_data)

    # generate the colors 
    variant_data$color <- "#B8B8B8"
    v.na <- variant_data[is.na(variant_data[,ld_col]),]
    v.notna <- variant_data[!is.na(variant_data[,ld_col]),]
    v.notna[v.notna[,ld_col] >= 0, "color"] <- "#357EBD" 
    v.notna[v.notna[,ld_col] >= 0.2, "color"] <- "#46B8DA" 
    v.notna[v.notna[,ld_col] >= 0.4, "color"] <- "#5CB85C" 
    v.notna[v.notna[,ld_col] >= 0.6, "color"] <- "#EEA236" 
    v.notna[v.notna[,ld_col] >= 0.8, "color"] <- "#D43F3A"
    
    variant_data <- rbind(v.na, v.notna)
    variant_data[variant_data[,marker_column] == ld_ref, "color"] <- "#9632B8"
    
    colors <- c("#B8B8B8", "#357EBD", "#46B8DA", "#5CB85C", "#EEA236", "#D43F3A", "#9632B8")
    
    # get y bounds for plot
    ylims <- c(0, max(variant_data$pvalue))

    # generate the variant track list, to combine later
    variant.tracks <- list()
    for (val in colors){
      cur_data <- variant_data[variant_data$color == val,]
      if (nrow(cur_data) == 0) next
      variant.tracks[[length(variant.tracks) + 1]] <- DataTrack(
        data = cur_data[,y_column], 
        start = cur_data[,pos_column], 
        end = cur_data[,pos_column], 
        chromosome = cur_data[,chr_column], 
        genome = genome_build, 
        name = title, 
        background.title = background_color,
        col.frame = background_color,
        frame = background_frame,
        baseline = horizontal_value,
        col.baseline = horizontal_color,
        lty.baseline = horizontal_style,
        col = cur_data$color[1],
        ylim = ylims)
    }

    variant.track <- OverlayTrack(variant.tracks,
      background.title = background_color,
      col.frame = background_color,
      frame = background_frame,
    )
    
  } else {
    print("No LD information found.")
    # generate the variant track
    variant.track <- DataTrack(
      data = variant_data[,y_column], 
      start = variant_data[,pos_column], 
      end = variant_data[,pos_column], 
      chromosome = variant_data[,chr_column], 
      genome = genome_build, 
      name = title, 
      background.title = background_color,
      col.frame = background_color,
      frame = background_frame,
      baseline = horizontal_value,
      col.baseline = horizontal_color,
      lty.baseline = horizontal_style,
      col = point_color)
  }
  return(variant.track)
}

makePlot <- function(chr, start, end, variant_data, variant_chr_column, variant_pos_column, variant_y_column,
                        variant_marker_column = NULL,
                        variant_ld_data = NULL,
                        variant_ld_ref = NULL,
                        genome_build = "hg19",
                        variant_should_log = T,
                        variant_horizontal_value = 5e-8,
                        variant_horizontal_color = "red",
                        variant_horizontal_style = "dashed",
                        variant_background_color = "#e6550d",
                        variant_background_frame = T,
                        variant_point_color = "#252525",
                        variant_title = "",
                        gene_highlight = NULL,
                        gene_title = NULL,
                        gene_background_color = "#31a354",
                        gene_frame = T,
                        bed_data = NULL, 
                        bed_titles = NULL,
                        bed_background_colors = NULL,
                        bed_frame = T){
  
  # make the data track for variant statistics
  print(min(variant_data[,variant_y_column]))
  variant.track <- makeVariantTrack(variant_data, 
                                    variant_chr_column, 
                                    variant_pos_column, 
                                    variant_y_column, 
                                    variant_marker_column,
                                    variant_ld_data,
                                    variant_ld_ref,
                                    genome_build, 
                                    variant_should_log, 
                                    variant_horizontal_value,
                                    variant_horizontal_color,
                                    variant_horizontal_style,
                                    variant_background_color,
                                    variant_background_frame,
                                    variant_point_color,
                                    variant_title)

  # make a mart object for gene coordinates
  mart <- makeMart(genome_build)
  
  # make gene track 
  gene.track <- makeGeneTrack(chr, 
                              start, 
                              end, 
                              mart, 
                              genome_build, 
                              gene_highlight, 
                              gene_title)

  # make bed tracks
  if (!is.null(bed_data)){
    bed.tracks <- list()
    for (bed_index in 1:length(bed_data)){
      bed.tracks[[bed_index]] <- makeBedTrack(bed_data[[bed_index]],
                                        bed_titles[bed_index], 
                                        genome_build, 
                                        bed_background_colors[bed_index],
                                        bed_frame)
    }
  }
  
  # make scale track
  scale.track <- GenomeAxisTrack(genome = genome_build, 
                              chromosome = chr)
  
  # make idiogram track
  idiogram.track <- IdeogramTrack(genome = genome_build, 
                              chromosome = chr)
  
  # generate the plot
  if (!is.null(bed_data)){
    plotTracks(
      c(list(idiogram.track, scale.track, variant.track, gene.track), bed.tracks),
      showTitle = TRUE,
      sizes=c(1, 2, 5, 4, rep(2, length(bed.tracks))),
      from = start,
      to = end
    )  
  } else {
    plotTracks(
      list(idiogram.track, scale.track, variant.track, gene.track),
      showTitle = TRUE,
      sizes=c(1, 2, 5, 4),
      from = start,
      to = end
    )
  }
}


# testing inputs
data.dir <- "/Users/tmajaria/Documents/projects/WGSAssociationVisualization/"
variant.file <- paste0(data.dir, "1kg-t2d.all.assoc.aug12.txt.gz")
ld.file <- paste0(data.dir, "1kg-t2d.chr20_60.9M-61.1M.ld.csv")
bed.file <- paste0(data.dir, "Islets.chromatinStates.bed")

chr <- 20
start <- 60900000
end <- 61100000
variant_data <- fread(variant.file, data.table = F, stringsAsFactors = F)
variant_data <- variant_data[!is.na(variant_data$chr),]
variant_data$MarkerName <- sub("chr", "", variant_data$MarkerName)
variant_chr_column = "chr"
variant_pos_column = "pos"
variant_y_column = "pvalue"
variant_marker_column = "MarkerName"
variant_ld_data <- fread(ld.file, data.table = F, stringsAsFactors = F)
variant_ld_ref = "20-61000005-A-G"
genome_build = "hg19"
variant_should_log = T
variant_horizontal_value = 5e-8
variant_horizontal_color = "red"
variant_horizontal_style = "dashed"
variant_background_color = "#e6550d"
variant_background_frame = T
variant_point_color = "#252525"
variant_title = "Waste hip ratio"
gene_highlight = "LYPLAL1"
gene_title = "Genomic Context"
gene_background_color = "#31a354"
gene_frame = T
bed_data <- list(fread(bed.file, data.table = F, stringsAsFactors = F))
bed_data[[1]][,1] <- sub("chr", "", bed_data[[1]][,1])
bed_titles = "Islet HMM"
bed_background_colors = c("#3182bd")
bed_frame = T

makePlot(chr, start, end, variant_data, variant_chr_column, variant_pos_column, variant_y_column,
         variant_marker_column,
         variant_ld_data,
         variant_ld_ref,
         genome_build,
         variant_should_log,
         variant_horizontal_value,
         variant_horizontal_color,
         variant_horizontal_style,
         variant_background_color,
         variant_background_frame,
         variant_point_color,
         variant_title,
         gene_highlight,
         gene_title,
         gene_background_color,
         gene_frame,
         bed_data, 
         bed_titles,
         bed_background_colors,
         bed_frame)