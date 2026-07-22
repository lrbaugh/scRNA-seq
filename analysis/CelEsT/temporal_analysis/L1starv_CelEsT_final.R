#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest
packages <- c("stringr", 
              "ggplot2",
              "openxlsx",
              "decoupleR",
              "ipkg",
              "tidyr",
              "ggrepel",
              "patchwork",
              "dplyr"
)

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Here define packages which need to be loaded through biocmanager

biocmanager_packages <- c("DESeq2")

bioc_package.check <- lapply(
  biocmanager_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
      }
      
      BiocManager::install(x, dependencies = TRUE)
      
      library(x, character.only = TRUE)
      
    }
  }
)



github_packages <- c("LBMC/wormRef")

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")}

if (!requireNamespace("remotes", quietly = TRUE)){
  install.packages("remotes")}

if (!requireNamespace("pryr", quietly = TRUE)){
  remotes::install_version("pryr",
                           version = "0.1.6",
                           repos   = "https://cloud.r-project.org",
                           upgrade = "never")
}

gh_package.check <- lapply(
  github_packages,
  FUN = function(x) {
    if (!require(str_remove(x, ".*\\/"), character.only = TRUE)) {
      
      if (!requireNamespace("devtools", quietly = TRUE)){
        install.packages("devtools")}
      
      devtools::install_github(x, build_vignettes = TRUE)
      
      library(str_remove(x, ".*\\/"), character.only = TRUE)
      
    }
  }
)

#### FUNCTIONS ####

# KUDOS to Kamil from biostars [https://www.biostars.org/p/171766/]
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  temp_counts <- counts[idx,]
  temp_effLen <- effLen[idx,]
  temp_featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(temp_counts), function(i) {
    rate = log(temp_counts[,i]) - log(temp_effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(temp_counts)
  rownames(tpm) <- rownames(temp_counts)
  
  return(tpm)
  
}

map2color <- function(x, 
                      pal,
                      limits = NULL){
  
  if(is.null(limits)){
    limits = range(x)
  }
  
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
  
}

do.volcano.plots <- function(decouple_in = L1starv_decouple){

lapply(unique(lookup_sheet$timepoint)[unique(lookup_sheet$timepoint) != "minus2hr"], function(tp){
  # tp <- unique(lookup_sheet$timepoint)[unique(lookup_sheet$timepoint) != "minus2hr"][5]
  temp_decouple <- decouple_in[decouple_in$condition == paste0("X", tp), ]
  
  temp_decouple_wide <- pivot_wider(
    data = temp_decouple[, c("source", "condition", "score")], 
    names_from = "source", 
    values_from = "score"
  )
  temp_decouple_wide <- data.frame(t(temp_decouple_wide))
  
  colnames(temp_decouple_wide) <- temp_decouple_wide[1, ]
  temp_decouple_wide <- temp_decouple_wide[-1, , drop = FALSE]
  
  temp_decouple_wide[] <- lapply(temp_decouple_wide, as.numeric)
  
  temp_decouple_wide_p <- as.data.frame(
    pivot_wider(data = temp_decouple[, c("source", "condition", "p_value")], 
                names_from = "source", values_from = "p_value")
  )
  row.names(temp_decouple_wide_p) <- unlist(temp_decouple_wide_p[, 1])
  temp_decouple_wide_p <- temp_decouple_wide_p[, -1, drop = FALSE]
  
  temp_decouple_wide_p <- data.frame(t(temp_decouple_wide_p))
  
  temp_decouple_wide_p[] <- lapply(temp_decouple_wide_p, as.numeric)
  
  DEdata_forbubbleplot <- data.frame(
    "mean_z" = unlist(temp_decouple_wide[row.names(temp_decouple_wide_p), ]),
    "p_val" = -log10(unlist(temp_decouple_wide_p)),
    "label" = toupper(unlist(wormRef::Cel_genes[match(row.names(temp_decouple_wide_p), wormRef::Cel_genes$sequence_name), "public_name"]))
  )
  
  DEdata_forbubbleplot[DEdata_forbubbleplot$p_val < 2.2, "label"] <- ""
  
  tempcolours_for_plot <- map2color(DEdata_forbubbleplot$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))
  
  return(
    ggplot(DEdata_forbubbleplot, aes(x = mean_z, y = p_val, size = 10^p_val, label = toupper(label))) + 
          geom_point(col = tempcolours_for_plot, alpha = 0.7) + 
          theme_classic() + 
          geom_label_repel(size = 2, max.overlaps = 40) +
          geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
          theme(
            legend.position = "none",
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10)
          ) + 
          ggtitle(tp) + 
          ylab(substitute("-log"[10]~"(p-value)")) + 
          xlab("TF activity") 
  )
  
})

}

do.time.plots <- function(decouple_in = L1starv_decouple,
                          logscale = FALSE,
                          replicates = FALSE,
                          smooth_line = FALSE) {
  
  tempplotlist <- lapply(unique(decouple_in$source), function(x){

    temp_decouple <- decouple_in[decouple_in$source == x, ]
    
    timepoints <- c(0, 2, 4, 6, 9, 12, 24, 48, 96, 192, 288)
    
    if(any(temp_decouple$time < 0)) {
      
      timepoints <- c(-2, timepoints)
      
    }
    
    tempplot <- ggplot(temp_decouple, aes(x = time, y = score)) + 
      theme_classic()
    
    if(any(temp_decouple$time < 0)) {
      
      tempplot <- tempplot + 
        annotate(
          "rect",
          xmin = -Inf,
          xmax = 0,
          ymin = -Inf,
          ymax = Inf,
          fill = "grey80",
          alpha = 0.5
        )
      
    }
    
      tempplot <- tempplot +
      theme(
        legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)
      ) + 
      ggtitle(paste(toupper(wormRef::Cel_genes[match(x, wormRef::Cel_genes$sequence_name), "public_name"]), x, sep = " - ")) + 
      ylab(substitute("TF activity score (AU)")) + 
      xlab("Time (h)") + 
      coord_cartesian(xlim = c(min(timepoints), 300))
    
    if(isTRUE(replicates)){
      
      newtemp <- temp_decouple %>% group_by(time) %>% summarise (mean = mean(score),
                                                                 sd = sd(score))
      
      newtemp[, "se"] <- newtemp$sd / sqrt(4)
      
      newtemp[, "ninetyfiveci_up"] <- newtemp$mean + (1.96 * newtemp$se)
      newtemp[, "ninetyfiveci_down"] <- newtemp$mean - (1.96 * newtemp$se)
      
      tempplot <- ggplot(data = newtemp, aes(x = time, y = mean)) + 
        annotate(
          "rect",
          xmin = -Inf,
          xmax = 0,
          ymin = -Inf,
          ymax = Inf,
          fill = "grey80",
          alpha = 0.5
        ) +
        geom_point() + 
        geom_errorbar(aes(ymax = ninetyfiveci_up,
                      ymin = ninetyfiveci_down)) + 
        theme_classic() + 
        theme(
          legend.position = "none",
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)
        ) + 
        ggtitle(paste(toupper(wormRef::Cel_genes[match(x, wormRef::Cel_genes$sequence_name), "public_name"]), x, sep = " - ")) + 
        ylab(substitute("TF activity score (AU)")) + 
        xlab("Time (h)") + 
        coord_cartesian(xlim = c(min(timepoints), 300))

      
    } else {

        tempplot <- tempplot + geom_point(aes(size = -log10(p_value)),
                   colour = "red")
      
    }
    
    if(isTRUE(logscale)){
      
      tempplot <- tempplot + 
        scale_x_continuous(trans = scales::pseudo_log_trans(base = 2),
                           breaks = timepoints,
                           labels = timepoints) 
      
    }
    
    if(isTRUE(smooth_line)){
    
    tempplot <- tempplot +
      geom_smooth(se = FALSE,
                  col = "grey45",
                  linetype = "dashed") 
    
  } else {
    
    tempplot <- tempplot +
      geom_line(col = "grey45",
                linetype = "dashed",
                alpha = 1)
    
  }
    
    return(tempplot)
    
  })
  
  names(tempplotlist) <- unique(decouple_in$source)
  
  tempplotlist
  
}

#### set working directory ####

dir.create("~/L1starv_CelEsT", showWarnings = FALSE)

setwd("~/L1starv_CelEsT")

dir.create("input", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("graphics", showWarnings = FALSE)
dir.create("plotdata", showWarnings = FALSE)

#### load input counts ####

ipkg::download_file(url = "https://github.com/IBMB-MFP/CelEsT-app/blob/main/www/CelEsTv1pt1_GRN.txt",
              destfile = "input/CelEsTv1pt1_GRN.txt")

CelEsT_v1pt1 <- read.table("input/CelEsTv1pt1_GRN.txt",
                           header = TRUE,
                           sep = "\t",
                           fill = TRUE)

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE173656&format=file&file=GSE173656%5FTimeSeries%5FRNAseq%5Foutput%2Exlsx",
              dest = "input/GSE173656_TimeSeries_RNAseq_output.xlsx")

raw_counts <- openxlsx::read.xlsx("input/GSE173656_TimeSeries_RNAseq_output.xlsx",
                                  sheet = 3)

row.names(raw_counts) <- raw_counts[, 1]

raw_counts <- raw_counts[, 2:ncol(raw_counts)]

lookup_sheet <- data.frame(row.names = colnames(raw_counts),
                           timepoint = str_extract(colnames(raw_counts),
                                                   "^.*hr"),
                           replicate = str_extract(colnames(raw_counts),
                                                   "rep[0-9]$"))

lookup_sheet[is.na(lookup_sheet$timepoint), "timepoint"] <- str_extract(colnames(raw_counts),
                                                             "^.*d")[!is.na(str_extract(colnames(raw_counts),
                                                                                        "^.*d"))]

lookup_sheet[, "timepoint_num"] <- as.numeric(str_extract(lookup_sheet$timepoint, 
                                               "[0-9]{1,2}"))

lookup_sheet[str_detect(lookup_sheet$timepoint, "minus"), "timepoint_num"] <- -lookup_sheet[str_detect(lookup_sheet$timepoint, "minus"), "timepoint_num"]

lookup_sheet[str_detect(lookup_sheet$timepoint, "d"), "timepoint_num"] <- lookup_sheet[str_detect(lookup_sheet$timepoint, "d"), "timepoint_num"] * 24

#### filter out unexpressed TFs ####

my_temp_lengths <- sapply(row.names(raw_counts), function(x){

    mean(wormRef::Cel_genes[wormRef::Cel_genes$sequence_name == x, "transcript_length"])

    })

  my_temp_lengths <- my_temp_lengths[!is.na(my_temp_lengths)]

  gene_there <- base::intersect(names(my_temp_lengths), row.names(raw_counts))

  counts_for_DE_TPM <- counts_to_tpm(counts = raw_counts[gene_there, ],
                                    featureLength = my_temp_lengths[gene_there],
                                    meanFragmentLength = rep(100, times = ncol(raw_counts)))

TF_TPMs <- counts_for_DE_TPM[row.names(counts_for_DE_TPM) %in% unique(CelEsT_v1pt1$source), ]

# remove TFs with TPM > 1 in fewer than 4 samples
TFs_to_censor <- row.names(TF_TPMs)[!apply(TF_TPMs, 1, function(x){sum(x > 1) > 3})]

write.xlsx(data.frame("seq_name" = TFs_to_censor,
                       "public_name" = toupper(wormRef::Cel_genes[match(TFs_to_censor, Cel_genes$sequence_name), "public_name"])),
            "output/censored_TFs.xlsx")

CelEsT_cens <- CelEsT_v1pt1[!CelEsT_v1pt1$source %in% TFs_to_censor, ]

#### variance stabilising transformation ####

L1starv_dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                              colData = lookup_sheet,
                                              design = ~ replicate + timepoint)

L1starv_vst <- vst(L1starv_dds)

L1starv_vst_df <- assay(L1starv_vst)

#### do decoupling ####

L1starvVST_decouple <- decoupleR::decouple(
  mat = L1starv_vst_df, 
  network = CelEsT_cens,
  .source = "source",
  .target = "target",
  statistics = "mlm",
  args = list(mlm = list(.mor = "weight")),
  consensus_score = FALSE
)

L1starvVST_decouple[, "TF"] <- toupper(wormRef::Cel_genes[match(L1starvVST_decouple$source, Cel_genes$sequence_name), "public_name"])
L1starvVST_decouple[, "time"] <- lookup_sheet[match(str_remove(L1starvVST_decouple$condition, "_rep.*$"), paste0(lookup_sheet$timepoint)), "timepoint_num"]
L1starvVST_decouple[, "rep"] <- str_extract(L1starvVST_decouple$condition, "rep.*$")

write.xlsx(list(L1starvVST_decouple),
           "output/L1starv_TFactivity.xlsx")

##### plot dynamics ####

# for VST plots, 4 replicates plotted separately

VST_logscale_plots <- do.time.plots(decouple_in = L1starvVST_decouple,
                                            logscale = TRUE,
                                    replicates = TRUE,
                                    smooth_line = FALSE)

pdf("graphics/L1starv_VST_TFplots_LOGSCALE.pdf")

VST_logscale_plots

dev.off()
