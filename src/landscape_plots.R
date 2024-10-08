
getCNstatefunc <- function(segments, dlpbins, ncores = 1){
  dlpbins_g <- plyranges::as_granges(dlpbins %>% dplyr::rename(seqnames = chr))
  segments_g <- plyranges::as_granges(segments %>% dplyr::rename(seqnames = chr), keep_mcols = TRUE)
  
  bins <- plyranges::find_overlaps(dlpbins_g, segments_g) %>%
    as.data.frame() %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::select(-width, -strand) %>%
    dplyr::select(cell_id, chr, start, end, dplyr::everything()) %>%
    signals::orderdf()
  
  return(bins)
}

#' @export
getCNstate <- function(segments, ncores = 1){
  
  dlpbins <- signals::getBins(binsize = 0.5e6)
  
  if (ncores == 1) {
    df <- data.table::rbindlist(lapply(unique(segments$cell_id),
                                       function(x){
                                         segments %>% dplyr::filter(cell_id == x) %>% getCNstatefunc(., dlpbins = dlpbins)
                                       })) %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(cell_id))
  } else{
    df <- data.table::rbindlist(parallel::mclapply(unique(segments$cell_id),
                                                   function(x){
                                                     segments %>% dplyr::filter(cell_id == x) %>%
                                                       getCNstatefunc(., dlpbins = dlpbins)
                                                   }, mc.cores = ncores)) %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(cell_id))
  }
  
  return(df)
}

remove_centromeres <- function(cndat, binsize = 0.5e6, extend = 6e6){
  ideogram_dat <- cytoband_map[["hg19"]]
  names(ideogram_dat) <- c("chr", "start", "end", "band", "colval")
  
  df_other <- data.frame(chr = "7", start = 56000001, end = 65500001, "1", "other")
  
  ideogram_dat <- ideogram_dat %>% 
    bind_rows(df_other) %>% 
    filter(colval %in% c("acen", "stalk", "gvar", "other")) %>% 
    mutate(start = start - extend) %>% 
    mutate(end =  end + extend) %>% 
    dplyr::mutate(chr = stringr::str_remove(chr, "chr")) %>% 
    dplyr::mutate(start = round(start / binsize) * binsize + 1, 
                  end = round(end / binsize) * binsize + 1) %>% 
    mutate(filt = TRUE) %>% 
    select(chr, start, end, filt) %>% 
    mutate(cell_id = "1") %>% 
    getCNstate() %>% 
    select(-cell_id)
  
  cndat <- left_join(cndat, ideogram_dat) %>% 
    replace_na(list(filt = FALSE)) %>% 
    filter(filt == FALSE)
  
  return(cndat)
}

plottinglist_ <- function(CNbins){
  #arrange segments in order, generate segment index and reorder CN state factor
  binsize <- CNbins$end[1] - CNbins$start[1] + 1
  bins <- signals::getBins(binsize = binsize) %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::mutate(idx = 1:dplyr::n())
  
  
  CNbins <- dplyr::full_join(bins, CNbins) %>%
    dplyr::mutate(idxs = forcats::fct_reorder(factor(idx), idx)) %>% 
    filter(!is.na(chr))
  
  #get breaks - first index of each chromosome
  chrbreaks <- bins %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(idx)
  
  chrticks <- bins %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(idx = round(median(idx))) %>%
    dplyr::pull(idx)
  
  chrlabels <- gtools::mixedsort(unique(CNbins$chr))
  chrlables <- chrlabels[!is.na(chrlabels)]
  minidx <- min(bins$idx)
  maxidx <- max(bins$idx)
  
  return(list(CNbins = CNbins, chrbreaks = chrbreaks, chrticks = chrticks, chrlabels = chrlabels, minidx = minidx, maxidx = maxidx))
}

genplot <- function(pl, tcga, mytitle = NULL, breaks = c(0.0, 0.5, 1.0), alphaval = 1.0,
                    limits = c(0, 1), ylabel = "", pcawg = NULL,pcawglinesize = 0.3,
                    yline = NULL,
                    chromspace_size = 0.18,
                    losslabels = c("0", "0.5", "1.0"), returnlist = FALSE,
                    loss_col = "#28536C",
                    gain_col = "#550000"){
  
  textsize <- 2
  filltype <- "stack"
  
  datgain <- pl$CNbins %>%
    mutate(frequency = gain) %>%
    arrange(chr, start) %>%
    mutate(x = rleid(frequency)) %>%
    add_count(x) %>%
    mutate(frequency = ifelse(is.na(frequency) & nn >10, 0, frequency)) %>%
    fill(frequency, .direction = "updown")
  
  cgains <- datgain %>%
    ggplot(aes(x = idx, y = frequency)) +
    geom_area(alpha = alphaval, position = filltype, aes(fill = "type")) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx)) + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, size = chromspace_size,
                        col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = c(gain_col),
                      labels = c("Gain")) + 
    scale_y_continuous(expand = c(0,0),  breaks = breaks, limits = limits) +
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(vjust=-1.2),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  if (!is.null(yline)){
    cgains <- cgains + geom_hline(yintercept = yline, lty = 2)
  }
  
  datloss <- pl$CNbins %>%
    mutate(frequency = loss) %>%
    arrange(chr, start) %>%
    mutate(x = rleid(frequency)) %>%
    add_count(x) %>%
    mutate(frequency = ifelse(is.na(frequency) & nn >10, 0, frequency)) %>%
    fill(frequency, .direction = "updown")
  
  closs <- datloss %>%
    ggplot(aes(x = idx, y = -frequency)) +
    geom_area(alpha = alphaval, position = filltype, aes(fill = "type")) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx), position = "top") + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, 
                        size = chromspace_size,
                        col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = c(loss_col),
                      labels = c("Loss")) + 
    scale_y_continuous(expand = c(0, 0), breaks = -breaks, 
                       labels = losslabels, limits = -rev(limits)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  if (!is.null(yline)){
    closs <- closs + geom_hline(yintercept = -yline, lty = 2)
  }
  
  if (!is.null(pcawg)){
    cgains <- cgains + geom_line(data = pl_tcga$CNbins, 
                                 aes(x = idx, y = gain), 
                                 col = gain_col, size = pcawglinesize)
    closs <- closs + geom_line(data = pcawg$CNbins, 
                               aes(x = idx, y = -loss), 
                               col = loss_col, size = pcawglinesize)
  }
  
  if (!is.null(mytitle)){
    title <- ggdraw() +
      draw_label(
        mytitle,
        size = 8,
        #fontface = 'bold',
        x = 0.0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(4, 2, 2, 7)
      )
    
    myplot <- cowplot::plot_grid(
      title,
      cgains + theme(legend.position = "none")+ ylab(ylabel), 
      closs + theme(legend.position = "none"), 
      align = "v", axis = "lr", ncol = 1, rel_heights = c(0.1, 1, 0.8))
  } else{
    myplot <- cowplot::plot_grid(
      cgains + theme(legend.position = "none")+ ylab(ylabel), 
      closs + theme(legend.position = "none"), 
      align = "v", axis = "lr", ncol = 1, rel_heights = c( 1, 0.8))
  }
  
  if (returnlist == TRUE){
    myplot <- list(gain = cgains, loss = closs, datgain = datgain, datloss = datloss)
  }
  
  return(myplot)
}

genplot_groups <- function(pl, 
                           pl_tcga = NULL,
                           group = "genotype", 
                           plotfrequency = FALSE,
                           scaling_tcga = 1,
                           persamplecutoff = 2, 
                           mytitle = NULL, 
                           breaks = c(0.0, 0.5, 1.0), 
                           alphaval = 1.0,
                           limits = c(0, 1), 
                           ylabel = "", 
                           pcawg = NULL,
                           tcgalinesize = 0.3,
                           losslabels = c("0", "0.5", "1.0"), 
                           returnlist = FALSE, 
                           clrs = unlist(config$clrs$genotype)){
  
  textsize <- 2
  filltype <- "stack"
  if (group == "all"){
    pl$CNbins["group"] <- "all"
  } else{
    pl$CNbins["group"] <- pl$CNbins[[group]]
  }
  
  nsamples <- unique(pl$CNbins$sample)
  nsamples <- length(nsamples[!is.na(nsamples)])
  
  datgain <- pl$CNbins %>%
    group_by(chr, start, end, group, idx) %>% 
    summarize(frequency = sum(gain >= persamplecutoff)) %>% 
    ungroup() %>% 
    filter(!is.na(group)) %>% 
    arrange(chr, start) %>%
    mutate(x = rleid(frequency)) %>%
    add_count(x) %>%
    mutate(frequency = ifelse(is.na(frequency) & n >10, 0, frequency)) %>%
    fill(frequency, .direction = "updown")
  if (plotfrequency == TRUE){
    datgain$frequency <- datgain$frequency / nsamples
  }
  if (group == "all"){
    datgain$group <- "Gain"
  }
  
  cgains <- datgain %>%
    na.omit() %>% 
    ggplot(aes(x = idx, y = frequency)) +
    geom_area(alpha = alphaval, position = filltype, aes(fill = group), stat = "identity") +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx)) + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = clrs) + 
    scale_y_continuous(expand = c(0,0),  breaks = breaks, limits = limits) +
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(vjust=-1.2),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  datloss <- pl$CNbins %>%
    group_by(chr, start, end, group, idx) %>% 
    summarize(frequency = sum(loss >= persamplecutoff)) %>% 
    ungroup() %>% 
    filter(!is.na(group)) %>% 
    arrange(chr, start) %>%
    mutate(x = rleid(frequency)) %>%
    add_count(x) %>%
    mutate(frequency = ifelse(is.na(frequency) & n >10, 0, frequency)) %>%
    fill(frequency, .direction = "updown")
  if (group == "all"){
    datloss$group <- "Loss"
  }
  if (plotfrequency == TRUE){
    datloss$frequency <- datloss$frequency / nsamples
  }
  
  closs <- datloss %>%
    ggplot(aes(x = idx, y = -frequency)) +
    geom_area(alpha = alphaval, position = filltype, aes(fill = group), stat = "identity") +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx), position = "top") + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = clrs) + 
    scale_y_continuous(expand = c(0, 0), breaks = -breaks, 
                       labels = losslabels, limits = -rev(limits)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  #scaling_gains <- max(pl_tcga$CNbins$Gain, na.rm = T)
  
  if (!is.null(pl_tcga)){
    cgains <- cgains + geom_line(data = pl_tcga$CNbins, 
                                 aes(x = idx, y = scaling_tcga * Gain), 
                                 col = "#550000", size = tcgalinesize)
    closs <- closs + geom_line(data = pl_tcga$CNbins, 
                               aes(x = idx, y = -scaling_tcga * Loss), 
                               col = "#28536C", size = tcgalinesize)
  }
  
  if (group == "all"){
    legpos <- "none"
  } else{
    legpos <- c(0.1, 0.1)
  }
  
  if (!is.null(mytitle)){
    title <- ggdraw() +
      draw_label(
        mytitle,
        fontface = 'bold',
        x = 0.05,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(4, 2, 2, 7)
      )
    
    myplot <- cowplot::plot_grid(
      title,
      cgains + theme(legend.position = "none")+ ylab(ylabel), 
      closs + theme(legend.position = legpos), 
      align = "v", axis = "lr", ncol = 1, rel_heights = c(0.1, 1, 0.8))
  } else{
    myplot <- cowplot::plot_grid(
      cgains + theme(legend.position = "none")+ ylab(ylabel), 
      closs + theme(legend.position = legpos), 
      align = "v", axis = "lr", ncol = 1, rel_heights = c( 1, 0.8))
  }
  
  if (returnlist == TRUE){
    myplot <- list(gain = cgains, loss = closs, datgain = datgain, datloss = datloss)
  }
  
  return(myplot)
}



genplotwithsamples <- function(plsample, mytitle = NULL, breaks = c(0.0, 0.5, 1.0), alphaval = 0.5,
                    limits = c(0, 1), breakssample = c(0,5, 10), limitssample = c(0, 10),
                    ylabel = "", pcawg = NULL,pcawglinesize = 0.3,
                    losslabels = c("0", "0.5", "1.0"),
                    losslabelssample = c("0", "5", "10"), returnlist = FALSE,
                    loss_col = "#28536C",
                    gain_col = "#550000"){
  
  textsize <- 2
  filltype <- "stack"
  
  datgain <- pl$CNbins %>%
    mutate(frequency = gain) %>%
    arrange(chr, start) %>%
    mutate(x = rleid(frequency)) %>%
    add_count(x) %>%
    mutate(frequency = ifelse(is.na(frequency) & n >10, 0, frequency)) %>%
    fill(frequency, .direction = "updown")
  
  cgains <- datgain %>%
    ggplot(aes(x = idx, y = frequency)) +
    geom_area(alpha = alphaval, position = filltype, aes(fill = "type")) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx)) + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = c(gain_col),
                      labels = c("Gain")) + 
    scale_y_continuous(expand = c(0,0),  breaks = breaks, limits = limits) +
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(vjust=-1.2),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  
  datgainsample <- plpersample$CNbins %>%
    mutate(gain = as.numeric(gain > 4)) %>% 
    mutate(frequency = gain) #%>%
    # arrange(chr, start) %>%
    # mutate(x = rleid(frequency)) %>%
    # add_count(x) %>%
    # mutate(frequency = ifelse(is.na(frequency) & n >10, 0, frequency)) %>%
    # fill(frequency, .direction = "updown")
  
  cgainssample <- datgainsample %>%
    ggplot(aes(x = idx, y = frequency)) +
    geom_area(alpha = 1.0, position = filltype, aes(fill = genotype)) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx)) + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = c("deepskyblue4", "plum4", "grey70", "red")) + 
    scale_y_continuous(expand = c(0,0),  breaks = breakssample, limits = limitssample) +
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(vjust=-1.2),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  datloss <- pl$CNbins %>%
    mutate(frequency = loss) %>%
    arrange(chr, start) %>%
    mutate(x = rleid(frequency)) %>%
    add_count(x) %>%
    mutate(frequency = ifelse(is.na(frequency) & n >10, 0, frequency)) %>%
    fill(frequency, .direction = "updown")
  
  closs <- datloss %>%
    ggplot(aes(x = idx, y = -frequency)) +
    geom_area(alpha = alphaval, position = filltype, aes(fill = "type")) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx), position = "top") + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    scale_fill_manual(values = c(loss_col),
                      labels = c("Loss")) + 
    scale_y_continuous(expand = c(0, 0), breaks = -breaks, 
                       labels = losslabels, limits = -rev(limits)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  
  datlosssample <- plpersample$CNbins %>%
    mutate(loss = as.numeric(loss> 4)) %>% 
    mutate(frequency = loss)
  
  closssample <- datlosssample %>%
    ggplot(aes(x = idx, y = -frequency)) +
    geom_area(alpha = 1.0, position = filltype, aes(fill = genotype)) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, 
                                expand = c(0, 0), 
                                guide = guide_axis(check.overlap = TRUE),
                                limits = c(pl$minidx, pl$maxidx), position = "top") + 
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("") +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    scale_fill_manual(values =  c("deepskyblue4", "plum4", "grey70", "red")) + 
    scale_y_continuous(expand = c(0, 0), breaks = -breakssample, 
                       labels = losslabelssample, limits = -rev(limitssample)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.title=element_blank()) +
    guides(guide = ggplot2::guide_axis(check.overlap = TRUE))
  
  if (!is.null(pcawg)){
    cgains <- cgains + geom_line(data = pcawg$CNbins, 
                                 aes(x = idx, y = gain), 
                                 col = gain_col, size = pcawglinesize)
    closs <- closs + geom_line(data = pcawg$CNbins, 
                               aes(x = idx, y = -loss), 
                               col = loss_col, size = pcawglinesize)
  }
  
  if (!is.null(mytitle)){
    title <- ggdraw() +
      draw_label(
        mytitle,
        fontface = 'bold',
        x = 0.05,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(4, 2, 2, 7)
      )
    
    myplot <- cowplot::plot_grid(
      title,
      cgains + theme(legend.position = "none")+ ylab(ylabel), 
      closs + theme(legend.position = "none"), 
      align = "v", axis = "lr", ncol = 1, rel_heights = c(0.1, 1, 0.8))
  } else{
    myplot <- cowplot::plot_grid(
      cgains + theme(legend.position = "none")+ ylab(ylabel), 
      closs + theme(legend.position = "none"), 
      align = "v", axis = "lr", ncol = 1, rel_heights = c( 1, 0.8))
  }
  
  if (returnlist == TRUE){
    myplot <- list(gain = cgains, loss = closs, datgain = datgain, datloss = datloss,
                   gainsample = cgainssample, losssample = closssample, datgainsample = datgainsample, datlosssample = datlosssample)
  }
  
  return(myplot)
}

