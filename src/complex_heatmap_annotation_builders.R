# heatmap
## col metadata annotation
### arbitrary tx anno (column) df function
### returns a list
### out$df is the input for `HeatmapAnnotation` `df`
### out$colors is the input for `HeatmapAnnotation` `col`
tx_col_df_maker <- function(meta){
  df <- gene_exon_tpm_meta %>% 
    select(transcript_id, recount_exon_id) %>% 
    unique()
  # collapse regions to list column
  df_l <- df %>% group_by(recount_exon_id) %>% summarise(tx = list(transcript_id))
  for (i in unique(df$transcript_id)){
    df_l[,i] <- grepl(i, df_l$tx)
  }
  out <- list()
  out$df <- df_l %>% select(-tx, -recount_exon_id) %>% data.frame()
  out$df[out$df == TRUE] <- 'TRUE'
  out$df[out$df == FALSE] <- 'FALSE'
  colors <- list()
  for (i in colnames(out$df)){
    colors[[i]] <- c('TRUE' = 'black', 'FALSE' = 'white')
  }
  out$colors <- colors
  out
}
# heatmap
## row metadata annotation
## same output as above
tx_row_df_maker <- function(meta){
  df <- data.frame(Condition = t(gene_exon_tpm) %>% as_tibble(rownames = 'name') %>% 
                     left_join(meta, by = 'name') %>% pull(Condition))
  colors <- c(pals::alphabet(n=26), pals::alphabet2(n=26), pals::glasbey(n=32))
  colors <- rep(colors, 100)
  names(colors) <- unique(df[,1])
  colors <- colors[!is.na(names(colors))]
  out <- list()
  out[['df']] <- df
  out[['colors']] <- colors
  out
  
  
}
