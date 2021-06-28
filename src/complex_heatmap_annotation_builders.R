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
  for (i in unique(df$transcript_id)){
    df[,i] <- grepl(i, df$transcript_id)
  }
  out <- list()
  out$df <- df %>% select(-transcript_id, -recount_exon_id) %>% data.frame()
  out$df[out$df == TRUE] <- 'TRUE'
  out$df[out$df == FALSE] <- 'FALSE'
  colors <- list()
  for (i in colnames(df)){
    colors[[i]] <- c('TRUE' = 'black', 'FALSE' = 'white')
  }
  out$colors <- colors
  out
}
# heatmap
## row metadata annotation
## same output as above
tx_row_df_maker <- function(meta){
  df <- data.frame(Genotype = t(gene_exon_tpm) %>% as_tibble(rownames = 'name') %>% 
                     left_join(meta, by = 'name') %>% pull(Genotype))
  colors <- c(pals::alphabet(n=26), pals::alphabet2(n=26), pals::glasbey(n=32))
  colors <- rep(colors, 100)
  names(colors) <- unique(df[,1])
  colors <- colors[!is.na(names(colors))]
  out <- list()
  out[['df']] <- df
  out[['colors']] <- colors
  out
  
  
}