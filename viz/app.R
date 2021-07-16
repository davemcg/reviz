#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(fst)
library(DBI)
library(RSQLite)
suppressMessages(library(recount3))
suppressMessages(library(cowplot))
#setwd('/data/swamyvs/reviz/viz) #setwd to app 
### run these to make fst's for app, change file paths accordingly

## Format for the new long sqlite file is `snaptron_id,sample_id,counts` (snaptron_id is a custom gene identifier)
## both snaptron_id and sample_id will need to be mapped to gene_name and SRA_accession accordingly
## theoretically this could be done upstream during the db (with little channge in run time) but right now this is where we're at
make_EiaD_st <- function(){
    
    EiaD_st <-  read_tsv('/data/swamyvs/EiaD_build/sampleTable_2019-05-25_tissues.tab') %>% 
        rename(sample_accession = `#sample_accession`)
    snaptron_st <- data.table::fread('/data/swamyvs/reviz/recount3_sra3h_samples.tsv') %>% as_tibble 
    snaptron_gtex_st <- data.table::fread('/data/swamyvs/reviz/recount3_gtex_samples.tsv') %>% as_tibble
    all_snaptron <- bind_rows(snaptron_st%>% select(rail_id, run_accession=external_id),
                              snaptron_gtex_st %>% select(rail_id, run_accession = run_acc))
    eiad_in_snaptron <- EiaD_st %>% inner_join(all_snaptron) %>% mutate(rail_id = as.character(rail_id))
    write_fst(eiad_in_snaptron, 'app_data/eiad_snaptron_sample_table.fst')
    ALL_STUDY_ACCESSIONS <- unique(eiad_in_snaptron$study_accession)
    ALL_SRS_ACCESSIONS <- unique(eiad_in_snaptron$sample_accession)
    save(ALL_STUDY_ACCESSIONS, ALL_SRS_ACCESSIONS, file = 'app_data/samples_and_studies.Rdata')
}
#make_EiaD_st()
cache_Snaptron_proj_info <- function(){
    r3_projects <- available_projects(organism = 'human' ) %>% as_tibble()
    write_fst(r3_projects, 'app_data/recount3_proj_info.fst')
    
}
## run this to pull raw snaptron sqlite files and reformat to long indexed sql files 
## https://github.com/vinay-swamy/reviz/blob/demo/reformat_snaptron_db.sh

## pull snaptron_id to gene_name from original snaptron sqlite file
## sqlite3 sql_files/sra3vh_genes.sqlite "SELECT snaptron_id, right_annotated FROM intron" > viz/app_data/snaptronid2gene_name.txt    
## sqlite3 sql_files/sra3vh_exons.sqlite "SELECT snaptron_id, right_annotated FROM intron" > viz/app_data/snaptronid2exon_name.txt

make_snaptron_id_converters <- function(){
    gene_mapper <- read_delim('app_data/snaptronid2gene_name.txt', delim = '|', col_names = c('snaptron_id', 'raw')) %>% 
        mutate(gene_id = str_split(raw, ':') %>% sapply(function(x) x[1]), 
               gene_name = str_split(raw, ':') %>% sapply(function(x) x[2]), 
               gene_type = str_split(raw, ':') %>% sapply(function(x) x[3])) %>% 
        select(-raw)
    ## I'm not totally sure how to interpret the way they annotate exons, should revisit this at a later date.
    exon_mapper <- read_delim('app_data/snaptronid2exon_name.txt', delim = '|', col_names = c('snaptron_id', 'raw')) %>% 
        mutate(left = str_split(raw, ':') %>% sapply(function(x) x[1]), 
               right = str_split(raw, ':') %>% sapply(function(x) x[2]),
               gene_id = left
               ) %>% 
        left_join(gene_mapper %>% select(-snaptron_id)) %>% 
        mutate(gene_name = replace(gene_name, is.na(gene_name), gene_id[is.na(gene_name)]),
               gene_type = replace(gene_name, is.na(gene_name), 'missing gene_name/unclear annotation') ) %>% 
        select(-raw, -left, -right)
    gene_mapper <- gene_mapper %>% mutate(snaptron_id = as.character(snaptron_id))
    exon_mapper <- exon_mapper %>% mutate(snaptron_id = as.character(snaptron_id))
    write_fst(gene_mapper, 'app_data/snaptron_id_2_gene_name.fst')
    write_fst(exon_mapper, 'app_data/snaptron_id_2_exon_name.fst')
    ALL_GENE_NAMES <- unique(gene_mapper$gene_name)
    ALL_EXON_NAMES <- unique(exon_mapper$gene_name)
    save(ALL_GENE_NAMES, ALL_EXON_NAMES, file = 'app_data/valid_gene_exon_names.Rdata')
}
#make_snaptron_id_converters()

pull_data_from_db <- function(metadata, mapper, table_name, gene, srs_acc=NA, study_acc=NA){
    # srs_acc = head(metadata$sample_accession)
    # study_acc= NA
    # gene= c('PAX6', 'VSX1')
    qdata <- tbl(CON, table_name)
    ## All the ID lookups could be be sped up with a Hashtable, but might add to load times
    t_sample_ids <- c()
    if(any(!is.na(srs_acc))){
        t_sample_ids <- c( t_sample_ids, filter(metadata, sample_accession %in% srs_acc) %>% pull(rail_id) )
    }
    if(any(!is.na(study_acc))){
        t_sample_ids <- c(t_sample_ids,  filter(metadata, study_accession %in% study_acc) %>% pull(rail_id) )
    }
    t_snaptron_id <- mapper %>% filter(gene_name %in% gene) %>% pull(snaptron_id)# haven't checked if snaptron_id<>gene_name is 1:1
    count_data <- qdata %>% filter(snaptron_id %in% t_snaptron_id, sample_id %in% t_sample_ids) %>% collect()
    count_data %>% 
        dplyr::rename(rail_id = sample_id) %>%  
        inner_join(mapper) %>% 
        left_join(metadata) %>% 
        mutate(counts = as.numeric(counts),
               counts = replace_na(counts, 0)) 

}
# plot_data <- pull_data_from_db(CORE_METADATA, GENE_NAME_MAPPER, 'gene_counts_long', c('PAX6', 'VSX1'), 
#                                srs_acc =head(ALL_SRS_ACCESSIONS), study_acc = 'SRP070148')
#plot_data <- pull_data_from_db(CORE_METADATA, EXON_NAME_MAPPER, 'exon_counts_long', 'PAX6', srs_acc =head(ALL_SRS_ACCESSIONS))

## make example user-proivided metadata
# CORE_METADATA %>% arrange(desc(Age_Days)) %>% head(100) %>% 
#     mutate(Age = as.numeric(Age_Days) %>% replace_na(-1)) %>% 
#     mutate(Developmental_time = case_when(
#         Age >= 90 ~ 'Oldest',
#         Age  >= 60 & Age < 90 ~ 'Oldish',
#         Age < 60 ~ 'Not Old'
#     )) %>% select(sample_accession, Developmental_time) %>% 
#     write_tsv('app_data/example_user_provided_metadata.tsv')



A <- Sys.time()
CORE_METADATA <- read_fst('app_data/eiad_snaptron_sample_table.fst')
RECOUNT_PROJ_INFO <- read_fst('app_data/recount3_proj_info.fst')
GENE_NAME_MAPPER <- read_fst('app_data/snaptron_id_2_gene_name.fst')
EXON_NAME_MAPPER <- read_fst('app_data/snaptron_id_2_exon_name.fst')
load('app_data/valid_gene_exon_names.Rdata')# loads all_gene_names, all_exon_names
load('app_data/samples_and_studies.Rdata')# all_study_accessions, all_srs_accessions 
DB_FILE <- 'app_data/EiaD_in_snaptron_reformatted.db'
METACOLS <- c( 'Tissue', 'Sub_Tissue', 'Origin')
B <- Sys.time()
#print(B-A) # Time difference of 0.7401347 secs
CON<- dbConnect(SQLite(), dbname = DB_FILE)



ui <- fluidPage(
    titlePanel('rEvIZ'),
    column(4,
        fileInput('user_metadata', 'Enter a csv with meta data'),
        actionButton('upload_user_data', 'Press to Upload Metadata'),
        selectizeInput('proj_ids', 'Select a SRA Project ID', choices = NULL, multiple = T),
        
        selectizeInput('acc_ids', 'Enter SRA Accession IDs', choices = NULL, multiple = T)
    ),
    column(8,   
        tabsetPanel(type = 'tabs',
                        tabPanel('Gene Expression Violin Plot',
                                 selectizeInput('vp_gene_input', 'Select a gene(s)', choices =NULL, multiple=T),
                                 selectizeInput('vp_color_choice', 'Select Metadata to color on',
                                                choices = METACOLS,
                                                selected = 'Sub_Tissue',
                                                multiple=F),
                                 actionButton('draw_vp', 'Draw Violin Plots'),
                                 plotOutput('vp_plot')
                        ),
                        tabPanel('Exon Heatmap',
                                 selectizeInput('hm_gene_input', 'Select one gene', choices = NULL, multiple =F),
                                 actionButton('draw_hm', 'Draw Heatmap'),
                                 InteractiveComplexHeatmapOutput("exon_hm")
                        )
        )
    
    
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    #input <- list(proj_ids = 'SRP163355', vp_gene_input  = "PAX6$",hm_gene_input= "PAX6", vp_color_choice = 'Sub_Tissue')
    updateSelectizeInput(session,inputId = 'proj_ids', choices = ALL_STUDY_ACCESSIONS, server = T, selected = 'SRP070148' )
    updateSelectizeInput(session,inputId = 'acc_ids', choices = ALL_SRS_ACCESSIONS, server = T, selected = ALL_SRS_ACCESSIONS[1:10])
    updateSelectizeInput(session,inputId = 'vp_gene_input', choices = ALL_GENE_NAMES, server = T , selected = c('PAX6', 'VSX1') )
    updateSelectizeInput(session,inputId = 'hm_gene_input', choices = ALL_EXON_NAMES, server = T , selected = 'PAX6')
    
    
    
    observeEvent(input$upload_user_data, {
        file <- isolate(input$user_metadata)
        new_meta <- read_tsv(file$datapath)
        CORE_METADATA <- CORE_METADATA %>% left_join(new_meta) 
        CORE_METADATA[is.na(CORE_METADATA)] <- '.'
        .GlobalEnv$CORE_METADATA <- CORE_METADATA
        ALL_SRS_ACCESSIONS <- unique(CORE_METADATA$sample_accession)
        updateSelectizeInput(session,inputId = 'acc_ids', choices = ALL_SRS_ACCESSIONS, server = T, selected = new_meta$sample_accession)
        METACOLS <- c(METACOLS, colnames(new_meta)[-1])
        updateSelectizeInput(session,inputId = 'vp_color_choice', choices = METACOLS, server = T )
    })
    
    observeEvent(input$draw_vp, {
        
        output$vp_plot <- renderPlot({
            if(length(input$user_metadata) >0){
                file <- isolate(input$user_metadata)
                new_meta <- read_tsv(file$datapath)
                CORE_METADATA <- CORE_METADATA %>% inner_join(new_meta) 
                cat('inside\n', file= stderr())
                cat(colnames(CORE_METADATA), file= stderr())
            }
            
            plot_data <- pull_data_from_db(CORE_METADATA, GENE_NAME_MAPPER, 'gene_counts_long',input$vp_gene_input, input$acc_ids, input$proj_ids)
            
            # for(i in input$user_metadata){
            #     cat(i, file = stderr())
            # }
            
            cat('outside plot data \n', file= stderr())
            cat(colnames(plot_data), file= stderr())
            cat('outside core_tight data \n', file= stderr())
            cat(colnames(CORE_METADATA), file= stderr())
            
            ggplot(plot_data, aes(x=!!as.symbol(input$vp_color_choice), 
                                  y=log2(counts+1), 
                                  fill = !!as.symbol(input$vp_color_choice)
                                  )) +
                geom_violin()+
                geom_boxplot(width = .1) +
                xlab('') + 
                ylab('log2(count+1)') +
                facet_wrap(~gene_name) +
                theme_cowplot()

        })
    })
    observeEvent(input$draw_hm, {
        id_cols <- colnames(plot_data)[!colnames(plot_data)%in% c('counts', 'rail_id', 'sample_accession', 'run_accession')]
        hm_plot_data <-  plot_data %>% pivot_wider(names_from = sample_accession, values_from = counts, id_cols = id_cols ) 
        mat <-hm_plot_data %>% select(-id_cols) %>% as.matrix %>% {log2(. +1)}
        rownames(mat) <- hm_plot_data$snaptron_id
        
        ht <- Heatmap(mat,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = F, 
                      name = 'log2(xxon counts+1)',
                      col = viridis::viridis(n=100)
        )
        ht <- draw(ht)
        makeInteractiveComplexHeatmap(input, output, session, ht)
    })
}    

# Run the application 
shinyApp(ui = ui, server = server)





