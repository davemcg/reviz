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
library(shinyjs)
library(cowplot)
library(glue)
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

load_exon_fst_data <- function(metadata, mapper, gene, srs_acc=NA, study_acc=NA){
    # srs_acc = head(metadata$sample_accession)
    # study_acc= NA
    # gene= c('PAX6', 'VSX1')
    ## All the ID lookups could be be sped up with a Hashtable, but might add to load times
    t_sample_ids <- c()
    if(any(!is.na(srs_acc))){
        t_sample_ids <- c( t_sample_ids, filter(metadata, sample_accession %in% srs_acc) %>% pull(rail_id) )
    }
    if(any(!is.na(study_acc))){
        t_sample_ids <- c(t_sample_ids,  filter(metadata, study_accession %in% study_acc) %>% pull(rail_id) )
    }
    t_snaptron_id <- mapper %>% filter(gene_name %in% gene) %>% pull(snaptron_id)# haven't checked if snaptron_id<>gene_name is 1:1
    fst_file <- glue('app_data/exon_fst/{gene}_.fst')
    count_data <- read_fst(fst_file) %>% filter(snaptron_id %in% t_snaptron_id, sample_id %in% t_sample_ids)
    count_data %>% 
        dplyr::rename(rail_id = sample_id) %>%  
        inner_join(mapper) %>% 
        left_join(metadata) %>% 
        mutate(counts = as.numeric(counts),
               counts = replace_na(counts, 0)) 
    
}
# plot_data <- load_exon_fst_data(CORE_METADATA, EXON_NAME_MAPPER,  'VSX1', 
#                                 srs_acc =head(ALL_SRS_ACCESSIONS))



## make example user-proivided metadata
# CORE_METADATA %>% arrange(desc(Age_Days)) %>% head(100) %>% 
#     mutate(Age = as.numeric(Age_Days) %>% replace_na(-1)) %>% 
#     mutate(Developmental_time = case_when(
#         Age >= 90 ~ 'Oldest',
#         Age  >= 60 & Age < 90 ~ 'Oldish',
#         Age < 60 ~ 'Not Old'
#     )) %>% select(sample_accession, Developmental_time) %>% 
#     write_tsv('app_data/example_user_provided_metadata.tsv'

#*****************************
#**CHANGE TO APPROPRIATE DIR**
#*****************************
setwd('/data/swamyvs/reviz/viz/app_data/example_data/')

A <- Sys.time()
CORE_METADATA <- read_fst('app_data/eiad_snaptron_sample_table.fst')
#RECOUNT_PROJ_INFO <- read_fst('app_data/recount3_proj_info.fst')
GENE_NAME_MAPPER <- read_fst('app_data/snaptron_id_2_gene_name.fst')
EXON_NAME_MAPPER <- read_fst('app_data/snaptron_id_2_exon_name.fst')
load('app_data/valid_gene_exon_names.Rdata')# loads all_gene_names, all_exon_names
load('app_data/samples_and_studies.Rdata')# all_study_accessions, all_srs_accessions 
GENE_DB_FILE <- 'sql_files/recount3All_gene_reformatted.db'
# ^this must match output of `reformat_snaptron_db-setup_appdata.sh`
METACOLS <- c( 'Tissue', 'Sub_Tissue', 'Origin')

B <- Sys.time()
#print(B-A) # Time difference of 0.7401347 secs
CON<- dbConnect(SQLite(), dbname = GENE_DB_FILE)



ui <- fluidPage(
    useShinyjs(),
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
                                 conditionalPanel("input.samplecheck == '_ERROR_'",{
                                     p(id = "element", "You have requested accessions which are not part of our database, but we will add them. Please try again in a few minutes")
                                 }),
                                 conditionalPanel("input.samplecheck == '_GOOD_'",{
                                     InteractiveComplexHeatmapOutput("exon_hm")
                                 }),
                                 hidden(
                                     selectInput('samplecheck', '',choices = NULL, multiple = F)    
                                 )
                        )
                                 
        )
    
    
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    cat(getwd(),file = stderr(), sep = '\n' )
    #input <- list(proj_ids = 'SRP163355', vp_gene_input  = "PAX6$",hm_gene_input= "PAX6", vp_color_choice = 'Sub_Tissue')
    updateSelectizeInput(session,inputId = 'proj_ids', choices = ALL_STUDY_ACCESSIONS, server = T, selected = 'SRP070148' )
    updateSelectizeInput(session,inputId = 'acc_ids', choices = ALL_SRS_ACCESSIONS, server = T, selected = ALL_SRS_ACCESSIONS[1:10])
    updateSelectizeInput(session,inputId = 'vp_gene_input', choices = ALL_GENE_NAMES, server = T , selected = c('PAX6', 'VSX1') )
    updateSelectizeInput(session,inputId = 'hm_gene_input', choices = ALL_EXON_NAMES, server = T , selected = 'PAX6')
    updateSelectizeInput(session,inputId = 'samplecheck', choices = c('_ERROR_','_GOOD_' ), selected = '_GOOD_', server = T )
    
    
    
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
        if(!all(new_meta[,1]%in% ALL_SRS_ACCESSIONS)){# if we are missing samples 
            updateSelectizeInput(session,inputId = 'samplecheck', choices = c('_ERROR_','_GOOD_' ), selected = '_ERROR_', server = T )
        }
    })
    
    observeEvent(input$draw_vp, {
        output$vp_plot <- renderPlot({
            if(length(input$user_metadata) >0){
                file <- isolate(input$user_metadata)
                new_meta <- read_tsv(file$datapath)
                CORE_METADATA_plot <- CORE_METADATA %>% inner_join(new_meta) 
                
            }else{
                CORE_METADATA_plot <- CORE_METADATA
            }
            
            plot_data <- pull_data_from_db(CORE_METADATA_plot, GENE_NAME_MAPPER, 'gene_counts_long',input$vp_gene_input, input$acc_ids, input$proj_ids)
            

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
        if(length(input$user_metadata) >0){
            file <- isolate(input$user_metadata)
            new_meta <- read_tsv(file$datapath)
            CORE_METADATA_plot <- CORE_METADATA %>% inner_join(new_meta) 
        }else{
            CORE_METADATA_plot <- CORE_METADATA
        }
        cat(getwd(), file = stderr())
        plot_data <- load_exon_fst_data(CORE_METADATA_plot, EXON_NAME_MAPPER, input$hm_gene_input, input$acc_ids, input$proj_ids)
        
        id_cols <- colnames(plot_data)[!colnames(plot_data)%in% c('counts', 'rail_id', 'sample_accession', 'run_accession')]
        hm_plot_data <-  plot_data %>% pivot_wider(names_from = sample_accession, values_from = counts, id_cols = id_cols ) 
        mat <-hm_plot_data %>% select(-id_cols) %>% as.matrix %>% {log2(. +1)}
        rownames(mat) <- hm_plot_data$snaptron_id
        mat[is.na(mat)] <- 0
        
        ht <- Heatmap(mat,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = F, 
                      name = 'log2(exon counts+1)',
                      col = viridis::viridis(n=100)
        )
        ht <- draw(ht)
        makeInteractiveComplexHeatmap(input, output, session, ht)
    })
}    

# Run the application 
shinyApp(ui = ui, server = server) 





