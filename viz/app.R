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
suppressMessages(library(recount3))
suppressMessages(library(cowplot))
# For making example data
# library(snapcount)
# library(SummarizedExperiment)
# get_snapcount_data_long <- function(db, gene, qfunc) {
#     sb <- QueryBuilder(compilation=db, regions=gene)
#     snapcount_raw_data <- qfunc(sb)
#     snap_count_long_data <- assay(snapcount_raw_data) %>%
#         as.matrix() %>%
#         t() %>%
#         as_tibble(rownames = 'name') %>%
#         left_join(colData(snapcount_raw_data) %>% as_tibble(rownames = 'name')) #%>%
#         #dplyr::rename(Expression = `643`)
#     return(snap_count_long_data)
# }
# get_snapcount_data_long('gtex', 'PAX6', query_exon) %>% write_tsv('example_app_data_exon.tsv')
# k <- get_snapcount_data_long('gtex', 'PAX6', query_exon)
# 


# - boxplot (one gene only????)
# - user input for metadata file (tsv with run_accession / metadata col) ?????
# - heatmap exon view (violin_heatmap_demo.R) has some code bits for this
# - bonus: splice graph exon view??????
    
genes <-  c('PAX6$', 'CRYGA', 'CRYB', 'CRYA', 'HSF4','PROX1$', 'RARB','SIX3$', 'PRX$')
#data <- data.table::fread('/Users/swamyvs/NIH/reviz/example_app_data.tsv') %>% as_tibble 
# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel('rEvIZ'),
    column(4,
        fileInput('user_metadata', 'Enter a csv with meta data'),
        selectizeInput('proj_ids', 'Select a SRA Project ID', choices = c('SRP163355',"SRP151763" ), multiple = F, 
                       selected = 'SRP163355')
    ),
    column(8,   
        tabsetPanel(type = 'tabs',
                        tabPanel('Gene Expression Violin Plot',
                                 selectizeInput('vp_gene_input', 'Select a gene(s)', choices =genes, multiple=T, 
                                                selected = genes[1]),
                                 plotOutput('vp_plot')
                        ),
                        tabPanel('Exon Heatmap',
                                 selectizeInput('hm_gene_input', 'Select one gene', choices = genes, multiple =F,
                                                selected = genes[1]),
                                 InteractiveComplexHeatmapOutput("exon_hm")
                        )
        )
    
    
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    #input <- list(proj_ids = 'SRP163355', vp_gene_input  = "PAX6$")
    output$vp_plot <- renderPlot({
        file <- input$user_metadata
        
        meta <- read_tsv(file$datapath)
        colnames(meta) <- c('name', 'Condition')
        r3_projects <- available_projects(organism = 'human' )
        # grab project info 
        proj_info <- subset(
            r3_projects,
            project == input$proj_ids & project_type == "data_sources"
        )

        rse_study_gene <- create_rse(proj_info, type = 'gene')
        assay(rse_study_gene, "counts") <- transform_counts(rse_study_gene)
        assays(rse_study_gene)$TPM <- recount::getTPM(rse_study_gene, length_var = "score")
        counts <- data.frame(assays(rse_study_gene)$counts)
        gene_quant <- counts[rowRanges(rse_study_gene) %>% 
                                                        as_tibble(rownames = 'id') %>% 
                                                        filter(grepl(paste0(input$vp_gene_input, collapse = '|'), 
                                                                     gene_name, ignore.case= TRUE)) %>% 
                                                        pull(gene_id) %>% unique(), ]
        gene_quant %>%
            as_tibble(rownames = 'gene_id') %>% 
            pivot_longer(cols = contains('SRR')) %>% 
            left_join(rowRanges(rse_study_gene) %>% as_tibble(), by = 'gene_id') %>% 
            left_join(meta, by = 'name') %>% 
        ggplot(aes(x=Condition, y=log2(value+1), fill = Condition)) +
            #geom_point() + 
            geom_boxplot() +
            #scattermore::geom_scattermore(pointsize = 3, position = position_dodge(width = 1)) +
            #geom_boxplot(width = 0.3)+
            theme_cowplot() +
            guides(x =  guide_axis(angle = 90)) +
            scale_fill_manual(values = c(pals::brewer.set1(n=12), pals::brewer.set2(n=10)) %>% unname()) +
            xlab('Gene') + facet_wrap(~gene_name, ncol = 5, scales = 'free_x')
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
