library(shiny)
library(bslib)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(clustree)


rm(list = ls())
sample_info <- read.csv("data/metadata_full.csv", header = F) 
colnames(sample_info) <- c('sample', "species", "ed", "sex", "folder", 'ref', "mito")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"
markers <- read.csv("data/markers/markers.csv") %>% 
  filter(Use == "Yes")
RDatas <- list.files("data/seur_objs/integrated/", full.names = T, pattern = ".RData")
RDatas <- RDatas[!grepl("MIXED", RDatas)]
#RDatas <- "data/seur_objs/integrated/Mallard_MIXED_MIXED_integrated_seurat.RData"
rerun = T

if (rerun == T & exists("DONE") == FALSE){
  seur_objs <- list()
  for (x in RDatas){
    print(x)
    load(x)
    si_tmp <- filter(sample_info, sample %in% unique(seurat_integrated$sample)) %>% 
      dplyr::select(sample, species, ed, sex, stage)
    new_md <- merge(si_tmp, seurat_integrated@meta.data, by.x = 'sample', by.y = 'sample')
    seurat_integrated@meta.data$species <- new_md$species
    seurat_integrated@meta.data$sex <- new_md$sex
    seurat_integrated@meta.data$stage <- new_md$stage
    
    species <- si_tmp$species %>% unique()
    sex = si_tmp$sex %>% unique()
    stage = si_tmp$stage %>% unique()
    
    seurat_integrated@meta.data$celltype <- "Unknown"
    seurat_integrated@meta.data$cluster_res <- "integrated_snn_res.0.1"
    DefaultAssay(seurat_integrated) <- "integrated"
    seur_objs[[paste(c(species, sex, stage), collapse = "_")]] <- seurat_integrated
  }
  clustrees <- lapply(seur_objs, function(x)(clustree(x) %>% 
                                               plot() %>%
                                               .$data %>%
                                               group_by(integrated_snn_res.) %>%
                                               summarise(no_clusters = n_distinct(cluster)) %>%
                                               dplyr::rename(resolution = `integrated_snn_res.`) %>%
                                               ggplot(aes(x = resolution, y = no_clusters)) +
                                               geom_point() + theme(text=element_text(size=20))))
  
  
  #save(clustrees, file = "data/clustrees.RData")
  rerun = F
  DONE = "DONE"
} 

load("data/seur_obj_celltype_metadata.RData")

change_metadata <- T
if (change_metadata == TRUE){
  for (i in names(save_seurats)){
    seur_objs[[i]]@meta.data <- save_seurats[[i]]
  }
}


ui <- shinyUI(fluidPage(
  # App title ----
  titlePanel("Cell type selection"),
  # Sidebar panel for inputs ----
  sidebarLayout(position = 'left', 
    sidebarPanel( width =2,
      selectInput("seur_obj", "seurat object", choices = names(seur_objs)), 
      selectInput("cluster", "Select Cluster:", choices = NULL),
      selectInput("celltype", "Assign Cell Type:", choices = NULL),
      actionButton("update", "Update Cell Type"),
      
      selectInput("cluster_no", "Select Cluster:", choices = NULL),
      actionButton("cluster_no_update", "Update cluster resolution"),
      #actionButton("reset_clus", "Reset clusters", class = "btn-warning"),
      
      actionButton("save", "Save choices", class = "btn-primary btn-lg")
    ),  # End sidebarPanel
  mainPanel( 
            fluidRow(
              column(5,plotOutput(outputId="DotDim", height = '1500px', width = '1500px')),  
              column(1,dataTableOutput(outputId="markertable"), offset = 6,),
            )
  )
              
  )))
  


# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  sos <- reactiveValues(data = seur_objs)
  
  obj <- reactive({
    sos$data[[input$seur_obj]]
  })
  
  obj_name <- reactive({input$seur_obj})
  
  mks_tmp <- 
    reactive({
      markers %>% filter(substring(sex, 1, 1) %in% c("B", obj()$sex)) %>% 
        .[,c("marker", "celltype", "sex", obj()$species[1])] %>% 
        .[.[,4] %in% rownames(obj()@assays$RNA),]
    }) 
  
  mks <- reactive({
    mks_tmp() %>%  .[,obj()$species[1]] %>% .[!is.na(.)] %>% unique() %>% 
      intersect(rownames(obj()@assays$RNA))
  })
  
  cells <- reactive({
    c("Unknown", mks_tmp() %>%  .[,'celltype'] %>% unique())
  })
  
  
  observe({
    updateSelectInput(session, "celltype", choices = cells())
    updateSelectInput(session, "cluster", choices = levels(Idents(obj())))
    updateSelectInput(session, "cluster_no", choices =  
                        colnames(obj()@meta.data)[grepl("integrated_snn", colnames(obj()@meta.data))], 
                      selected = unique(c(obj()$cluster_res[1])))
    
  })
  
  
  output$DotDim <- renderPlot(
    res = 50,
    ggarrange(
      ggarrange(DimPlot(obj(), label = T, label.size = 8, group.by = obj()$cluster_res[1]), 
              DotPlot(obj(), features = mks(), assay = "RNA", group.by = obj()$cluster_res[1])+ coord_flip(), 
              clustrees[[obj_name()]],
              nrow = 1), 
      ggarrange(DimPlot(obj(), label = T, group.by = 'celltype', , label.size = 5),
                DotPlot(obj(), features = mks(), assay = "RNA", group.by = 'celltype') + coord_flip() +
                  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)), 
                FeaturePlot(obj(), features = 'nFeature_RNA'), nrow = 1), 
      ggarrange(DimPlot(obj(), label = T, group.by = "Phase"), FeaturePlot(obj(), feature = 'mitoRatio'), nrow = 1),
      nrow = 3
  )
  )
  output$markertable <- DT::renderDT(
    mks_tmp(), options = list(pageLength = nrow(mks_tmp()))
  )
 

  observeEvent(input$update, {
    seur_obj <- obj()
    seur_obj@meta.data$celltype[Idents(seur_obj) == input$cluster] <- input$celltype
    seur_objs[[input$seur_obj]] <<- seur_obj
    sos$data[[input$seur_obj]] <- seur_obj
  }) 
  
  observeEvent(input$cluster_no_update, {
    seur_obj <- obj()
    Idents(seur_obj) <- input$cluster_no
    seur_obj@meta.data$cluster_res <- input$cluster_no
    seur_obj@meta.data$Ident_group <- input$cluster_no
    seur_objs[[input$seur_obj]] <<- seur_obj
    sos$data[[input$seur_obj]] <- seur_obj
  })
  
  observeEvent(input$save, {
    save_seurats[[input$seur_obj]] <<- obj()@meta.data
    save(save_seurats, file = "data/seur_obj_celltype_metadata.RData")
  })
  
  # observeEvent(input$reset_clus{
  #   seur_obj <- obj()
  #   seur_obj@meta.data$celltype <- "Unknown"
  #   seur_objs[[input$seur_obj]] <<- seur_obj
  #   sos$data[[input$seur_obj]] <- seur_obj
  # })
}

shinyApp(ui = ui, server = server)
