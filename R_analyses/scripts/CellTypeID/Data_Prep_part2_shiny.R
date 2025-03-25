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

sample_info <- read.csv("data/metadata_full.csv", header = F) 
colnames(sample_info) <- c('sample', "species", "ed", "sex", "folder", 'ref', "mito")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"
markers <- read.csv("data/markers/markers.csv")

rerun = F
if (rerun == T){
  seur_objs <- list()
  parm <- read.csv("data/seur_objs/cluster_params.csv")
  RDatas <- list.files("data/seur_objs/integrated/", full.names = T, pattern = ".RData")
  for (x in RDatas){
    print(x)
    load(x)
    si_tmp <- filter(sample_info, sample %in% unique(seurat_integrated$sample))
    species <- si_tmp$species[1]
    sex = si_tmp$sex[1]
    stage = si_tmp$stage[1]
    seurat_integrated@meta.data$celltype <- "Unknown"
    
    seurat_integrated$sex <- sex
    seurat_integrated$stage <- stage
    seurat_integrated$species <- species
    DefaultAssay(seurat_integrated) <- "integrated"
    seurat_integrated <- RunPCA(seurat_integrated, verbose = F)
    parm_tmp <- parm[parm$sex == sex & parm$stage == stage & parm$species == species,]
    
    pcs <- parm_tmp$PCs
    seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:pcs, verbose = F)
    seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:pcs, verbose = F)
    
    print("MAPS RUN")
    ######### Determining resolution to use --------
    seurat_integrated <- FindNeighbors(object=seurat_integrated, dims=1:pcs, verbose = F)
    seurat_integrated <- FindClusters(object=seurat_integrated, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                               0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                                               2.5, 3), verbose = F)
    Idents(seurat_integrated) <- parm_tmp$res
    

    si_tmp <- filter(sample_info, sample %in% unique(seurat_integrated$sample))
    seur_objs[[paste(species, sex, stage, sep = "_")]] <- seurat_integrated
  }
}
if (rerun == "load"){
  load("data/seur_objs/clusted_object.RData")
}

ui <- page_sidebar(
  # App title ----
  title = "Cell type selection",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    # Input: Slider for the number of bins ----
    selectInput(
      inputId = "seur_obj",
      label = "seurat object",
      choices = names(seur_objs),
     # selected = names(seur_objs)[[1]]
    ), 
    selectInput("cluster", "Select Cluster:", choices = NULL),
    selectInput("celltype", "Assign Cell Type:", choices = NULL),
    actionButton("update", "Update Cell Type"),
    
    selectInput("cluster_no", "Select Cluster:", choices = NULL),
    actionButton("cluster_no_update", "Update cluster resolution"),
    
    actionButton("save", "Save choices")
    
  ),
  # Output: Histogram ----
  plotOutput(outputId = "DotDim1"),
  plotOutput(outputId = "DotDim2"), 
  tableOutput("markerTable"),
  
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  sos <- reactiveValues(data = seur_objs)
  
  obj <- reactive({
    sos$data[[input$seur_obj]]
  })
  
  mks_tmp <- 
    reactive({
      markers %>% filter(substring(sex, 1, 1) %in% c("B", obj()$sex))
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
                        colnames(obj()@meta.data)[grepl("integrated_snn", colnames(obj()@meta.data))]
    )
    
  })
  
 
  output$DotDim1 <- renderPlot(
    res = 50,
    ggarrange(DimPlot(obj(), label = T, label.size = 8), 
              DotPlot(obj(), features = mks(), assay = "RNA") + coord_flip(), 
              clustree(obj(), verbose = F) %>% 
                plot() %>% 
                .$data %>% 
                group_by(integrated_snn_res.) %>% 
                summarise(no_clusters = n_distinct(cluster)) %>% 
                dplyr::rename(resolution = `integrated_snn_res.`) %>%
                ggplot(aes(x = resolution, y = no_clusters)) +
                geom_point() + theme(text=element_text(size=20)),
              ggtexttable(mks_tmp()),
              nrow = 1)
  )
  output$DotDim2 <- renderPlot({
    ggarrange(DimPlot(obj(), label = T, group.by = 'celltype', , label.size = 8),
              DotPlot(obj(), features = mks(), assay = "RNA", group.by = 'celltype') + coord_flip() + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    )
  })
  
  observeEvent(input$update, {
    seur_obj <- obj()
    seur_obj@meta.data$celltype[Idents(seur_obj) == input$cluster] <- input$celltype
    seur_objs[[input$seur_obj]] <<- seur_obj
    sos$data[[input$seur_obj]] <- seur_obj
  }) 
  
  observeEvent(input$cluster_no_update, {
    seur_obj <- obj()
    Idents(seur_obj) <- input$cluster_no
    seur_objs[[input$seur_obj]] <<- seur_obj
    sos$data[[input$seur_obj]] <- seur_obj
  })
  
  
}
shinyApp(ui = ui, server = server)
