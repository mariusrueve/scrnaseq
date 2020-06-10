#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Seurat)
library(ggplot2)
library(readxl)

source("MyFeaturePlot.R")

param = list()
param$col = "palevioletred"

# read_excel_allsheets = function(filename, tibble = FALSE) {
#     # I prefer straight data.frames
#     # but if you like tidyverse tibbles (the default with read_excel)
#     # then just pass tibble = TRUE
#     sheets = readxl::excel_sheets(filename)
#     x = lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
#     if(!tibble) x = lapply(x, as.data.frame)
#     names(x) = sheets
#     x
# }
#sc_markers = read_excel_allsheets("/lager2/rcug_cd/scrnaseq/projects/pbmc/results/markers.xlsx")

sc = readRDS(file = "pbmc_2020-06-08.rds")
#seurat_genes = rownames(sc[[DefaultAssay(sc)]])
#gene_ensembl_combi = paste(seurat_genes," (",sc[["RNA"]][[]][,1],")", sep="")

features_names_ids = paste(rownames(sc[["RNA"]][[]]), " (", sc[["RNA"]][[]][,1], ")", sep = "")


# Define UI for application that draws a histogram
ui = fluidPage(
        titlePanel("Einzeldarstellungen von Genen"),

        sidebarPanel(
            fileInput("rds_file","Choose Seurat file:", accept = ".rds", buttonLabel = "Browse..."),
            selectInput("genes", "Gene:", features_names_ids, multiple = TRUE),
            downloadButton('download_plots')
        ),
        
        mainPanel(
            #splitLayout(cellWidths = c("50%","50%"),uiOutput('out_feature'), uiOutput('out_ridge'))
            tabsetPanel(type = "tabs",
                        tabPanel("Feature", uiOutput('ui_feature')),
                        #tabPanel("Feature", plotOutput('plot_feature_test')),
                        tabPanel("Ridge", uiOutput('ui_ridge')),
                        tabPanel("DotPlot", plotOutput('plot_dotplot'))
            )
        )
)

# Define server logic required to draw a histogram
server = function(input, output) {
  
  # observe({
  #   sc = readRDS(file = input$rds_file)
  #   seurat_genes = rownames(sc[[DefaultAssay(sc)]])
  #   updateSelectInput(
  #     session,
  #     "genes",
  #     choices=names(seurat_genes)
  #   )
  # })

    ########################################
    # gen plot containers
    output$ui_feature = renderUI({
        out_feature = list()

        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out_feature[[i]] =  plotOutput(outputId = paste0("plot_feature",i),
                                           width = "100%",
                                           height = "600px")
        }
        return(out_feature)
    })
    # render plots
    observe({
        for (i in 1:length(input$genes)){
            local({  #because expressions are evaluated at app init
                ii = i
                output[[paste0('plot_feature',ii)]] = renderPlot({
                  return(Seurat::FeaturePlot(sc, features=unlist(strsplit(input$genes[[ii]]," \\("))[c(T,F)], cols=c("lightgrey", param$col), combine=FALSE))
                })
            })
        }
    })

    ########################################
    # Test
    output$plot_feature_test = renderPlot({
      p = Seurat::FeaturePlot(sc, features=input$genes, cols=c("lightgrey", param$col), combine=FALSE)
      names(p) = input$genes
      for (i in input$genes) p[[i]] = PlotMystyle(p[[i]], title=i)
      patchwork::wrap_plots(p, ncol=1)
    })

    ########################################
    #Ridge Plot
    output$ui_ridge = renderUI({
        out_ridge = list()
        
        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out_ridge[[i]] =  plotOutput(outputId = paste0("plot_ridge",i),
                                         width = "100%",
                                         height = "600px")
        }  
        return(out_ridge) 
    })
    observe({  
        for (i in 1:length(input$genes)){  
            local({  #because expressions are evaluated at app init
                ii = i 
                output[[paste0('plot_ridge',ii)]] = renderPlot({ 
                    return(Seurat::RidgePlot(sc, features=unlist(strsplit(input$genes[[ii]]," \\("))[c(T,F)], combine=FALSE))
                })
            })
        }                                  
    })
    
    ########################################
    output$plot_dotplot = renderPlot({
      Seurat::DotPlot(sc, features=unlist(strsplit(input$genes," \\("))[c(T,F)], cols=c("lightgrey", param$col))
    })
    
    output$downloadPlot = downloadHandler(
        filename = function() { paste(input$dataset, '.png', sep='') },
        content = function(file) {
            ggsave(file, plot = plotInput(), device = "png")
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
