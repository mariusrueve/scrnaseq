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

read_excel_allsheets <- function(filename, tibble = TRUE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}
sc_markers = read_excel("/lager2/rcug_cd/scrnaseq/projects/pbmc/results/markers.xlsx")
seurat_genes = sc_markers[["gene"]]

# Define UI for application that draws a histogram
ui <- fluidPage(
        titlePanel("Einzeldarstellungen von Genen"),

        sidebarPanel(
            selectInput("genes", "Gene:", seurat_genes, multiple = TRUE),
            downloadButton('foo')
        ),
        
        mainPanel(
            splitLayout(cellWidths = c("50%","50%"),uiOutput('out_umap'), uiOutput('out_ridge'))
        )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    reactive({
        sc_markers = read_excel_allsheets("/lager2/rcug_cd/scrnaseq/projects/pbmc/results/markers.xlsx")
        seurat_genes = sc_markers[["gene"]]
    })
    
    output$out_umap = renderUI({
        out = list()
        
        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out[[i]] <-  plotOutput(outputId = paste0("plot_umap",i))
        }  
        return(out) 
    })
    observe({  
        for (i in 1:length(input$genes)){  
            local({  #because expressions are evaluated at app init
                ii <- i 
                output[[paste0('plot_umap',ii)]] <- renderPlot({ 
                        return(FeaturePlot(sc, features=input$genes[[ii]], cols=c("lightgrey", param$col), combine=FALSE))
                })
            })
        }                                  
    })

        
    #Ridge Plot
    output$out_ridge = renderUI({
        out = list()
        
        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out[[i]] <-  plotOutput(outputId = paste0("plot",i))
        }  
        return(out) 
    })
    observe({  
        for (i in 1:length(input$genes)){  
            local({  #because expressions are evaluated at app init
                ii <- i 
                output[[paste0('plot',ii)]] <- renderPlot({ 
                    return(RidgePlot(sc, features=input$genes[[ii]], combine=FALSE))
                })
            })
        }                                  
    })
    
    output$foo = downloadHandler(
        filename = 'test.png',
        content = function(file) {
            device <- function(..., width, height) {
                grDevices::png(..., width = width, height = height,
                               res = 300, units = "in")
            }
            ggsave(file, plot = output$out_ridge , device = device)
        })
}

# Run the application 
shinyApp(ui = ui, server = server)
