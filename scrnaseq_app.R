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

seurat_genes = sc.markers[["gene"]]

# Define UI for application that draws a histogram
ui <- fluidPage(
        titlePanel("Einzeldarstellungen von Genen"),

        sidebarPanel(
            selectInput("genes", "Gene:", seurat_genes, multiple = TRUE),
        ),
        
        mainPanel(
            splitLayout(cellWidths = c("50%","50%"),uiOutput('out_umap'), uiOutput('out_ridge'))
        )
)



# Define server logic required to draw a histogram
server <- function(input, output) {

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
}

# Run the application 
shinyApp(ui = ui, server = server)
