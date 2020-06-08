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
library(patchwork)
library(readxl)

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

#sc = readRDS(file = "pbmc_2020-05-27.rds")
seurat_genes = sc@misc[["gene_lists"]]


PlotMystyle = function(p, title=NULL, col=NULL, legend_title=NULL, legend_position=NULL) {
    p = p + theme_light() + theme(panel.border = element_blank())
    if (!is.null(title)) p = p + ggtitle(title) #+ theme(plot.title = element_text(hjust=0.5))
    if (length(col) > 0) p = p + scale_fill_manual(values=col)
    if (!is.null(legend_title)) {
        p = p + labs(color=legend_title, fill=legend_title)
    } else {
        p = p + theme(legend.title = element_blank()) 
    }
    if (!is.null(legend_position)) p = p + theme(legend.position=legend_position)
    return(p)
}

# Define UI for application that draws a histogram
ui = fluidPage(
        titlePanel("Einzeldarstellungen von Genen"),

        sidebarPanel(
            fileInput("rds_file","Choose Seurat file:", accept = ".rds", buttonLabel = "Browse..."),
            selectInput("plots", "Plots:", c("FeaturePlot", "RidgePlot"), multiple = TRUE),
            selectInput("genes", "Gene:", seurat_genes, multiple = TRUE),
            downloadButton('download_plots')
        ),
        
        mainPanel(
            #splitLayout(cellWidths = c("50%","50%"),uiOutput('out_feature'), uiOutput('out_ridge'))
            tabsetPanel(type = "tabs",
                        # tabPanel("Feature", uiOutput('ui_feature')),
                        tabPanel("Feature", plotOutput('plot_feature')),
                        tabPanel("Ridge", uiOutput('ui_ridge'))
            )
        )
)

# Define server logic required to draw a histogram
server = function(input, output) {
    
    # # gen plot containers
    # output$ui_feature = renderUI({
    #     out_feature = list()
    # 
    #     if (length(input$genes)==0){return(NULL)}
    #     for (i in 1:length(input$genes)){
    #         out_feature[[i]] =  plotOutput(outputId = paste0("plot_feature",i))
    #     }
    #     return(out_feature)
    # })
    # # render plots
    # observe({
    #     for (i in 1:length(input$genes)){
    #         local({  #because expressions are evaluated at app init
    #             ii = i
    #             output[[paste0('plot_feature',ii)]] = renderPlot({
    #                 return(Seurat::FeaturePlot(sc, features=input$genes[[ii]], cols=c("lightgrey", param$col), combine=FALSE))
    #             })
    #         })
    #     }
    # })

    out$plot_feature = renderPlot({
      p = Seurat::FeaturePlot(sc, features=input$genes, cols=c("lightgrey", param$col), combine=FALSE)
      names(p) = input$genes
      for (i in features) p[[i]] = PlotMystyle(p[[i]], title=i)
      patchwork::wrap_plots(p, ncol=2)
    })

    #Ridge Plot
    output$ui_ridge = renderUI({
        out_ridge = list()
        
        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out_ridge[[i]] =  plotOutput(outputId = paste0("plot_ridge",i))
        }  
        return(out_ridge) 
    })
    observe({  
        for (i in 1:length(input$genes)){  
            local({  #because expressions are evaluated at app init
                ii = i 
                output[[paste0('plot_ridge',ii)]] = renderPlot({ 
                    return(Seurat::RidgePlot(sc, features=input$genes[[ii]], combine=FALSE))
                })
            })
        }                                  
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
