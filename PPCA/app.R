#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load libraries
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinybusy)
library(shinyWidgets)

library(tidyverse)
library(plotly)
library(cowplot)

library(lubridate)

library(PoissonPCA)
# library(factoextra)

# options(shiny.maxRequestSize=30*1024^2) 
# set.seed(1961)
theme_set(theme_bw())
appalette<-c("#000000","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
             "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

ui <- dashboardPage(
    dashboardHeader(title = "Microbiome PPCA"),
    
    # Sidebar content
    dashboardSidebar(
        
        sidebarMenu(
            menuItem("Plotting", tabName = "plotting", icon = icon("th")),
            menuItem("Other widgets", tabName = "other", icon = icon("dashboard"))
        ),
        
        fileInput("raw_counts",
                  span("Upload ASV/OTU sample counts table",
                       tags$a("(Required Format .csv)",
                              href="#",
                              onclick = "window.open('asv_input_format_example.png', 'newwindow'); return false;")
                  ),
                  accept = c(".csv"), multiple = FALSE, placeholder = "No samples here yet"
        ),
        
        fileInput("metadata", 
                  span("Select metadata",
                       tags$a("(Required Format .tsv)",
                              href="#",
                              onclick = "window.open('meta_input_format_example.png', 'newwindow'); return false;")
                  ),
                  accept = c(".csv", ".txt", ".tsv"),multiple = FALSE, placeholder = "No metadata here yet"
        ),
        
        
        selectizeInput("disc_select",
                    label="Select categorical variable",
                    choices = "No metadata yet",
                    multiple = FALSE
        ),    
        
        selectizeInput("disc_filter_by",
                    label="Filter samples by",
                    choices = "Filter samples",
                    multiple = TRUE,
                    selected = NULL
        ),
        
        selectizeInput("cts_select",
                    label="Select continuous variable",
                    choices = "No metadata yet",
                    multiple = FALSE
        ),  
        # dateRangeInput("daterange",label="Select date range"),
        
        numericInput("cutoff",label="Sparsity threshold",value = 0.1),
        
        selectizeInput("disc_order",
                    label="Select categorical order",
                    choices = "No variable selected",
                    multiple = TRUE
        ),
        
        selectizeInput("ind_samples",
                    label="Select sample(s)",
                    choices = "No sample(s) selected",
                    multiple = TRUE),
        
        actionButton("ppca_scores", label = "Plot PPCA scores",lib = "font-awesome", icon = icon("table-cells")),
        
        actionButton("taxaplot", label="Plot taxa",lib = "font-awesome", icon = icon("chart-column")),
        
        downloadBttn(
            outputId = "downloadRes",
            label="Download results",
            style = "bordered",
            color = "primary"
        ),
        
        downloadBttn(
            outputId = "downloadPlot",
            label="Download plots",
            style = "bordered",
            color = "primary"
        )
        # downloadButton("downloadRes",label= "Download")
        
        
        
    ),#end dashboardSidebar
    
    dashboardBody(
        
        tabItems(
            tabItem(tabName = "plotting",
                    # fluidRow(column(12,dataTableOutput("test_table"))),
                    # fluidRow(box(textOutput("test"))),
                    # fluidRow(box(title="Company:",textOutput("meta.filter.by"))),
                    fluidRow(column(width=12,plotlyOutput("scores_plot")))
                    
                    
            )#end tabItem
            
        ),#end endtabItems
        
        tags$img(src = "logo_luminultra.png", style="height:129.5px; width:100%")
    )#end dashboardBody
    
    
    
    
    
    
    
    
    
    
)#end dashboardPage


server<-function(input, output, session){
    cdata <- session$clientData
    rvs<-reactiveValues()

    # Read the input CSV files
    observeEvent(input$raw_counts,{
        rvs$raw.asv<-read_csv(input$raw_counts$datapath)
        # rvs$asv<-rvs$raw.asv
    })
    
    
    observeEvent(input$taxonomy,{
        taxonomy<-read_csv(input$taxonomy$datapath)
        if(names(taxonomy)[1]=="otuId"){rename(taxonomy,asvId=otuId)}
        rvs$taxo<-mutate(taxonomy,full_name=paste(
            # Kingdom,
            # Phylum,
            # Class,
            Order,
            Family,
            Genus,
            Species
            ,"ASV:",str_sub(asvId,1,3)
            # ,"Confidence:",signif(as.numeric(Confidence),2)
        ),
        .before="Kingdom")%>%
            select(asvId,full_name)%>%
            # mutate(full_name = str_remove_all(full_name, paste(remove.taxa, collapse = "|")))%>%
            mutate(full_name=trimws(gsub("\\s+", " ",full_name)))%>%
            mutate(full_name=trimws(gsub("[A-Z]_[0-9]__", " ",full_name)))%>%
            mutate(full_name=trimws(gsub("[a-z]__", " ",full_name)))%>%
            mutate(full_name=trimws(gsub("[a-z]_", " ",full_name)))
        
        
        
        # Replace asvIds in counts/prop table with matching taxonomy full name
        # temp.asv<-raw.asv()
        # names(temp.asv)<-taxonomy$Full_name[match(names(temp.asv), taxonomy$asvId)]
        # temp.asv
        
    })#end taxonomy input
    
    observeEvent(input$metadata, {
        raw.metadata <- read_tsv(input$metadata$datapath)%>%
            mutate(date=as.Date(date,format="%d-%b-%y"))
        
        updateSelectInput(session, "disc_select", choices = names(raw.metadata))
        updateSelectInput(session, "cts_select", choices = names(raw.metadata))
        updateSelectizeInput(session, "ind_samples", choices = unique(raw.metadata$sampleId), server = TRUE)
        rvs$raw.metadata<-raw.metadata
        
    })#end metadata input
    
    
    #Select metadata discrete variable column
    observeEvent(input$disc_select,{
        rvs$disc_select<-input$disc_select
        if(!is.null(rvs$raw.metadata)){
            meta.init<-rvs$raw.metadata%>%pull(rvs$disc_select)%>%unique()%>%sort()
            updateSelectInput(session, "disc_filter_by", choices = meta.init)
            updateSelectInput(session, "disc_order", choices = meta.init)
        }else{
            updateSelectInput(session, "disc_filter_by",choices = NULL)
            updateSelectInput(session, "disc_order",choices = NULL)
        }
    })
    
    observeEvent(input$disc_order,{
        rvs$disc_order<-input$disc_order
        #TEST
        # output$test<-renderText(rvs$disc_order)
    })
    
    #Combine metadata and asvs and filter by
    observeEvent(input$disc_filter_by,{
        if(is.null(rvs$raw.metadata)){return(NULL)}
        rvs$disc_filter_by<-input$disc_filter_by
        output$disc.filter.by<-renderText(input$disc_filter_by)
    }) #end filter by
    
    #Select metadata continuous variable column
    observeEvent(input$cts_select,{
        rvs$cts_select<-input$cts_select
    })
    
    
    
    
    observeEvent(input$ppca_scores,{
        withProgress(message="Microbiome PCA modelling and plotting",detail= "This can take some time depending on sample size and number of ASVs/OTUs",{
            if(is.null(rvs$raw.asv)|is.null(rvs$raw.metadata)#|is.null(rvs$taxo)
               ){return(NULL)}else{
            asv<-rvs$raw.asv
            metadata<-rvs$raw.metadata%>%
                filter(
                    .data[[input$disc_select]]%in%c(input$disc_filter_by,input$ind_samples)
                    &sampleId%in%asv$sampleId
                    )%>%
                arrange(sampleId)%>%
                mutate(cts_var=if_else(.data[[input$cts_select]]%in%c("X","#DIV/0!"),NA,if_else(.data[[input$cts_select]]%in%c("<LOD"),0,as.numeric(.data[[input$cts_select]]))),.after = .data[[input$cts_select]])
            
            filt.asv<-asv%>%
                filter(sampleId%in%metadata$sampleId)%>%
                arrange(sampleId)%>%
                column_to_rownames("sampleId")%>%
                select(which(colSums(.)>quantile(colSums(.),input$cutoff)))
            # select(which(sapply(., var)>quantile(sapply(., var),0.8)))
            
            # Poisson PCA
            k<-10
            ppca<-Poisson_Corrected_PCA(filt.asv,
                                            k=k,
                                            # k=dim(filt.asv)[2],
                                            transformation="log",
                                            seqdepth="minvar" #for log transform
                                            # seqdepth="compositional"
                )
            rownames(ppca$scores)<-metadata$sampleId
            ppca_ve<-(ppca$sdev^2)/sum(ppca$sdev^2)*100 # variance explained
            # cumsum(ppca_ve) # cumulative variance explained
            Lambda<-ppca$means
            ppcl<-ppca$loadings

            ppcs_meta<-data.frame(ppca$scores)%>%rownames_to_column("sampleId")%>%
                left_join(select(metadata,sampleId,date,email,cts_var), by="sampleId")
            names(ppcs_meta)[2:(k+1)]<-paste0("PPC",1:k)
            rvs$ppcs_meta<-ppcs_meta
            # output$test_table<-renderDataTable({ppcs_meta})
            scores_plot<-ggplot(ppcs_meta,
                                aes(x=PPC1,y=PPC2, colour=cts_var,
                                    sampleId=sampleId,
                                    email=email,
                                    date=date
                                ))+
                geom_point(size=3,alpha=(2/4))+
                xlab(paste0("PPC1 (",round(ppca_ve[1],1),"%)"))+
                ylab(paste0("PPC2 (",round(ppca_ve[2],1),"%)"))+
                ggtitle(paste0(round(cumsum(ppca_ve)[2],1),"% variance explained"))
                theme( #legend.position = "top",
                    title = element_text(size = 16, face = "bold") ,
                    axis.title = element_text(size = 16, face = "bold"),
                    axis.text = element_text(size = 16, face = "bold"),
                    strip.text = element_text(size = 14, face = "bold"),
                    legend.title =element_text(size = 14, face = "bold"),
                    legend.text = element_text(size = 12, face = "bold")
                )

            rvs$scores_plot<-scores_plot
                
            output$scores_plot<-renderPlotly({ggplotly(scores_plot, tooltip = c("sampleId","email","date"))}) # use tooltip to specify the ggplot aes you want to include in hover text
            } #end if
        }) #end warning progress
    }) #end observeEvent
    
    
    
    output$downloadRes <- downloadHandler(
        filename = function() {
            paste0("ppca_results_", Sys.Date(), ".csv")
        },
        content = function(file) {
            write.csv(rvs$ppcs_meta, file)
        }
    )
    
    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste0("ppca_scores_plot_", Sys.Date(), ".pdf")
        },
        content = function(file) {
            ggsave(filename = file, plot=rvs$scores_plot, height = 12, width = 12)
        }
    )
    
    
    
} #end server

# Run the application 
shinyApp(ui = ui, server = server)
