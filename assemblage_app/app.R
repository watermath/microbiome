
# Load libraries
# options(repos = "https://cran.rstudio.com/")
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinybusy)

library(tidyverse)
library(plotly)
library(cowplot)
library(ggrepel)
library(ggdendro)
library(vegan)
library(lubridate)
library(psych)
library(corrplot)

library(BiocManager)
options(repos = BiocManager::repositories())
library(NMF)

# set.seed(1961)
theme_set(theme_bw())
appalette<-c("#000000","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
            "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
# Define UI

ui <- dashboardPage(
    dashboardHeader(title = "Microbiome Community Assembly",
                    titleWidth = 400
                    ),

# Sidebar content
    dashboardSidebar(
        width=300,
        sidebarMenu(
            menuItem("Assemblages", tabName = "assemblages", icon = icon("th")),
            menuItem("Other plots", tabName = "other_plots", icon = icon("dashboard"))
        ),
        
        fileInput("raw_asv_counts",
                  span("Upload ASV/OTU counts table CSV",
                       tags$a("(Required Format)",
                              href="#",
                              onclick = "window.open('asv_input_format_example.png', 'newwindow'); return false;")
                  ),
                  accept = c(".csv"), multiple = FALSE, placeholder = "ASV  table"
        ),
        
        actionButton("rel_abund", label = "Counts to relative abundance",lib = "font-awesome", icon = icon("chart-pie")),
        actionButton("counts", label = "Relative abundance to counts",lib = "font-awesome", icon = icon("chart-simple")),
        
        
        fileInput("taxonomy", 
                  span("Upload taxonomy CSV",
                       tags$a("(Required Format)",
                              href="#",
                              onclick = "window.open('taxonomy_input_format_example.png', 'newwindow'); return false;")
                  ),
                  accept = c(".csv"),multiple = FALSE, placeholder = "Taxonomy"
        ),
        
        fileInput("metadata",
                  span("Upload metadata CSV",
                       tags$a("(Required Format)",
                              href="#",
                              onclick = "window.open('meta_input_format_example.png', 'newwindow'); return false;")
                  ),
                  accept = c(".csv"),multiple = FALSE, placeholder = "Metadata"
        ),
    
    selectInput("meta_filter",
                label="Filter",
                choices = "Choose metadata category",
                multiple = FALSE,
                selected = NULL
    ),
    
    selectInput("meta_filter_by",
                label="Filter by",
                choices = "Choose environment",
                multiple = FALSE,
                selected = NULL
    ),
    
    selectInput("meta_select",
                label="Select",
                choices = "Select metadata variable",
                multiple = FALSE,
                selected = NULL
    ),
    
    selectInput("meta_facet",
                label="Facet by",
                choices = "Facet categorical variable",
                multiple = FALSE,
                selected = NULL
    ),
    
    # dateRangeInput("daterange",label="Select date range"),
    
    actionButton("abdiv", label = "Diversity",lib = "font-awesome", icon = icon("globe")),
     
    sliderInput(inputId = "k1",
                label = "Choose k assemblages",
                value = 5, min=2, max=8, step = 1, ticks = FALSE),
    
    numericInput("ndom",label="Select number of most predominant taxa",value = 10, step = 1),
    
    actionButton("mod_asm", label = "Model community assemblages",lib = "font-awesome", icon = icon("group-arrows-rotate")),
    
    actionButton("plot_asm", label = "Plot assemblages",lib = "font-awesome", icon = icon("chart-column")),
    
    
    sliderInput(inputId = "k2",
                label = "Compare second set of assemblages",
                value = 5, min=2, max=8, step = 1,ticks = FALSE),
    
    actionButton("compare_asm", label = "Compare", icon = icon("code-compare"))
    
    
    #ADD MORE SIDEBAR THINGS HERE
    
    ),# end dashboardSidebar
    
#Body content
    dashboardBody(
        tabItems(
            tabItem(tabName = "assemblages",
                    
                    # fluidRow(box(verbatimTextOutput("test1"))),
                    # fluidRow(box(verbatimTextOutput("test2"))),
                    fluidRow(box(title="Environment filter:",textOutput("meta.filter.by"))),
                    fluidRow(
                        box(
                            title = paste("Assemblages and environmental variables"),
                            plotOutput("asmbar",height = "600", width = "100%"))
                        ),
                        fluidRow(box(
                            title = "Assemblage beta-diversity",
                            plotOutput("asmbetadiv", width = "100%")
                        )),
                    fluidRow(div(style='overflow-x: scroll;overflow-y: scroll;',uiOutput("domtaxabar", width = "100%" )))
                    )
            ),
        
        tags$img(src = "cwrs_logo.png", style="height:200px; width:100%")        
        )# end dashboardBody
    )

#Functions for server logic




# Define server
server <- function(input, output, session){
    cdata <- session$clientData
    rvs<-reactiveValues()
    
    # Read the input CSV files
    observeEvent(input$raw_asv_counts,{
      # showModal(modalDialog(
      #   title = "ASV/OTU table input format required",
      #   HTML('<img src="asv_input_format_example.png" />'),
      #   easyClose = TRUE,
      #   footer = NULL
      # ))
        rvs$raw.asv<-read.csv(input$raw_asv_counts$datapath, header = TRUE, row.names = 1, check.names = F)
        rvs$asv<-rvs$raw.asv
        })
    
    observeEvent(input$rel_abund,{
        rvs$asv<-as.data.frame(t(apply(rvs$raw.asv, 1, function(row) row/sum(row))))#count to relative abundance proportion
        showModal(modalDialog(
          title = "ASV counts converted to relative abundances.",
          paste0("Taxa are now proportions of read depth."),
          fade=TRUE,
          easyClose = TRUE,
          footer = NULL
        ))
    })
    
    observeEvent(input$counts,{
      rvs$asv<-rvs$raw.asv
      showModal(modalDialog(
        title = "ASVs are now counts",
        fade=TRUE,
        easyClose = TRUE,
        footer = NULL
      ))
    })
    
    observeEvent(input$taxonomy,{
        taxonomy<-read.csv(input$taxonomy$datapath, header = TRUE,check.names = F)
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
        raw.metadata <- read.csv(input$metadata$datapath,header=T,row.names = 1, check.names = FALSE)
        raw.metadata$Date<-as.Date(raw.metadata$Date, format="%Y-%m-%d")
        updateSelectInput(session, "meta_filter", label = "Filter", choices = names(raw.metadata))
        updateSelectInput(session, "meta_select", label = "Select", choices = names(raw.metadata))
        updateSelectInput(session, "meta_facet", label = "Facet", choices = names(raw.metadata))
        rvs$raw.metadata<-raw.metadata
        
    })#end metadata input
    
    
    
   #Filter metadata sample rows
    observeEvent(input$meta_filter,{
        if(is.null(rvs$raw.metadata)){return(NULL)}
        updateSelectInput(session, "meta_filter_by", label = "Filter by", choices = select(rvs$raw.metadata,input$meta_filter)%>%pull()%>%unique())
        
    })
    
    #Select metadata variable columns
    observeEvent(input$meta_select,{
        rvs$meta.select<-input$meta_select
        # output$test1<-renderPrint(rvs$meta.select)
    })
    
    #Facet metadata variable columns
    observeEvent(input$meta_facet,{
        rvs$meta.facet<-input$meta_facet
        # output$test2<-renderPrint(rvs$meta.facet)
    })
    
    #Combine metadata and asvs and filter by
    observeEvent(input$meta_filter_by,{
        if(is.null(rvs$raw.metadata)){return(NULL)}
        output$meta.filter.by<-renderText(input$meta_filter_by)
        
    }) #end filter by
    
    
    
    observeEvent(input$mod_asm,{
        withProgress(message = "Modelling assemblages",detail = "This may take a while depending on sample size and number of taxa",
        {#progress expr
          if(is.null(rvs$asv)|is.null(rvs$raw.metadata)|is.null(rvs$taxo)){return(NULL)}else{
            
            rvs$meta.asv<-full_join(rvs$raw.metadata%>%
                                      rownames_to_column("sampleId"),
                                    rvs$asv%>%
                                      rownames_to_column("sampleId"),by="sampleId")%>%
              filter(.data[[input$meta_filter]]%in%input$meta_filter_by)
            
                meta.asv<-rvs$meta.asv
                asv<-meta.asv%>%select(sampleId,names(rvs$asv))%>%column_to_rownames("sampleId")%>%
                  select(which(colSums(.)!=0))%>%
                  select(which(!is.na(colSums(.))))
                
    
    # Unsupervised NMF assemblage models 
    
            
            # if(input$mod_asm!="NMF"){return(NULL)}else{ # if statement if there are multiple modelling methods
                
                # NMF k1 assmeblage decomp and get weight and type matrices
                # k1<-6;k2<-8
                unsup.nmf1<-nmf(asv,input$k1,'KL')
                w.nmf1<-as.data.frame(basis(unsup.nmf1))
                T.eg1<-as.data.frame(coef(unsup.nmf1))
                
                # Assemblage weights to proportions and labelling
                asv.asm1<-as.data.frame((apply(T.eg1, 1, function(row) row/sum(row))))
                names(asv.asm1)<-paste("NMF",1:input$k1,sep="")
                rvs$asv.asm1<-asv.asm1
                
                samp.asm1 <- data.frame(t(apply(w.nmf1, 1, function(row) row/sum(row))))
                names(samp.asm1)<-paste("NMF",1:input$k1,sep="")
                rvs$samp.asm1<-samp.asm1
                
                # Add metadata categories
                meta<-meta.asv%>%
                    # rownames_to_column("sampleId")%>%
                    select(sampleId,rvs$meta.select,rvs$meta.facet)
                
                meta.wide1<-samp.asm1%>%
                    rownames_to_column("sampleId")%>%
                    left_join(meta, by="sampleId")%>%column_to_rownames("sampleId")
                rvs$meta.wide1<-meta.wide1
                
                # Pivot long the sample-asm for plotting
                meta.long1<-pivot_longer(meta.wide1,cols = c(1:input$k1), names_to="asm",values_to="mix")
                rvs$meta.long1<-meta.long1
                
                
                
                
                
                # saveRDS(list(samp.asm=samp.asm1, asv.asm=asv.asm1,metawide.asm=meta.wide1,metalong.asm=meta.long1),paste("nmf_asm",input$meta_filter_by,".rds",sep = ""))
                
            # }#end NMF if statement
            
            
        # return(list(samp.asm=samp.asm1, asv.asm=asv.asm1,metawide.asm=meta.wide1,metalong.asm=meta.long1)) #use this if using eventReactive
          }#end if statement check for inputs asvs, taxonomy, metadata
    }) #end progress msg
})#end nmf modelling output
    
    observeEvent(input$plot_asm,{
      if(is.null(rvs$meta.long1)){return(NULL)}else{
            asmbar<-ggplot(rvs$meta.long1, aes(x=.data[[rvs$meta.select]], y=mix, fill=asm)) + 
                geom_bar(stat="identity")+
                facet_grid(.data[[rvs$meta.facet]]~.)+
                theme(axis.text.x = element_text(angle=45, hjust=1),
                      plot.title = element_text(size=12, face = "bold"),
                      axis.title = element_text(size = 12, face = "bold"),
                      strip.text = element_text(size = 12, face = "bold",color = "black"),
                      strip.background =element_rect(fill="white"),
                      legend.title = element_text(size=12,face = "bold"),
                      legend.text = element_text(size=12,face = "bold"), 
                      axis.text = element_text(size = 12, face = "bold")
                )+
                xlab(paste(rvs$meta.select))+ylab("Mixture Proportion")+
                guides(fill=guide_legend(title = "Assemblage"))+
                scale_fill_manual(name = "Assemblage", values = appalette)
            # geom_vline(xintercept = c(unique(floor_date(meta.long1$Date,unit = "year"))), lty="dashed",lwd=1, colour="grey50")
            if(is.Date(rvs$meta.select)){
                asmbar<-asmbar+scale_x_date(
                    # date_breaks = "1 months", date_labels = "%b-%y"
                )}
            output$asmbar<-renderPlot({asmbar})#end asm barplot1 output
        # output$testprint<-renderPrint({str(samp.asm1%>%rownames_to_column("sampleId"))})
        
        
        #filter predom asvs in asm by cutoff value   
        # dom.cut<-0.05
        # dom.asv<-asv.asm1%>%filter(if_any(everything(), ~ . > dom.cut))
        #or filter predom asvs in asm by top n
        dom.asv<-rvs$asv.asm1%>%rownames_to_column("asvId")%>%
            pivot_longer(cols = -asvId)%>%
            group_by(name)%>%
            slice_max(order_by = value,n=input$ndom)%>% #n predom asvs
            pivot_wider()
        
        dom.taxa<-left_join(dom.asv,rvs$taxo,by="asvId")%>%
            select(-asvId)%>% #run here to get wide table
            pivot_longer(cols=-c(full_name))%>%filter(!(full_name%in%c("   ","    ASV: 234 Confidence: 0.95","metagenome")))
        
        taxa.axis<-unique(dom.taxa$full_name)
        
        output$domtaxabar<-renderUI({
            
            renderPlot({
                
                ggplot(dom.taxa, aes(x=name, y=factor(full_name,levels = taxa.axis), fill = value)) +
                    geom_tile()+
                    # coord_flip()+
                    scale_fill_gradient2(low = "lightgrey", high = "black", na.value = "white")+
                    ylab("Taxon")+xlab("Assemblage")+labs(fill = "Mixing Weight \n Proportion")+
                scale_y_discrete(limits=rev)+
                theme_bw()+
                theme(axis.title = element_text(size = 12, face = "bold"),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
                
            })
            
        })
        
        meta.asv<-rvs$meta.asv
        asv<-meta.asv%>%select(sampleId,names(rvs$asv))%>%column_to_rownames("sampleId")
        asv.bc.mds = metaMDS(asv, distance = "bray", trymax = 200)#gives stress
        
        mds.scores = as.data.frame(scores(asv.bc.mds)$sites)
        
        samp.asm<-rvs$samp.asm1
        samp.asm$dom.asm<-colnames(samp.asm)[max.col(samp.asm)]
        samp.asm$dom.asm.wt<-do.call(pmax,samp.asm%>%select(-dom.asm))
        
        # output$testtable<-renderDataTable({head(samp.asm%>%rownames_to_column("sampleId"))})
        # output$testtable2<-renderDataTable({head(mds.scores%>%rownames_to_column("sampleId"))})
        
        mds.scores$dom.asm<-samp.asm$dom.asm
        mds.scores$dom.asm.wt<-samp.asm$dom.asm.wt
        
        
        output$asmbetadiv<-renderPlot({
            ggplot(data = mds.scores, aes(x = NMDS1, y = NMDS2)) +
                geom_point(data=mds.scores,aes(colour=factor(dom.asm)),shape=19, size = 3,alpha=0.5)+
                geom_text_repel(aes(label=ifelse(dom.asm.wt<=0.5,as.character(round(dom.asm.wt,digits=2)),'')), min.segment.length = 0,
                                # hjust=0.1,vjust=0, 
                                fontface="bold") +
                # geom_segment(data = en_coord_cont,aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), size=1, alpha = 0.5, colour = "grey30")+
                # geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(en_coord_cont)) +
                scale_colour_manual(
                    # values = c("gray75","grey40","black"),
                    values = appalette,
                    guide= guide_legend(override.aes = list(shape = 15)))+
                # scale_shape_manual(values=c(3,17,16)) +
                theme_light()+
                theme(plot.title = element_text(size=12),axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                      # panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
                      # axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
                      legend.title = element_blank(),
                      legend.text = element_text(size = 12, colour = "grey30"))
                # geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), shape = "diamond", size = 5, alpha = 0.6, colour = "navy") +
                # geom_text_repel(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04),label = row.names(en_coord_cat), fontface = "bold")+ggtitle(paste("Beta ","(stress:",asv.bc.mds$stress,")",sep=""))
        })#end betadiv plot
        
      }#end if statement check on mod_asm
    }) 
     
    
    
}


# Run the Shiny App
shinyApp(ui = ui, server = server)
