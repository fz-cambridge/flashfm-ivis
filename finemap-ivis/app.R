# finemap-ivis version_2022_03_10

####----Libraries----
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(dplyr)
library(ggplot2)
library(shiny)
library(shinydashboard) 
library(plotly)
library(DT)
library(gaston)
library(corrplot)
library(RColorBrewer)
#display.brewer.all(type = 'qual')
library(randomcoloR)
library(igraph)
library(networkD3)
library(stringr)

# if pre-load data
#load("data/R_data.RData")
####----pre-loaded data as examples----
#users will overwrite it if they upload new files
# finemap.config <- read.csv("data/finemap.config", sep="")
# finemap.ld <- read.table("data/finemap.ld", quote="\"", comment.char="")
# finemap.snp <- read.csv("data/finemap.snp", sep="")
# finemap.z <- read.csv("data/finemap.z", sep="")

####----Max size of input file----
options(shiny.maxRequestSize = 100*1024^2) 

#generate a vector of color based on snpGroups----
set.seed(1)
n <- 4
col_vector = distinctColorPalette(n)
pal <- c("gray", col_vector)
pal <- setNames(pal, c("0", "cs1@99.5%+@99%+@95%+@90%", "cs1@99.5%+@99%+@95%", "cs1@99.5%+@99%", "cs1@99.5%"))
# pal_2 <- c("gray", "dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue")
# pal_2 <- setNames(pal_2, c("0", "cs1@99.5%+@99%+@95%+@90%", "cs1@99.5%+@99%+@95%", "cs1@99.5%+@99%", "cs1@99.5%"))

#Step 1: function from finemap-example.R
format.finemap <- function(fconfig="finemap1.config",fsnp="finemap1.snp",fz="finemap1.z",fld="finemap.ld") {
  # constuct modPP, a data.frame with models and their PPs, equivalent to mpp.pp$PP[[1]][,1]
  fmresults <- read.table(fconfig, header = TRUE, as.is = TRUE, sep = " ")
  modPP <- data.frame( str = fmresults$config, PP = fmresults$prob, stringsAsFactors = FALSE)     
  modsnps <- strsplit(modPP$str, ",")
  snpmods <- sapply(modsnps, function(x) paste(x, collapse = "%"))
  modPP$str <- snpmods
  # gwas summary statistics and snp MPPs
  snpPP <- read.table(fsnp, header = TRUE, as.is = TRUE, sep = " ")        
  MPP <- snpPP[,c("rsid","prob")]
  # LD, defined as cor^2
  r <- read.table(fld, header = FALSE, as.is = TRUE, sep = " ")   
  gwas <- read.table(fz,header = TRUE, as.is = TRUE, sep = " ")   # needed to get snp names for ld matrix  
  #gwas$p <- 2*(1-pnorm(abs(gwas$beta/gwas$se)))
  gwas$p <- 2*(pnorm(abs(gwas$beta/gwas$se), lower.tail = F))
  r2 <- r^2
  colnames(r2) <- gwas$rsid
  rownames(r2) <- gwas$rsid
  return(list(PP=modPP,MPP=MPP,LD=r2, gwas=gwas))
}
#fstub='finemap'
#tt <- format.finemap(fconfig=paste0(fstub,".config"),fsnp=paste0(fstub,".snp"),fz=paste0(fstub,".z"),fld=paste0(fstub,".ld"))  

# Step 2: define credible sets based on a specific significant level (i.e. 0.99)----
credset <- function(modPP,cred=0.99) {
  tmp <- modPP[order(modPP, decreasing = TRUE)]
  cpp <- cumsum(tmp)
  wh <- which(cpp <= cred)
  if (!length(wh)) wh <- 1
  wh <- c(wh, max(wh) + 1)
  keepmodPP <- tmp[wh]
  mods <- names(keepmodPP)
  usnps <- unique(unlist(strsplit(mods,"%")))
  return(usnps)
}
M <- 1 # number of traits 
# #cs1_selected <- csM_selected <- vector("list",M)
# cs1_selected <- vector("list",M)
# for(i in 1:M){
#   #tmp <- mpp.pp$PP[[i]][,1] # trait i, single-trait fine-mapping
#   tmp <- tt$PP[,-1]
#   names(tmp) <- tt$PP[,1] 
#   cs1_selected[[i]] <- credset(tmp,0.99)
#   #tmp <- mpp.pp$PP[[i]][,2] # trait i, multi-trait fine-mapping
#   #csM_selected[[i]] <- credset(tmp,0.99)
# }



####----Shiny user interface----
ui <- dashboardPage(
  #Dashboard header----
  #Dashboard header----
  dashboardHeader(
    title = "finemap-ivis R package",
    # Drop-down menu for messages
    dropdownMenu(type = "messages", badgeStatus = "success",
                 headerText = 'Contact info:',
                 messageItem("Contact 1: Feng",
                             "feng.zhou[AT]mrc-bsu.cam.ac.uk"
                 ),
                 messageItem("Contact 2: Jenn",
                             "jennifer.asimit[AT]mrc-bsu.cam.ac.uk"
                 )
    ),
    
    # Drop-down menu for notifications
    dropdownMenu(type = "notifications", #badgeStatus = "warning",
                 headerText = 'Update and version info:',
                 notificationItem(icon = icon("users"), status = "info",
                                  "Update: this is Version 0.3"
                 ),
                 notificationItem(icon = icon("users"), status = "info",
                                  "Update: modified on 10-Mar-2022"
                 )
    ),
    
    # Drop-down menu for tasks, with progress bar
    dropdownMenu(type = "tasks", badgeStatus = "danger",
                 headerText = 'Source code:',
                 taskItem(value = 100, color = "aqua",
                          text = "Source code #1: flashfm", 
                          href = 'https://github.com/jennasimit/flashfm'
                 ),
                 taskItem(value = 100, color = "green",
                          text = "Source code #2: flashfm-ivis", 
                          href = 'https://github.com/fz-cambridge/flashfm-ivis'
                 ),
                 taskItem(value = 100, color = "yellow",
                          text = "Source code #3: finemap-ivis", 
                          href = 'http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/'
                 )
    )
  ),
  
  #Dashboard sidebar----
  dashboardSidebar(
    #Data inputs----
    radioButtons(inputId = "control_widgets_0", 
                 label = h4("Select data source:"),
                 choices = list("Pre-loaded example" = 0, 
                                "Upload 4 key FINEMAP outputs" = 1), 
                 selected = 0),
    
    fileInput("inFiles",
              label="Upload FINEMAP results",
              multiple = TRUE),
    
    #Control widgets----
    radioButtons(inputId = "control_widgets_1", 
                 label = h4("Credible sets:"),
                 choices = list("Show all SNPs (default)" = 0, 
                                "Show cs1@99.5%" = 1,
                                "Show cs1@99%" = 2, 
                                "Show cs1@95%" = 3, 
                                "Show cs1@90%" = 4, 
                                "Show not in cs1@99%" = 5), 
                 selected = 0),
    
    sliderInput(inputId = "control_widgets_2", 
                label = h4("Define MPP value"), 
                min = 0, max = 1, value = c(0, 1)),   
    
    #Dashboards----
    sidebarMenu(
        menuItem("Dashboard", icon = icon("dashboard"), startExpanded = TRUE,
              menuSubItem("Input values", tabName = "ch_0", 
                          icon = icon("angle-right")),   
              menuSubItem("Output plots", tabName = "ch_1", 
                          icon = icon("angle-right")),
              menuSubItem("README", tabName = "README", 
                          icon = icon("angle-right")),
              menuSubItem(text = "Go back to flashfm-ivis", 
                          href = "http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/", 
                          icon = icon("angle-right"))
        )
    )
  ),
  
  #Dashboard body----
  dashboardBody(  
    tabItems(
      #Dashboard_ch_0----
      tabItem(tabName = "ch_0",
              #check FINEMAP input files
              h2("Input values"),
              h3("Check files in the current temporary server"),
              h4("--- (at least 4 files with ext: .config; .ld; .snp; .z)---"),
              fluidRow(column(12, verbatimTextOutput("value"))),
              br(),
              h4("--- 4 key files' location in the server ---"),
              fluidRow(column(12, DT::dataTableOutput("file_locations"))),
              br(),
              h4("--- tt$PP ---"),
              fluidRow(column(12, DT::dataTableOutput("contents_tt_PP"))),
              br(),
              h4("--- tt$MPP ---"),
              fluidRow(column(12, DT::dataTableOutput("contents_tt_MPP"))),
              br(),
              h4("--- tt$LD ---"),
              fluidRow(column(12, DT::dataTableOutput("contents_tt_LD"))),
              br(),
              h4("--- tt$gwas ---"),
              fluidRow(column(12, DT::dataTableOutput("contents_tt_gwas"))),
              br(),
              h4("--- credible sets at different confidence intervals ---"),              
              fluidRow(column(12, verbatimTextOutput("credible_set"))),
      ),
      
      #Dashboard_ch_1----
      tabItem(tabName = "ch_1",
              h2("Output plots"),
              fluidRow(
                column(# width should be between 1 and 12
                  width=12,
                  tabBox(
                    #title = "Result_1 charts",
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "tabset1", height = "1200px",
                    tabPanel("Regional_Association_Plot_1", 
                             plotlyOutput("manhattan_1a", height = "700px")
                             ),
                    tabPanel("Regional_Association_Plot_2", 
                             plotlyOutput("manhattan_2", height = "700px"),
                             sliderInput(inputId = "control_widgets_4", 
                                         label = h4("Define user's own credible set"), 
                                         min = 0.90, max = 1, step =0.001, value = 0.99),
                             sliderInput(inputId = "control_widgets_6", 
                                         label = h4("Define MPP value"), 
                                         min = 0, max = 1, step =0.005, value = c(0, 1)),  
                             radioButtons(inputId = "control_widgets_5", 
                                          label = h4("Define the plot color:"),
                                          inline=T,
                                          choices = list("Blue" = 0, 
                                                         "Red" = 1,
                                                         "Pink" = 2, 
                                                         "Green" = 3, 
                                                         "Orange" = 4, 
                                                         "Black" = 5), 
                                          selected = 0)
                    ),
                    tabPanel("LD_matrix", 
                             box(plotlyOutput("Plot_LD1"),
                                 title = "LD Heatmap",
                                 #status = "warning",
                                 solidHeader = TRUE, collapsible = TRUE,
                                 width=NULL),
                             box(plotOutput("Plot_LD2"),
                                 title = "LD Plot", 
                                 solidHeader = TRUE, collapsible = TRUE,
                                 width=NULL)
                             ),
                    tabPanel("Venn_diagram", 
                             plotlyOutput("fig_ch4_tab2_venn_cs1"),
                             #tableOutput("table_ch4_tab2_venn_cs1"),
                             box(div(style = 'overflow-y: scroll; overflow-x: scroll; ', 
                                     DT::dataTableOutput('table_ch4_tab2_venn_cs1_table')),
                                 title = "Table of credible sets",
                                 #status = "warning", 
                                 solidHeader = TRUE,
                                 collapsible = TRUE, width=NULL)
                             ),
                    tabPanel("SNP_network", 
                             radioButtons(inputId = "control_widgets_3b_cs", 
                                          label = h4("Credible sets:"),
                                          choices = list("Show all SNPs" = 0, 
                                                         "Show cs1@99.5%  (default)" = 1,
                                                         "Show cs1@99%" = 2,
                                                         "Show cs1@95%" = 3,
                                                         "Show cs1@90%" = 4), 
                                          selected = 1),
                             sliderInput(inputId = "control_widgets_3b_finemap",
                                         label = h4("Define mpp.pp$PP value (finemap)"),
                                         min = 0, max = 1, step =0.001, value = c(0.05, 1),
                                         width = '800px'),
                             #textOutput("control_widgets_3b_output"),
                             forceNetworkOutput(outputId = "ch1_network_snp_finemap"),
                             downloadButton('d', 'Download network as html')
                    ),
                    width=NULL
                  )
                )
              )
      ),

      #Dashboard_ch_3----
      tabItem(tabName = "README",
              h2("README"),
              br(),
              h3("About this package:"),
              br(),
              h4("GitHub: "),
              br(),
              h4("flashfm R package: https://github.com/jennasimit/flashfm"),
              h4("flashfm-ivis R package: https://github.com/fz-cambridge/flashfm-ivis"),
              br(),
              h4("YouTube: "),
              #br(),
              #HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/videoseries?list=PLcy5X5WM9r0AqRrEb5vUTfRkcBRUoOc8T" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
              fluidRow(
                fluidRow(
                  column(6, 
                         box(
                           width = NULL, 
                           title = "flashfm-ivis and finemap-ivis overview", 
                           HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/hDbV9qNzsZo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                         )),
                  column(6, 
                         box(
                           width = NULL, 
                           title = "finemap-ivis user data input", 
                           HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/dTRdfsMumyY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                         )),
                )
              )
              
      )
    )
  )
) #End of user interface


####----Shiny server----
server <- function(input, output) {
  
  output$value <- renderPrint({
    if (input$control_widgets_0 == 0){
      print("Using pre-loaded example...")
    }else{
      str(input$inFiles)
    }
    })
  
  tt <- reactive({
    if (input$control_widgets_0 == 0){
      fstub='data/finemap'
      #format.finemap(fconfig=get(paste0(fstub,".config")),fsnp=get(paste0(fstub,".snp")),fz=get(paste0(fstub,".z")),fld=get(paste0(fstub,".ld")))  
      format.finemap(fconfig=paste0(fstub,".config"),fsnp=paste0(fstub,".snp"),fz=paste0(fstub,".z"),fld=paste0(fstub,".ld"))  
    }else{
    file <- input$inFiles
    req(file)
    for (i in 1:length(file$datapath)){
      ext_check <- tools::file_ext(file$datapath[i])
      if (ext_check == "config"){fconfig = file$datapath[i]}
      if (ext_check == "ld"){fld = file$datapath[i]}
      if (ext_check == "snp"){fsnp = file$datapath[i]}
      if (ext_check == "z"){fz = file$datapath[i]}
    }
    format.finemap(fconfig=fconfig,fsnp=fsnp,fz=fz,fld=fld) 
    }
  })
  
  output$file_locations <- DT::renderDataTable({
    if (input$control_widgets_0 == 0){
        fconfig = "Using pre-loaded example..."
        fld = "Using pre-loaded example..."
        fsnp = "Using pre-loaded example..."
        fz = "Using pre-loaded example..."
        location = rbind(fconfig, fld, fsnp, fz)
        colnames(location) = c('location')
      DT::datatable(location, options = list(dom = 't'))
    }else{
    file <- input$inFiles
    req(file)
    for (i in 1:length(file$datapath)){
      ext_check <- tools::file_ext(file$datapath[i])
      if (ext_check == "config"){fconfig = file$datapath[i]}
      if (ext_check == "ld"){fld = file$datapath[i]}
      if (ext_check == "snp"){fsnp = file$datapath[i]}
      if (ext_check == "z"){fz = file$datapath[i]}
    }
    location = rbind(fconfig, fld, fsnp, fz)
    colnames(location) = c('location')
    DT::datatable(location, options = list(dom = 't'))
    }
  })
  
  output$contents_tt_PP <- DT::renderDataTable({
    DT::datatable(data.frame(tt()$PP), 
                  extensions = 'Buttons',
                  options = list(
                    paging = TRUE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    scrollX=TRUE,
                    ordering = TRUE,
                    dom = 'lftiprBRSPQ',
                    buttons = list( 
                      list(extend = 'csv',   filename =  "tt_PP"),
                      list(extend = 'excel', filename =  "tt_PP"))
                  ),
                  class = "display")
  })
  
  output$contents_tt_MPP <- DT::renderDataTable({
    DT::datatable(data.frame(tt()$MPP), 
                  extensions = 'Buttons',
                  options = list(
                    paging = TRUE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    scrollX=TRUE,
                    ordering = TRUE,
                    dom = 'lftiprBRSPQ',
                    buttons = list( 
                      list(extend = 'csv',   filename =  "tt_PP"),
                      list(extend = 'excel', filename =  "tt_PP"))
                  ),
                  class = "display")
  })
  
  output$contents_tt_LD <- DT::renderDataTable({
    DT::datatable(data.frame(round(tt()$LD[1:5,1:5], 6)),
                  options = list(dom = 't'))
  })
  
  GWAS_all_final <- reactive({
    df_mpp = data.frame(tt()$MPP)
    colnames(df_mpp) = c("rsid","mpp")
    merge(data.frame(tt()$gwas),df_mpp, by.x=c("rsid"), by.y=c("rsid"))
  })
  
  # output$contents_tt_gwas <- DT::renderDataTable({
  #   DT::datatable(GWAS_all_final(),              
  #                 extensions = 'Buttons',
  #                 options = list(
  #                   paging = TRUE,
  #                   searching = TRUE,
  #                   fixedColumns = TRUE,
  #                   autoWidth = TRUE,
  #                   scrollX=TRUE,
  #                   ordering = TRUE,
  #                   dom = 'lftiprBRSPQ',
  #                   buttons = list( 
  #                     list(extend = 'csv',   filename =  "tt_PP"),
  #                     list(extend = 'excel', filename =  "tt_PP"))
  #                 ),
  #                 class = "display")
  # })

  output$credible_set <- renderPrint({
    cl = c(0.995, 0.99, 0.95, 0.90)
    cs1_selected <- vector("list",length(cl))
    for(i in 1:length(cl)){
      tmp <- tt()$PP[,-1]
      names(tmp) <- tt()$PP[,1]
      cs1_selected[[i]] <- credset(tmp,cl[i])
    }
    names(cs1_selected) = c("cs1 @ 0.995","cs1 @ 0.99","cs1 @ 0.95","cs1 @ 0.90")
    print(cs1_selected)
  })
  
  GWAS_all_final_cs <- reactive({
    cl = c(0.995, 0.99, 0.95, 0.90)
    cs1_selected <- vector("list",length(cl))
    for(i in 1:length(cl)){
      tmp <- tt()$PP[,-1]
      names(tmp) <- tt()$PP[,1]
      cs1_selected[[i]] <- credset(tmp,cl[i])
    }
    GWAS_all_final_cs_0 = GWAS_all_final()
    GWAS_all_final_cs_0$p_log10 = -log10(GWAS_all_final_cs_0$p)
    for (i in 1:length(GWAS_all_final_cs_0$rsid)){
      if (GWAS_all_final_cs_0$rsid[i] %in% cs1_selected[[1]]){
        GWAS_all_final_cs_0$'cs1@99.5%'[i] = 1
      } else {GWAS_all_final_cs_0$'cs1@99.5%'[i] = 0}
      if (GWAS_all_final_cs_0$rsid[i] %in% cs1_selected[[2]]){
        GWAS_all_final_cs_0$'cs1@99%'[i] = 1
      } else {GWAS_all_final_cs_0$'cs1@99%'[i] = 0}
      if (GWAS_all_final_cs_0$rsid[i] %in% cs1_selected[[3]]){
        GWAS_all_final_cs_0$'cs1@95%'[i] = 1
      } else {GWAS_all_final_cs_0$'cs1@95%'[i] = 0}
      if (GWAS_all_final_cs_0$rsid[i] %in% cs1_selected[[4]]){
        GWAS_all_final_cs_0$'cs1@90%'[i] = 1
      } else {GWAS_all_final_cs_0$'cs1@90%'[i] = 0}      
    }
    
    for (i in 1:length(GWAS_all_final_cs_0$rsid)){
      if (GWAS_all_final_cs_0$'cs1@90%'[i] == 1){
        GWAS_all_final_cs_0$cs1_group[i]  = "cs1@99.5%+@99%+@95%+@90%"
      } else if (GWAS_all_final_cs_0$'cs1@95%'[i] == 1){
        GWAS_all_final_cs_0$cs1_group[i]  = "cs1@99.5%+@99%+@95%"
      } else if (GWAS_all_final_cs_0$'cs1@99%'[i] == 1){
        GWAS_all_final_cs_0$cs1_group[i]  = "cs1@99.5%+@99%"
      } else if (GWAS_all_final_cs_0$'cs1@99.5%'[i] == 1){
        GWAS_all_final_cs_0$cs1_group[i]  = "cs1@99.5%"
      } else {GWAS_all_final_cs_0$cs1_group[i]  = 0}
    }
    
    GWAS_all_final_cs_0 = GWAS_all_final_cs_0[GWAS_all_final_cs_0$p != 0, ] ##CHECK!!!
    GWAS_all_final_cs_0
  })
  
  output$contents_tt_gwas <- DT::renderDataTable({
    DT::datatable(GWAS_all_final_cs(),              
                  extensions = 'Buttons',
                  options = list(
                    paging = TRUE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    scrollX=TRUE,
                    ordering = TRUE,
                    dom = 'lftiprBRSPQ',
                    buttons = list( 
                      list(extend = 'csv',   filename =  "tt_PP"),
                      list(extend = 'excel', filename =  "tt_PP"))
                  ),
                  class = "display")
  })
  
  ##CHECK!!!-done
  # control and show credible sets in plots
  GWAS_select_t1_a_cs <- reactive({
    if (input$control_widgets_1 == 1){
      GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99.5%'==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99%'==1, ]
    }else if (input$control_widgets_1 == 3){
      GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@95%'==1, ]
    }else if (input$control_widgets_1 == 4){
      GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@90%'==1, ]
    }else if (input$control_widgets_1 == 5){
      GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99%'==0, ]
    }else{
      GWAS_all_final_cs()
    }
  })
  GWAS_select_t1_a <- reactive({
    GWAS_select_t1_a_cs()[GWAS_select_t1_a_cs()$mpp>=input$control_widgets_2[1] & GWAS_select_t1_a_cs()$mpp<=input$control_widgets_2[2], ]
  })
  # # 
  output$manhattan_1a <- renderPlotly({
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}

    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}

    fig_1a <- plot_ly(GWAS_select_t1_a(), x = ~position, y = ~p_log10,
                      type = 'scatter', mode = 'markers',
                      size = ~mpp, color = ~cs1_group,
                      colors = pal,
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP ID:', rsid,
                                    '<br>Allele 2:', allele2,
                                    '<br>Allele 1:', allele1,
                                    '<br>maf:', signif(maf,4),
                                    '<br>ps:', position,
                                    '<br>pval:', signif(p,4),
                                    '<br>-log10(pval):', signif(p_log10,4),
                                    '<br>MPP:', signif(mpp,4),
                                    '<br>Credible_set:', cs1_group))
    fig_1a <- fig_1a %>% layout(title = 'Interactive Regional Association Plots',
                                xaxis = list(title = paste0('Chromosome_', GWAS_select_t1_a()$chromosome[1], '_position'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t1_a()$position),
                                                       max(GWAS_select_t1_a()$position)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t1_a()$p_log10)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    fig_1a <- fig_1a %>% layout(shapes = list(vline(GWAS_select_t1_a()$position[which(GWAS_select_t1_a()$p_log10==max(GWAS_select_t1_a()$p_log10))]), hline(7.3)))
    fig_1a <- fig_1a %>% layout(annotations = list(
      list(x = 0.01 , y = 0.01, align = 'left',
           text = "NOTE: Size = MPP, Color = Credible_sets",
           showarrow = F, xref='paper', yref='paper')) )
    fig_1a <- fig_1a %>% layout(legend = list(orientation = 'h'))
    fig_1a
  })
  
  
  # pal_2 <- c("gray", "dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue")
  # pal_2 <- setNames(pal_2, c("0", "cs1@99.5%+@99%+@95%+@90%", "cs1@99.5%+@99%+@95%", "cs1@99.5%+@99%", "cs1@99.5%"))
  pal_defined <- reactive({
    if (input$control_widgets_5 == 1){
      pal_2 = c("gray", "red")
    }else if (input$control_widgets_5 == 2){
      pal_2 = c("gray", "pink")
    }else if (input$control_widgets_5 == 3){
      pal_2 = c("gray", "darkgreen")
    }else if (input$control_widgets_5 == 4){
      pal_2 = c("gray", "orange")
    }else if (input$control_widgets_5 == 5){
      pal_2 = c("gray", "black")
    }else{
      pal_2 = c("gray", "dodgerblue")
    }
    pal_2 <- setNames(pal_2, c("0", "user defined cs"))
  })
  
  GWAS_all_final_cs_2 <- reactive({
    cl = input$control_widgets_4[1]
    cs1_selected <- vector("list",length(cl))
    for(i in 1:length(cl)){
      tmp <- tt()$PP[,-1]
      names(tmp) <- tt()$PP[,1]
      cs1_selected[[i]] <- credset(tmp,cl[i])
    }
    GWAS_all_final_cs_0 = GWAS_all_final()
    GWAS_all_final_cs_0$p_log10 = -log10(GWAS_all_final_cs_0$p)
    for (i in 1:length(GWAS_all_final_cs_0$rsid)){
      if (GWAS_all_final_cs_0$rsid[i] %in% cs1_selected[[1]]){
        GWAS_all_final_cs_0$'user defined cs'[i] = 1
      } else {GWAS_all_final_cs_0$'user defined cs'[i] = 0}
    }
    
    for (i in 1:length(GWAS_all_final_cs_0$rsid)){
      if (GWAS_all_final_cs_0$'user defined cs'[i] == 1){
        GWAS_all_final_cs_0$cs1_group_2[i]  = "user defined cs"
      } else {GWAS_all_final_cs_0$cs1_group_2[i]  = 0}
    }
    GWAS_all_final_cs_0 = GWAS_all_final_cs_0[GWAS_all_final_cs_0$p != 0, ] ##CHECK!!!
    GWAS_all_final_cs_0
  })
  
  GWAS_select_t2_a_cs <- reactive({
    GWAS_all_final_cs_2()
  })
  GWAS_select_t2_a <- reactive({
    GWAS_select_t2_a_cs()[GWAS_select_t2_a_cs()$mpp>=input$control_widgets_6[1] & GWAS_select_t2_a_cs()$mpp<=input$control_widgets_6[2], ]
  })
  
  output$manhattan_2 <- renderPlotly({
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_2 <- plot_ly(GWAS_select_t2_a(), x = ~position, y = ~p_log10,
                      type = 'scatter', mode = 'markers',
                      size = ~mpp, color = ~cs1_group_2,
                      colors = pal_defined(),
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP ID:', rsid,
                                    '<br>Allele 2:', allele2,
                                    '<br>Allele 1:', allele1,
                                    '<br>maf:', signif(maf,4),
                                    '<br>ps:', position,
                                    '<br>pval:', signif(p,4),
                                    '<br>-log10(pval):', signif(p_log10,4),
                                    '<br>MPP:', signif(mpp,4),
                                    '<br>Credible_set:', cs1_group_2))
    fig_2 <- fig_2 %>% layout(title = 'Interactive Regional Association Plots',
                                xaxis = list(title = paste0('Chromosome_', GWAS_select_t1_a()$chromosome[1], '_position'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t2_a()$position),
                                                       max(GWAS_select_t2_a()$position)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t2_a()$p_log10)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    fig_2<- fig_2 %>% layout(shapes = list(vline(GWAS_select_t2_a()$position[which(GWAS_select_t2_a()$p_log10==max(GWAS_select_t2_a()$p_log10))]), hline(7.3)))
    fig_2 <- fig_2 %>% layout(annotations = list(
      list(x = 0.01 , y = 0.01, align = 'left',
           text = "NOTE: Size = MPP, Color = Credible_sets",
           showarrow = F, xref='paper', yref='paper')) )
    fig_2 <- fig_2 %>% layout(legend = list(orientation = 'h'))
    fig_2
  })

  output$Plot_LD1 <- renderPlotly({
    LD_selected_0 = as.matrix(tt()$LD)
    LD_selected = LD_selected_0[ ,GWAS_select_t1_a()$rsid]
    LD_selected = LD_selected[GWAS_select_t1_a()$rsid, ]
    fig_ld <- plot_ly(x = colnames(LD_selected), y = rownames(LD_selected), z = LD_selected, colors = "Greys", type = "heatmap")
    fig_ld
  })
  
  output$Plot_LD2 <- renderPlot({
    LD_selected_0 = as.matrix(tt()$LD)
    LD_selected = LD_selected_0[ ,GWAS_select_t1_a()$rsid]
    LD_selected = LD_selected[GWAS_select_t1_a()$rsid, ]
    if (input$control_widgets_1 == 0){
      hist(LD_selected)
    }else if (input$control_widgets_1 == 5){
      hist(LD_selected)
    }else{
      # LD_selected_0 = as.matrix(tt()$LD)
      # LD_selected = LD_selected_0[ ,GWAS_select_t1_a()$rsid]
      # LD_selected = LD_selected[GWAS_select_t1_a()$rsid, ]
      LD.plot(LD_selected, snp.positions = GWAS_select_t1_a()$position)
    }
  })
  
  #Venn diagram cs_1----
  fig_ch4_tab2_venn_cs1_data <- reactive({x_cs1_1 = c();
                                          x_cs1_2 = c();
                                          x_cs1_3 = c();
                                          x_cs1_4 = c();
  for (i in 1:length(GWAS_all_final_cs()$rsid)){
    if (GWAS_all_final_cs()$'cs1@99.5%'[i]==1){
      x_cs1_1 = c(x_cs1_1, GWAS_all_final_cs()$rsid[i])
    }
    if (GWAS_all_final_cs()$'cs1@99%'[i]==1){
      x_cs1_2 = c(x_cs1_2, GWAS_all_final_cs()$rsid[i])
    }
    if (GWAS_all_final_cs()$'cs1@95%'[i]==1){
      x_cs1_3 = c(x_cs1_3, GWAS_all_final_cs()$rsid[i])
    }
    if (GWAS_all_final_cs()$'cs1@90%'[i]==1){
      x_cs1_4 = c(x_cs1_4, GWAS_all_final_cs()$rsid[i])
    }
  }
  x_cs1_all <- list('cs1@99.5%' = x_cs1_1, 
                    'cs1@99%' = x_cs1_2,
                    'cs1@95%' = x_cs1_3,
                    'cs1@90%' = x_cs1_4)
  })
  
  output$fig_ch4_tab2_venn_cs1 <- renderPlotly({     
    ggVennDiagram(fig_ch4_tab2_venn_cs1_data(), show_intersect = TRUE)
  })
  
  output$table_ch4_tab2_venn_cs1_table <- DT::renderDataTable({
    check_x_cs1_table = process_region_data(Venn(fig_ch4_tab2_venn_cs1_data()))[,c(1,2,5,4,3)]
    DT::datatable(data.frame(check_x_cs1_table), 
                  #options = list(dom = 't')
                  extensions = 'Buttons',
                  
                  options = list(
                    paging = FALSE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    ordering = TRUE,
                    dom = 'tB',
                    buttons = list( 
                      list(extend = 'csv',   filename =  "cs1_table"),
                      list(extend = 'excel', filename =  "cs1_table"))
                  ),
                  class = "display"
    )
  })
  
  
  
  
  ##CHECK!!!
  pp_1_rownames_cs1 <- reactive({ 
    if (input$control_widgets_3b_cs == 1){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99.5%'==1, ]$rsid
    }else if (input$control_widgets_3b_cs == 2){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99%'==1, ]$rsid
    }else if (input$control_widgets_3b_cs == 3){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@95%'==1, ]$rsid
    }else if (input$control_widgets_3b_cs == 4){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@90%'==1, ]$rsid
    }else{
      cs1_nk = GWAS_all_final_cs()$rsid
    }
    pp_1_rownames = tt()$PP$str
    pp_1_rownames = str_split(pp_1_rownames, "%", simplify = TRUE)
    ##cs1_selected
    fun_check_cs <- function(x,cs){
      if (x %in% cs){
        return(x)
      }else{
        return("NA")
      }
    }
    pp_1_rownames_cs1 = apply(pp_1_rownames, 1:2, fun_check_cs, cs = cs1_nk)
    pp_1_rownames_cs1 = pp_1_rownames_cs1[rowSums(pp_1_rownames_cs1 == "NA")<dim(pp_1_rownames_cs1)[2]-1,]
    pp_1_rownames_cs1
  })
    
  pp_cs1_link <- reactive({     
    pp_1_rownames_cs1_link = c()
    for (i in 1:dim(pp_1_rownames_cs1())[1]){
      tmp = combn(pp_1_rownames_cs1()[i,][pp_1_rownames_cs1()[i,] != 'NA'], 2)
      tmp = as.data.frame(t(tmp))
      tmp$level = i
      tmp$value = tt()$PP[i,2]
      tmp$trait_linecolor = "gray"
      tmp$color = "chocolate"
      pp_1_rownames_cs1_link = rbind(pp_1_rownames_cs1_link, tmp)
    }
    colnames(pp_1_rownames_cs1_link)[1:2] = c('source', 'target')
    #pp_cs1_link = pp_1_rownames_cs1_link
    pp_1_rownames_cs1_link
  })
  
  pp_cs1_node <- reactive({     
    ## pp_nodes
    if (input$control_widgets_3b_cs == 1){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99.5%'==1, ]$rsid
    }else if (input$control_widgets_3b_cs == 2){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@99%'==1, ]$rsid
    }else if (input$control_widgets_3b_cs == 3){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@95%'==1, ]$rsid
    }else if (input$control_widgets_3b_cs == 4){
      cs1_nk = GWAS_all_final_cs()[GWAS_all_final_cs()$'cs1@90%'==1, ]$rsid
    }else{
      cs1_nk = GWAS_all_final_cs()$rsid
    }
    pp_cs1_node = cs1_nk
    pp_cs1_node = as.data.frame(unique(pp_cs1_node))
    colnames(pp_cs1_node) <- 'name'
    pp_cs1_node$id = c(0:(length(pp_cs1_node[,1])-1))
    pp_cs1_node$role_size_1 = 0
    for (i in 1:length(pp_cs1_node[,1])){
      pp_cs1_node$role_size_1[i] = sum(pp_cs1_node$name[i] == pp_1_rownames_cs1())+1
    }
    colnames(pp_cs1_node) <- c('name', 'id', 'role_size_total')
    pp_cs1_node$group_color = "blue"
    pp_cs1_node
  })
  
  plt_nodes_snp <- reactive({  
    plt_nodes_snp <- data.frame(
      name =pp_cs1_node()$name,
      n_color=pp_cs1_node()$group_color,
      n_size =pp_cs1_node()$role_size_total/10
    )
    plt_nodes_snp
  })
    
  plt_links_snp <- reactive({  
    plt_links_snp <- data.frame(
      source=pp_cs1_link()$source,
      target=pp_cs1_link()$target,
      l_color=pp_cs1_link()$color,
      l_size=pp_cs1_link()$value*10
    )
    for (i in 1:length(plt_links_snp$source)){
      plt_links_snp$source_id[i] = pp_cs1_node()$id[which(plt_links_snp$source[i]==pp_cs1_node()$name)]
      plt_links_snp$target_id[i] = pp_cs1_node()$id[which(plt_links_snp$target[i]==pp_cs1_node()$name)]
    }
    plt_links_snp  = plt_links_snp[plt_links_snp$l_size>=input$control_widgets_3b_finemap[1]*10 & plt_links_snp$l_size<=input$control_widgets_3b_finemap[2]*10, ]
    plt_links_snp 
  })
  
  #CHECK
  # plt_nodes_snp_udpated <- reactive({
  #   plt_links_snp_updated = unique(c(plt_links_snp()$source,plt_links_snp()$target))
  #   plt_nodes_snp_updated = plt_nodes_snp()[plt_nodes_snp()$name %in% plt_links_snp_updated,]
  #   plt_nodes_snp_updated
  # })

    # #
  network <- reactive({
                forceNetwork(Links = plt_links_snp(),
                             #Nodes = plt_nodes_snp_udpated(),
                             Nodes = plt_nodes_snp(),
                             Source = 'source_id',
                             Target = 'target_id',
                             Value = 'l_size',
                             NodeID = 'name',
                             Nodesize = 'n_size',
                             Group = 'n_color',
                             #height = 500,
                             #width = 800,
                             fontSize = 10,
                             fontFamily = "serif",
                             linkDistance = 150,
                             charge = -30, #numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value).
                             linkColour = plt_links_snp()$l_color,
                             opacity = 0.8,
                             zoom = T,
                             legend = F,
                             bounded = T,
                             opacityNoHover = 1)
  })
  
  output$ch1_network_snp_finemap <- renderForceNetwork({
    network()
  })
  
  output$d <- downloadHandler(
    filename = function() {
      'snp_network.html'
    },
    content = function(file) {
      saveNetwork(network(), file)
    }
  )
  # 
  
} #End of Shiny server


####----Shiny app----
shinyApp(ui, server)