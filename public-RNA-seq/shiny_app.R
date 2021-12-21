library(shiny)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(tidyr)
library(dplyr)
library(shinythemes)
library(plotly)
library(DESeq2)
library(heatmaply)
library(shinyWidgets)
library(RColorBrewer)

# load data

# loading meta_data
meta_all <- read.csv("data/meta_all.csv")
meta_rep <- read.csv("data/meta_rep.csv")

# counts and abundance (TPM)
txi.all <- readRDS("data/txi.all.RData")
tpm <- as.data.frame(txi.all$abundance)

# Variance stabalised log2 (VST) 
vst.all <- readRDS("data/vst.all.RData")
vst <- as.data.frame(assay(vst.all))

# z-score of TPM
tpm_z <- as.data.frame(t(scale(t(tpm), center = T, scale = T)))

# circadian expression

circ <- read.csv("data/LM_circadian.csv", header = TRUE, check.names = FALSE, row.names = 1)

# scRNA-seq

scRNA_mo <- readRDS("data/scRNA_monocle")


# merging/collapsing biological replicates 
# TPM
tpm_m <- tpm
colnames(tpm_m) <- meta_all$condition.study
tpm_m <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(tpm_m)), # for each unique column name
         function(col) rowMeans(tpm_m[names(tpm_m) == col]) # calculate row means
  ))

# VST
vst_m <- vst
colnames(vst_m) <- meta_all$condition.study
colnames(vst) <- meta_all$condition.study
vst_m <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(vst_m)), # for each unique column name
         function(col) rowMeans(vst_m[names(vst_m) == col]) # calculate row means
  ))

# z-score TPM
tpm_z_m <- tpm_z
colnames(tpm_z_m) <- meta_all$condition.study
tpm_z_m <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(tpm_z_m)), # for each unique column name
         function(col) rowMeans(tpm_z_m[names(tpm_z_m) == col]) # calculate row means
  ))


# Define UI 
ui <- navbarPage("Marpo-seq", theme = shinytheme("flatly"),
                 
  tabPanel("Sample overview", 
           
           p("Marpo-seq is a interactive data visualisation tool for transcriptomic data of ", em("Marchantia polymoprha"), "."),
              
           h3("Meta data and sample overview"),
           
           fluidPage(
             
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(checkboxGroupInput("studies", 
                                               label="Studies", inline=T, 
                                               choices = as.character(unique(meta_all$study_accession)), 
                                               selected = as.character(unique(meta_all$study_accession))), 
                            selectInput("PCA_colour", 
                                        label = "Colour by",
                                        choices = c("Study", "Tissue","Genotype", "Age"),
                                        selected = "Study"), width = 2
                            ),  
               
               mainPanel(plotlyOutput(outputId = "PCA", height = "800px")) 
               )
             ), 
                 
           
           
           p("Data was collected from:"),
           p(a("PRJDB4420", href="https://doi.org/10.1093/pcp/pcw005", target="_blank" ), "- Higo A et al. 2016, Plant Cell Physiol: Reproductive tissues"), 
           p(a("PRJDB5890", href="https://doi.org/10.1007/s10265-018-1044-7", target="_blank" ), "- Kubo H et al. 2018, J Plant Res: Regulation of anthocyanins (R2R3-MYBs)"),                
           p(a("PRJDB6579", href="https://doi.org/10.1016/j.cub.2017.12.053", target="_blank" ), "- Yamaoka S et al. 2018, Curr Biol: Gametangium development (MpBNB)"),
           p(a("PRJDB7023", href="https://doi.org/10.1016/j.cell.2017.09.030", target="_blank" ), "- Bowman JL et al. 2017, Cell: Various tissues"),
           p(a("PRJNA218052", href="https://doi.org/10.1186/1471-2164-14-915", target="_blank" ), "- Sharma N et al. 2013, BMC Genomics: Reproductive tissues"),
           p(a("PRJNA251267", href="https://doi.org/10.1016/j.cell.2017.09.030", target="_blank" ), "- Bowman JL et al. 2017, Cell: Various tissues"),
           p(a("PRJNA265205", href="https://doi.org/10.1093/molbev/msu303", target="_blank" ), "- Frank MH et al. 2015, Mol Biol Evol: Sporophyte/apical cell"),
           p(a("PRJNA350270", href="https://doi.org/10.1016/j.cell.2017.09.030", target="_blank" ), "- Bowman JL et al. 2017, Cell: Spore germination timecourse"),
           p(a("PRJNA397394", href="https://doi.org/10.7554/eLife.33399", target="_blank" ), "- Mutte SK et al. 2018, eLife: Auxin treatment"),              
           p(a("PRJNA433456", href="https://doi.org/10.1111/nph.15090", target="_blank" ), "- Flores-Sandoval E et al. 2018, New Phytologist: MpARF3, MpMIR160 knockout"), 
           p(a("PRJDB5847", href="https://doi.org/10.1104/pp.18.00761", target="_blank" ), "- Jahan A et al. 2019, Plant Physiology: ABA response"), 
           p(a("PRJDB8103", href="https://doi.org/10.1093/pcp/pcy161", target="_blank" ), "- Ikeda Y et al. 2018, Plant and Cell Physiology: Loss of methylation (mpmet)"), 
           p(a("PRJNA488718", href="https://doi.org/10.1007/s10265-019-01095-w", target="_blank" ), "- Arai H et al. 2019, J Plant Res: MpBHLH12 mutants"), 
           
           p("This tool was developed by Marius Rebmann (Haseloff lab, Department of Plant Sciences, Universty of Cambridge), for feedback, questions, bugs and suggestions please contact: mr728@cam.ac.uk")
           
  ), 
  

  tabPanel("Search Genes",
  fluidPage(
  
    selectizeInput("gene",label="Select genes", 
                   choices = rownames(tpm), 
                   selected="Mapoly0024s0116", multiple=T,
                   options = list(maxOptions = 5, create=T)), 
    

 
    
    # Sidebar panel for inputs ----
      
      # Input: Slider for the number of bins ----
  
      
      
      

    # Main panel for displaying outputs ----
      
      tabsetPanel(
        tabPanel("Bulk Expression", fluid=TRUE, 
                 sidebarLayout(
                 sidebarPanel(
                   selectInput("normalisation", 
                               label = "Select normalisation",
                               choices = c("VST", "TPM", "Z-score (TPM)"),
                               selected = "VST"),
                   
                   checkboxGroupInput("studies_g", 
                                      label="Studies", inline=T, 
                                      choices = as.character(unique(meta_all$study_accession)), 
                                      selected = as.character(unique(meta_all$study_accession))),
                   
                   sliderInput("rep", "Minimum replicates",
                               min = 1, max = 5, value = 2),
                   
                   selectInput("sort1", 
                               label = "Sort first",
                               choices = c("Study", "Tissue","Genotype", "Age"),
                               selected = "Study"),  
                   
                   selectInput("sort2", 
                               label = "Sort second",
                               choices = c("Study", "Tissue","Genotype", "Age"),
                               selected = "Age"), 
                   
                   selectInput("cluster", 
                               label = "Clustering",
                               choices = c("No", "Samples","Genes", "Both"),
                               selected = "No"), 
                   width=2
                 
              ), 
                 
                 mainPanel(
                   plotlyOutput(outputId = "heatmap_gene", height = "800px"), 
                   plotlyOutput(outputId = "rel_expr", height = "400px") )
                 
                 
                 
                   ) 
                ),
        
        tabPanel("Circadian",  fluid=TRUE, 
                 sidebarLayout(
                   sidebarPanel(
                     
                selectInput("circ_data", 
                                 label = "Select dataset",
                                 choices = c("Light-Light", "Light-Dark", "All"),
                                 selected = "Light-Light"),
                
                 shinyWidgets::sliderTextInput("fdr","FDR",
                                               choices=c(0.0001, 0.001, 0.01, 0.05, 0.1),
                                               selected=0.01, grid = T), 
                 width=2
                 
                   ), 
                 mainPanel(plotlyOutput(outputId = "circadian_gene", height = "800px"), 
                           plotlyOutput(outputId = "circadian_overview", height = "400px") 
                           )
                 )
                 
                 ), 
        tabPanel("scRNA-seq",  fluid=TRUE, 
                 sidebarLayout(
                   sidebarPanel(
                     
                     width=2
                     
                   ), 
                   mainPanel(plotlyOutput(outputId = "UMAP", height = "800px"))
                   )
                 )
                 
        )
                 
        
      )
        
      
    
  ),
  tabPanel("Differential expression", 
           
           fluidPage(
                     
                     # App title ----
                     #titlePanel("M. polymorpha transcriptomics"),
                     
                     # Sidebar layout with input and output definitions ----
                     sidebarLayout(
                       
                       # Sidebar panel for inputs ----
                       sidebarPanel(
           
           radioButtons("studies_DE", 
                              label="Select study", inline=F, 
                              choices = c("PRJDB4420" = "PRJDB4420",
                                          "PRJDB5890"="PRJDB5890", 
                                          "PRJDB6579"="PRJDB6579",
                                          "PRJDB7023"="PRJDB7023",
                                          "PRJNA265205"="PRJNA265205", 
                                          "PRJNA350270"="PRJNA350270", 
                                          "PRJNA397394"="PRJNA397394", 
                                          "PRJNA433456"="PRJNA433456", 
                                          "PRJDB5847"="PRJDB5847", 
                                          "PRJDB8103"="PRJDB8103", 
                                          "PRJNA488718"="PRJNA488718"), 
                              selected = c("PRJDB4420")), 
           
           
            uiOutput("selectCondA"), 
           
            uiOutput("selectCondB"),
           
           shinyWidgets::sliderTextInput("padj","adjusted p-value:",
                                         choices=c(0.0001, 0.001, 0.01, 0.05, 0.1),
                                         selected=0.01, grid = T),
           
            #sliderInput("padj", "adjusted p-value",
               #        min = 0.0001, max = 0.1, value = 0.01),
           
           width = 2
           
           ),
           
           mainPanel( 
             tabsetPanel(
               tabPanel("Vulcan", plotlyOutput(outputId = "vulcan", height = "800px")),
               tabPanel("MA_plot", plotlyOutput(outputId = "MA_plot", height = "800px")),
               tabPanel("TopDE_genes", plotlyOutput(outputId = "topHM", height = "800px"), 
                        fluidRow(
                          column(6,
                         selectInput("cluster2", 
                                    label = "Clustering",
                                    choices = c("No", "Genes"),
                                    selected = "Genes")), 
                         column(6,
                        sliderInput("n_DE_genes", "nGenes to display",
                                    min = 5, max = 100, value = 20))))
          
          )
        )
      )
    )
  )
)


# Define server logic 
server <- function(input, output) {
  
  # PCA plot
  
  output$PCA <- renderPlotly({
    
    sample_sub <- filter(meta_all, study_accession %in% input$studies)$run_accession
    
    vst.sub <- vst.all[,sample_sub]
    
    rv <- rowVars(assay(vst.sub))
    select <- order(rv, decreasing = TRUE)[seq_len(min(1000, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(vst.sub)[select, ]))
    percentVar <- round(100*(pca$sdev^2/sum(pca$sdev^2)), digits = 1)
    
    pcaData <- inner_join(data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[,3], 
                                     run_accession = colnames(vst.sub)), 
                          filter(meta_all, study_accession %in% input$studies)
                          ,by = "run_accession")
    
    PCA_colour <- switch(input$PCA_colour, 
                         "Study"=c("study_accession"),
                         "Tissue"=c("tissue"),
                         "Genotype"=c("genotype"), 
                         "Age"=c("age"))
    
    plot_ly(pcaData, type = "scatter3d", x = ~PC1, y= ~PC2, z =~PC3,
            color=~get(PCA_colour), mode = "markers", marker = list(size=5),
            text = ~paste("Study", study_accession,
                          "</br> Sample: ", condition, 
                          "</br> Replicate: ", rep)) %>% layout(
                            scene = list(
                              xaxis = list(title = paste0("PC1: ",percentVar[1],"% variance")),
                              yaxis = list(title = paste0("PC2: ",percentVar[2],"% variance")),
                              zaxis = list(title = paste0("PC3: ",percentVar[3],"% variance"))
                            ))
  })
  
  
  output$selectCondA <- renderUI({
    radioButtons("cond_a", 
                 label="Condition A", inline=F, 
                 choices =  filter(meta_rep, study_accession %in% input$studies_DE)$condition, 
                 selected = filter(meta_rep, study_accession %in% input$studies_DE)$condition[1])
  })
  
  output$selectCondB <- renderUI({
    radioButtons("cond_b", 
                 label="Condition B", inline=F, 
                 choices =  filter(meta_rep, study_accession %in% input$studies_DE)$condition, 
                 selected = filter(meta_rep, study_accession %in% input$studies_DE)$condition[2])
  })
  
  
  
  # Bulk_expression heatmap of gene specified by user
  
  output$heatmap_gene <- renderPlotly({
    
    df_bulk <- switch(input$normalisation, 
                      "VST"=vst_m, 
                      "TPM" = tpm_m, 
                      "Z-score (TPM)"=tpm_z_m)
    
    
    # filter by number of replicates
    # these expressions need to be reactive or direktly depend on a input
    
    # sorting
    sort_1 <- switch(input$sort1, 
                     "Study"=c("study_accession"),
                     "Tissue"=c("tissue"),
                     "Genotype"=c("genotype"), 
                     "Age"=c("age"))
    
    sort_2 <- switch(input$sort2, 
                     "Study"=c("study_accession"),
                     "Tissue"=c("tissue"),
                     "Genotype"=c("genotype"), 
                     "Age"=c("age"))
    
    
    df_h <- reactive(df_bulk[,match(dplyr::arrange_(filter(meta_rep, replicates >= input$rep & 
                                                      study_accession %in% input$studies_g), sort_1, sort_2)$condition.study, colnames(df_bulk), nomatch=0)])
    
    
    # get study name for heatmap annotation
    
    annot_study <- as.data.frame(arrange_(filter(meta_rep, replicates >= input$rep & study_accession %in% input$studies_g), sort_1, sort_2)$study_accession)
    colnames(annot_study) <- c("study")
    rownames(annot_study) <- colnames(df_h())

    
    clustering <- switch(input$cluster, 
                        "No"="none",
                        "Samples"="row",
                        "Genes"="column", 
                        "Both"="both")
     
     heatmaply(t(df_h()[input$gene,]),
               dendrogram = clustering, row_side_colors = annot_study)
    
  })
  
  # relative expression rank 
  
  output$rel_expr <- renderPlotly({
    
    rel_expr <- as.data.frame(rowSums(tpm_m)/ncol(tpm_m))
    
    colnames(rel_expr) <- c("TPM")
    
    rel_expr$Rank <- rank(-rel_expr$TPM)
    
    rel_expr$TPM <- ifelse(rel_expr$TPM < 1e-3, 1e-3, rel_expr$TPM) # cap minimum expression for visualisation
  
    goi_expr <- as.data.frame(melt(as.matrix(t(tpm_m[input$gene,]))))
    
    colnames(goi_expr) <- c("Sample", "Mapoly_ID", "TPM")
    
    goi_rel_expr <- rel_expr[input$gene,]
    
    goi_rel_expr$Mapoly_ID <- rownames(goi_rel_expr)
    
    goi_expr <- left_join(goi_expr, goi_rel_expr[2:3], by="Mapoly_ID")
    
    
    plot_ly() %>%
      add_trace(data = rel_expr, x=~Rank, y=~TPM, marker = list(color='black'), 
                hovertemplate = paste("<b>",rownames(rel_expr),"</b>", 
                                      "<br>%{yaxis.title.text}: %{y}",
                                      "<br>%{xaxis.title.text}: %{x}"), showlegend=FALSE) %>%
      add_trace(data = goi_expr, x=~Rank, y=~TPM, color=~Mapoly_ID, hovertemplate = paste("<b>",goi_expr$Mapoly_ID,"</b>",
                                                                        "<br>",goi_expr$Sample,
                                                                        "<br>%{yaxis.title.text}: %{y}", 
                                                                        "<br>%{xaxis.title.text}: %{x}")) %>%  
      layout(yaxis = list(type = "log"))
  
  
    })
  
  # circadian expression
  
  circ_sig <- reactive({
    
    circ[circ$fdr_BH_LL_4884p222834 < input$fdr, ]
    
  })
  
  datapoints <- reactive({
    
    switch(input$circ_data, 
           "Light-Light"=7:19,
           "Light-Dark"=1:6,
           "All"=1:19)
  })
  
  circ_sig_goi <- reactive({
  
  melt(t(circ_sig()[datapoints()][input$gene,]/(rowSums(circ_sig()[datapoints()][input$gene,])/length(circ_sig()[datapoints()][input$gene,])))) # melt dataframe to enable plotting with ggplot
  
  })
  
  
  
  output$circadian_overview <- renderPlotly({
    
    ggplotly(ggplot(circ_sig(), aes(x=phase_LL, y=log10(amplitude_LL))) + geom_point(alpha=0.2) + 
               geom_rug(alpha=0.2) + geom_point(data = circ_sig()[input$gene,], aes(x=phase_LL, y=log10(amplitude_LL)), color=brewer.pal(dim(circ_sig()[input$gene,])[1], "Set2")[1:dim(circ_sig()[input$gene,])[1]]) + 
               scale_x_continuous(breaks = seq(0, 32, 4), limits = c(0,32)) + theme_light() )
  })
  
  output$circadian_gene <- renderPlotly({
    
    circ_sig_goi <- circ_sig_goi()
    
    colnames(circ_sig_goi) <- c("Timepoint", "Mapoly_ID", "normalised_expression")
    
    plot_ly(data = circ_sig_goi, x=~Timepoint, y=~normalised_expression, color=~Mapoly_ID, mode='lines')
    
  })
  
  
  # scRNA-seq
  
  output$UMAP <- renderPlotly({
    
    d1 <- data.frame(pData(scRNA_mo), t(scRNA_mo@reducedDimS))
    d1$goi <- as.numeric(exprs(scRNA_mo[paste(input$gene, ".v3.1", sep = ""),]))
    
    plot_ly(data = d1, x=~X1, y=~X2, color=~log10(goi+1))
  
  })
  # DE-genes
  
  de.sub <- reactive({
                      DESeq(DESeqDataSetFromTximport(lapply(txi.all, 
                            function(x) if(is.matrix(x)) return(x[,as.character(filter(meta_all, 
                                                                  study_accession %in% input$studies_DE)$X)]) else return(x)), 
                                           filter(meta_all, study_accession %in% input$studies_DE), ~ condition))
                    })
  res <- reactive({results(de.sub(), contrast = c("condition",input$cond_a, 
                                        input$cond_b))
              })
  
  output$vulcan <- renderPlotly({
    
    p <- ggplot(as.data.frame(res()),
           aes(log2FoldChange, -log10(padj), color=padj<input$padj)) + 
      geom_point(aes(text=rownames(res())), size=0.5) + scale_color_manual(values=c("grey50", "red3")) +
      annotate("text", x=0,y=-10, label = paste("DE genes at padj = ", as.character(input$padj), " : ",
                                                length(which(as.data.frame(res())$padj<input$padj)))) +
      theme_light() + ggtitle(paste(input$cond_a, 
                                    input$cond_b, sep="-vs-")) 
    
    ggplotly(p)
    
  })
  
  output$MA_plot <- renderPlotly({
    
    p <- ggplot(as.data.frame(res()), 
           aes(baseMean,log2FoldChange, color=padj<input$padj)) + 
      geom_point(aes(text=rownames(res())), size=0.5) + scale_x_log10() + scale_color_manual(values=c("grey50", "red3")) + 
      theme_light() + ggtitle(paste(input$cond_a,input$cond_b, sep="-vs-"))
    
    ggplotly(p)
    
  })
  
  
  output$topHM <- renderPlotly({
    
    
  clustering2 <- switch(input$cluster2, 
                         "No"="none",
                         "Genes"="row")

  top_genes <- head(arrange(tibble::rownames_to_column(as.data.frame(res()),var="Phytozome.ID"), padj),input$n_DE_genes)$Phytozome.ID

  #pheatmap(vst[top_genes,paste(c(input$cond_a,input$cond_b),input$studies_DE, sep = ".")], cluster_rows=T, show_rownames=T,
  #         cluster_cols=F)
  
  heatmaply(vst[top_genes,paste(c(input$cond_a,input$cond_b),input$studies_DE, sep = ".")], dendrogram = clustering2)
  
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
