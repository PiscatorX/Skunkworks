library(tidyverse)
library(plotly)
library(dplyr)
library(caret)
library(feather)

reticulate::virtualenv_create("py3k", python = "python3")

reticulate::virtualenv_install(c('pandas',
                                 'scikit-learn==0.24.0'),
                       envname = "py3k",
                       ignore_installed = TRUE)

reticulate::use_virtualenv("py3k")
reticulate::source_python("reformat_zeo.py")
reticulate::source_python("random_forest.py")

ui <- fluidPage(
    titlePanel("Zeolite RF prediction engine"
      
    ),
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          h3("Adsorbent properties"),
          
          column(4, selectInput("Adsorbent",
                                "Adsorbent:",
                                choices = c('CuAgY', 'CuCeY', 'NiCeY', 'AgY', 'clinoptilolite', 'CeY', 'CuY', 'CuX', 'NiY', 'CsY', 'NaY', 'MCM-22', 'AgCeY', 'AgX', 'CuHY'),
                                multiple = FALSE)),
          column(4,numericInput("SA",
                                "Surface area:",
                                min = 141.4,
                                max = 720.0,
                                value = 0)),
          column(4,numericInput("Vmicro",
                                "Micropore volume:",
                                min = 0.08,
                                max = 0.666,
                                value = 0)),
        ),
        fluidRow(
          column(4,numericInput("Vmeso",
                                "Mesopore volume:",
                                min = 0.04,
                                max = 0.18,
                                value = 0)),
          
          column(4,  numericInput("pore_size",
                                  "Pore size:",
                                  min = 0.83,
                                  max = 3.1729,
                                  value = 0)),
          
          column(4, numericInput("Si_Al",
                                 "Si/Al ratio:",
                                 min = 1.031,
                                 max = 46.0,
                                 value = 0)),
        ),
        h2("Metal properties"),
        fluidRow(
          column(4,   numericInput("Na.",
                                   "Na+:",
                                   min = 0.02,
                                   max = 1.042538354,
                                   value = 0),
          ),
          column(4,   numericInput("Ag.",
                                   "Ag+:",
                                   min = 0.297,
                                   max = 1.18845,
                                   value = 0)
          ),
          column(4,   numericInput("Cu.",
                                   "Cu+:",
                                   min = 0.169,
                                   max = 0.467,
                                   value = 0)
          )),
          fluidRow(
          column(4,  numericInput("Ce.4",
                                  "Ce+4:",
                                  min = 0.21,
                                  max = 2.29489,
                                  value = 0)),
          column(4, numericInput("Cs.2",
                                 "Cs+2:",
                                 min = 0.666666667,
                                 max = 0.666666667,
                                 value = 0)),
          column(4, numericInput("Ni.2",
                                 "Ni+2:",
                                 min = 0.024232082,
                                 max = 0.3125,
                                 value = 0)),
        ),
        h3("Adsorbate properties"),
        fluidRow(
        column(4, selectInput("adsorbate", 
                             "Adsorbate:", 
                             choices = c('TP', 'BT', 'DBT'),
                             multiple = FALSE)
               ),
        column(4, numericInput("dipole_moment",
                              "Dipole moment:",
                              min = 0.51,
                              max = 0.79,
                              value = 0)
        ),
        column(4,   numericInput("chemical_hardness",
                                 "Chemical hardness:",
                                 min = 3.0401,
                                 max = 5.602,
                                 value = 0))
               
        ),
        fluidRow(
          column(4, numericInput("kinetic_diameter",
                                 "Kinetic diameter:",
                                 min = 0.77,
                                 max = 0.91,
                                 value = 0)
                 )
        ),
        h3("Adsorption conditions"),
        fluidRow(
         column(4, numericInput("C_0",
                                "Initial concentration:",
                                min = 8.0,
                                max = 3405.0,
                                value = 0)),

         column(4, selectInput("solvent",
                               "Solvent:",
                               choices = c('cyclohexane', '1-octane', 'n-octane', 'n-heptane', 'iso-octane', 'n-Heptane', 'n-Octane', 'ether', 'hexadecane'),
                               multiple = FALSE)),

         column(4, numericInput("oil_adsorbent_ratio",
                                "Oil adsorbent ratio:",
                                min = 14.0,
                                max = 260.0,
                                value = 0)),


        ),
        fluidRow(
          column(4, numericInput("Temp",
                                 "Temperature:",
                                 min = 20.0,
                                 max = 80.0,
                                 value = 0)),
        ),
        fluidRow(
          column(8, fileInput(
                    "zeolite_file",
                    "Zeolite file input",
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
          )),
          column(8, downloadButton("downloadExample", label = "Example data")),
       
),

                  


          
                 ),
      
      mainPanel(
        
        tabsetPanel(type = "tabs",
                    tabPanel("Prediction results",
                            h3("Random Forests prediction results"),
                            
                            fluidRow(
                              column(3, actionButton("RunRandomForest","Predict adsoption capacity")),
                              column(3, downloadButton("downloadOutput", label = "Download")),
                            ),
                             fluidRow(
                               column(12, plotlyOutput('RF_gauge')),
                             ),
                            fluidRow(
                              column(6, dataTableOutput('zeo_ml_out')),
                            ),
                            ),
                    tabPanel("Plots", 
                               h3("Zeolite reference database PCA"),
                               fluidRow(
                                 column(3, selectInput("discrete",
                                           "Colour by:",
                                           choices = c("Adsorbent", "adsorbate", "solvent"),
                                           multiple = FALSE)),
                                  column(3, textInput("Label", 
                                                       "Label:", 
                                                       "MyZeolite"))
                                 
                                 ),
                             fluidRow(
                               column(3, actionButton("LoadZeo1","Load zeolite PCA plot")),
                               column(3, actionButton("LoadZeo2","Load zeolite PCA 3d plot"))),
                             plotlyOutput('PCAplot1'),
                             plotlyOutput('PCAplot2')
                             
                             ),
                    tabPanel("Data Entered",
                             dataTableOutput('entry'),
                    ),
                    tabPanel("About",
                             br(),
                             h5('Written by Andrew Ndhlovu and Liberty L. Mguni'),
                    )
        )
        
      

             
      )
  )
)
  
server <- function(input, output) {
  

  zeo_ref <- read_feather("zeolites_DBX.feather") 
  
  
  zeolite_encoded <- read_feather("ZeoX_Final_encoded_V2x.feather") 
  
  
  prop_cols <- c("x_Na.", "x_Ag.", "x_Cu.", "x_Ce.4", "x_Cs.2", 
                 "x_Ni.2", "R_Na.", "R_Ag.", "R_Cu.", "R_Ce.4", "R_Cs.2", "R_Ni.2")
  
  metal_properties <- cbind(zeo_ref["Adsorbent"], zeo_ref[prop_cols]) %>%
                      distinct()
  
  get_prop <- function(adsorbent){
    
    return  (metal_properties[metal_properties$Adsorbent == adsorbent,] %>%
               mutate(across(everything(), ~replace_na(.x, 0))) %>%
               data.frame())
  }
  
  
  get_zeo <- reactive({entry <- data.frame(
                                   Adsorbent=input$Adsorbent,
                                   SA=input$SA,
                                   Vmicro=input$Vmicro,
                                   Vmeso=input$Vmeso,
                                   pore_size=input$pore_size,
                                   Si_Al=input$Si_Al,
                                   `Na+`=input$`Na.`,
                                   `Ag+`=input$`Ag.`,
                                   `Cu+`=input$`Cu.`,
                                   `Ce+4`=input$`Ce.4`,
                                   `Cs+2`=input$`Cs.2`,
                                   `Ni+2`=input$`Ni.2`,
                                   adsorbate=input$adsorbate,
                                   dipole_moment=input$dipole_moment,
                                   chemical_hardness=input$chemical_hardness,
                                   kinetic_diameter=input$kinetic_diameter,
                                   C_0=input$C_0,
                                   solvent=input$solvent,
                                   oil_adsorbent_ratio=input$oil_adsorbent_ratio,
                                   Temp=input$Temp)
  
                          return(entry)
  })

    
    get_prop <- function(adsorbent){
      
      return  (metal_properties[metal_properties$Adsorbent == adsorbent,] %>%
                 mutate(across(everything(), ~replace_na(.x, 0))) %>%
                 data.frame())
      
    }
    
    
    
    numeric_cols  <- zeo_ref %>%
                     select_if(is.numeric) %>%
                     colnames()
    
    get_entryx <- reactive({
      
      data_in <- get_zeo()
      metal_prop <- get_prop(input$Adsorbent)
      
      return(cbind(data_in[c(1:12)], metal_prop ,data_in[c(13:20)]))
    })
    

    
    
    get_data <- reactive({
          
          zeo_ml <- list() 
      
          inFile <- input$zeolite_file
          
          if (is.null(inFile)){
            zeo_ml$zeo_in <- get_entryx()
            
          }
          else{ 
            data_in <- read.table(inFile$datapath,
                                    sep = "\t",
                                    header = T,
                                  stringsAsFactors = F)
             zeo_ml$zeo_in <- data_in %>% select(-c(Label))
             zeo_ml$Label <- data_in$Label   
          
          }
          
          zeo <- FormatZeo(zeo_ml$zeo_in)
          zeo_ml$entryx <- zeo$get_entry()
          
          stopifnot(ncol(zeolite_encoded) == ncol(zeo_ml$entryx))
          colnames(zeo_ml$entryx) <- colnames(zeolite_encoded)
          
          zeo <- ZeoRandomForest(zeo_ml$entryx)
          zeo$RF_predict()                    
          zeo_ml$y_pred <- zeo$RF_predict()
          
       
       return(zeo_ml)
    })
    

    
observeEvent(input$RunRandomForest,{
  
    zeom_ml <- get_data()
    
    if (is.null(input$zeolite_file)){ 
      zeom_ml$Label <- input$Label
    }
    
    model_out <- data.frame(Label=zeom_ml$Label, Capacity = zeom_ml$y_pred)
    
    
    output$downloadOutput <- downloadHandler(
      filename = function() {
         return("zeolite_ml.tsv")
      },
      content = function(file) {
        write.table(model_out, file, sep = "\t", row.names = FALSE, quote = F)
      }
    )
  
    
    
    
    
    output$zeo_ml_out <- renderDataTable({ model_out })
    
    output$RF_gauge <- renderPlotly({
      
    CAPACITY_MAX =  60.8
    CAPACITY_MIN =  0.5
    ypred_mean = mean(zeom_ml$y_pred)
    
    gauge_min <- if(ypred_mean < CAPACITY_MIN) ypred_mean else CAPACITY_MIN
    gauge_max <- if(ypred_mean > CAPACITY_MAX) ypred_mean else CAPACITY_MAX
    
    plot_ly(
        domain = list(x = c(0, 1), y = c(0, 1)),
        value = ypred_mean,
        title = list(text = "Predicted Capacity (mg S per g)"),
        type = "indicator",
        mode = "gauge+number", width = 500,
        height = 500)  %>%
        layout(margin = list(l=gauge_min, r= gauge_max))

})
   
})   
    
    
    

    pca_data <- reactive({

            zeo_ml <- get_data()
            PCA <- list()
            zeolite_encodedx <- rbind(zeolite_encoded, zeo_ml$entryx)
            
            Scaler <- preProcess(zeolite_encodedx, method = list(center = numeric_cols, scale = numeric_cols))
            zeolite <- predict(Scaler, zeolite_encodedx)
            zeoDB_pca <- prcomp(zeolite)
            PCA$data <- data.frame(zeoDB_pca$x, stringsAsFactors = F)
            
            color_by <- rbind(zeo_ref[input$discrete], zeo_ml$zeo_in[input$discrete])
            PCA$data <- cbind(PCA$data, color_by)
        
            
            PCA$n_i <- nrow(zeolite_encoded)
            PCA$n_j <- nrow(PCA$data)
            PCA$dada$Label <- character()
            
            PCA$data$Label[c((PCA$n_i+1):PCA$n_j)] <- if (is.null(input$zeolite_file)) input$Label else zeo_ml$Label 
            
            return(PCA)

    })

        

    
       
observeEvent(input$LoadZeo1,{
     
output$PCAplot1 <- renderPlotly({

            PCA <- pca_data()
            
            ggplotly(ggplot(data = PCA$data[c(1:PCA$n_i),], aes_string(x = "PC1" , y = "PC2", color = input$discrete )) +
                    xlab(label = "PC1") +
                    ylab(label = "PC2") +
                    geom_point(size = 1) +
                    geom_point(data= filter(PCA$data, row_number() > PCA$n_i), aes_string(x = "PC1" , y = "PC2", color = input$discrete), size = 1.75, shape = 23) +
                    geom_text(data= filter(PCA$data, row_number() > PCA$n_i), aes(label = Label, colour = "black"), nudge_y = 0.25) +
                    theme(panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          legend.key = element_rect(fill = "white"),
                          legend.position="top",
                          axis.text = element_text(size=18, colour = "black"),
                          axis.title = element_text(size=18, colour = "black"),
                          panel.border = element_rect(colour = "black", fill=NA)), 
                    width = 750,
                    height = 750) %>%
                            layout(legend = list(orientation = "h", x = 0.4, y = -0.2))})

})    


observeEvent(input$LoadZeo2,{
  
  color_tag <- unlist(pca_data()[input$discrete])
  
  output$PCAplot1 <-  renderPlotly(plot_ly(data = pca_data()$data, 
                                            x = ~PC1, 
                                            y = ~PC2, 
                                            z = ~PC3, 
                                            color = ~unlist(pca_data()$data[input$discrete]), width = 750,
                                           height = 750)  %>%
                                      add_markers(size = 5))     
                      
})




get_prop <- function(adsorbent){
  
  return  (metal_properties[metal_properties$Adsorbent == adsorbent,] %>%
             mutate(across(everything(), ~replace_na(.x, 0))) %>%
             data.frame() %>%
             select(-c(Adsorbent)))
  
}


output$entry <- renderDataTable({

      inFile <- input$zeolite_file
      
      if (is.null(inFile)){
        
        return(get_zeo())
        
      }
      else{
        
        return(read.table(inFile$datapath,
                 sep = "\t",
                 header = T))
      }


})



output$downloadExample <- downloadHandler(
  filename = function() {
    return("zeo_sample.tsv")
  },
  content = function(file) {
    file.copy("zeo_sample.tsv", file)
  }
)






}    
shinyApp(ui, server)