library(tidyverse)
library(plotly)
library(dplyr)
library(caret)
library(feather)

reticulate::virtualenv_create("py3k")
reticulate::virtualenv_install("numpy", "scipy", "sklearn", "pandas",  ignore_installed = T)
reticulate::use_virtualenv("py3k")
reticulate::py_install("pandas", pip = T)
reticulate::source_python("reformat_zeo.py")

ui <- fluidPage(
    titlePanel("Zeolie RF prediction engine"
      # app title/description
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



          
                 ),
      
      mainPanel(
        
        tabsetPanel(type = "tabs",
                    tabPanel("Plots", 
                             fluidRow(
                               h3("Zeolite reference database PCA"),
                               selectInput("discrete",
                                           "Colour categorical variable:",
                                           choices = c("Adsorbent", "adsorbate", "solvent"),
                                           multiple = FALSE)),
                             fluidRow(
                               column(4, actionButton("LoadZeo1","Load zeolite PCA plot")),
                               column(4, actionButton("LoadZeo2","Load zeolite PCA 3d plot"))),
                             plotlyOutput('PCAplot1'),
                             plotlyOutput('PCAplot2')
                             
                             ),
                    tabPanel("Data Entered",
                             dataTableOutput('entry')
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
  
  
  get_zeo <- reactive({entry <- data.frame(Adsorbent=input$Adsorbent,
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
    
    
    output$entry <- renderDataTable(get_entryx())

    pca_data <- reactive({
      
            zeo_in <- get_entryx()
            zeo <- FormatZeo(zeo_in)
            entryx <- zeo$get_entry()
            stopifnot(ncol(zeolite_encoded) == ncol(entryx))
            colnames(entryx) <- colnames(zeolite_encoded)
            zeolite_encodedx <- rbind(zeolite_encoded, entryx)
            Scaler <- preProcess(zeolite_encodedx, method = list(center = numeric_cols, scale = numeric_cols))
            zeolite <- predict(Scaler, zeolite_encodedx)
            zeoDB_pca <- prcomp(zeolite)
            pca_df <- data.frame(zeoDB_pca$x)
            label_col <- rbind(zeo_ref[input$discrete], zeo_in[input$discrete])

            pca_df <- cbind(pca_df, label_col)
            pca_df

    })
    
   
observeEvent(input$LoadZeo1,{
     
output$PCAplot1 <- renderPlotly({

            nx <- nrow(zeolite_encoded)
            pca_df <- pca_data() 
            
            ggplotly(ggplot(data = pca_df[c(1:nx),], aes_string(x = "PC1" , y = "PC2", color = input$discrete )) +
                    xlab(label = "PC1") +
                    ylab(label = "PC2") +
                    geom_point(size = 1) +
                    geom_point(data= filter(pca_df, row_number() > nx), aes_string(x = "PC1" , y = "PC2", color = input$discrete), size = 2, shape = 23) +
                    theme(panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          legend.key = element_rect(fill = "white"),
                          legend.position="top",
                          axis.text = element_text(size=18, colour = "black"),
                          axis.title = element_text(size=18, colour = "black"),
                          panel.border = element_rect(colour = "black", fill=NA))) %>%
                            layout(legend = list(orientation = "h", x = 0.4, y = -0.2), width = 750,
                                   height = 750)})

})    


observeEvent(input$LoadZeo2,{
  
  color_tag <- unlist(pca_data()[input$discrete])
  
  output$PCAplot1 <-  renderPlotly(plot_ly(data = pca_data(), 
                                            x = ~PC1, 
                                            y = ~PC2, 
                                            z = ~PC3, 
                                            color = ~unlist(pca_data()[input$discrete]), width = 750,
                                           height = 750)  %>%
                                      add_markers(size = 12))     
                      
})




get_prop <- function(adsorbent){
  
  return  (metal_properties[metal_properties$Adsorbent == adsorbent,] %>%
             mutate(across(everything(), ~replace_na(.x, 0))) %>%
             data.frame() %>%
             select(-c(Adsorbent)))
  
}



}    
shinyApp(ui, server)