library(tidyverse)
library(plotly)
library(caret)
library(shiny)




ui <- shinyUI(fluidPage(

  selectInput("descrete", 
              "Categorical variable:", 
              choices = c("Adsorbent", "adsorbate", "solvent"),
              multiple = FALSE), 
  mainPanel(plotOutput("plot2")),
  
  actionButton("LoadZeo","Load zeolite database")
  
))

server <- shinyServer(function(input, output) {
  
  output$category <-  renderText({ input$descrete })
  
   pca_data <- reactive({
    
    zeo_ref <- read.table("zeolites_db.txt", header = T, sep = "\t") 
    zeolite_encoded <- read.table("ZeoX_Final_encoded_V2x.tsv", header = T, sep = "\t")
    numeric_cols  <- zeo_ref %>%
      select_if(is.numeric) %>%
      colnames()

    Scaler <- preProcess(zeolite_encoded, method = list(center = numeric_cols, scale = numeric_cols))
    zeolite <- predict(Scaler, zeolite_encoded)

    zeoDB_pca <- prcomp(zeolite)
    pca_df <- data.frame(zeoDB_pca$x)
    pca_df <- cbind(pca_df, zeo_ref[input$descrete])
    pca_df
})
  
  

observeEvent(input$LoadZeo,{

   output$plot2 <- renderPlot(ggplot(pca_data(), aes_string(x="PC1",y="PC2", color = input$descrete )) + geom_point())
       
})
 
   
   

})
  
  

shinyApp(ui, server)
























# output$static <- renderTable(head(mtcars))
# observeEvent(input$regen,{
#   rv$m <- data.frame(x=rnorm(n),y=rnorm(n))
# })
# 
# 
# 
# 
# n <- 100  
# rv <- reactiveValues(m=data.frame(x=rnorm(n),y=rnorm(n)))
# 
# observeEvent(input$regen,{
#   rv$m <- data.frame(x=rnorm(n),y=rnorm(n))
# })
# 
# output$myPlot <- renderPlotly({
#   plot_ly() %>%  add_markers(data=rv$m,x=~x,y=~y  )
# })  
# 
# 
# })
# 
# zeolite_encoded <- read.table("ZeoX_Final_encoded_V2x.tsv", header = T, sep = "\t")
# numeric_cols  <- zeo_ref %>%
#   select_if(is.numeric) %>%
#   colnames()
# 
# Scaler <- preProcess(zeolite_encoded, method = list(center = numeric_cols, scale = numeric_cols))
# zeolite <- predict(Scaler, zeolite_encoded)
# 
# zeoDB_pca <- prcomp(zeolite)
# pca_df <- data.frame(zeoDB_pca$x)
# COL <- "Adsorbent"
# pca_df <- cbind(pca_df, zeo_ref[COL])
# pca_df