library(shiny)


ui <- fluidPage(
 
  

  selectInput("Adsorbent", 
              "Adsorbent:", 
              choices = c('CuAgY', 'CuCeY', 'NiCeY', 'AgY', 'clinoptilolite', 'CeY', 'CuY', 'CuX', 'NiY', 'CsY', 'NaY', 'MCM-22', 'AgCeY', 'AgX', 'CuHY'),
              multiple = FALSE),
  
  selectInput("Adsorbate", 
              "Adsorbate:", 
              choices = c('TP', 'BT', 'DBT'),
              multiple = FALSE),
  
  selectInput("Solvent", 
              "Solvent:", 
              choices = c('cyclohexane', '1-octane', 'n-octane', 'n-heptane', 'iso-octane', 'n-Heptane', 'n-Octane', 'ether', 'hexadecane'),
              multiple = FALSE),
  
  numericInput("Sa",
               "Sa",
               min = 141.4,
               max = 720.0,
               value = 0),
  
  numericInput("Vmicro",
               "Vmicro",
               min = 0.08,
               max = 0.666,
               value = 0),
  
  # numericInput("Vmeso",
  #              "Vmeso",
  #              min = 0.04,
  #              max = 0.18,
  #              value = 0),
  # 
  # numericInput("Pore_Size",
  #              "Pore_Size",
  #              min = 0.83,
  #              max = 3.1729,
  #              value = 0),
  # 
  # numericInput("Si_Al",
  #              "Si_Al",
  #              min = 1.031,
  #              max = 46.0,
  #              value = 0),
  # 
  # numericInput("Na+",
  #              "Na+",
  #              min = 0.02,
  #              max = 1.042538354,
  #              value = 0),
  # 
  # numericInput("Ag+",
  #              "Ag+",
  #              min = 0.297,
  #              max = 1.18845,
  #              value = 0),
  # 
  # numericInput("Cu+",
  #              "Cu+",
  #              min = 0.169,
  #              max = 0.467,
  #              value = 0),
  # 
  # numericInput("Ce+4",
  #              "Ce+4",
  #              min = 0.21,
  #              max = 2.29489,
  #              value = 0),
  # 
  # numericInput("Cs+2",
  #              "Cs+2",
  #              min = 0.666666667,
  #              max = 0.666666667,
  #              value = 0),
  # 
  # numericInput("Ni+2",
  #              "Ni+2",
  #              min = 0.024232082,
  #              max = 0.3125,
  #              value = 0),
  # 
  # numericInput("Dipole_Moment",
  #              "Dipole_Moment",
  #              min = 0.51,
  #              max = 0.79,
  #              value = 0),
  # 
  # numericInput("Chemical_Hardness",
  #              "Chemical_Hardness",
  #              min = 3.0401,
  #              max = 5.602,
  #              value = 0),
  # 
  # numericInput("Kinetic_Diameter ",
  #              "Kinetic_Diameter ",
  #              min = 0.77,
  #              max = 0.91,
  #              value = 0),
  
  tableOutput("tablex")
  
)


  
  
               #      get_zeo <- reactive(data.frame(SA=input$SA,
               #                          Vmicro=input$Vmicro,
               #                          Vmeso=input$Vmeso,
               #                          pore_size=input$pore_size,
               #                          Si_Al=input$Si_Al,
               #                          `Na+`=input$`Na+`,
               #                          `Ag+`=input$`Ag+`,
               #                          `Cu+`=input$`Cu+`,
               #                          `Ce+4`=input$`Ce+4`,
               #                          `Cs+2`=input$`Cs+2`,
               #                          `Ni+2`=input$`Ni+2`,
               #                          `x_Na+`=input$`x_Na+`,
               #                          `x_Ag+`=input$`x_Ag+`,
               #                          `x_Cu+`=input$`x_Cu+`,
               #                          `x_Ce+4`=input$`x_Ce+4`,
               #                          `x_Cs+2`=input$`x_Cs+2`,
               #                          `x_Ni+2`=input$`x_Ni+2`,
               #                          `R_Na+`=input$`R_Na+`,
               #                          `R_Ag+`=input$`R_Ag+`,
               #                          `R_Cu+`=input$`R_Cu+`,
               #                          `R_Ce+4`=input$`R_Ce+4`,
               #                          `R_Cs+2`=input$`R_Cs+2`,
               #                          `R_Ni+2`=input$`R_Ni+2`,
               #                          adsorbate=input$adsorbate,
               #                          dipole_moment=input$dipole_moment,
               #                          chemical_hardness=input$chemical_hardness,
               #                          kinetic_diameter=input$kinetic_diameter,
               #                          C_0=input$C_0,
               #                          solvent=input$solvent,
               #                          oil_adsorbent_ratio=input$oil_adsorbent_ratio,
               #                          Temp=input$Temp,
               #                          Capacity=input$Capacity))
               # 
               # metal_props <- reactive(c('x_Na+'=26.2,
               #                           'x_Ag+'=14.6, 
               #                           'x_Cu+'=14.0, 
               #                           'x_Ce+4'=3.038, 
               #                           'x_Cs+2'=0.666666667, 
               #                           'x_Ni+2'=27.2, 
               #                           'R_Na+'=1.02, 
               #                           'R_Ag+'=1.15, 
               #                           'R_Cu+'=0.77, 
               #                           'R_Ce+4'=0.87,
               #                           'R_Cs+2'=1.67, 
               #                           'R_Ni+2'=0.7))  
server <- function(input, output) {


  

}

shinyApp(ui, server)


#     
