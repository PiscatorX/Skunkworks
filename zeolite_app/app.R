library(shiny)

ui <- fluidPage(

  selectInput("Metals", 
              "Metals:", 
              choices = c('Na+', 'Ag+', 'Cu+', 'Ce+4', 'Cs+2', 'Ni+2'),
              multiple = TRUE),
  
  sliderInput("Sa",
              "Sa",
              141.4,
              720.0,
              492.55379213483155),
  
  sliderInput("Vmicro",
              "Vmicro",
              0.08,
              0.6659999999999999,
              0.2569763636363636),
  
  sliderInput("Vmeso",
              "Vmeso",
              0.04,
              0.18,
              0.0808421052631579),
  
  sliderInput("Pore_Size",
              "Pore_Size",
              0.83,
              3.1729,
              1.3090391509433963),
  
  sliderInput("Si_Al",
              "Si_Al",
              1.031,
              46.0,
              7.295000006952248),
  
  sliderInput("Na+",
              "Na+",
              0.02,
              1.0425383540000002,
              0.49715446279710146),
  
  sliderInput("Ag+",
              "Ag+",
              0.297,
              1.18845,
              0.8477998749909909),
  
  sliderInput("Cu+",
              "Cu+",
              0.16899999999999998,
              0.467,
              0.2289083723975904),
  
  sliderInput("Ce+4",
              "Ce+4",
              0.21,
              2.29489,
              0.9387571653543307),
  
  sliderInput("Cs+2",
              "Cs+2",
              0.666666667,
              0.666666667,
              0.666666667),
  
  sliderInput("Ni+2",
              "Ni+2",
              0.024232082000000002,
              0.3125,
              0.10820202342857144),
  
  sliderInput("Dipole_Moment",
              "Dipole_Moment",
              0.51,
              0.79,
              0.6308426966292134),
  
  sliderInput("Chemical_Hardness",
              "Chemical_Hardness",
              3.0401,
              5.602,
              4.074383707865168),
  
  sliderInput("Kinetic_Diameter ",
              "Kinetic_Diameter ",
              0.77,
              0.91,
              0.8232584269662919),
  
  

)

server <- function(input, output, session) {
  
}

shinyApp(ui, server)