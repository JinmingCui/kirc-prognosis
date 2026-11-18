library(shiny)
library(rms)

fit <- readRDS("model.rds")
surv_fun <- Survival(fit)

ui <- fluidPage(
  titlePanel("AQP1 Survival Calculator (TCGA-KIRC)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("aqp1", "AQP1 expression", choices = c("Low", "High")),
      sliderInput("age", "Age (years)", min = 26, max = 90, value = 60),
      selectInput("stage", "Pathologic stage",
                  choices = c("Stage I", "Stage II", "Stage III", "Stage IV"))
    ),
    mainPanel(
      h4("Predicted overall survival"),
      tableOutput("result"),
      p(em("For research and educational purposes only. Not for clinical use."))
    )
  )
)

server <- function(input, output) {
  output$result <- renderTable({
    patient <- data.frame(
      AQP1_group = factor(input$aqp1, levels = c("Low", "High")),
      age = input$age,
      stage = factor(input$stage,
                     levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
    )
    lp <- predict(fit, newdata = patient)
    s <- surv_fun(c(1, 3, 5), lp)
    data.frame(Time = c("1 year", "3 years", "5 years"),
               Survival = paste0(round(100 * s, 1), "%"))
  })
}

shinyApp(ui, server)
