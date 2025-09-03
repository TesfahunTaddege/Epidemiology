############################################################
# Author: Tesfahun Taddege
# Email: ttaddege@gmail.com
# Project: Shiny App: Applied Calculus for Public Health and Epidemiology
# Date: Aug 01, 2025
############################################################

library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Applied Calculus for Public Health and Epidemiology"),
 fluidRow(
    column(12,
           h5("Author: Your Name Here"),
           h5("Email: your.email@domain.com"),
           h5("Organization: Your Organization Here"),
           hr()
    )
  ),
  sidebarLayout(
    sidebarPanel(
      h4("Choose Scenario"),
      selectInput("scenario", "Problem:", 
                  choices = c("Exponential Growth (P2)",
                              "Logistic Vaccination Uptake (P2.2)",
                              "Dose-Response (P3)",
                              "Cumulative Incidence (P4)",
                              "Screening Uptake ODE (P6)",
                              "Resource Allocation (P8)",
                              "Multivariable Surface (P7)",
                              "Extrema of Multivariate Function (P9)")),
      hr(),
      # P2
      conditionalPanel(
        condition = "input.scenario == 'Exponential Growth (P2)'",
        sliderInput("I0", "Initial incidence", 10, 200, 50),
        sliderInput("doubling", "Doubling time (days)", 2, 14, 7),
        sliderInput("tmax", "Time horizon (days)", 10, 100, 30)
      ),
      # P2.2
      conditionalPanel(
        condition = "input.scenario == 'Logistic Vaccination Uptake (P2.2)'",
        sliderInput("k", "Growth rate k", 0.05, 1, 0.3, step=0.05),
        sliderInput("t0", "Midpoint t0", 1, 50, 14),
        sliderInput("tmax_log", "Time horizon (days)", 10, 100, 40)
      ),
      # P3
      conditionalPanel(
        condition = "input.scenario == 'Dose-Response (P3)'",
        sliderInput("beta", "Beta (risk per dose unit)", 0.0005, 0.01, 0.002, step=0.0005),
        sliderInput("dosemax", "Max dose", 10, 500, 100)
      ),
      # P4
      conditionalPanel(
        condition = "input.scenario == 'Cumulative Incidence (P4)'",
        sliderInput("I0_cum", "Initial incidence", 50, 500, 200),
        sliderInput("decay", "Decay rate", 0.01, 0.2, 0.05, step=0.01),
        sliderInput("tmax_cum", "Time horizon (days)", 10, 100, 30)
      ),
      # P6
      conditionalPanel(
        condition = "input.scenario == 'Screening Uptake ODE (P6)'",
        sliderInput("kappa", "Recruitment rate", 0.01, 0.2, 0.08, step=0.01),
        sliderInput("mu", "Attrition rate", 0.01, 0.2, 0.02, step=0.01),
        sliderInput("y0", "Initial uptake fraction", 0, 1, 0.1, step=0.05),
        sliderInput("tmax_screen", "Time horizon (days)", 10, 200, 60)
      ),
      # P8
      conditionalPanel(
        condition = "input.scenario == 'Resource Allocation (P8)'",
        sliderInput("B", "Total staff-hours", 10, 200, 100),
        sliderInput("a1", "Throughput coefficient site 1", 1, 10, 5),
        sliderInput("a2", "Throughput coefficient site 2", 1, 10, 3)
      ),
      # P7
      conditionalPanel(
        condition = "input.scenario == 'Multivariable Surface (P7)'",
        sliderInput("alpha", "Coefficient α", -2, 2, 1, step=0.1),
        sliderInput("beta2", "Coefficient β", -2, 2, 0.5, step=0.1),
        sliderInput("gamma", "Interaction γ", -1, 1, 0.2, step=0.1)
      ),
      # P9
      conditionalPanel(
        condition = "input.scenario == 'Extrema of Multivariate Function (P9)'",
        sliderInput("a", "Coefficient a (x²)", -2, 2, -1, step=0.1),
        sliderInput("b", "Coefficient b (y²)", -2, 2, -1, step=0.1),
        sliderInput("c", "Interaction c (xy)", -2, 2, 0.5, step=0.1)
      )
    ),
    mainPanel(
      plotOutput("plot"),
      verbatimTextOutput("text")
    )
  )
)

server <- function(input, output) {
  
  output$plot <- renderPlot({
    scenario <- input$scenario
    
    if (scenario == "Exponential Growth (P2)") {
      t <- 0:input$tmax
      k <- log(2)/input$doubling
      I <- input$I0*exp(k*t)
      df <- data.frame(t,I)
      ggplot(df, aes(t,I)) + geom_line(color="blue") +
        labs(title="Exponential Incidence Growth",
             x="Time (days)", y="Cases/day")
    }
    
    else if (scenario == "Logistic Vaccination Uptake (P2.2)") {
      t <- 0:input$tmax_log
      V <- 1/(1+exp(-input$k*(t-input$t0)))
      Vprime <- (input$k*exp(-input$k*(t-input$t0)))/((1+exp(-input$k*(t-input$t0)))^2)
      df <- data.frame(t,V,Vprime)
      ggplot(df, aes(t)) +
        geom_line(aes(y=V, color="Coverage")) +
        geom_line(aes(y=Vprime, color="Daily uptake rate")) +
        labs(title="Vaccination Uptake (Logistic)", y="Value") +
        scale_color_manual(values=c("Coverage"="green", "Daily uptake rate"="red"))
    }
    
    else if (scenario == "Dose-Response (P3)") {
      dose <- 1:input$dosemax
      risk <- 1 - exp(-input$beta*dose)
      df <- data.frame(dose,risk)
      ggplot(df, aes(dose,risk)) + geom_line(color="blue") +
        labs(title="Dose-Response: Infection Risk", 
             x="Dose", y="Probability of infection")
    }
    
    else if (scenario == "Cumulative Incidence (P4)") {
      t <- 0:input$tmax_cum
      I <- input$I0_cum*exp(-input$decay*t)
      C <- cumsum(I)
      df <- data.frame(t,I,C)
      ggplot(df, aes(t)) +
        geom_line(aes(y=I, color="Daily incidence")) +
        geom_line(aes(y=C, color="Cumulative cases")) +
        labs(title="Incidence and Cumulative Cases") +
        scale_color_manual(values=c("Daily incidence"="red","Cumulative cases"="blue"))
    }
    
    else if (scenario == "Screening Uptake ODE (P6)") {
      t <- 0:input$tmax_screen
      ystar <- input$kappa/(input$kappa+input$mu)
      y <- ystar + (input$y0 - ystar)*exp(-(input$kappa+input$mu)*t)
      df <- data.frame(t,y)
      ggplot(df, aes(t,y)) + geom_line(color="blue") +
        labs(title="Screening Uptake Dynamics", 
             x="Time (days)", y="Fraction screened")
    }
    
    else if (scenario == "Resource Allocation (P8)") {
      x <- 0:input$B
      y <- input$B - x
      Tval <- input$a1*sqrt(x) + input$a2*sqrt(y)
      df <- data.frame(x,y,Tval)
      ggplot(df, aes(x,Tval)) + geom_line(color="purple") +
        labs(title="Total Tests vs Staff Allocation",
             x="Staff-hours at site 1", y="Total tests")
    }
    
    else if (scenario == "Multivariable Surface (P7)") {
      x <- seq(-5,5,0.2)
      y <- seq(-5,5,0.2)
      grid <- expand.grid(x=x, y=y)
      grid$z <- input$alpha*grid$x^2 + input$beta2*grid$y^2 + input$gamma*grid$x*grid$y
      ggplot(grid, aes(x,y,z=z)) +
        geom_contour_filled() +
        labs(title="Prevalence Surface", x="Exposure", y="SES")
    }
    
    else if (scenario == "Extrema of Multivariate Function (P9)") {
      x <- seq(-5,5,0.2)
      y <- seq(-5,5,0.2)
      grid <- expand.grid(x=x, y=y)
      grid$z <- input$a*grid$x^2 + input$b*grid$y^2 + input$c*grid$x*grid$y
      ggplot(grid, aes(x,y,z=z)) +
        geom_contour_filled() +
        labs(title="Extrema of Multivariate Function",
             x="x", y="y")
    }
  })
  
  output$text <- renderPrint({
    if (input$scenario == "Resource Allocation (P8)") {
      x <- 0:input$B
      y <- input$B - x
      Tval <- input$a1*sqrt(x) + input$a2*sqrt(y)
      opt <- which.max(Tval)
      list(optimal_x = x[opt], optimal_y = y[opt], max_tests = Tval[opt])
    }
    else if (input$scenario == "Extrema of Multivariate Function (P9)") {
      # Hessian test
      H <- matrix(c(2*input$a, input$c,
                    input$c, 2*input$b), 2, 2)
      eig <- eigen(H)$values
      if (all(eig > 0)) result <- "Local minimum"
      else if (all(eig < 0)) result <- "Local maximum"
      else result <- "Saddle point"
      list(Hessian=H, Eigenvalues=eig, Classification=result)
    }
  })
}

shinyApp(ui, server)


