library(shiny)
library(shinythemes)
library(plotly)
library(DT)
library(readr)
library(stringr)
library(tidyr)
library(shinyjs)
library(sortable)
library(tibble)
library(rclipboard)
library(dplyr)
library(limma)

options(shiny.maxRequestSize = 10000*1024^2)


ui <- shinyUI(navbarPage(title = "Linear Modeler",
                         theme=shinytheme('flatly'),
                         sidebarLayout(
                           sidebarPanel(
                             fileInput("datfile", label = "Locate RDS object"),
                             actionButton('load_results','Load limma Session'),
                             br(),
                             br(),
                             selectInput('itiss','Select Tissue',
                                         choices=c('heart','brain',
                                                   'liver', 'kidney'), 
                                         selected = 'heart'),
                             br(),
                             selectInput('ngenes','Choose N Genes to Plot',
                                         choices=c(50,100,200,300,500,1000),
                                         selected=100),
                             br(),
                             actionButton('run_contr','Run Contrast'),
                            width = 2),
                           mainPanel(
                                  fluidRow(
                                  
                                     column(5,
                                      
                                      checkboxGroupInput('group_in','Select Members of Group',
                                                         choices=c("SummerActive","30C","25C","20C","12C","4C"),
                                                        inline=T),
                                      br(),
                                      radioButtons('within_rel','Choose Within Group Relation',
                                                   choices = c('Difference', 'Average'),
                                                   selected = 'Difference'),
                                      br(),
                                      actionButton('add_group','Add Group'),
                                      br(),
                                      fluidRow(
                                      column(3,uiOutput('group_select'),uiOutput('group2_select')),
                                      column(3,br(),br(),radioButtons('tsign','t sign', choices=c('+','-'),selected = '+'))
                                      ),
                                      br(),
                                      DT::dataTableOutput("groups"),
                                      br()
                                      ),
                                     
                                  column(7,
                                         plotOutput("model_graph"),
                                         br(),
                                         DT::dataTableOutput("contrast_results")
                                        )
                                  ),
                                  
                                  uiOutput('datable')
                           )
                                
                                  
                        )
))