#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#' @import shiny Biostrings DT GenomicFeatures GenomicRanges TnT dplyr magrittr plotly readr reshape2 rtracklayer shinythemes stringr
#' @importFrom shinyjs hidden
NULL

# manually load Tnt (should use import once TnT be available from bioconductor 3.8)
require(TnT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme=shinythemes::shinytheme("yeti"),

  navbarPage("repeatCraft",
    # 1st tab panel
    tabPanel("General",
             shinyjs::useShinyjs(),
      sidebarPanel(width=3,
        selectInput("sampleset","Sample Data:",c(" "),selected = " ", multiple = FALSE),
        hr(),
        fileInput("rmgffu","Upload .gff", multiple = FALSE, accept = c(".gff")),
        checkboxInput("rcgffb","GFF from RepeatCraft pipeline?",value = FALSE),
        fileInput("rmoutu","Upload repeatmasker .out", multiple = FALSE, accept = c(".out")),
        hr(),
        h4("Optional files:"),
        fileInput("tegffu","Upload TEclass .gff", multiple = FALSE, accept = c(".gff")),
        fileInput("ltrgffu","Upload LTR_FINER .gff", multiple = FALSE, accept = c(".gff")),
        fileInput("sgffu","Upload SINE_SCAN .gff", multiple = FALSE, accept = c(".gff")),
        shinyjs::hidden(
          div(id="mainoption",
              numericInput("shortsize","Short TEs size",value=100, min=1),
              numericInput("gapsize","TE gap size",value=150, min=1),
              checkboxInput("loosemerge","Loose mode",value = FALSE))
        ),
        shinyjs::hidden(
          div(id="ltroption",
              numericInput("maxltrsize","Max. LTR size",value=10000,min=1),
              numericInput("ltrflanksize","LTR flank size",value=200,min=1),
              checkboxInput("mergeltr","merged gff",value=FALSE))
        ),
        shinyjs::hidden(
          div(id="rcmerge",
              fileInput("rcmerge","Merged GFF from RepeatCraft",multiple = FALSE,accept = c(".gff")))
        ),
        hr(),
        actionButton("upload","Upload",icon=icon("upload")),
        shinyjs::hidden(
          div(id="run",
              actionButton("run","Run TEmerger",icon=icon("angle-double-right")))
        ),
        shinyjs::hidden(
         div(id="downgff",downloadButton("downloadgff","Download TEmerger.gff"))
        )

      ),
        mainPanel(
          tabsetPanel(
            tabPanel("Pie chart",
              fluidRow(checkboxInput("pieUnknown", label = "Include Unknown", value = TRUE)),
              fluidRow(plotly::plotlyOutput("rePie"))
            ),
            tabPanel("Bar chart",plotly::plotlyOutput("reBar",height=700))
          )
        )
    ),

    tabPanel("Family age",
             fluidRow(
               column(4,uiOutput("ageClass"))
             ),
             fluidRow(
               column(4,DT::DTOutput("ageDT")),
               column(8,plotly::plotlyOutput("ageBar",height = 500))
             )
             ),

    # 2nd tab panel
    tabPanel("DataTable",

      # Seletors
      fluidRow(
        column(2,
               uiOutput("chromselector")
               ),
        column(2,
               uiOutput("typeselector")
               ),
        column(2,
               numericInput("minsize","Min.",value=1,min = 1)
               ),
        column(2,
               uiOutput("maxsize")
               ),
        column(2,
               checkboxInput("shortte","show short TEs",value = TRUE),
               checkboxInput("tegroupattr","TEgroup", value = FALSE),
               uiOutput("ltrgroupattr")
               )
      ),
        # Table
        fluidRow(
          column(
           DT::DTOutput("table"), width=12
          )
        )
      ),

    # 3rd tabpanel
    tabPanel("Range Plot",
      fluidRow(), # gff upload
      fluidRow(
        sidebarPanel(width=3,
          uiOutput("plotchrom"),
          uiOutput("plotstart"),
          uiOutput("plotend"),
          hr(),
          checkboxInput("showunknown", label = "show Unknown", value = TRUE),
          uiOutput("radio"),
          helpText("Note:",
                    "set range larger than 100000 will be slow."),
          actionButton("plotgo","Plot")
          # plot option
        ),
        mainPanel(width=9,
                  TnT::TnTOutput("tntplot",width = "auto",height = "90%")
        )
          )
    ),

    # 4th tabpanel
    tabPanel("Repeat plot",
             fluidRow(
               column(2,uiOutput("plot2class"))
             ),
             fluidRow(
               column(2,DT::DTOutput("rptable")),
               column(4,DT::DTOutput("membertable")),
               column(6,
                      tabsetPanel(type="tabs",
                                  tabPanel("Viewer",
                                           fluidRow(fileInput("genegff","Gene Annotation (.gff)",multiple=FALSE,accept = c(".gff"))),
                                           fluidRow(TnT::TnTOutput("rpplot"))),
                                  tabPanel("ORFs",
                                           fluidRow(column(6,fileInput("orfFa","Genome (.fa/.fasta)",multiple = FALSE,accept = c(".fa",".fasta"))),
                                                    column(4,numericInput("orfn","Max. number of ORFs to return",value=10,min=1))
                                                    ),
                                           fluidRow(column(4,actionButton("orfGo","Find ORFs",icon=icon("angle-double-right"))),
                                                    column(4,downloadButton("downloadorf","Download Result (.tsv)")),
                                                    column(4,downloadButton("downloadfa","Fasta"))),
                                           #fluidRow(DT::dataTableOutput("orfdf")),
                                           fluidRow(TnT::TnTOutput("orftnt"))
                                           )

                      ))
             )
    )
  )
))
