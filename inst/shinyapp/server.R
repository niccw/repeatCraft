#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

require(magrittr)
require(ggplot2)
# manually load Tnt (should use import once TnT be available from bioconductor 3.8)
require(TnT)

shinyServer(function(input, output,session) {


#' @import rtracklayer
#' @importFrom shinyjs show
NULL

    # Compute these first -----------------------------------------------------

    # Sample
    checksample <- reactive({
        if(input$sampleset == "Hydra magnipappilata"){
        # read gff to GRanges
          #tegff <<- rtracklayer::import.gff(tegff_p)
          #ltrgff <<- rtracklayer::import.gff(ltrgff_p)
          #sgff <<- rtracklayer::import.gff(sgff_p)
          #rmgff <<- rtracklayer::import.gff(rmgff_p)
          #rmout <<- readr::read_table2(rmout_p,col_names = FALSE,skip = 3)
          #mergermgff <<- import.gff3(mergermgff_p)
          #LTRlabel(TRUE)
          #fapath <- "hsym.fa"
        }
    })

    # Upload
    checkupload <- reactive({

      input$go
      input$upload
      if(input$upload == 0){
        return()
      }
      if(input$run == 0){
        return()
      }
 
      isolate({
        if((!is.null(input$rmgffu)&(!is.null(input$rmoutu)))){
          print("inside checkupload")
          # Run temerger in background
          # Check if LTR.gff is uploaded
          if(!is.null(input$ltrgffu)){
            ltrgffp <- NULL
          }else{
            ltrgffp <- input$ltrgffu$datapath
          }

          temergewrap(input$rmgffu$datapath,input$rmoutu$datapath,ltrgffp,input$shortsize,input$gapsize,input$maxltrsize,input$ltrflanksize,input$mergeltr,input$loosemerge)

          rmout <<- readr::read_table2(input$rmoutu$datapath,col_names = FALSE,skip = 3)
         
           if(!is.null(input$rcmerge)){
            mergermgff <<- rtracklayer::import.gff3(input$rcmerge$datapath)
          }else if (!runrc){
            tntshowmerge <<- FALSE
          }

          # Optional file
          if (!is.null(input$tegffu)){
            tegff <<- import.gff(input$tegffu$datapath)
          }else{
            tegff <<- NULL
          }
          if (!is.null(input$ltrgffu)){
            ltrgff <<- import.gff(input$ltrgffu$datapath)
          }else{
            ltrgff <<- NULL
          }
          if (!is.null(input$sgffu)){
            sgff <<- import.gff(input$sgffu$datapath)
          }else{
            sgff <<- NULL
          }
        }
      })
    })

    # Prepare data

    parsefile <- reactive({

      # Check input source
      checksample()
      checkupload()

      if(!exists("rmgff")){
        return()
      }

      print("before parsing rmgff")
      # Parse rmgff -------------------------------------------------------------------
      if("LTRgroup" %in% colnames(rtracklayer::mcols(rmgff))){
        rmgff$LTRgroup[sapply(rmgff$LTRgroup, function(x) length(x)==0)] <- NA # todo: find a better way to flat compressed list object
        # parse row with more than one label, prevent unconsistant length after unlist
        parsemultiple <-function(r){
          if (length(r) > 1){
            return(paste(r,collapse = "|"))
          }
          else{
            return(r)
          }
        }
        rmgff$LTRgroup <- lapply(rmgff$LTRgroup,parsemultiple)
        rmgff$LTRgroup <- unlist(rmgff$LTRgroup)
        LTRlabel(TRUE)
      }else{
        LTRlabel(FALSE)
      }

      # rmgff to dataframe for datatable output
      gff <<- as.data.frame(rmgff)
      gff <- dplyr::mutate(gff,size = end-start)
      if(LTRlabel()){
        gffdf <<- gff[,c("seqnames","type","start","end","strand","ID","TEgroup","LTRgroup","shortTE","size")]
      }else{
        gffdf <<- gff[,c("seqnames","type","start","end","strand","ID","TEgroup","shortTE","size")]
      }
      })

      # TnT track ---------------------------------------------------------------
    tntT <- reactive({
      if(is.null(rmgff)){
        return()
      }

      rmftrack <<- TnT::FeatureTrack(rmgff,tooltip = as.data.frame(rmgff), name=paste(rmgff$ID,":",rmgff$type))

      if(!is.null(tegff)){
        teftrack <<- TnT::FeatureTrack(tegff,tooltip = as.data.frame(tegff), name=paste(tegff$ID,":",tegff$type))
        teftrack$color <<- "green"
      }

      if(!is.null(ltrgff)){
        ltrftrack <<- TnT::FeatureTrack(ltrgff,tooltip = as.data.frame(ltrgff), name=ltrgff$type)
        ltrftrack$color <<- "red"
      }
      if(!is.null(sgff)){
        sineftrack <<- TnT::FeatureTrack(sgff,tooltip = as.data.frame(sgff))
        sineftrack$color <<- "orange"
      }

      # tnt rmftrack color

      teuniqV <- unique(gff$TEgroup)
      temapc <- mapc(teuniqV)
      tecolor <<- temapc[match(gff$TEgroup,teuniqV)]

      if(LTRlabel()){
        ltruniqV <- unique(gff$LTRgroup)
        ltrmapc <- mapc(ltruniqV)
        ltrcolor <<- ltrmapc[match(gff$LTRgroup,ltruniqV)]
      }

      # mergermgfftrack color
      if(tntshowmerge){
        mergeV <- unique(mergermgff$TEgroup)
        mergec <- mapc(mergeV)
        mergecolor <<- mergec[match(mergermgff$TEgroup,mergeV)]
      }

  })


    # Selector object for dynamic UI ------------------------------------------

    # dynamic inputs used in general
    output$downloadgff <- renderUI({
      if(dgff()){
        downloadButton("downloadgff","Download TEmerger .gff",icon=icon("download"))
      }
    })

    # dynamic inputs used in datatable
    output$chromselector <- renderUI({
      validate(need(exists("rmgff"),""))
      selectizeInput("chrom","Chrom",c("All",as.character(unique(gff$seqnames))),selected="All",multiple=TRUE)
    })
    output$typeselector <- renderUI({
      validate(need(exists("rmgff"),""))
      selectizeInput("repeattype","Repeat Class:",c("All",as.character(unique(gff$type))),selected="All",multiple=TRUE)
    })
    output$maxsize <- renderUI({
      validate(need(exists("rmgff"),""))
      numericInput("maxsize","Max.",value = max(gffdf$size), max=max(gffdf$size))
    })

    output$ltrgroupattr <- renderUI({
      if(LTRlabel()){
        checkboxInput("ltrgroupattr","LTRgroup", value = FALSE)
      }
    })

    # dynamic inputs used in the plot
    output$plotchrom <- renderUI({
      validate(need(exists("rmgff"),""))
      selectInput("pchrom","Chrom",unique(gff$seqnames),multiple = FALSE)
    })
    # start and end depend on chrom
    chromsize <- reactive({
      validate(need(exists("rmgff"),""))
      tail(end(rmgff[GenomicRanges::seqnames(rmgff)==input$pchrom]),1)
    })



    output$plotstart <- renderUI({
      validate(need(exists("rmgff"),""))
      numericInput("pstart","Start",value=1, max=chromsize())
    })

    output$plotend <- renderUI({
      validate(need(exists("rmgff"),""))
      numericInput("pend","End",value = chromsize(),max=chromsize())
    })

    output$radio <- renderUI({
      if(LTRlabel()){
        radioButtons("radio",h4("Group label"),choices = list("LTR group"=1,"TE group"=2,"None"=3), selected = 3)
      }else{
        radioButtons("radio",h4("Group label"),choices = list("TE group"=2,"None"=3), selected = 3)
      }
    })

    # dynamic inputs used in repeatplot
    output$plot2class <- renderUI({
      validate(need(exists("rmgff"),""))
      selectInput("plot2class","Repeat Class",unique(gff$t), multiple = FALSE)
    })

    # dynamic inputs used in age plot
    output$ageClass <- renderUI({
      validate(need(exists("rmgff"),""))
      selectInput("ageClass","Repeat Class",unique(gff$t), multiple = FALSE)
    })


    # Reactive ----------------------------------------------------------------

    classdfre <- reactive({

      validate(need(exists("rmgff"),""))

      parsefile()

      classdf <<- gff[gff$type == input$plot2class,"ID"] %>% table %>% as.data.frame()
      try({colnames(classdf) <<- c("Family","Count")}, silent = TRUE)
      try({classdf <<- classdf[base::rev(base::order(classdf$Count)),]},silent = TRUE)
    })

    subdfre <- reactive({
      classdfre()
      family <- classdf[input$rptable_rows_selected,1]
      subdf <- gff[gff$ID == family,c("seqnames","start","end","Tstart","Tend")]
      subdf <- dplyr::mutate(subdf,csize = as.integer(Tend)-as.integer(Tstart))
      subdf <<- subdf[base::rev(base::order(subdf$csize)),c("seqnames","start","end","csize")]
    })

    ageclassre <- reactive({

      parsefile()

      ageclassdf <- gff[gff$type == input$ageClass,"ID"] %>% table %>% as.data.frame()
      colnames(ageclassdf) <- c("Family","Count")
      ageclassdf <<- ageclassdf[base::rev(base::order(ageclassdf$Count)),]
    })

    orf <- reactive({
      # Make sure user have selected a repeat
      if(length(input$membertable_rows_selected) == 0 ){
        return()
      }

      # Input: repeat coordinate
      rfamily <- classdf[input$rptable_rows_selected,1]
      pchrom <- subdf$seqnames[input$membertable_rows_selected]
      pstart <- as.integer(subdf$start[input$membertable_rows_selected])
      pend <- as.integer(subdf$end[input$membertable_rows_selected])

      # Store the location for file name
      filechrom(pchrom)
      filestart(pstart)
      fileend(pend)

      # Extract seq by get1seq.py
      getseq(input$orfFa$datapath,pchrom,pstart,pend,rfamily)

      # Find ORF
      orfdf <<- as.data.frame(systemPipeR::predORF(targetDNAstring, n=input$orfn, mode="orf", longest_disjoint=TRUE, strand="both"))
    })



    # Observe -----------------------------------------------------------------

    # Run TEmerger button
    observe({
      req(input$upload)
      print("after req upload")
      if(!input$rcgffb){
        shinyjs::show("mainoption")
        if(!is.null(input$ltrgffu)){
          shinyjs::show("ltroption")
        }
        shinyjs::show("run")
      }
    })
    
    # Upload output
    observe({
      req(input$upload)
      print("upload output")
      if(input$rcgffb){
        shinyjs::show("rcmerge")
        shinyjs::show("run")
      }
    })


    # Download button
    observe({
      if(dgff() == TRUE){
        shinyjs::show("downgff")
      }
    })



    # Outputs -----------------------------------------------------------------

    # outputs in general
    output$rePie <- plotly::renderPlotly({

      parsefile()

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      if(input$pieUnknown == FALSE){
        t2 <- gff$t[gff$t != "Unknown"]
        t2 <- gsub("$",";",t2)
      }else{
        t2 <- gsub("$",";",gff$t)
      }

      typetable <- stringr::str_extract(t2,".*?(?=/|;)") %>% table() %>% as.data.frame()
      colnames(typetable) <- c("Class","Count")

      plotly::plot_ly(data = typetable, labels = ~`Class`, values= ~Count,type="pie",textposition = "inside",
                   textinfo = "label+percent" )

    })

    output$reBar <- plotly::renderPlotly({

      parsefile()
      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      subrmout <- rmout[,c("X2","X11")]

      subrmout$X11 <- gsub("$",";",subrmout$X11)
      subrmout$X11 <- stringr::str_extract(subrmout$X11,".*?(?=/|;)")
      plotly::plot_ly(subrmout,x=~X2,type="histogram",alpha=0.7,color=~X11, xbins=list(start=0,end=round(max(subrmout$X2)), size =0.5)) %>% plotly::layout(barmode="stack",xaxis=list(title="percentage of div."), yaxis =list(title="count"))

    })

    # outputs in datatable
    output$table <- DT::renderDataTable({

      parsefile()
      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      # Filter df based on input
      # if tableattr = 0, all not selected; 1: only ltr; 3: only te; 4: both
      tableattr <- 0  # ltr:1  te:3

      print(input$tegroupattr)
      if(input$tegroupattr){
        tableattr <- tableattr + 3
      }

      # Make sure input$ltrgoupattr exist (if no LTRgrou in rmgff, checkbox would not be made)
      if(LTRlabel()){
        if(input$ltrgroupattr){
          tableattr <- tableattr + 1
        }
      }

      # Filter step by step
      if(input$chrom != "All"){
        gffdf <- dplyr::filter(gffdf, seqnames == input$chrom)
      }
      if(input$repeattype != "All"){
        gffdf <- dplyr::filter(gffdf,type == input$repeattype)
      }


      gffdf <- dplyr::filter(gffdf,size >= input$minsize)
      gffdf <- dplyr::filter(gffdf,size <= input$maxsize)

      if(input$shortte == FALSE){
        gffdf <- dplyr::filter(gffdf,shortTE == "F")
      }
      if(tableattr == 1){
        gffdf <- dplyr::filter(gffdf,!is.na(LTRgroup))
      }
      if(tableattr == 3){
        gffdf <- dplyr::filter(gffdf,!is.na(TEgroup))
      }
      if(tableattr == 4){
        gffdf <- dplyr::filter(gffdf, ((!is.na(LTRgroup)) & (!is.na(TEgroup))))
      }

    # print datatable
      DT::datatable(gffdf)
    })

    # output in regionplot
    output$tntplot <- TnT::renderTnT({


      validate(need(exists("rmgff"),"Please select sample data / upload files"))
      parsefile()
      tntT()

      input$plotgo

      if(input$plotgo == 0){
        return()
      }

      isolate({ # Only execute when plotgo is clicked

        # Filter rmgff (show unknown and label)
        if(input$showunknown == FALSE){
          rmgff_nounknown <- rmgff[rmgff$type != "Unknown"]
          # labeling
          if (as.integer(input$radio) == 3){
            rmftrack_nounknown <- TnT::FeatureTrack(rmgff_nounknown,tooltip = as.data.frame(rmgff_nounknown), name=paste(rmgff_nounknown$ID,":",rmgff_nounknown$type))
          }
          else if(as.integer(input$radio) == 1){
            rmftrack_nounknown <- TnT::FeatureTrack(rmgff_nounknown,tooltip = as.data.frame(rmgff_nounknown), name=paste(rmgff$ID,":",rmgff$type), color =ltrcolor)
          }
          else if (as.integer(input$radio) == 2){
            rmftrack_nounknown <- TnT::FeatureTrack(rmgff_nounknown,tooltip = as.data.frame(rmgff_nounknown), name=paste(rmgff$ID,":",rmgff$type), color = tecolor)
          }
        }
        else{
          if (as.integer(input$radio) == 3){
            rmftrack <- TnT::FeatureTrack(rmgff,tooltip = as.data.frame(rmgff), name=paste(rmgff$ID,":",rmgff$type))
          }
          else if(as.integer(input$radio) == 1){
            rmftrack <- TnT::FeatureTrack(rmgff,tooltip = as.data.frame(rmgff), name=paste(rmgff$ID,":",rmgff$type), color =ltrcolor)
          }
          else if (as.integer(input$radio) == 2){
            rmftrack <- TnT::FeatureTrack(rmgff,tooltip = as.data.frame(rmgff), name=paste(rmgff$ID,":",rmgff$type), color = tecolor)
          }
        }

        # Filter mergermgff 
        if(tntshowmerge){
          if(input$showunknown == FALSE){
            mergermgff_nounknown <- mergermgff[mergermgff$type != "Unknown"]
            mergermgfftrack_nounknown <-  TnT::FeatureTrack(mergermgff_nounknown,tooltip = as.data.frame(mergermgff_nounknown), name=paste(mergermgff_nounknown$ID,":",mergermgff_nounknown$type), color = mergecolor)
          }else{
            mergermgfftrack <- TnT::FeatureTrack(mergermgff,tooltip = as.data.frame(mergermgff), name=paste(mergermgff$ID,":",mergermgff$type), color = mergecolor)
          }
        }


        # Only show track with TE
        trackls <<- c()
        checkrecord <- function(gff,tk){
          if (length(gff[GenomicRanges::seqnames(gff)==input$pchrom & rtracklayer::start(gff) >= input$pstart & rtracklayer::end(gff)<= input$pend]) >= 1){
            trackls <<- c(trackls,tk)
          }
        }

        if(input$showunknown == FALSE){
           checkrecord(rmgff_nounknown,rmftrack_nounknown)
          if(tntshowmerge){
            checkrecord(mergermgff_nounknown,mergermgfftrack_nounknown)
          }
         }else{
           checkrecord(rmgff,rmftrack)
           if(tntshowmerge){
           checkrecord(mergermgff,mergermgfftrack)
           }
         }

        if(!is.null(tegff)){
          checkrecord(tegff,teftrack)
        }
        if(!is.null(ltrgff)){
          checkrecord(ltrgff,ltrftrack)
        }
        if(!is.null(sgff)){
          checkrecord(sgff,sineftrack)
        }

        # todo: use validate to show message
        if(is.null(trackls)){
          print("empty trackls")
          return()
        }

        TnT::TnTGenome(trackls,view.range = GenomicRanges::GRanges(input$pchrom, IRanges::IRanges((input$pstart)+200,(input$pend)+200)) )
      })
    })

    # output in repeatplot
    output$rptable <- DT::renderDataTable({

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      classdfre()
      DT::datatable(classdf,selection = "single",rownames=FALSE, extensions = "Scroller",options= list(scrollY=500,scroller=TRUE))
    })

    output$membertable <- DT::renderDataTable({

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      classdfre()
      subdfre()
      DT::datatable(subdf,selection = "single",extensions = "Scroller",rownames=FALSE,options = list(searching = FALSE,scrollY = 500,scroller = TRUE))
    })

    output$rpplot <- TnT::renderTnT({

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      classdfre()
      subdfre()
      tntT()

      pchrom <- subdf$seqnames[input$membertable_rows_selected]
      pstart <- subdf$start[input$membertable_rows_selected]
      pend <- subdf$end[input$membertable_rows_selected]
      TnT::TnTGenome(rmftrack,view.range = GenomicRanges::GRanges(pchrom, IRanges::IRanges((pstart)-1000,(pend)+1000)) )
    })

    output$orfdf<- DT::renderDataTable({

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      input$orfGo
      if(input$orfGo == 0){
        return()
      }

      if(is.null(input$orfFa)){
        return()
      }

      isolate({
        classdfre()
        subdfre()

        # Make sure user have selected a repeat
        orf()
        DT::datatable(orfdf)
      })
    })

    output$orftnt <- TnT::renderTnT({
      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      input$orfGo
      if(input$orfGo == 0){
        return()
      }

      if(is.null(input$orfFa)){
        return()
      }

      isolate({
        classdfre()
        subdfre()
        orf()

        # Reverse dataframe to GRanges (6 frames)
        pchrom <- subdf$seqnames[input$membertable_rows_selected]
        pstart <- as.integer(subdf$start[input$membertable_rows_selected])
        pend <- as.integer(subdf$end[input$membertable_rows_selected])

        orfdf$seqnames <- pchrom
        orfdf$start <- pstart + orfdf$start
        orfdf$end <- pstart + orfdf$end

        porfdf1 <- orfdf[(orfdf$strand == "+") & (orfdf$inframe2end==1),]
        porfdf2 <- orfdf[(orfdf$strand == "+") & (orfdf$inframe2end==2),]
        porfdf3 <- orfdf[(orfdf$strand == "+") & (orfdf$inframe2end==3),]
        norfdf1 <- orfdf[(orfdf$strand == "-") & (orfdf$inframe2end==1),]
        norfdf2 <- orfdf[(orfdf$strand == "-") & (orfdf$inframe2end==2),]
        norfdf3 <- orfdf[(orfdf$strand == "-") & (orfdf$inframe2end==3),]

        p1 <- GenomicRanges::makeGRangesFromDataFrame(porfdf1,keep.extra.columns = FALSE)
        p2 <- GenomicRanges::makeGRangesFromDataFrame(porfdf2,keep.extra.columns = FALSE)
        p3 <- GenomicRanges::makeGRangesFromDataFrame(porfdf3,keep.extra.columns = FALSE)
        n1 <- GenomicRanges::makeGRangesFromDataFrame(norfdf1,keep.extra.columns = FALSE)
        n2 <- GenomicRanges::makeGRangesFromDataFrame(norfdf2,keep.extra.columns = FALSE)
        n3 <- GenomicRanges::makeGRangesFromDataFrame(norfdf3,keep.extra.columns = FALSE)

        p1t <- TnT::FeatureTrack(p1,tooltip = as.data.frame(p1),color="red")
        p2t <- TnT::FeatureTrack(p2,tooltip = as.data.frame(p2),color="orange")
        p3t <- TnT::FeatureTrack(p3,tooltip = as.data.frame(p3),color="yellow")
        n1t <- TnT::FeatureTrack(n1,tooltip = as.data.frame(n1),color="blue")
        n2t <- TnT::FeatureTrack(n2,tooltip = as.data.frame(n2),color="purple")
        n3t <- TnT::FeatureTrack(n3,tooltip = as.data.frame(n3),color="grey")


        dfls <- list(porfdf1,porfdf2,porfdf3,norfdf1,norfdf2,norfdf3)
        trackls2check <- list(p1t,p2t,p3t,n1t,n2t,n3t)

        trackls <- c(rmftrack)

        for (i in 1:length(dfls)){
          if(nrow(dfls[[i]])>0){
            trackls <- c(trackls,trackls2check[i])
          }
        }

        TnT::TnTGenome(trackls,view.range = GRanges(pchrom, IRanges::IRanges((pstart-200),(pend+200))))
      })
    })

    # outputs for ageplot
    output$ageDT <- DT::renderDataTable({

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      ageclassre()
      # ageclassdf <- gff[gff$t == input$ageClass,"family"] %>% table %>% as.data.frame()
      # colnames(ageclassdf) <- c("Family","Count")
      # ageclassdf <- ageclassdf[base::rev(base::order(ageclassdf$Count)),]
      DT::datatable(ageclassdf,selection = "single", rownames = FALSE,extensions = "Scroller",options = list(scrollY = 500,scroller = TRUE))
    })

    output$ageBar <- plotly::renderPlotly({

      validate(need(exists("rmgff"),"Please select sample data / upload files"))

      if(length(input$ageDT_rows_selected)==0){
        return()
      }
      ageclassre()
      f <- ageclassdf$Family[input$ageDT_rows_selected]
      testdf <- rmout[rmout$X10 == f ,"X2"]
      p <- ggplot2::ggplot(testdf,aes(X2)) + ggplot2::geom_histogram(color="#337175",fill="#00AFBB",bins = round(max(testdf)*2)) +
        ggplot2::scale_x_continuous(breaks =  pretty(testdf$X2,n=max(testdf$X2/2))) + ggplot2::xlab("percentage of div.") + ggplot2::ggtitle(as.character(f))+ ggplot2::theme_classic()
      plotly::ggplotly()
    })




    # Export ------------------------------------------------------------------


    # Export reformated gff (after TEmerger)
    output$downloadgff <- downloadHandler(
      filename = "TEmerger.gff",
      content = function(file){
        rtracklayer::export.gff3(rmgff,file,"")
      }
    )

    # Export ORFs table
    output$downloadorf <- downloadHandler(
      filename =  function(){
        paste0(filechrom(),"_",filestart(),"_",fileend(),".tsv")
      },
      content = function(file){
        write.table(orfdf,file,sep="\t",row.names = FALSE)
      }
    )

    # Export selected TE.fa
    output$downloadfa <- downloadHandler(
      filename =  function(){
        paste0(filechrom(),"_",filestart(),"_",fileend(),".fa")
      },
      content = function(file){
        tigger::writeFasta(targetDNAstring,file,width=60)
      }
    )


    # Cleaning ----------------------------------------------------------------

    session$onSessionEnded(
      function(){
        clean1 <- paste0("rm -r"," ",parsedgff)
        clean2 <- paste0("rm -r"," ",parsedgff_b)
        clean3 <- paste0("rm -r"," ",fuseltrgff)
        clean4 <- paste0("rm -r"," ",fuseltrgff_b)
        clean5 <- paste0("rm -r"," ",tobesort_label)
        clean6 <- paste0("rm -r"," ",tobesort_merge)
        system(clean1)
        system(clean2)
        system(clean3)
        system(clean4)
        system(clean5)
        system(clean6)
      }
    )

})
