#' @import shiny Biostrings DT GenomicFeatures GenomicRanges TnT dplyr magrittr plotly readr reshape2 rtracklayer shinythemes stringr
#' @importFrom shinyjs hidden
#' @importFrom systemPipeR predORF
NULL

# change upload limit (1gb)
options(shiny.maxRequestSize=1000*1024^2)

# Flags
id <- NULL
runrc <- FALSE
tntshowmerge <- TRUE
LTRlabel <- reactiveVal(TRUE)
dgff <- reactiveVal(FALSE)
mainoptionf <- reactiveVal(FALSE)
ltroptionf <- reactiveVal(FALSE)

# ReactiveVal
filechrom <- reactiveVal()
filestart <- reactiveVal()
fileend <- reactiveVal()

# todo: add the switch that if some file is missing, skip the downstream processing (add flag?)
# check the rm.gff format (if it is the original .gff, parse it using refomatRmGff.py

# Sample input of hsym
# sampledir <- system.file("shinyapp/exampleData", package = "repeatcraft")
# tegff_p <- paste0(sampledir,"/hsym_teclass.gff")
# ltrgff_p <- paste0(sampledir,"/hsym_onlyFullLTR.fuse_label.gff")  # with group label
# sgff_p <- paste0(sampledir,"/hsym_sine.gff")
# rmgff_p <- paste0(sampledir,"/hsym.fa.reformat.combine.gff")
# rmout_p <- paste0(sampledir,"/hsym.fa.out")
#fa <- paste0(sampledir,"/hsym.fa")


# Functions ---------------------------------------------------------------
# Color for group track
mapc <-function(uniV){
  mappedC <- list()
  for (i in uniV){
    if (is.na(i)| i == " "){
      mappedC <- c(mappedC,"#000000")
    }else{
      mappedC <- c(mappedC,"#4286f4")
    }
  }
  return(unlist(mappedC))
}


# Run Temerger  ---------------------------------------------
scriptdir <- paste0(system.file("pythonScript", package = "repeatcraft"),"/")

# tmp files path
parsedgff <<- tempfile()
parsedgff_b <<- tempfile()
fuseltrgff <<- tempfile()
fuseltrgff_b <- tempfile()
tobesort_label <- tempfile()
tobesort_merge <- tempfile()


temergewrap <- function(inputpath,outinputpath,ltrinputpath,shortsize,gapsize,maxltrsize,ltrflanksize,mergeltr,loosemerge){
  
  gffsample <- readr::read_lines(inputpath,n_max=5)
  print(gffsample)
  # Check repeatmasker.gff
  if(!(startsWith(gffsample[2],"##repeatcraft"))){
    showNotification("Raw .gff from repeatmasker. Running TEmerger ...",type="message",duration=2)
    runrc <<- TRUE
  }
  if(startsWith(gffsample[2],"##repeatcraft")){
    showNotification("RepeatCraft.gff, reading...", type="message", duration=5)
    rmgff <<- rtracklayer::import.gff3(inputpath)
    return()
  }

  # make sure python3 is installed
  showNotification("Raw .gff from repeatmasker. Running TEmerger ...",type="message",duration=5)
  if(length(system("which python3"))==0){
    showNotification("Please install python3",type="error",closeButton=FALSE)
  }else{ # start
    showNotification("Step 1: reformat gff...",type="message",duration=2)
    temerge1 <- paste0("python3"," ",scriptdir,"reformatRMgff.py"," ",inputpath," ",outinputpath,">",parsedgff)
    system(temerge1)

    showNotification("Step 2: add short read label using default value...",type="message",duration=2)
    temerge2 <- paste0("python3"," ",scriptdir,"filterShortTE.py"," ","-r"," ",parsedgff," ","-s"," ",shortsize," > ",parsedgff_b)
    print(temerge2)
    system(temerge2)

    if(!exists("ltrgff")){
      showNotification("LTR_FINDER gff is missing, skip TEmerge step 3...",type="warning",duration=10)
      missltr <- TRUE
    }else{
      showNotification("Step 3: add LTRgroup label using default value...",type="message",duration=5)

      # check if the gff has been fused or not
      if(!mergeltr){
      # filter and fuse the ltrgff
        ltrfuse1 <- paste0("grep 'LTR_retrotransposon'"," ",ltrinputpath,">",fuseltrgff)
        system(ltrfuse1)
        ltrfuse2 <- paste0("cut -f1-8"," ",fuseltrgff ,"|sort|uniq >",fuseltrgff_b)
        system(ltrfuse2)
        ltrfuse3 <- paste0("sort -k1,1 -k4,4n -k5,5n"," ",fuseltrgff_b,">",fuseltrgff)
        system(ltrfuse3)
        ltrfuse4 <- paste0("python3"," ",scriptdir,"combineGFFoverlap.py"," ",fuseltrgff," ","LTR_retrotransposon"," ","T >",fuseltrgff_b)
        temerge3 <-  paste0("python3"," ",scriptdir,"fuseLtrFinderResult.py"," ","-r"," ",parsedgff_b," ","-l"," ",fuseltrgff_b," ","-m"," ",maxltrsize," ","-f"," ",ltrflanksize," >",parsedgff)
      }else{
        temerge3 <-  paste0("python3"," ",scriptdir,"fuseLtrFinderResult.py"," ","-r"," ",parsedgff_b," ","-l"," ",input$ltrgff$datapath," ","-m"," ",maxltrsize," ","-f"," ",ltrflanksize," >",parsedgff)
      }
       missltr <- FALSE
    }
    
    if (!loosemerge){ # strict merge (add label)
      showNotification("Step 4: add TEgroup label  (strict merge)...",type="message",duration=5)
      if (missltr){
        temerge4 <- paste0("python3"," ",scriptdir,"combineTE.py"," ",parsedgff_b," ",gapsize," ",">",parsedgff)
        system(temerge4)
        rmgff <<- rtracklayer::import.gff(parsedgff)
        dgff(TRUE)
      }else{
        temerge4 <- paste0("python3"," ",scriptdir,"combineTE.py"," ",parsedgff," ",gapsize," ",">",parsedgff_b)
        system(temerge4)
        rmgff <<- rtracklayer::import.gff(parsedgff_b)
        dgff(TRUE)
      }
      # true merge
      showNotification("Step 5: Merge gff row according to TEgroup label", type ="message", duration = 3)
      if(missltr){
        temerge5 <- paste0("python3"," ",scriptdir,"trueMerge.py"," ",parsedgff," ",">",parsedgff_b)
        system(temerge5)
        mergermgff <<- rtracklayer::import.gff(parsedgff_b)
      }else{
        temerge5 <- paste0("python3"," ",scriptdir,"trueMerge.py"," ",parsedgff_b," ",">",parsedgff )
        system(temerge5)
        mergermgff <<- rtracklayer::import.gff(parsedgff)
      }
      
    }else{ # loose merge (add label)
      showNotification("Step 4: add TEgroup label (loose merge)...",type="message",duration=5)
      if (missltr){
        temerge4 <- paste0("python3"," ",scriptdir,"extraTEmergep_pall.py"," ",parsedgff_b," ",gapsize," ",">",tobesort_label)
        system(temerge4)
        sortc <- paste0("sort -k1,1 -k4,4n -k5,5n"," ",tobesort_label,">",parsedgff)
        system(sortc)
        rmgff <<- rtracklayer::import.gff(parsedgff)
        dgff(TRUE)
      }else{
        temerge4 <- paste0("python3"," ",scriptdir,"extraTEmergep_pall.py"," ",parsedgff," ",gapsize," ",">",tobesort_label)
        system(temerge4)
        sortc <- paste0("sort -k1,1 -k4,4n -k5,5n"," ",tobesort_label,">",parsedgff_b)
        system(sortc)
        rmgff <<- rtracklayer::import.gff(parsedgff_b)
        dgff(TRUE)
      }
      # true merge
      showNotification("Step 5: Merge gff row according to TEgroup label", type ="message", duration = 3)
      if(missltr){
        temerge5 <- paste0("python3"," ",scriptdir,"extraTEtruemerge.py"," ",parsedgff," ",">",tobesort_merge)
        system(temerge5)
        sortc <- paste0("sort -k1,1 -k4,4n -k5,5n"," ",tobesort_merge,">",parsedgff_b)
        system(sortc)
        mergermgff <<- rtracklayer::import.gff(parsedgff_b)
      }else{
        temerge5 <- paste0("python3"," ",scriptdir,"extraTEtruemerge.py"," ",parsedgff_b," ",">",tobesort_merge )
        system(temerge5)
        sortc <- paste0("sort -k1,1 -k4,4n -k5,5n"," ",tobesort_merge,">",parsedgff_b)
        system(sortc)
        mergermgff <<- rtracklayer::import.gff(parsedgff)
      }
    }

  }
}


# Find ORF ----------------------------------------------------------------


getseq <- function(fa,chromName,rstart,rend,repeatfamily){
  getseq1 <- paste0("python3"," ",scriptdir,"get1seq.py"," ",fa," ",chromName)
  targetseq <- system(getseq1,intern = TRUE)

  targetseq <- paste(targetseq[2:length(targetseq)], collapse = "")
  repeatseq <- substr(targetseq,rstart,rend)
  targetDNAstring <<- Biostrings::DNAStringSet(repeatseq)
  names(targetDNAstring) <<- repeatfamily
}

findorf <- function(fa,chromName,rstart,rend,repeatfamily,orfN){
  getseq1 <- paste0("python3"," ",scriptdir,"get1seq.py"," ",fa," ",chromName)
  targetseq <- system(getseq1,intern = TRUE)

  targetseq <- paste(targetseq[2:length(targetseq)], collapse = "")
  repeatseq <- Biostrings::substr(targetseq,rstart,rend)
  targetDNAstring <- Biostrings::DNAStringSet(repeatseq)
  names(targetDNAstring) <<- repeatfamily

  orfResult <- as.data.frame(systemPipeR::predORF(targetDNAstring, n=orfN, mode="orf", longest_disjoint=TRUE, strand="both"))
  return(orfResult)
}

