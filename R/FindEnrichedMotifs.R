#' Find conserved sequences that are enriched specifically in your promoters of interest vs. control promoters
#'
#' This function takes in an input list of promoters of interest, control promoters, and a custom database of promoters.
#' It splits up promoters of interest into "bits" and BLASTs each bit against all other promoters.
#' Promoters that are found more in your promoters of interest list, than in the control promoters list are considered enriched.
#' Several parameters can be adjusted to increase or decrease precision of BLASTs, and to reduce running time.
#'
#' @param promoters_input Promoters of interest input file
#' @param c_promoters_input Control promoters input file
#' @return A hitsframe containing enrichment of each specific bit
#' @export
FindEnrichedMotifs <- function(Job_name,
                               comp_filter = TRUE,
                               com_threshold = 0.9,
                               bit_size = 16,
                               bit_step = 5,
                               promoter_size = 1000,
                               min_bitratio = 1.5,
                               num_threads = 3,
                               min_match_size = 4,
                               gap_open = 0,
                               gap_extend = 2,
                               penalty_mismatch = -2,
                               match_reward = 1,
                               e_value = 10,
                               promoters_input = 'input.txt',
                               c_promoters_input = 'control.txt') {

  #User inputs:
  #Job_name <- 'top20induced'
  #comp_filter <- TRUE     #Determines if low complexity bits should be removed (FALSE by default)
  #com_threshold <- 0.9    #If true, determines complexity threshold (AT and GC content higher than this will be filtered)
  #bit_size <- 16          #Set the size of each bit to be BLASTed
  #bit_step <- 5           #Set the number of nucleotides to be stepped in generating bits
  #promoter_size <- 1000   #Set the promoter size (max = 2000)
  #min_bitratio <- 1.5     #Sets the minimum ratio of enrichment for hits between promoters of interest and control promoters

  #BLAST parameters
  #num_threads <- 3
  #min_match_size <- 4     #Default: 4
  #gap_open <- 0           #Default: 0
  #gap_extend <- 2         #Default: 2
  #penalty_mismatch <- -2  #Default: -2
  #match_reward <- 1       #Default: 1
  #e_value <- 10           #Default: 10

  promoters_of_interest <- readLines(promoters_input)
  control_promoters_input <- readLines(c_promoters_input)

  #Promoter inputs
  All_Tomato_promoters <- read.csv("All Tomato 2k promoters.csv", sep = ";")
  Selected_promoters_df <- All_Tomato_promoters[All_Tomato_promoters$LocusID %in% promoters_of_interest,]
  Selected_promoters <- substr(Selected_promoters_df$Sequence, 1, promoter_size)
  Selected_promoters_df$Sequence <- Selected_promoters
  Control_promoters_df <- All_Tomato_promoters[All_Tomato_promoters$LocusID %in% control_promoters_input,]
  Control_promoters_cut <- substr(Control_promoters_df$Sequence, 1, promoter_size)
  Control_promoters_df$Sequence <- Control_promoters_cut
  Control_promoters_df <- Control_promoters_df[1:nrow(Selected_promoters_df),]

  #Generate bitstable
  print("Generating bits...")
  iterations <- round((promoter_size - bit_size) / bit_step)
  current_iteration <- NULL
  bitsframe <- data.frame(bitID = character(),
                          bitSeq = character(),
                          stringsAsFactors = F)

  #For each promoter, split into bits, and generate bits ID for each bit
  for (i in 1:nrow(Selected_promoters_df)) {

    step <- 0
    for (j in 1:iterations) {
      current_iteration <- (i-1)*iterations + j
      bitsframe[current_iteration,1] <- paste0('bit_', Selected_promoters_df[i,1], '_', j) #Generate bit identifier
      bitsframe[current_iteration,2] <- substr(Selected_promoters[i], (j+step), (j+step+bit_size-1)) #Extract bit sequence
      step <- step + bit_step
    }

  }

  #Fix bitsframe, remove NNs, empty and incomplete bits
  bitsframe <- bitsframe[rowSums(is.na(bitsframe)) == 0,]
  bitsframe <- subset(bitsframe, nchar(as.character(bitsframe$bitSeq)) >= (bit_size-1))
  bitsframe <- bitsframe[!grepl('N',bitsframe$bitSeq),]
  bitsframe$ATcontent <- (stringr::str_count(bitsframe$bitSeq, 'A') + stringr::str_count(bitsframe$bitSeq, 'T')) / (bit_size + 1)
  bitsframe$GCcontent <- (stringr::str_count(bitsframe$bitSeq, 'G') + stringr::str_count(bitsframe$bitSeq, 'C')) / (bit_size + 1)

  if (comp_filter == TRUE) {
    bitsframe <- subset(bitsframe, bitsframe$ATcontent <= com_threshold)
    bitsframe <- subset(bitsframe, bitsframe$GCcontent <= com_threshold)
  }

  print("Done.")

  #Generate fasta db for reference sequences
  print("Generating BLAST database...")

  db <- character(nrow(Selected_promoters_df)*2)
  db[c(TRUE,FALSE)] <- paste0(">", Selected_promoters_df$LocusID)
  db[c(FALSE,TRUE)] <- Selected_promoters_df$Sequence
  writeLines(db, paste0("C:\\db\\", 'experimental', ".fa"))
  dbc <- character(nrow(Control_promoters_df)*2)
  dbc[c(TRUE,FALSE)] <- paste0(">", Control_promoters_df$LocusID)
  dbc[c(FALSE,TRUE)] <- Control_promoters_df$Sequence
  writeLines(dbc, paste0("C:\\db\\", 'control', ".fa"))

  #Create BLAST databases
  rBLAST::makeblastdb('C:\\db\\experimental.fa', dbtype = "nucl")
  rBLAST::makeblastdb('C:\\db\\control.fa', dbtype = "nucl")
  blex <- rBLAST::blast(db="C:\\db\\experimental.fa")
  blco <- rBLAST::blast(db="C:\\db\\control.fa")

  print("Done.")

  #convert bitsframe to DNA string
  print("Aligning bits...")

  dna <- Biostrings::DNAStringSet(bitsframe$bitSeq)
  names(dna) <- bitsframe$bitID

  #Run BLASTS and store hits
  Hitsframe <- data.frame(bitID = character(),
                          bithits = numeric(),
                          bitconts = numeric(),
                          bitratio = numeric(),
                          stringsAsFactors = F)


  resex <- predict(blex, dna, BLAST_args = paste0("-task \"blastn-short\"",
                                                     " -num_threads ", num_threads,
                                                     " -word_size ", min_match_size,
                                                     " -gapopen ", gap_open,
                                                     " -gapextend ", gap_extend,
                                                     " -penalty ", penalty_mismatch,
                                                     " -evalue ", e_value,
                                                     " -reward ", match_reward))

  resco <- predict(blco, dna, BLAST_args = paste0("-task \"blastn-short\"",
                                                     " -num_threads ", num_threads,
                                                     " -word_size ", min_match_size,
                                                     " -gapopen ", gap_open,
                                                     " -gapextend ", gap_extend,
                                                     " -penalty ", penalty_mismatch,
                                                     " -evalue ", e_value,
                                                     " -reward ", match_reward))
  print("Finished.")

  print("Estimating ratios...")
  Sys.time.start <- Sys.time()
  p100 <- nrow(bitsframe)
  for (z in 1:p100) {

    current_bit <- bitsframe$bitID[z]
    Hitsframe[z,1] <- current_bit

    selected_resex <- resex[resex$QueryID %in% current_bit,]
    selected_resco <- resco[resco$QueryID %in% current_bit,]

    Hitsframe[z,2] <- nrow(selected_resex)
    Hitsframe[z,3] <- nrow(selected_resco)
    Hitsframe[z,4] <- Hitsframe[z,2] / Hitsframe[z,3]

    if (z %% 100 == 0) {
      percentage <- round((z / p100) * 100, digits = 2)
      print(paste0('Finished: ', percentage, '%'))
    }

  }

  Hitsframe$bitSeq <- bitsframe$bitSeq

  #Finishing up
  Sys.time.current <- Sys.time()
  diff.time <- Sys.time.current - Sys.time.start
  cat("Finished in:", diff.time / 60, "minutes.\n")

  return(Hitsframe)

}
