#' This function returns PWMs and seqlogos of enriched queries.
#' @param input A hitsframe containing enriched queries
#' @return A results dataframe containing PWMs and plots for each input query
#' @export
PWMFromHits <- function(min_bitratio = 1.5,
                        num_threads = 3,
                        min_match_size = 4,
                        gap_open = 0,
                        gap_extend = 2,
                        penalty_mismatch = -2,
                        match_reward = 1,
                        e_value = 10,
                        promoter_data = 'All Tomato 2k promoters.csv',
                        input = input){

  #User inputs
  #min_bitratio <- 1.5
  #num_threads <- 3
  #min_match_size <- 4
  #gap_open <- 0
  #gap_extend <- 2
  #penalty_mismatch <- -2
  #match_reward <- 1
  #e_value <- 10
  #promoter_data <- 'All Tomato 2k promoters.csv'

  All_Tomato_promoters <- read.csv(promoter_data, sep = ";")

  #Make BLAST objects
  rBLAST::makeblastdb('C:\\db\\experimental.fa', dbtype = "nucl")
  rBLAST::makeblastdb('C:\\db\\control.fa', dbtype = "nucl")
  blex <- rBLAST::blast(db="C:\\db\\experimental.fa")
  blco <- rBLAST::blast(db="C:\\db\\control.fa")

  #Select enriched bits for further analysis
  Hitsframe <- input
  Hitsframe <- subset(Hitsframe, Hitsframe$bitratio >= min_bitratio)
  Hitsframe <- Hitsframe[order(Hitsframe$bitratio),]
  dna_q <- Biostrings::DNAStringSet(Hitsframe$bitSeq)
  names(dna_q) <- Hitsframe$bitID

  #BLAST all queries
  blast_results <- predict(blex, dna_q, BLAST_args = paste0("-task \"blastn-short\"",
                                                            " -num_threads ", num_threads,
                                                            " -word_size ", min_match_size,
                                                            " -gapopen ", gap_open,
                                                            " -gapextend ", gap_extend,
                                                            " -penalty ", penalty_mismatch,
                                                            " -evalue ", e_value,
                                                            " -reward ", match_reward))

  #Extract unique queries
  query_ids <- unique(blast_results$QueryID)

  #Initialize results dataframe
  results <- data.frame(QueryID = character(), PWM = list(), PWMLogo = list(), stringsAsFactors = FALSE)
  for (i in 1:length(query_ids)) {

    #Extract BLAST result for the current ID
    query_results <- blast_results[blast_results$QueryID == query_ids[i],]

    #Initialize empty matrix to hold aligned sequences
    max_length <- nchar(dna_q[i])
    num_queries <- nrow(query_results)
    aligned_seqs <- matrix("", nrow = num_queries, ncol = max_length)

    # Loop over each query ID and extract the aligned sequences
    for (j in 1:nrow(query_results)) {
      # Extract the blast results for the current query ID
      seq_id <- query_results[j,'SubjectID']
      seq <- All_Tomato_promoters[All_Tomato_promoters$LocusID %in% seq_id,'Sequence']
      s_start <- query_results[j, 'S.start']
      s_end <- query_results[j, 'S.end']
      q_start <- query_results[j, 'Q.start']
      q_end <- query_results[j, 'Q.end']
      q_size <- nchar(dna_q[i])

      if (s_start < s_end) {
        #If the alignment is in sense: Calculate subject start and end
        alignment_start <- s_start - (q_start - 1)
        alignment_end <- alignment_start + q_size - 1
        aligned_seq <- substr(seq, alignment_start, alignment_end)

        if (nchar(aligned_seq) < q_size) {
          next
        }

      } else {
        alignment_start <- s_end - (q_size-q_end)
        alignment_end <- alignment_start + q_size - 1
        aligned_seq <- Biostrings::reverseComplement(DNAString(substr(seq, alignment_start, alignment_end)))
        aligned_seq <- as.character(aligned_seq)

        if (nchar(aligned_seq) < q_size) {
          next
        }

      }

      aligned_seqs[j,] <- strsplit(aligned_seq, "")[[1]]

    }

    #Add pwm and pwmlogo to results dataframe
    pwm <- Biostrings::consensusMatrix(aligned_seqs, as.prob = FALSE)
    pwm <- as.matrix(pwm)
    pwmlogo <- ggplot2::ggplot() + ggseqlogo::geom_logo(pwm, method = 'bits') + ggseqlogo::theme_logo()
    results[i, 'QueryID'] <- query_ids[i]
    results$PWM[[i]] <- list(pwm)
    results$PWMLogo[i] <- list(pwmlogo)

    percentage <- (i / length(query_ids)) * 100

    if (round(percentage) %% 10 == 0) {
      cat("Finished: ", percentage, '%\n')
    }

  }

  results$bithits <- Hitsframe$bithits
  results$bitconts <- Hitsframe$bitconts
  results$bitratio <- Hitsframe$bitratio
  results$bitSeq <- Hitsframe$bitSeq

  return(results)
}
