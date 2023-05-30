#'Align and highlight enriched bits on promoters
#'Be sure to have performed a thorough preselection before using this function
#'@export
ShowEnrichedLocations <- function(Job_name = 'Test-y',
                                  min_bitratio = 1.5,
                                  num_threads = 3,
                                  min_match_size = 4,
                                  gap_open = 0,
                                  gap_extend = 2,
                                  penalty_mismatch = -2,
                                  match_reward = 1,
                                  e_value = 10,
                                  promoter_data = 'All Tomato 2k promoters.csv',
                                  data = Hitsframe,
                                  promoter_size = 1000,
                                  Comp_Filter = TRUE,
                                  filter_strength = 0.8,
                                  return_f = 'txt'){

  #Make directory
  if (dir.exists(paste0('./', Job_name, '/promoterAlignments/')) == FALSE) {
    dir.create(paste0('./', Job_name, '/promoterAlignments/'))
  }

  #Make BLASTable object
  rBLAST::makeblastdb('C:\\db\\experimental.fa', dbtype = "nucl")
  blex <- rBLAST::blast(db="C:\\db\\experimental.fa")

  #Prepare promoters
  All_Tomato_promoters <- read.csv(promoter_data, sep = ";")
  sequence_length <- nchar(All_Tomato_promoters$Sequence)
  start_point <- sequence_length - promoter_size + 1
  end_point <- sequence_length
  Sized_promoters <- substr(All_Tomato_promoters$Sequence, start_point, end_point)
  All_Tomato_promoters$Sequence <- Sized_promoters

  #Prepare bits for sequence extraction
  Hitsframe_l <- data[data$bitratio >= min_bitratio,]
  dna_q <- Biostrings::DNAStringSet(Hitsframe_l$bitSeq)
  names(dna_q) <- Hitsframe_l$bitID

  #Blast all bits
  blast_results <- predict(blex, dna_q, BLAST_args = paste0("-task \"blastn-short\"",
                                                            " -num_threads ", num_threads,
                                                            " -word_size ", min_match_size,
                                                            " -gapopen ", gap_open,
                                                            " -gapextend ", gap_extend,
                                                            " -penalty ", penalty_mismatch,
                                                            " -evalue ", e_value,
                                                            " -reward ", match_reward))

  #Extract unique promoters
  promoters <- unique(blast_results$SubjectID)

  #For every unique promoter:
  for (i in 1:length(promoters)) {
    #Find all bits belonging to promoter
    current_promoter <- promoters[i]
    current_promoter_sequence <- All_Tomato_promoters[All_Tomato_promoters$LocusID %in% current_promoter,'Sequence']
    current_bits_alignment <- blast_results[blast_results$SubjectID == current_promoter,]

    #Add column to current_hits_alignment for promoter and bit matching
    current_bits_alignment$PromMatch <- substr(current_bits_alignment$QueryID, 5, 18)

    #Initialize matrix
    row_number <- nrow(current_bits_alignment)
    col_number <- nchar(current_promoter_sequence) + 3
    alignment_matrix <- matrix(nrow = row_number, ncol = col_number)

    #Assign promoter to alignment matrix
    current_promoter_sequence <- strsplit(current_promoter_sequence, "")[[1]]
    alignment_matrix[1,2] <- paste0(current_promoter, '_promoter')
    for (sequ in 1:length(current_promoter_sequence)) {
      alignment_matrix[1,(2 + sequ)] <- current_promoter_sequence[sequ]
    }

    #For all bits, align sequence
    for (j in 1:nrow(current_bits_alignment)) {
      current_alignment <- current_bits_alignment[j,]
      alignment_row <- 1

      #Bit is own bit
      if (current_alignment$PromMatch == current_promoter && current_alignment$Perc.Ident == 100 && current_alignment$Alignment.Length == nchar(Hitsframe_l[Hitsframe_l$bitID == current_alignment[1,1],'bitSeq'])) {

        alignment_sequence <- strsplit(Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitSeq'], "")[[1]]
        alignment_start <- current_alignment$S.start
        alignment_end <- current_alignment$S.end

        if (alignment_start < 3 || alignment_end > promoter_size) {
          next
        }

        if (Comp_Filter == TRUE) {
          str_alignment <- Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitSeq']
          alignment_count_A <- stringr::str_count(str_alignment, pattern = 'A')
          alignment_count_T <- stringr::str_count(str_alignment, pattern = 'T')
          alignment_count_C <- stringr::str_count(str_alignment, pattern = 'C')
          alignment_count_G <- stringr::str_count(str_alignment, pattern = 'G')

          if ((alignment_count_A + alignment_count_T) / nchar(str_alignment) > filter_strength || (alignment_count_G + alignment_count_C) / nchar(str_alignment) > filter_strength ) {
            next
          }

        }

        #Check if row is empty between S.start and S.end
        alignment_allotment <- alignment_matrix[alignment_row, (alignment_start - 3):alignment_end]
        while (!all(is.na(alignment_allotment))) {
          alignment_row <- alignment_row + 1
          alignment_allotment <- alignment_matrix[alignment_row, (alignment_start - 3):alignment_end]
        }

        #Align sequence to alignment_matrix
        alignment_matrix[alignment_row, alignment_start - 1] <- 'OB'
        alignment_matrix[alignment_row, alignment_start] <- current_alignment$QueryID[1]
        alignment_matrix[alignment_row, alignment_start + 1] <- Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitratio']
        for (seqa in 1:length(alignment_sequence)) {
          alignment_matrix[alignment_row,(alignment_start + seqa + 1)] <- alignment_sequence[seqa]
        }

        #Bit attaching to this promoter
      } else {

        #If alignment is sense
        if (current_alignment$S.start < current_alignment$S.end) {

          alignment_start <- current_alignment$S.start
          alignment_end <- current_alignment$S.end
          query_start <- current_alignment$Q.start
          query_end <- current_alignment$Q.end

          if (alignment_start < 3 || alignment_end > promoter_size) {
            next
          }

          if (Comp_Filter == TRUE) {
            str_alignment <- substr(Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitSeq'], query_start, query_end)
            alignment_count_A <- stringr::str_count(str_alignment, pattern = 'A')
            alignment_count_T <- stringr::str_count(str_alignment, pattern = 'T')
            alignment_count_C <- stringr::str_count(str_alignment, pattern = 'C')
            alignment_count_G <- stringr::str_count(str_alignment, pattern = 'G')

            if ((alignment_count_A + alignment_count_T) / nchar(str_alignment) > filter_strength || (alignment_count_G + alignment_count_C) / nchar(str_alignment) > filter_strength ) {
              next
            }

          }

          alignment_allotment <- alignment_matrix[alignment_row, (alignment_start - 3):alignment_end]
          alignment_sequence <- strsplit(substr(Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID, 'bitSeq'], query_start, query_end), "")[[1]]

          while (!all(is.na(alignment_allotment))) {
            alignment_row <- alignment_row + 1
            alignment_allotment <- alignment_matrix[alignment_row, (alignment_start - 3):alignment_end]
          }

          #Align sequence to alignment matrix
          alignment_matrix[alignment_row, alignment_start] <- current_alignment$QueryID[1]
          alignment_matrix[alignment_row, alignment_start + 1] <- Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitratio']
          for (seqa in 1:length(alignment_sequence)) {
            alignment_matrix[alignment_row,(alignment_start + seqa + 1)] <- alignment_sequence[seqa]
          }


          #If alignment is antisense
        } else {

          alignment_start <- current_alignment$S.end
          alignment_end <- current_alignment$S.start
          query_start <- current_alignment$Q.start
          query_end <- current_alignment$Q.end

          if (alignment_start < 3 || alignment_end > promoter_size) {
            next
          }

          if (Comp_Filter == TRUE) {
            str_alignment <- substr(Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitSeq'], query_start, query_end)
            alignment_count_A <- stringr::str_count(str_alignment, pattern = 'A')
            alignment_count_T <- stringr::str_count(str_alignment, pattern = 'T')
            alignment_count_C <- stringr::str_count(str_alignment, pattern = 'C')
            alignment_count_G <- stringr::str_count(str_alignment, pattern = 'G')

            if ((alignment_count_A + alignment_count_T) / nchar(str_alignment) > filter_strength || (alignment_count_G + alignment_count_C) / nchar(str_alignment) > filter_strength ) {
              next
            }

          }

          alignment_allotment <- alignment_matrix[alignment_row, (alignment_start - 3):alignment_end]
          alignment_sequence <- Biostrings::DNAString(substr(Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID, 'bitSeq'], query_start, query_end))
          alignment_sequence <- strsplit(as.character(Biostrings::reverseComplement(alignment_sequence)), "")[[1]]

          while (!all(is.na(alignment_allotment))) {
            alignment_row <- alignment_row + 1
            alignment_allotment <- alignment_matrix[alignment_row, (alignment_start - 3):alignment_end]
          }

          #Align sequence to alignment matrix
          alignment_matrix[alignment_row, alignment_start] <- current_alignment$QueryID[1]
          alignment_matrix[alignment_row, alignment_start + 1] <- Hitsframe_l[Hitsframe_l$bitID == current_alignment$QueryID,'bitratio']
          for (seqa in 1:length(alignment_sequence)) {
            alignment_matrix[alignment_row, (alignment_start + seqa + 1)] <- alignment_sequence[seqa]
          }

        }
      }
    }

    alignment_matrix <- alignment_matrix[rowSums(is.na(alignment_matrix)) != ncol(alignment_matrix),]

    if (return_f == 'csv') {
      write.csv(alignment_matrix, file = paste0('./', Job_name, '/promoterAlignments/', current_promoter, '.csv'), na = "")
    } else if (return_f == 'txt') {
      alignment_matrix[nchar(alignment_matrix) > 1] <- NA
      write.table(alignment_matrix, file = paste0('./', Job_name, '/promoterAlignments/', current_promoter, '.txt'), na = " ", sep = "", quote = FALSE, col.names = FALSE, row.names = FALSE)
    } else {
      cat("Please provide valid output function \'txt\'or \'csv\'.")
    }

  }

}

