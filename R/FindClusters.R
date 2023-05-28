#' FindClusters
#'
#' @description
#' This function performs hierarchical clustering on sequences based on sequence alignment scores.
#'
#' @param df A data frame containing columns 'bitSeq' with sequences to be clustered and 'bitID' with sequence names.
#'
#' @return
#' An object of class 'hclust' representing the hierarchical clustering of sequences.
#'
#' @examples
#' # clustering <- FindClusters(Hitsframe)
#'
#' @export
FindClusters <- function(df) {
  # Check if the required columns 'bitSeq' and 'bitID' exist in the data frame
  if (!all(c("bitSeq", "bitID") %in% names(df))) {
    stop("The input data frame must contain columns named 'bitSeq' and 'bitID'.")
  }
  
  # Extract bit sequences and IDs from the input data frame
  bitSeqs <- DNAStringSet(df$bitSeq)
  bitIDs <- df$bitID
  
  # Initialize a similarity matrix with row and column names
  similarity_matrix <- matrix(nrow = length(bitSeqs), ncol = length(bitSeqs))
  rownames(similarity_matrix) <- bitIDs
  colnames(similarity_matrix) <- bitIDs
  
  iterations <- length(bitSeqs)

  for(i in 1:iterations) {
    for(j in i:length(bitSeqs)) {
      alignment <- Biostrings::pairwiseAlignment(pattern = bitSeqs[i], subject = bitSeqs[j], substitutionMatrix=NULL)
      score <- Biostrings::score(alignment)
      similarity_matrix[i, j] <- score
      similarity_matrix[j, i] <- score
    }
    cat(i, '/', iterations, ' iterations.\n')
  }
  
  normalized_similarity_matrix <- (similarity_matrix - min(similarity_matrix)) / (max(similarity_matrix) - min(similarity_matrix))
  distance_matrix <- as.dist(1 - normalized_similarity_matrix)

  # Perform hierarchical clustering
  clustering <- stats::hclust(distance_matrix, method = "average")
  
  return(clustering)
}
