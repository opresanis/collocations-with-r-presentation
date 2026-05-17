collocates <-
function(corpus, node, left = 5, right =5,
                       to_lower = FALSE, corrected = FALSE,
                       measures = "all") {

  # Measures implemented based on Brezina (2018), p.72
  # measures <- c("MU", "MI", "MI2", "MI3", "LL_short", "LL_long", "Z-score", "T-score",
  #               "Dice", "LogDice", "LogRatio", "MS", "DeltaP", "all" )


  ##################################
  ## ADD EVENT HANDLER
  ##################################

  ###########
  # Step 1  #
  ###########

  # optional: all tokens to lower case
  if (to_lower) {
    corpus <- lapply(corpus, tolower)
    node <- tolower(node)

  }

  # Flatten token lists into one long vector
  all_tokens <- unlist(corpus, use.names = FALSE)



  ###########
  # Step 2  #
  ###########
  # Observed frequencies

  # Calculate word frequencies
  total_freq_table <- table(all_tokens)

  # Total number of tokens in the corpus ('N')
  # cf. Evert (2009), Brezina (2018) p. 70
  N <- sum(total_freq_table)

  # Calculate 'R1' without window size correction
  R1_unadj <- total_freq_table[[node]]
  # shouldn't be the case
  if (is.null(R1_unadj)) R1_unadj <- 0

  # Extract all collocates in window around each node
  collocates <- character()
  for (tokens in corpus) {

    # optional: all tokens to lower case
    if (to_lower) {
      tokens <- tolower(tokens)
    }

    # get all the node word positions
    positions <- which(tokens == node)

    for (pos in positions) {
      left_indices <- max(1, pos - left):(pos - 1)
      right_indices <- (pos + 1):min(length(tokens), pos + right)
      # Note: it may be the case that an individual token may be
      # counted more than once as a collocate - if the windows between
      # two notes overlap
      collocates <- c(collocates, tokens[left_indices], tokens[right_indices])

      # the below alternative code would avoid double counting
      # but since we are interested in collocation, overlapping counting might
      # be preferred

      # all_indices <- unique(c(left_indices, right_indices))
      # collocates <- c(collocates, tokens[all_indices])
    }

  }

  # Result: list of all collocates in the defined window

  # Count collocates
  # This produces the 'O11' values
  colloc_freq_table <- table(collocates)

  # onwards, we work with a data frame
  colloc_df <- data.frame(
    token = names(colloc_freq_table),
    O11 = as.integer(colloc_freq_table),
    stringsAsFactors = FALSE
  )

  # Add total frequency of the collocates from full corpus
  # This are the 'C1' values
  colloc_df$C1 <- as.integer(total_freq_table[colloc_df$token])

  # Replace NA with 0 in total_freq if token not found (should not happen)
  colloc_df$C1[is.na(colloc_df$C1)] <- 0

  # calculate the other values of the contingency table
  # Note: we need to use pmax so that the O21 values
  # can't become negative in case of O11>C1
  colloc_df$O21 <- pmax((colloc_df$C1 - colloc_df$O11), 0)

  # allow for corrected and uncorrected use
  if (corrected) {
    # Using R1 as a variable rather than a column in the data frame,
    # is computationally more efficient when calculating the measures
    R1 <- R1_unadj * (left + right)
    colloc_df$R1 <-  R1
    colloc_df$O12 <- R1 - colloc_df$O11
  } else {
    # R1 can stay as it is
    # maybe we exclude this column later
    R1 <- R1_unadj
    colloc_df$R1 <-  R1_unadj
    colloc_df$O12 <- R1_unadj - colloc_df$O11
  }

  # Calculate 'R2'
  colloc_df$R2 <- N - colloc_df$R1
  R2 <- N - R1

  # calculate the remaining observed frequencies
  # and margins
  colloc_df$C2 <- N - colloc_df$C1
  colloc_df$O22 <- colloc_df$C2 - colloc_df$O12

  # rearrange a bit
  colloc_df <- colloc_df[, c("token", "O11", "O21", "C1", "O12", "O22", "C2", "R1", "R2")]

  # change all int variables to double precision num, to avoid any issues with
  # integer overflow in the calculations of the measures
  colloc_df[] <- lapply(colloc_df, function(x) if (is.integer(x)) as.numeric(x) else x)

  ###########
  # Step 3  #
  ###########
  # calculate expected values
  colloc_df$E11 <- colloc_df$R1 * colloc_df$C1 / N
  colloc_df$E12 <- colloc_df$R1 * colloc_df$C2 / N
  colloc_df$E21 <- R2 * colloc_df$C1 / N
  colloc_df$E22 <- R2 * colloc_df$C2 / N

  ###########
  # Step 4  #
  ###########
  # Calculate measures

  # Measure 2: MU (Evert 2009)
  if ("MU" %in% measures | "all" %in% measures ) {
    colloc_df$MU <- colloc_df$O11/colloc_df$E11
  }

  # Measure 3: MI using formula in Brezina (2018)
  # Formula is also in Evert (2002)
  if ("MI" %in% measures | "all" %in% measures ) {
    colloc_df$MI <- log2(
      colloc_df$O11/colloc_df$E11
    )
  }

  # Measure 4: MI2
  if ("MI2" %in% measures | "all" %in% measures ) {
    colloc_df$MI2 <- log2(
      colloc_df$O11^2/colloc_df$E11
    )
  }

  # Measure 5: MI3
  if ("MI3" %in% measures | "all" %in% measures ) {
    colloc_df$MI3 <- log2(
      colloc_df$O11^3/colloc_df$E11
    )
  }

  # Measure 6a: Log Likelihood short

  # Note: for a few cases, with very few collocation counts, we
  # may get O21 as being 0 (or even negative; but this has been handled earlier in the code)
  # for words which are exclusively appear in the context of the node
  # it feels 'wrong' to correct for 'Inf/-Inf' log values
  if ("LL_short" %in% measures | "all" %in% measures ) {
    colloc_df$LL_short <- 2 * ((colloc_df$O11 * log(colloc_df$O11 / colloc_df$E11)) +
                                 (colloc_df$O21 * log(colloc_df$O21 / colloc_df$E21)))
  }

  # Measure 6b: Log Likelihood long
  if ("LL_long" %in% measures | "all" %in% measures ) {
    colloc_df$LL_long <-
      2 * (colloc_df$O11 * log(colloc_df$O11 / colloc_df$E11) +
             colloc_df$O21 * log(colloc_df$O21 / colloc_df$E21) +
             colloc_df$O12 * log(colloc_df$O12 / colloc_df$E12) +
             colloc_df$O22 * log(colloc_df$O22 / colloc_df$E22)
      )
  }

  # Measure 7: Z-score
  if ("Z-score" %in% measures | "all" %in% measures ) {
    colloc_df$Z_score <- (colloc_df$O11 - colloc_df$E11) / sqrt(colloc_df$E11)
  }

  # Measure 8: T-score
  if ("T-score" %in% measures | "all" %in% measures ) {
    colloc_df$T_score <- (colloc_df$O11 - colloc_df$E11) / sqrt(colloc_df$O11)
  }

  # Measure 9: Dice
  if ("Dice" %in% measures | "all" %in% measures ) {
    colloc_df$Dice <- (2 * colloc_df$O11) / (R1 + colloc_df$C1)
  }

  # Measure 10: Log Dice
  if ("LogDice" %in% measures | "all" %in% measures ) {
    colloc_df$LogDice <- 14 + log2((2 * colloc_df$O11) / (R1 + colloc_df$C1))
  }

  # Measure 11: Log ratio
  if ("LogRatio" %in% measures | "all" %in% measures ) {
    colloc_df$LogRatio <- log2((colloc_df$O11 * R2) / (colloc_df$O21 * R1))
  }

  # Measure 12: Minimum sensitivity
  if ("MS" %in% measures | "all" %in% measures ) {
    colloc_df$MS <- pmin(colloc_df$O11 / colloc_df$C1, colloc_df$O11 / R1)
  }

  # Measure 13: Delta P
  if ("DeltaP" %in% measures | "all" %in% measures ) {
    colloc_df$Delta_p1 <- colloc_df$O11 / R1 - colloc_df$O21 / R2
    colloc_df$Delta_p2 <- colloc_df$O11 / colloc_df$C1 - colloc_df$O12 / colloc_df$C2
  }

  # Measure 14: Cohen's D
  # Will be calculated using a dedicated function


  # Prettify the returned df:
  # rename the first three columns
  # to make them similar to the wording in LancsBox X
  colnames(colloc_df)[colnames(colloc_df) == "token"] <- "collocate"
  colnames(colloc_df)[colnames(colloc_df) == "O11"] <- "freq_collocation"
  colnames(colloc_df)[colnames(colloc_df) == "C1"] <- 'freq_corpus'

  # Reduce output df to not include (most) observed frequencies
  # and expected frequencies. We keep only the frequencies and association measures
  colloc_df <- subset(colloc_df, select = -c(O21, O12, O22, C2, R1, R2,
                                             E11, E12, E21, E22))


  

  return(colloc_df)
}
