graph_to_tibble <- function(g_input, col) {

  # This function converts an input to a symmetric matrix, then reconverts
  #   it to a long tibble so that both of the following rows exist:
  #
  # X, Y, a
  # Y, X, a

  require(igraph)

  g <- graph_from_data_frame(g_input[, 1:2], directed = FALSE)
  E(g)$weights <- unlist(g_input[col])

  gmatrix <- as_adjacency_matrix(g, attr = "weights", type = "both",
                                 sparse = FALSE) %>%
    as.matrix()

  gtibble <- gmatrix %>%
    as_tibble(rownames = "ROI1") %>%
    pivot_longer(-ROI1, names_to = "ROI2", values_to = "value") %>%
    mutate(
      value = if_else(ROI1 == ROI2, NA, value)
    )

  return(gtibble)

}

summarize_tibble_by_net <- function(tibble, networks,
                                    drop = c("OAN", "vMN", "pMN")) {

  # This function takes the tibble from `graph_to_tibble()`, and assigns each
  #   ROI to a network in `networks`, then summarizes all net-to-net counts.
  #
  # `tibble` contains a logical column `sig` with predetermined significance
  #   values.

  long <- tibble %>%
    mutate(
      net1 = networks$net_short[match(ROI1, networks$roi)],
      net2 = networks$net_short[match(ROI2, networks$roi)],
    )

  result <- long %>%
    group_by(net1, net2) %>%
    summarize(
      n = n(),
      sig = sum(sig, na.rm = TRUE),
      .groups = "keep"
    ) %>%
    mutate(
      prop = sig / n
    )

  if (!all(is.na(drop))) {

    # Drop any networks in the drop list, if provided

    result <- result %>%
      filter(
        !(net1 %in% drop),
        !(net2 %in% drop)
      )

  }

  return(result)

}

# Bootstrap ====

shuffle_matrix <- function(mat, tibble = TRUE) {

  # Give a matrix of pvalues (square, symmetric, logical), return a shuffled
  #   symmetric matrix with the same number of significant edges.
  # If `tibble = TRUE`, reformat like the long format from `graph_to_tibble()`

  total_n <- sum(mat, na.rm = TRUE) / 2

  # Flatten the matrix, assigning the upper triangle to be NA
  mat_lower <- mat
  mat_lower[upper.tri(mat_lower, diag = TRUE)] <- NA
  v <- as.vector(mat_lower)

  # Replace all non-NA values with FALSE
  v2 <- rep(NA, length(v))
  v2[!is.na(v)] <- FALSE

  # Replaces some of those FALSE values with TRUEs at random
  new_random_lower <- sample(which(!is.na(v)), size = total_n)
  v2[new_random_lower] <- TRUE

  # Turn that new matrix into a matrix, and add it as its own upper triangle,
  #   creating the symmetry
  mat_final <- matrix(v2, nrow = nrow(mat_lower), ncol = ncol(mat_lower))
  mat_final[upper.tri(mat_final)] <- t(mat_final)[upper.tri(mat_final)]
  mat_final[diag(mat_final)] <- FALSE

  # Replace names from original matrix
  colnames(mat_final) <- colnames(mat)
  rownames(mat_final) <- rownames(mat)

  if (tibble) {

    # Pivot longer, if requested

    mat_final <- mat_final %>%
      as_tibble(rownames = "ROI1") %>%
      pivot_longer(-ROI1, names_to = "ROI2", values_to = "sig")

  }

  return(mat_final)

}

bootstrap_pmatrix <- function(pvalues, nets, n = 10000,
                              drop = c("OAN", "vMN", "pMN")) {

  require(infer)

  # Given a long-matrix-style list of pvalues, calculate the empirical p-value
  #   for the number of nodes in each network pair

  est_time <- round(4 / 10000 * n, 1)

  message(paste("Total estimated time:", est_time, "minutes"))

  # Calculate the ground truth
  pvalues_summary <- summarize_tibble_by_net(pvalues, nets, drop = drop)

  pmatrix <- pvalues %>%
    select(ROI1, ROI2, sig) %>%
    pivot_wider(names_from = ROI2, values_from = sig) %>%
    column_to_rownames(var = "ROI1") %>%
    as.matrix()

  message("Creating shuffled matrices (can take a minute or two) ...")

  # Create n shuffled matrices
  bootstrap1 <- tibble(rep = 1:n) %>%
    mutate(
      mat = map(rep, ~shuffle_matrix(pmatrix), .progress = TRUE),
    )

  message("Summarizing reps ... ")

  # Summarize each shuffled matrix by network and join true values
  bootstrap2 <- bootstrap1 %>%
    mutate(
      result = map(mat, ~summarize_tibble_by_net(.x, nets_ordered, drop = drop),
                   .progress = TRUE)
    ) %>%
    select(-mat) %>%
    unnest(result) %>%
    group_by(net1, net2, n) %>%
    nest() %>%
    left_join(
      pvalues_summary,
      by = join_by(net1, net2, n)
    )

  message("Calculating p values ...")

  # Based on the bootstrapped null empirical distribution, estimate the p
  #   value for the actual number of significant connections
  bootstrap2p <- bootstrap2 %>%
    mutate(

      pvalue = map2_dbl(data, sig,
                        ~get_p_value(tibble(sig = .x$sig), .y,
                                     direction = "greater")[[1]],
                        .progress = TRUE)

    )

  # Return the p-values without the massive data frames
  bootstrap2p %>%
    select(-data) %>%
    return()

}

remove_tri <- function(x, upper = FALSE) {

  # Many tibbles I'm working with have the rows:
  #   X, Y, a,
  #   Y, X, a
  # This deletes the second row, mostly for visualization purposes

  id_cols <- x[, 1:2] %>%
    mutate(
      x = TRUE
    )

  wide <- id_cols %>%
    pivot_wider(names_from = colnames(id_cols)[2], values_from = x) %>%
    column_to_rownames(colnames(.)[1])

  if (upper)
    wide[upper.tri(wide)] <- FALSE
  else
    wide[lower.tri(wide)] <- FALSE

  rows_to_drop <- wide %>%
    rownames_to_column(colnames(id_cols)[1]) %>%
    pivot_longer(-colnames(id_cols)[1], names_to = colnames(id_cols)[2],
                 values_to = "keep")

  result <- left_join(x, rows_to_drop) %>%
    filter(
      keep
    ) %>%
    select(-keep)

  return(result)

}
