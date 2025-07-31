# Load necessary libraries
library(dplyr)
library(purrr)
library(tibble)
library(valr)

# # --------------------------------------------------------------------------- #
# # Main Function to Create Sliding Windows Across Exons
# # --------------------------------------------------------------------------- #

# #' Create sliding windows across multiple exons as if they were concatenated.
# #'
# #' @param exons A data frame with exon coordinates (chrom, start, end).
# #' @param window_width The width of each sliding window (numeric).
# #' @param window_step The step size for sliding the window (numeric).
# #' @return A tibble with the genomic coordinates for each window. For windows
# #'   that span exon junctions, it returns a row for each exon segment. It also
# #'   includes columns indicating if a window spans a junction and the full
# #'   genomic start/end of the window.
# make_slidingwindow <- function(exons, window_width, window_step) {

#   # Prepare exon data with 0-based cumulative coordinates
#   # This creates a "concatenated" coordinate system.
#   exons_processed <- exons %>%
#     arrange(chrom, start) %>% # Ensure correct order before creating cumulative coords
#     mutate(
#       length = end - start
#     ) %>%
#     mutate(
#       cumulative_end = cumsum(length),
#       cumulative_start = cumulative_end - length
#     )

#   # Calculate the total length of the concatenated sequence
#   total_length <- sum(exons_processed$length)

#   # 5. Generate window start positions on the 0-based concatenated sequence
#   # For a 0-based system, the last possible start is total_length - window_width
#   window_starts <- seq(from = 0, to = total_length - window_width, by = window_step)

#   # 6. Create a tibble for windows in the cumulative coordinate system
#   windows_cumulative <- tibble(
#     window_id = seq_along(window_starts),
#     win_start_cum = window_starts,
#     win_end_cum = window_starts + window_width
#   )

#   # 7. Define a nested helper function to map a window to its genomic segments
#   map_window_to_genome <- function(win_id, win_start_c, win_end_c, exon_data) {

#     # Find the exon(s) that this window overlaps with
#     # An interval [A,B) overlaps [C,D) if A < D and C < B
#     overlapping_exons <- exon_data %>%
#       filter(
#         win_start_c < cumulative_end & cumulative_start < win_end_c
#       )

#     # Check if the window spans a junction (i.e., touches more than one exon)
#     spans_junction <- nrow(overlapping_exons) > 1

#     # Calculate the overall genomic start and end for the entire window
#     first_exon <- head(overlapping_exons, 1)
#     last_exon <- tail(overlapping_exons, 1)
#     full_genomic_start <- first_exon$start + (win_start_c - first_exon$cumulative_start)
#     full_genomic_end <- last_exon$start + (win_end_c - last_exon$cumulative_start)

#     # For each overlapping exon, calculate the genomic coordinates of the segment
#     result <- overlapping_exons %>%
#       dplyr::mutate(
#         # Calculate the genomic start for this segment of the window
#         segment_start = start + pmax(0, win_start_c - cumulative_start),
#         # Calculate the genomic end for this segment of the window
#         segment_end = start + pmin(length, win_end_c - cumulative_start),
#         # Add the other required columns
#         window_id = win_id,
#         spans_junction = spans_junction,
#         genomic_start = full_genomic_start,
#         genomic_end = full_genomic_end
#       ) %>%
#       mutate(
#         unique_window_id = paste(chrom, genomic_start, genomic_end, strand, collapse = ":")
#       ) %>%
#       # Select and reorder the final columns
#       dplyr::select(
#         unique_window_id, chrom, start = segment_start, end = segment_end, strand,
#         spans_junction, genomic_start, genomic_end
#       )

#     return(result)
#   }


#   # 8. Apply the mapping function to every window.
#   # pmap_dfr iterates over each row of `windows_cumulative` and row-binds the results.
#   final_windows <- pmap_dfr(
#     .l = list(
#       win_id = windows_cumulative$window_id,
#       win_start_c = windows_cumulative$win_start_cum,
#       win_end_c = windows_cumulative$win_end_cum
#     ),
#     .f = map_window_to_genome,
#     exon_data = exons_processed
#   )

#   return(final_windows)
# }


# --------------------------------------------------------------------------- #
# Optimized Fast Function
# --------------------------------------------------------------------------- #

#' (Fast Version) Create sliding windows across exons using a vectorized approach.
#'
#' @param exons A data frame with exon coordinates (chrom, start, end, strand).
#' @param window_width The width of each sliding window (numeric).
#' @param window_step The step size for sliding the window (numeric).
#' @return A tibble with the genomic coordinates for each window segment.
#'
#' @note This function uses `valr::bed_intersect` for a fast, single-pass
#' interval intersection, which is significantly more performant.
make_slidingwindow <- function(exons, window_width, window_step) {

  # 1. Prepare exon data with cumulative coordinates
  exons_processed <- exons %>%
    arrange(start) %>%
    mutate(
      length = end - start,
      cumulative_end = cumsum(length),
      cumulative_start = cumulative_end - length
    )

  # 2. Generate all windows in the "cumulative" coordinate space
  total_length <- sum(exons_processed$length)
  
  # Return NULL if the total exon length is smaller than a window, preventing errors.
  if (total_length < window_width) {
    return(NULL)
  }
  
  window_starts <- seq(from = 0, to = total_length - window_width, by = window_step)
  
  # Create a valr-compatible tibble for cumulative windows
  windows_cumulative <- tibble(
    chrom = "cumulative", # Use a dummy chromosome for the intersection
    start = window_starts,
    end = window_starts + window_width,
    window_id = seq_along(window_starts)
  )

  # 3. Create a valr-compatible tibble for cumulative exons
  exons_cumulative <- exons_processed %>%
    mutate(chrom_for_join = "cumulative") %>% # Create the dummy chromosome name
    select(
      chrom = chrom_for_join, # Use dummy chrom for the intersection
      start = cumulative_start,
      end = cumulative_end,
      orig_chrom = chrom, # Keep original chrom name
      orig_start = start,
      orig_end = end,
      orig_strand = strand
    )

  # 4. Perform a single, fast interval intersection to find all overlaps
  # This is the key performance improvement.
  joined_data <- bed_intersect(
    x = windows_cumulative,
    y = exons_cumulative
  ) %>%
    filter(.overlap > 0)

  # 5. Vectorized calculation of final coordinates
  result <- joined_data %>%
    # Calculate segment coordinates based on the overlap
    mutate(
      segment_start = orig_start.y + pmax(0, start.x - start.y),
      segment_end = orig_start.y + pmin(end.y - start.y, end.x - start.y)
    ) %>%
    # Group by each window to calculate full coordinates and junction status
    group_by(window_id.x) %>%
    mutate(
      spans_junction = n() > 1,
      # Get the first exon's start for the full window start calculation
      first_exon_start = first(orig_start.y),
      first_exon_cum_start = first(start.y),
      # Get the last exon's start for the full window end calculation
      last_exon_start = last(orig_start.y),
      last_exon_cum_start = last(start.y)
    ) %>%
    ungroup() %>%
    # Calculate the full start/end for the entire window
    mutate(
      genomic_start = first_exon_start + (start.x - first_exon_cum_start),
      genomic_end = last_exon_start + (end.x - last_exon_cum_start)
    ) %>%
    unite(col = "unique_window_id", c("orig_chrom.y", "genomic_start", "genomic_end", "orig_strand.y"), sep = ":", remove = FALSE) %>%
    # Select and rename final columns
    select(
      unique_window_id,
      chrom = orig_chrom.y,
      start = segment_start,
      end = segment_end,
      strand = orig_strand.y,
      spans_junction,
      genomic_start,
      genomic_end
    )

  return(result)
}

# # --------------------------------------------------------------------------- #
# # Example Usage and Verification
# # --------------------------------------------------------------------------- #

# # 1. Your example exon data
# exons_df <- tibble::tribble(
#   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#   "chr1", 100,    200,   "A",    ".",     "+", # length = 100
#   "chr1", 500,    550,   "B",    ".",     "+"  # length = 50
# )

# # 2. Define window parameters
# WINDOW_WIDTH <- 52
# WINDOW_STEP <- 1

# # 3. Call the function to generate the windows
# final_windows <- make_slidingwindow(exons = exons_df,
#                                     window_width = WINDOW_WIDTH,
#                                     window_step = WINDOW_STEP)


# # 4. Display the results
# print("--- Final Genomic Windows (Segments) ---")
# # View the first few windows, which are entirely within the first exon
# print(head(final_windows, 5))

# print("--- Windows Spanning the Exon A -> Exon B Junction ---")
# # Exon A has length 100. Junction-spanning windows will have a cumulative start
# # before 100 and a cumulative end after 100.
# # The first window to span the junction starts at cumulative pos 48 (window_id 49).
# print(filter(final_windows, window_id > 48 & window_id < 52))

# # 5. Verification
# # Verify a non-spanning window has the correct width and coordinates
# print("--- Verification: A non-spanning window ---")
# non_spanning_window <- filter(final_windows, window_id == 10)
# print(non_spanning_window)
# print(paste("Segment width (end - start):", non_spanning_window$end - non_spanning_window$start))
# print(paste("Full window width (genomic_end - genomic_start):", non_spanning_window$genomic_end - non_spanning_window$genomic_start))


# # Verify a spanning window has a TRUE flag and correct total width from its segments
# print("--- Verification: A spanning window ---")
# spanning_window_segments <- filter(final_windows, window_id == 50)
# print(spanning_window_segments)
# total_segment_length <- sum(spanning_window_segments$end - spanning_window_segments$start)
# print(paste("Sum of segment lengths:", total_segment_length))
# print(paste("Genomic distance (genomic_end - genomic_start):", spanning_window_segments$genomic_end[1] - spanning_window_segments$genomic_start[1]))
