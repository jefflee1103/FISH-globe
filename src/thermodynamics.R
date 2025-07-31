###############################
## Functions for calculating thermodynamics parameters for target RNA sequence
## and attaching these into the dataframe
## 
## Functions are from my HCRv3 probe designer 
###############################

sourceCpp("./src/thermo_calc.cpp")
# sourceCpp("./src/find_optimal_probe_set_dp_cpp.cpp")

get_thermodynamic_parameters <- function(candidate_probes, temperature, Na, oligo_conc) {
  
  # Using a standard lapply loop.
  all_results <- lapply(candidate_probes$target_sequence, function(target_seq) {
    
    # Calculate thermodynamics for the full probe and its two halves
    full_thermo <- calculate_thermodynamics_cpp(target_seq, temperature, Na, oligo_conc)
    
    first_half_seq <- substr(target_seq, 1, nchar(target_seq) / 2 - 1)
    first_half_thermo <- calculate_thermodynamics_cpp(first_half_seq, temperature, Na, oligo_conc)
    
    second_half_seq <- substr(target_seq, nchar(target_seq) / 2 + 2, nchar(target_seq))
    second_half_thermo <- calculate_thermodynamics_cpp(second_half_seq, temperature, Na, oligo_conc)
    
    # Return a list containing all results for the probe
    list(
      Full = full_thermo,
      First = first_half_thermo,
      Second = second_half_thermo
    )
  })
  
  # Restructure the list of results into the format expected by the next function
  output <- list(
    Full = purrr::map(all_results, "Full"),
    First = purrr::map(all_results, "First"),
    Second = purrr::map(all_results, "Second")
  )
  
  return(output)
}


passed_a_comp <- function(theProbeSeq){

  A_composition <- str_count(theProbeSeq, pattern = "A") / nchar(theProbeSeq)
  theVerdict <- if_else(A_composition < 0.28, TRUE, FALSE)
  return(theVerdict)
}

passed_a_stack <- function(theProbeSeq){

  theVerdict <- if_else(
    str_detect(theProbeSeq, "AAAA"), FALSE, TRUE
  )
  return(theVerdict)
}

passed_c_comp <- function(theProbeSeq){

  C_composition <- str_count(theProbeSeq, pattern = "C") / nchar(theProbeSeq)
  theVerdict <- if_else(
    C_composition > 0.22 & C_composition < 0.28, TRUE, FALSE)
  return(theVerdict)
}

passed_c_stack <- function(theProbeSeq){

  theVerdict <- if_else(
    str_detect(theProbeSeq, "CCCC"), FALSE, TRUE
  )
  return(theVerdict)
}

passed_c_spec_stack <- function(theProbeSeq){
  max_c_comp <- map_dbl(1:7, ~ {
    substring <- str_sub(theProbeSeq, start = .x, end = .x + 5)
    C_composition <- str_count(substring, pattern = "C") / nchar(substring)
  }) %>% max()
  theVerdict <- if_else(max_c_comp <= 0.5, TRUE, FALSE)
  return(theVerdict)
} 


annotate_probes <- function(candidate_probes, thermodynamics){
  # Takes candidate_probes df and adds thermodynamic and nucleotide composition parameters
  
  candidate_probes %>%
    mutate(rev_comp = str_revcomp(target_sequence)) %>%
    mutate(
      GC_content = str_count(target_sequence, pattern = "G|C") / nchar(target_sequence),
      A_content  = str_count(target_sequence, pattern = "A") / nchar(target_sequence),
      C_content  = str_count(target_sequence, pattern = "C") / nchar(target_sequence)
    ) %>%
    bind_cols(dG = map_dbl(thermodynamics$Full, function(x) x[["dG"]])) %>%
    bind_cols(Tm = map_dbl(thermodynamics$Full, function(x) x[["Tm"]])) %>%
    bind_cols(dG_1st_half = map_dbl(thermodynamics$First, function(x) x[["dG"]])) %>%
    bind_cols(Tm_1st_half = map_dbl(thermodynamics$First, function(x) x[["Tm"]])) %>%
    bind_cols(dG_2nd_half = map_dbl(thermodynamics$Second, function(x) x[["dG"]])) %>%
    bind_cols(Tm_2nd_half = map_dbl(thermodynamics$Second, function(x) x[["Tm"]])) %>%
    mutate(
      passed_a_comp       = passed_a_comp(rev_comp),
      passed_c_comp       = passed_c_comp(rev_comp),
      passed_a_stack      = passed_a_stack(rev_comp),
      passed_c_stack      = passed_c_stack(rev_comp),
      passed_c_spec_stack = map_lgl(rev_comp, passed_c_spec_stack)
    ) %>%
    return()
}