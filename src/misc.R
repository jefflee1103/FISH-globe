str_revcomp <- function(seq){
  # For reverse complementation of character string
  chartr("ATGC", "TACG", seq) %>% stringi::stri_reverse()
}