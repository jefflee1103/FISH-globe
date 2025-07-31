str_revcomp <- function(seq){
  # For reverse complementation of character string
  chartr("ATGC", "TACG", seq) %>% stringi::stri_reverse()
}

#' Import and convert a GTF/GFF file into a valr compatible bed tbl format
#'
#' @description This function will output a tibble with the
#' required chrom, start, and end columns, as well as other columns depending
#' on content in GTF/GFF file.
#'
#' @param path path to gtf or gff file
#' @param zero_based if TRUE, convert to zero based
#'
#' @examples
#'
#' gtf <- read_gtf(valr_example("hg19.gencode.gtf.gz"))
#' head(gtf)
#'
#' @importFrom rtracklayer import
#' @export
import_gtf <- function(path, zero_based = TRUE) {
  gtf <- rtracklayer::import(path)
  gtf <- as.data.frame(gtf)
  gtf <- dplyr::mutate_if(gtf, is.factor, as.character)
  res <- dplyr::rename(gtf, chrom = seqnames)

  if (zero_based) {
    res <- dplyr::mutate(res, start = start - 1L)
  }

  tibble::as_tibble(res)
}
