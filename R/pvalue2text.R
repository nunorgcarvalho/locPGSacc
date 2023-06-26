#' @title pvalue2text
#' @description Convert an inputted p-value into a string to be formatted by [plotmath()]. Uses base10 notation
#' @inherit locPGSacc author
#' 
#' @param p-value numerical p-value to be converted
#' @param round_digits Integer of total digits to round p-value to
#' @param alpha level of significance used so that close values do not get mistaken for being on the border. E.g. p=0.0499999 will never be displayed as p=0.05
#' 
#' @returns Returns a string to be parsed through [plotmath()] to show up as base10 scientific notation (e.g. 1.23 x 10^-4).
#' 
#' @export

pvalue2text <- function(p_value, round_digits=3, alpha=0.05) {
  # extends shown digits in cases where rounded p == 0.05
  p_text <- round(p_value,round_digits)
  while (p_text==alpha | !exists("p_text") | p_value != alpha) {
    p_text <- round(p_value,round_digits)
    round_digits <- round_digits + 1
  }
  # uses scientific notation when < 0.001
  if (p_value==0) {p_text <- "0"
  } else if (abs(p_value) < 0.001) {
    p_text <- formatC(p_value,format="E", digits=round_digits-1)
    p_split <- strsplit(p_text,"E")[[1]]
    #p_text_stem <- as.numeric(substr(p_text,1,4))
    #p_text_exp <- as.numeric(substr(p_text,6,10))
    p_text <- paste0(p_split[1],"%*%10^",p_split[2])
  }
  
  # returns p-value text and significance color 
  return(p_text)
}
