# helper function for converting a p-value to legible text when parsed

pvalue_to_text <- function(p_value, round_digits=3, alpha=0.05) {
  # extends shown digits in cases where rounded p == 0.05
  p_text <- round(p_value,round_digits)
  while (p_text==alpha | !exists("p_text")) {
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
