# Author: Eric Kalosa-Kenyon
#
# Analysis script for power method c code

## @knitr READ

# subroutines
html_table_width <- function(kable_output, width){
  width_html <- paste0(paste0('<col width="', width, '">'), collapse = "\n")
  sub("<table>", paste0("<table>\n", width_html), kable_output)
}

## @knitr GRAPH
a = 0
