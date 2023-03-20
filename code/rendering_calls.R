library(rmarkdown)
library(knitr)

# Render RMDs to .html and .md and convert .RMD to .R

## Preparations

render(input = "code/Preparations.Rmd",
       envir = new.env())

purl(input="code/Preparations.Rmd",
     output="code/Preparations.R",
     documentation = 2)
