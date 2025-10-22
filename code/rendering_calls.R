library(rmarkdown)
library(knitr)

# Render RMDs to .html and .md and convert .RMD to .R

## Preparations

render(input = "code/Preparations_ESS.Rmd",
       envir = new.env())

purl(input="code/Preparations_ESS.Rmd",
     output="code/Preparations_ESS.R",
     documentation = 2)


render(input = "code/Preparations_UN.Rmd",
       envir = new.env())

purl(input="code/Preparations_UN.Rmd",
     output="code/Preparations_UN.R",
     documentation = 2)


render(input = "code/Preparations_GGGI.Rmd",
       envir = new.env())

purl(input="code/Preparations_GGGI.Rmd",
     output="code/Preparations_GGGI.R",
     documentation = 2)

render(input = "code/Preparations_GDP.Rmd",
       envir = new.env())

purl(input="code/Preparations_GDP.Rmd",
     output="code/Preparations_GDP.R",
     documentation = 2)

## Analysis

render(input = "code/Analysis.Rmd",
       envir = new.env())

purl(input="code/Analysis.Rmd",
     output="code/Analysis.R",
     documentation = 2)

## Analysis youth

render(input = "code/Analysis_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_youth.Rmd",
     output="code/Analysis_youth.R",
     documentation = 2)

