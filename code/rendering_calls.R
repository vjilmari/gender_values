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

## Invariance examinations

render(input = "code/Invariance.Rmd",
       envir = new.env())

purl(input="code/Invariance.Rmd",
     output="code/Invariance.R",
     documentation = 2)

## Analysis

render(input = "code/Analysis.Rmd",
       envir = new.env())

purl(input="code/Analysis.Rmd",
     output="code/Analysis.R",
     documentation = 2)

## Analysis_with_years (instead of essrounds)

render(input = "code/Analysis_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_with_years.Rmd",
     output="code/Analysis_with_years.R",
     documentation = 2)

## Analysis_power_with_years (instead of essrounds)

render(input = "code/Analysis_power_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_power_with_years.Rmd",
     output="code/Analysis_power_with_years.R",
     documentation = 2)

## Analysis_benevolence_with_years (instead of essrounds)

render(input = "code/Analysis_benevolence_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_benevolence_with_years.Rmd",
     output="code/Analysis_benevolence_with_years.R",
     documentation = 2)

## Analysis_achievement_with_years (instead of essrounds)

render(input = "code/Analysis_achievement_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_achievement_with_years.Rmd",
     output="code/Analysis_achievement_with_years.R",
     documentation = 2)

## Analysis_universalism_with_years (instead of essrounds)

render(input = "code/Analysis_universalism_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_universalism_with_years.Rmd",
     output="code/Analysis_universalism_with_years.R",
     documentation = 2)

## Analysis_stimulation_with_years (instead of essrounds)

render(input = "code/Analysis_stimulation_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_stimulation_with_years.Rmd",
     output="code/Analysis_stimulation_with_years.R",
     documentation = 2)

## Analysis_conformity_with_years (instead of essrounds)

render(input = "code/Analysis_conformity_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_conformity_with_years.Rmd",
     output="code/Analysis_conformity_with_years.R",
     documentation = 2)

## Analysis_tradition_with_years (instead of essrounds)

render(input = "code/Analysis_tradition_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_tradition_with_years.Rmd",
     output="code/Analysis_tradition_with_years.R",
     documentation = 2)

## Analysis_security_with_years (instead of essrounds)

render(input = "code/Analysis_security_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_security_with_years.Rmd",
     output="code/Analysis_security_with_years.R",
     documentation = 2)

## Analysis_self_direction_with_years (instead of essrounds)

render(input = "code/Analysis_self_direction_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_self_direction_with_years.Rmd",
     output="code/Analysis_self_direction_with_years.R",
     documentation = 2)

## Analysis_hedonism_with_years (instead of essrounds)

render(input = "code/Analysis_hedonism_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_hedonism_with_years.Rmd",
     output="code/Analysis_hedonism_with_years.R",
     documentation = 2)


## Analysis_with_years (instead of essrounds), exclude RQ3 because models run so slowly

render(input = "code/Analysis_with_years_no_RQ3.Rmd",
       envir = new.env())

purl(input="code/Analysis_with_years_no_RQ3.Rmd",
     output="code/Analysis_with_years_no_RQ3.R",
     documentation = 2)

## Analysis_youth_with_years (instead of essrounds)

render(input = "code/Analysis_youth_with_years.Rmd",
       envir = new.env())

purl(input="code/Analysis_youth_with_years.Rmd",
     output="code/Analysis_youth_with_years.R",
     documentation = 2)


## Analysis_youth

render(input = "code/Analysis_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_youth.Rmd",
     output="code/Analysis_youth.R",
     documentation = 2)

## Analysis_universalism

render(input = "code/Analysis_universalism.Rmd",
       envir = new.env())

purl(input="code/Analysis_universalism.Rmd",
     output="code/Analysis_universalism.R",
     documentation = 2)

## Analysis_achievement

render(input = "code/Analysis_achievement.Rmd",
       envir = new.env())

purl(input="code/Analysis_achievement.Rmd",
     output="code/Analysis_achievement.R",
     documentation = 2)

## Analysis_benevolence

render(input = "code/Analysis_benevolence.Rmd",
       envir = new.env())

purl(input="code/Analysis_benevolence.Rmd",
     output="code/Analysis_benevolence.R",
     documentation = 2)

## Analysis_power

render(input = "code/Analysis_power.Rmd",
       envir = new.env())

purl(input="code/Analysis_power.Rmd",
     output="code/Analysis_power.R",
     documentation = 2)

## Analysis_stimulation

render(input = "code/Analysis_stimulation.Rmd",
       envir = new.env())

purl(input="code/Analysis_stimulation.Rmd",
     output="code/Analysis_stimulation.R",
     documentation = 2)

## Analysis_power_youth

render(input = "code/Analysis_power_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_power_youth.Rmd",
     output="code/Analysis_power_youth.R",
     documentation = 2)

## Analysis_universalism_youth

render(input = "code/Analysis_universalism_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_universalism_youth.Rmd",
     output="code/Analysis_universalism_youth.R",
     documentation = 2)

## Analysis_benevolence_youth

render(input = "code/Analysis_benevolence_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_benevolence_youth.Rmd",
     output="code/Analysis_benevolence_youth.R",
     documentation = 2)

## Analysis_achievement_youth

render(input = "code/Analysis_achievement_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_achievement_youth.Rmd",
     output="code/Analysis_achievement_youth.R",
     documentation = 2)


## Analysis_stimulation_youth

render(input = "code/Analysis_stimulation_youth.Rmd",
       envir = new.env())

purl(input="code/Analysis_stimulation_youth.Rmd",
     output="code/Analysis_stimulation_youth.R",
     documentation = 2)

## VBMT_reliability

render(input = "code/VBMT_reliability.Rmd",
       envir = new.env())

purl(input="code/VBMT_reliability.Rmd",
     output="code/VBMT_reliability.R",
     documentation = 2)

## cl_interaction_and_simple_slope_power

render(input = "code/cl_interaction_and_simple_slope_power.Rmd",
       envir = new.env())

purl(input="code/cl_interaction_and_simple_slope_power.Rmd",
     output="code/cl_interaction_and_simple_slope_power.R",
     documentation = 2)