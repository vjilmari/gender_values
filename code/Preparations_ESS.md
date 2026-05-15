---
title: "Data preparations for European Social Survey"
output: 
  html_document: 
    toc: true
    keep_md: true
date: "2026-05-15"
---



# Packages


``` r
library(rio)
library(dplyr)
```

# Read Data

This reads directly from the zip


``` r
zip_path <- "../data/ESS-data_subset.zip"
tmp_dir  <- tempdir()

unzip(zip_path, files = "Datafile-subset.sav", exdir = tmp_dir, overwrite = TRUE)

d <- import(file.path(tmp_dir, "Datafile-subset.sav"))

table(d$essround,useNA="always")
```

```
## 
##     1     2     3     4     5     6     7     8     9    10    11  <NA> 
## 42359 47537 43000 56752 52458 54673 40185 44387 49519 59685 50116     0
```

``` r
table(d$cntry,useNA="always")
```

```
## 
##    AL    AT    BE    BG    CH    CY    CZ    DE    DK    EE    ES    FI    FR    GB    GR    HR    HU 
##  1201 17579 19045 15479 18309  6750 20090 36845 12408 18149 21296 21095 20809 22663 15315  8098 18760 
##    IE    IL    IS    IT    LT    LU    LV    ME    MK    NL    NO    PL    PT    RO    RS    RU    SE 
## 24250 17124  4817 13043 13017  3187  5173  4087  1429 20024 17402 19131 19254  2146  5111 12458 19446 
##    SI    SK    TR    UA    XK  <NA> 
## 14732 12734  4272 12648  1295     0
```

# Calculate the number of waves per each country


``` r
gd<-data.frame(cntry=unique(d$cntry))
gd$waves<-NA

for (i in 1:nrow(gd)){
  gd[i,"waves"]<-length(unique(d[d$cntry==unique(d$cntry)[i],
                                 "essround"]))
}
gd  
```

```
##    cntry waves
## 1     AT     8
## 2     BE    11
## 3     CH    11
## 4     CZ     9
## 5     DE    11
## 6     DK     8
## 7     ES    11
## 8     FI    11
## 9     FR    11
## 10    GB    11
## 11    GR     6
## 12    HU    11
## 13    IE    11
## 14    IL     8
## 15    IT     6
## 16    LU     2
## 17    NL    11
## 18    NO    11
## 19    PL    11
## 20    PT    11
## 21    SE    11
## 22    SI    11
## 23    EE    10
## 24    IS     6
## 25    SK     8
## 26    TR     2
## 27    UA     6
## 28    BG     7
## 29    CY     7
## 30    RU     5
## 31    HR     5
## 32    LV     4
## 33    RO     1
## 34    LT     7
## 35    AL     1
## 36    XK     1
## 37    ME     3
## 38    RS     3
## 39    MK     1
```

``` r
# attach number of rounds variable to data

d<-left_join(
  x=d,
  y=gd,
  by="cntry")

table(d$waves,useNA="always")
```

```
## 
##      1      2      3      4      5      6      7      8      9     10     11   <NA> 
##   6071   7459   9198   5173  20556  45823  35246  59845  20090  18149 313061      0
```

# Gender


``` r
# code binary gender variable
attributes(d$gndr)
```

```
## $label
## [1] "Gender"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##      Male    Female No answer 
##         1         2         9
```

``` r
table(d$gndr,useNA="always")
```

```
## 
##      1      2   <NA> 
## 249642 290062    967
```

``` r
d$gndr.bin<-ifelse(d$gndr==2,0,d$gndr)
table(d$gndr.bin,useNA="always")
```

```
## 
##      0      1   <NA> 
## 290062 249642    967
```

``` r
# code effect coded gender variable
d$gndr.c<-d$gndr.bin-0.5
table(d$gndr.c,useNA="always")
```

```
## 
##   -0.5    0.5   <NA> 
## 290062 249642    967
```

# Childlessness


``` r
# first based only on the single variable
# this variable has a large number of missing values
# because it is a follow up question
# but the prior question is only present in ESS1-8
attributes(d$chldhhe)
```

```
## $label
## [1] "Ever had children living in household"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##            Yes             No Not applicable        Refusal     Don't know      No answer 
##              1              2              6              7              8              9
```

``` r
table(d$chldhhe,useNA="always")
```

```
## 
##      1      2   <NA> 
## 169277 166856 204538
```

``` r
table(d$essround,d$chldhhe,useNA="always")
```

```
##       
##            1     2  <NA>
##   1    11678 12850 17831
##   2    13101 14563 19873
##   3    13095 13203 16702
##   4    15844 17679 23229
##   5    16666 15584 20208
##   6    17504 15863 21306
##   7    13347 12299 14539
##   8    14928 13858 15601
##   9    17087 15005 17427
##   10   18509 19971 21205
##   11   17518 15981 16617
##   <NA>     0     0     0
```

``` r
prop.table(table(d$essround,d$chldhhe,useNA="always"),1)
```

```
##       
##                1         2      <NA>
##   1    0.2756911 0.3033594 0.4209495
##   2    0.2755959 0.3063508 0.4180533
##   3    0.3045349 0.3070465 0.3884186
##   4    0.2791796 0.3115133 0.4093072
##   5    0.3177018 0.2970758 0.3852225
##   6    0.3201580 0.2901432 0.3896988
##   7    0.3321389 0.3060595 0.3618017
##   8    0.3363147 0.3122085 0.3514768
##   9    0.3450595 0.3030150 0.3519255
##   10   0.3101114 0.3346067 0.3552819
##   11   0.3495490 0.3188802 0.3315708
##   <NA>
```

``` r
d$childless1<-
  case_when(d$chldhhe==2~1,
            d$chldhhe==1~0,
            TRUE~NA_integer_)
table(d$childless1,useNA="always")
```

```
## 
##      0      1   <NA> 
## 169277 166856 204538
```

``` r
prop.table(table(d$childless1,useNA="always"))
```

```
## 
##         0         1      <NA> 
## 0.3130869 0.3086091 0.3783040
```

``` r
prop.table(table(d$childless1))
```

```
## 
##         0         1 
## 0.5036013 0.4963987
```

``` r
# then based on the two variables
# this latter includes those who live with a kid currently

attributes(d$chldhm)
```

```
## $label
## [1] "Children living at home or not"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
## Respondent lives with children at household grid                                         Does not 
##                                                1                                                2 
##                                    Not available 
##                                                9
```

``` r
table(d$chldhhe,d$chldhm,useNA="always")
```

```
##       
##             1      2   <NA>
##   1      2125 113831  53321
##   2       412 115174  51270
##   <NA> 142325   6475  55738
```

``` r
d$childless2<-
  case_when(d$chldhhe==2 & d$chldhm==2~1, # childless if ever and currently childless
            d$chldhhe==1 | d$chldhm==1 ~0, # non-childless if ever OR currently has kids in the HH
            TRUE~NA_integer_)
table(d$childless2,useNA="always")
```

```
## 
##      0      1   <NA> 
## 312014 115174 113483
```

``` r
prop.table(table(d$childless2,useNA="always"))
```

```
## 
##         0         1      <NA> 
## 0.5770866 0.2130205 0.2098929
```

``` r
prop.table(table(d$childless2))
```

```
## 
##         0         1 
## 0.7303904 0.2696096
```

``` r
table(d$essround,d$childless2,useNA="always")
```

```
##       
##            0     1  <NA>
##   1    27938 12612  1809
##   2    31917 14439  1181
##   3    29200 13098   702
##   4    37050 17498  2204
##   5    36346 15554   558
##   6    38473 15841   359
##   7    27651 12295   239
##   8    30325 13837   225
##   9    17087     0 32432
##   10   18509     0 41176
##   11   17518     0 32598
##   <NA>     0     0     0
```

``` r
# check the responses in rshipa
attributes(d$rship10)
```

```
## $label
## [1] "Tenth person in household: relationship to respondent"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##      Husband/wife/partner Son/daughter/step/adopted      Parent/parent-in-law            Other relative 
##                         1                         2                         3                         4 
##        Other non-relative            Not applicable                   Refusal                Don't know 
##                         5                         6                         7                         8 
##                 No answer 
##                         9
```

``` r
attributes(d$rshipa20)
```

```
## $label
## [1] "Twentieth person in household: Relationship to respondent"
## 
## $format.spss
## [1] "F2.0"
## 
## $labels
##               Husband/wife/partner   Son/daughter/step/adopted/foster               Parent/parent-in-law 
##                                  1                                  2                                  3 
## Brother/sister/step/adopted/foster                     Other relative                 Other non-relative 
##                                  4                                  5                                  6 
##                     Not applicable                            Refusal                         Don't know 
##                                 66                                 77                                 88 
##                          No answer 
##                                 99
```

``` r
hh_rshipa=paste0("rshipa",2:24)
hh_rship=paste0("rship",2:9)

d$hh_child1<-apply(d[,hh_rship], 1, function(row) as.integer(any(row == 2)))
table(d$hh_child1,useNA="always")
```

```
## 
##      0      1   <NA> 
##     64  16951 523656
```

``` r
d$hh_child2<-apply(d[,hh_rshipa], 1, function(row) as.integer(any(row == 2)))
table(d$hh_child2,useNA="always")
```

```
## 
##      1   <NA> 
## 179347 361324
```

``` r
d$hh_child<-ifelse(d$hh_child1==1 | d$hh_child2==1,1,0)
table(d$hh_child,useNA="always")
```

```
## 
##      1   <NA> 
## 196298 344373
```

``` r
# code NA as 0
d$hh_child[is.na(d$hh_child)]<-0
table(d$hh_child,useNA="always")
```

```
## 
##      0      1   <NA> 
## 344373 196298      0
```

``` r
table(d$chldhm,d$hh_child,useNA="always")
```

```
##       
##             0      1   <NA>
##   1       962 143900      0
##   2    235463     17      0
##   <NA> 107948  52381      0
```

``` r
# try to apply this same system for all waves
# if there is a child in the HH, then the person in not childless
# if there is not a child in the HH, then the second question is used
# to define if the person is childless


d$childless3<-
  case_when(d$hh_child==1 | d$chldhhe==1~0, # non-childless if kid at home OR ever had
            d$hh_child==0 & d$chldhhe==2~1, # childless if no kid at household AND never was
            TRUE~NA_integer_)
table(d$childless3,useNA="always")
```

```
## 
##      0      1   <NA> 
## 363111 166323  11237
```

``` r
prop.table(table(d$childless3,useNA="always"))
```

```
## 
##          0          1       <NA> 
## 0.67159326 0.30762331 0.02078343
```

``` r
prop.table(table(d$childless3))
```

```
## 
##         0         1 
## 0.6858475 0.3141525
```

# essround


``` r
# center around midpoint (first round + last round)/2 = (1+11)/2 = 6
d$essround.c<-d$essround-6
table(d$essround.c)
```

```
## 
##    -5    -4    -3    -2    -1     0     1     2     3     4     5 
## 42359 47537 43000 56752 52458 54673 40185 44387 49519 59685 50116
```

# Age


``` r
attributes(d$agea)
```

```
## $label
## [1] "Age of respondent, calculated"
## 
## $format.spss
## [1] "F3.0"
## 
## $labels
## Not available 
##           999
```

``` r
d$age.c<-scale(d$agea,center = T,scale=T)

# belong to the focal age range of childlessness examination

d$age_included<-
  case_when(d$agea>17 & d$agea<46~1,
            TRUE~0)
table(d$age_included,useNA="always")
```

```
## 
##      0      1   <NA> 
## 317390 223281      0
```

``` r
table(d$essround,
      d$age_included,useNA="always")
```

```
##       
##            0     1  <NA>
##   1    22684 19675     0
##   2    25745 21792     0
##   3    24078 18922     0
##   4    31575 25177     0
##   5    30348 22110     0
##   6    31684 22989     0
##   7    24080 16105     0
##   8    26360 18027     0
##   9    31373 18146     0
##   10   37453 22232     0
##   11   32010 18106     0
##   <NA>     0     0     0
```

# Education


``` r
attributes(d$eduyrs)
```

```
## $label
## [1] "Years of full-time education completed"
## 
## $format.spss
## [1] "F4.2"
## 
## $labels
##    Refusal Don't know  No answer 
##         77         88         99
```

``` r
table(d$eduyrs,useNA="always")
```

```
## 
##     0     1     2     3     4     5     6     7     8     9    10    11  11.5    12    13    14  14.5 
##  3171   857  1454  3239 10022  6304 10788 10587 31589 30549 34377 54172     4 89982 49998 37361     1 
##    15    16    17    18  18.5    19    20    21    22    23    24    25    26    27    28    29    30 
## 37720 36630 27389 22316     1 10956 10607  3603  2972  1731  1187  1168   350   219   158    89   307 
##    31    32    33    34    35    36    37    38    39    40    41    42    43    44    45    47    48 
##    54    36    28    13    37    20    10    14     9   102     2     8     4     6    16     3     6 
##    50    51    54    55    56    60    65    67    69    76  <NA> 
##    17     4     2     2     2     3     1     1     1     3  8409
```

``` r
d$eduyrs.c<-scale(d$eduyrs,center = T,scale=T)
```

# Religiosity


``` r
attributes(d$rlgdgr)
```

```
## $label
## [1] "How religious are you"
## 
## $format.spss
## [1] "F2.0"
## 
## $labels
## Not at all religious                    1                    2                    3                    4 
##                    0                    1                    2                    3                    4 
##                    5                    6                    7                    8                    9 
##                    5                    6                    7                    8                    9 
##       Very religious              Refusal           Don't know            No answer 
##                   10                   77                   88                   99
```

``` r
table(d$rlgdgr,useNA="always")
```

```
## 
##     0     1     2     3     4     5     6     7     8     9    10  <NA> 
## 82328 30323 37924 41349 33566 86913 52115 60586 53626 22827 33752  5362
```

``` r
d$rlgdgr.c<-scale(d$rlgdgr,center=T,scale=T)

# country mean-level religiosity

reld<-d %>%
  group_by(cntry) %>%
  summarise(rlgdgr.cm=mean(rlgdgr.c,na.rm=T))
reld
```

```
## # A tibble: 39 × 2
##    cntry rlgdgr.cm
##    <chr>     <dbl>
##  1 AL       0.843 
##  2 AT       0.0319
##  3 BE      -0.0265
##  4 BG      -0.117 
##  5 CH       0.112 
##  6 CY       0.661 
##  7 CZ      -0.736 
##  8 DE      -0.262 
##  9 DK      -0.180 
## 10 EE      -0.401 
## # ℹ 29 more rows
```

``` r
d<-
  left_join(
    x=d,
    y=reld,
    by="cntry")
```

# Country x time datasets


``` r
d$cntry_time<-paste0(d$cntry,"_",d$essround)
table(d$cntry_time,useNA="always")
```

```
## 
##  AL_6  AT_1 AT_10 AT_11  AT_2  AT_3  AT_7  AT_8  AT_9  BE_1 BE_10 BE_11  BE_2  BE_3  BE_4  BE_5  BE_6 
##  1201  2257  2003  2354  2256  2405  1795  2010  2499  1899  1341  1594  1778  1798  1760  1704  1869 
##  BE_7  BE_8  BE_9 BG_10 BG_11  BG_3  BG_4  BG_5  BG_6  BG_9  CH_1 CH_10 CH_11  CH_2  CH_3  CH_4  CH_5 
##  1769  1766  1767  2718  2239  1400  2230  2434  2260  2198  2040  1523  1384  2141  1804  1819  1506 
##  CH_6  CH_7  CH_8  CH_9 CY_10 CY_11  CY_3  CY_4  CY_5  CY_6  CY_9  CZ_1 CZ_10  CZ_2  CZ_4  CZ_5  CZ_6 
##  1493  1532  1525  1542   875   685   995  1215  1083  1116   781  1360  2476  3026  2018  2386  2009 
##  CZ_7  CZ_8  CZ_9  DE_1 DE_10 DE_11  DE_2  DE_3  DE_4  DE_5  DE_6  DE_7  DE_8  DE_9  DK_1  DK_2  DK_3 
##  2148  2269  2398  2919  8725  2420  2870  2916  2751  3031  2958  3045  2852  2358  1506  1487  1505 
##  DK_4  DK_5  DK_6  DK_7  DK_9 EE_10 EE_11  EE_2  EE_3  EE_4  EE_5  EE_6  EE_7  EE_8  EE_9  ES_1 ES_10 
##  1610  1576  1650  1502  1572  1542  1293  1989  1517  1661  1793  2380  2051  2019  1904  1729  2283 
## ES_11  ES_2  ES_3  ES_4  ES_5  ES_6  ES_7  ES_8  ES_9  FI_1 FI_10 FI_11  FI_2  FI_3  FI_4  FI_5  FI_6 
##  1844  1663  1876  2576  1885  1889  1925  1958  1668  2000  1577  1563  2022  1896  2195  1878  2197 
##  FI_7  FI_8  FI_9  FR_1 FR_10 FR_11  FR_2  FR_3  FR_4  FR_5  FR_6  FR_7  FR_8  FR_9  GB_1 GB_10 GB_11 
##  2087  1925  1755  1503  1977  1771  1806  1986  2073  1728  1968  1917  2070  2010  2052  1149  1684 
##  GB_2  GB_3  GB_4  GB_5  GB_6  GB_7  GB_8  GB_9  GR_1 GR_10 GR_11  GR_2  GR_4  GR_5 HR_10 HR_11  HR_4 
##  1897  2394  2352  2422  2286  2264  1959  2204  2566  2799  2757  2406  2072  2715  1592  1563  1484 
##  HR_5  HR_9  HU_1 HU_10 HU_11  HU_2  HU_3  HU_4  HU_5  HU_6  HU_7  HU_8  HU_9  IE_1 IE_10 IE_11  IE_2 
##  1649  1810  1685  1849  2118  1498  1518  1544  1561  2014  1698  1614  1661  2046  1770  2017  2286 
##  IE_3  IE_4  IE_5  IE_6  IE_7  IE_8  IE_9  IL_1 IL_10 IL_11  IL_4  IL_5  IL_6  IL_7  IL_8 IS_10 IS_11 
##  1800  1764  2576  2628  2390  2757  2216  2499  1308   906  2490  2294  2508  2562  2557   903   842 
##  IS_2  IS_6  IS_8  IS_9  IT_1 IT_10 IT_11  IT_6  IT_8  IT_9 LT_10 LT_11  LT_5  LT_6  LT_7  LT_8  LT_9 
##   579   752   880   861  1207  2640  2865   960  2626  2745  1659  1365  1677  2109  2250  2122  1835 
##  LU_1  LU_2 LV_10 LV_11  LV_4  LV_9 ME_10 ME_11  ME_9 MK_10  NL_1 NL_10 NL_11  NL_2  NL_3  NL_4  NL_5 
##  1552  1635  1023  1252  1980   918  1278  1609  1200  1429  2364  1470  1695  1881  1889  1778  1829 
##  NL_6  NL_7  NL_8  NL_9  NO_1 NO_10 NO_11  NO_2  NO_3  NO_4  NO_5  NO_6  NO_7  NO_8  NO_9  PL_1 PL_10 
##  1845  1919  1681  1673  2036  1411  1337  1760  1750  1549  1548  1624  1436  1545  1406  2110  2065 
## PL_11  PL_2  PL_3  PL_4  PL_5  PL_6  PL_7  PL_8  PL_9  PT_1 PT_10 PT_11  PT_2  PT_3  PT_4  PT_5  PT_6 
##  1442  1716  1721  1619  1751  1898  1615  1694  1500  1511  1838  1373  2052  2222  2367  2150  2151 
##  PT_7  PT_8  PT_9  RO_4 RS_10 RS_11  RS_9  RU_3  RU_4  RU_5  RU_6  RU_8  SE_1 SE_10 SE_11  SE_2  SE_3 
##  1265  1270  1055  2146  1505  1563  2043  2437  2512  2595  2484  2430  1999  2287  1230  1948  1927 
##  SE_4  SE_5  SE_6  SE_7  SE_8  SE_9  SI_1 SI_10 SI_11  SI_2  SI_3  SI_4  SI_5  SI_6  SI_7  SI_8  SI_9 
##  1830  1497  1847  1791  1551  1539  1519  1252  1248  1442  1476  1286  1403  1257  1224  1307  1318 
## SK_10 SK_11  SK_2  SK_3  SK_4  SK_5  SK_6  SK_9  TR_2  TR_4 UA_11  UA_2  UA_3  UA_4  UA_5  UA_6  XK_6 
##  1418  1442  1512  1766  1810  1856  1847  1083  1856  2416  2661  2031  2002  1845  1931  2178  1295 
##  <NA> 
##     0
```

# Same-gender partnership


``` r
# must be done separately for ESS1

# use a function that checks if there are 
# same-sex partners living in the household

has_samegndr_partner<-function(data,hh_gndr,hh_rshipa){
  
  gndr_self<-data[,"gndr"]
  test_mat<-matrix(NA,nrow=nrow(data),ncol=length(hh_gndr))
  
  for (i in 1:length(hh_gndr)){
    test_mat[,i]<-
      (gndr_self == data[,hh_gndr[i]] & data[,hh_rshipa[i]]==1)
  }
  
  result<-rowSums(test_mat,na.rm=T)
  return(result)
  
}


hh_gndr=paste0("gndr",2:24)
hh_rshipa=paste0("rshipa",2:24)

d$same_gndr210<-has_samegndr_partner(data=d,
                        hh_gndr=hh_gndr,
                        hh_rshipa=hh_rshipa)
table(d$same_gndr210,useNA="always")
```

```
## 
##      0      1      2      3      4   <NA> 
## 536227   4416     22      5      1      0
```

``` r
d$same_gndr1<-has_samegndr_partner(
  data=d,
  hh_gndr=paste0("gndr",2:9),
  hh_rshipa=paste0("rship",2:9)
)

table(d$same_gndr1,useNA="always")
```

```
## 
##      0      1      2   <NA> 
## 540361    308      2      0
```

``` r
# combine

d$same_gndr_partner<-
  case_when(d$same_gndr1>0 | d$same_gndr210>0~1,
            TRUE~0)
table(d$same_gndr_partner,useNA="always")
```

```
## 
##      0      1   <NA> 
## 535917   4754      0
```

# Values

For rounds 1-10 the variable names are the same. For round 11, the variables have an additional "a" at the end. First merge them to same variables.


``` r
table(d$essround,is.na(d$impdiff)) 
```

```
##     
##      FALSE  TRUE
##   1  37546  4813
##   2  44361  3176
##   3  41370  1630
##   4  55045  1707
##   5  51158  1300
##   6  53763   910
##   7  39103  1082
##   8  43562   825
##   9  48486  1033
##   10 37113 22572
##   11     0 50116
```

``` r
table(d$essround,is.na(d$impdiffa))
```

```
##     
##      FALSE  TRUE
##   1      0 42359
##   2      0 47537
##   3      0 43000
##   4      0 56752
##   5      0 52458
##   6      0 54673
##   7      0 40185
##   8      0 44387
##   9      0 49519
##   10     0 59685
##   11 49251   865
```

``` r
table(!is.na(d$impdiff),!is.na(d$impdiffa))
```

```
##        
##          FALSE   TRUE
##   FALSE  39913  49251
##   TRUE  451507      0
```

``` r
value_items <- c(
  "impdiff",
  "impenv",
  "impfree",
  "impfun",
  "imprich",
  "impsafe",
  "imptrad",
  "ipadvnt",
  "ipbhprp",
  "ipcrtiv",
  "ipeqopt",
  "ipfrule",
  "ipgdtim",
  "iphlppl",
  "iplylfr",
  "ipmodst",
  "iprspot",
  "ipshabt",
  "ipstrgv",
  "ipsuces",
  "ipudrst"
)
```


``` r
for (x in value_items) {
  var_a   <- paste0(x, "a")
  var_out <- paste0(x, "_comb")
  d[[var_out]] <- ifelse(!is.na(d[[var_a]]), d[[var_a]], d[[x]])
}

table(d$essround,is.na(d$impdiff_comb))
```

```
##     
##      FALSE  TRUE
##   1  37546  4813
##   2  44361  3176
##   3  41370  1630
##   4  55045  1707
##   5  51158  1300
##   6  53763   910
##   7  39103  1082
##   8  43562   825
##   9  48486  1033
##   10 37113 22572
##   11 49251   865
```




``` r
## reverse code all value items

d$impdiff.R<-7-d$impdiff_comb
d$impenv.R<-7-d$impenv_comb
d$impfree.R<-7-d$impfree_comb
d$impfun.R<-7-d$impfun_comb
d$imprich.R<-7-d$imprich_comb
d$impsafe.R<-7-d$impsafe_comb
d$imptrad.R<-7-d$imptrad_comb
d$ipadvnt.R<-7-d$ipadvnt_comb
d$ipbhprp.R<-7-d$ipbhprp_comb
d$ipcrtiv.R<-7-d$ipcrtiv_comb
d$ipeqopt.R<-7-d$ipeqopt_comb
d$ipfrule.R<-7-d$ipfrule_comb
d$ipgdtim.R<-7-d$ipgdtim_comb
d$iphlppl.R<-7-d$iphlppl_comb
d$iplylfr.R<-7-d$iplylfr_comb
d$ipmodst.R<-7-d$ipmodst_comb
d$iprspot.R<-7-d$iprspot_comb
d$ipshabt.R<-7-d$ipshabt_comb
d$ipstrgv.R<-7-d$ipstrgv_comb
d$ipsuces.R<-7-d$ipsuces_comb
d$ipudrst.R<-7-d$ipudrst_comb
```

Calculate the ten value aggregates from the combined and reverse scored items


``` r
# conformity
attributes(d$ipfrule)
```

```
## $label
## [1] "Important to do what is told and follow rules"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$ipbhprp)
```

```
## $label
## [1] "Important to behave properly"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$con<-rowMeans(d[,c("ipfrule.R","ipbhprp.R")],na.rm=T)

# tradition
attributes(d$ipmodst)
```

```
## $label
## [1] "Important to be humble and modest, not draw attention"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$imptrad)
```

```
## $label
## [1] "Important to follow traditions and customs"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$tra<-rowMeans(d[,c("ipmodst.R","imptrad.R")],na.rm=T)

# benevolence

attributes(d$iphlppl)
```

```
## $label
## [1] "Important to help people and care for others well-being"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$iplylfr)
```

```
## $label
## [1] "Important to be loyal to friends and devote to people close"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$ben<-rowMeans(d[,c("iphlppl.R","iplylfr.R")],na.rm=T)

# universalism

attributes(d$ipeqopt)
```

```
## $label
## [1] "Important that people are treated equally and have equal opportunities"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$ipudrst)
```

```
## $label
## [1] "Important to understand different people"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$impenv)
```

```
## $label
## [1] "Important to care for nature and environment"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$uni<-rowMeans(d[,c("ipeqopt.R","ipudrst.R","impenv.R")],na.rm=T)

# self-diretion

attributes(d$ipcrtiv)
```

```
## $label
## [1] "Important to think new ideas and being creative"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$impfree)
```

```
## $label
## [1] "Important to make own decisions and be free"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$sdi<-rowMeans(d[,c("ipcrtiv.R","impfree.R")],na.rm=T)

# stimulation

attributes(d$impdiff)
```

```
## $label
## [1] "Important to try new and different things in life"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$ipadvnt)
```

```
## $label
## [1] "Important to seek adventures and have an exciting life"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$sti<-rowMeans(d[,c("impdiff.R","ipadvnt.R")],na.rm=T)

# hedonism

attributes(d$ipgdtim)
```

```
## $label
## [1] "Important to have a good time"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$impfun)
```

```
## $label
## [1] "Important to seek fun and things that give pleasure"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$hed<-rowMeans(d[,c("ipgdtim.R","impfun.R")],na.rm=T)

# achievement

attributes(d$ipshabt)
```

```
## $label
## [1] "Important to show abilities and be admired"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$ipsuces)
```

```
## $label
## [1] "Important to be successful and that people recognize achievements"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$ach<-rowMeans(d[,c("ipshabt.R","ipsuces.R")],na.rm=T)

# power

attributes(d$imprich)
```

```
## $label
## [1] "Important to be rich, have money and expensive things"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$iprspot)
```

```
## $label
## [1] "Important to get respect from others"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$pow<-rowMeans(d[,c("imprich.R","iprspot.R")],na.rm=T)

# security

attributes(d$impsafe)
```

```
## $label
## [1] "Important to live in secure and safe surroundings"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
attributes(d$ipstrgv)
```

```
## $label
## [1] "Important that government is strong and ensures safety"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##  Very much like me            Like me   Somewhat like me   A little like me        Not like me 
##                  1                  2                  3                  4                  5 
## Not like me at all            Refusal         Don't know          No answer 
##                  6                  7                  8                  9
```

``` r
d$sec<-rowMeans(d[,c("impsafe.R","ipstrgv.R")],na.rm=T)
```

## Within-country standardized values


``` r
# center individual scores around country mean-levels on values
cntry_values<-d %>%
  group_by(cntry) %>%
  summarise(con.cm=mean(con,na.rm=T),
            tra.cm=mean(tra,na.rm=T),
            ben.cm=mean(ben,na.rm=T),
            uni.cm=mean(uni,na.rm=T),
            sdi.cm=mean(sdi,na.rm=T),
            sti.cm=mean(sti,na.rm=T),
            hed.cm=mean(hed,na.rm=T),
            ach.cm=mean(ach,na.rm=T),
            pow.cm=mean(pow,na.rm=T),
            sec.cm=mean(sec,na.rm=T),
            con.csd=sd(con,na.rm=T),
            tra.csd=sd(tra,na.rm=T),
            ben.csd=sd(ben,na.rm=T),
            uni.csd=sd(uni,na.rm=T),
            sdi.csd=sd(sdi,na.rm=T),
            sti.csd=sd(sti,na.rm=T),
            hed.csd=sd(hed,na.rm=T),
            ach.csd=sd(ach,na.rm=T),
            pow.csd=sd(pow,na.rm=T),
            sec.csd=sd(sec,na.rm=T))

# add the means and sds to original data frame
d<-left_join(
  x=d,
  y=cntry_values,
  by="cntry"
)  


d$con.cz<-(d$con-d$con.cm)/d$con.csd
d$tra.cz<-(d$tra-d$tra.cm)/d$tra.csd
d$ben.cz<-(d$ben-d$ben.cm)/d$ben.csd
d$uni.cz<-(d$uni-d$uni.cm)/d$uni.csd
d$sdi.cz<-(d$sdi-d$sdi.cm)/d$sdi.csd
d$sti.cz<-(d$sti-d$sti.cm)/d$sti.csd
d$hed.cz<-(d$hed-d$hed.cm)/d$hed.csd
d$ach.cz<-(d$ach-d$ach.cm)/d$ach.csd
d$pow.cz<-(d$pow-d$pow.cm)/d$pow.csd
d$sec.cz<-(d$sec-d$sec.cm)/d$sec.csd
```

## Within-country standardized value items


``` r
cntry_items<-d %>%
  group_by(cntry) %>%
  summarise(
    impdiff.R.cm=mean(impdiff.R,na.rm=T),
    impenv.R.cm=mean(impenv.R,na.rm=T),
    impfree.R.cm=mean(impfree.R,na.rm=T),
    impfun.R.cm=mean(impfun.R,na.rm=T),
    imprich.R.cm=mean(imprich.R,na.rm=T),
    impsafe.R.cm=mean(impsafe.R,na.rm=T),
    imptrad.R.cm=mean(imptrad.R,na.rm=T),
    ipadvnt.R.cm=mean(ipadvnt.R,na.rm=T),
    ipbhprp.R.cm=mean(ipbhprp.R,na.rm=T),
    ipcrtiv.R.cm=mean(ipcrtiv.R,na.rm=T),
    ipeqopt.R.cm=mean(ipeqopt.R,na.rm=T),
    ipfrule.R.cm=mean(ipfrule.R,na.rm=T),
    ipgdtim.R.cm=mean(ipgdtim.R,na.rm=T),
    iphlppl.R.cm=mean(iphlppl.R,na.rm=T),
    iplylfr.R.cm=mean(iplylfr.R,na.rm=T),
    ipmodst.R.cm=mean(ipmodst.R,na.rm=T),
    iprspot.R.cm=mean(iprspot.R,na.rm=T),
    ipshabt.R.cm=mean(ipshabt.R,na.rm=T),
    ipstrgv.R.cm=mean(ipstrgv.R,na.rm=T),
    ipsuces.R.cm=mean(ipsuces.R,na.rm=T),
    ipudrst.R.cm=mean(ipudrst.R,na.rm=T),
    impdiff.R.csd=sd(impdiff.R,na.rm=T),
    impenv.R.csd=sd(impenv.R,na.rm=T),
    impfree.R.csd=sd(impfree.R,na.rm=T),
    impfun.R.csd=sd(impfun.R,na.rm=T),
    imprich.R.csd=sd(imprich.R,na.rm=T),
    impsafe.R.csd=sd(impsafe.R,na.rm=T),
    imptrad.R.csd=sd(imptrad.R,na.rm=T),
    ipadvnt.R.csd=sd(ipadvnt.R,na.rm=T),
    ipbhprp.R.csd=sd(ipbhprp.R,na.rm=T),
    ipcrtiv.R.csd=sd(ipcrtiv.R,na.rm=T),
    ipeqopt.R.csd=sd(ipeqopt.R,na.rm=T),
    ipfrule.R.csd=sd(ipfrule.R,na.rm=T),
    ipgdtim.R.csd=sd(ipgdtim.R,na.rm=T),
    iphlppl.R.csd=sd(iphlppl.R,na.rm=T),
    iplylfr.R.csd=sd(iplylfr.R,na.rm=T),
    ipmodst.R.csd=sd(ipmodst.R,na.rm=T),
    iprspot.R.csd=sd(iprspot.R,na.rm=T),
    ipshabt.R.csd=sd(ipshabt.R,na.rm=T),
    ipstrgv.R.csd=sd(ipstrgv.R,na.rm=T),
    ipsuces.R.csd=sd(ipsuces.R,na.rm=T),
    ipudrst.R.csd=sd(ipudrst.R,na.rm=T))

# add the means and sds to original data frame
d<-left_join(
  x=d,
  y=cntry_items,
  by="cntry"
)  

d$impdiff.R.cz<-(d$impdiff.R-d$impdiff.R.cm)/d$impdiff.R.csd
d$impenv.R.cz<-(d$impenv.R-d$impenv.R.cm)/d$impenv.R.csd
d$impfree.R.cz<-(d$impfree.R-d$impfree.R.cm)/d$impfree.R.csd
d$impfun.R.cz<-(d$impfun.R-d$impfun.R.cm)/d$impfun.R.csd
d$imprich.R.cz<-(d$imprich.R-d$imprich.R.cm)/d$imprich.R.csd
d$impsafe.R.cz<-(d$impsafe.R-d$impsafe.R.cm)/d$impsafe.R.csd
d$imptrad.R.cz<-(d$imptrad.R-d$imptrad.R.cm)/d$imptrad.R.csd
d$ipadvnt.R.cz<-(d$ipadvnt.R-d$ipadvnt.R.cm)/d$ipadvnt.R.csd
d$ipbhprp.R.cz<-(d$ipbhprp.R-d$ipbhprp.R.cm)/d$ipbhprp.R.csd
d$ipcrtiv.R.cz<-(d$ipcrtiv.R-d$ipcrtiv.R.cm)/d$ipcrtiv.R.csd
d$ipeqopt.R.cz<-(d$ipeqopt.R-d$ipeqopt.R.cm)/d$ipeqopt.R.csd
d$ipfrule.R.cz<-(d$ipfrule.R-d$ipfrule.R.cm)/d$ipfrule.R.csd
d$ipgdtim.R.cz<-(d$ipgdtim.R-d$ipgdtim.R.cm)/d$ipgdtim.R.csd
d$iphlppl.R.cz<-(d$iphlppl.R-d$iphlppl.R.cm)/d$iphlppl.R.csd
d$iplylfr.R.cz<-(d$iplylfr.R-d$iplylfr.R.cm)/d$iplylfr.R.csd
d$ipmodst.R.cz<-(d$ipmodst.R-d$ipmodst.R.cm)/d$ipmodst.R.csd
d$iprspot.R.cz<-(d$iprspot.R-d$iprspot.R.cm)/d$iprspot.R.csd
d$ipshabt.R.cz<-(d$ipshabt.R-d$ipshabt.R.cm)/d$ipshabt.R.csd
d$ipstrgv.R.cz<-(d$ipstrgv.R-d$ipstrgv.R.cm)/d$ipstrgv.R.csd
d$ipsuces.R.cz<-(d$ipsuces.R-d$ipsuces.R.cm)/d$ipsuces.R.csd
d$ipudrst.R.cz<-(d$ipudrst.R-d$ipudrst.R.cm)/d$ipudrst.R.csd
```

# Marriage


``` r
attributes(d$marital)
```

```
## $label
## [1] "Legal marital status"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##       Married     Separated      Divorced       Widowed Never married       Refusal    Don't know 
##             1             2             3             4             5             7             8 
##     No answer 
##             9
```

``` r
table(d$marital,useNA="always")
```

```
## 
##      1      2      3      4      5   <NA> 
##  46436   1319   5891   8114  24406 454505
```

``` r
attributes(d$maritala)
```

```
## $label
## [1] "Legal marital status"
## 
## $format.spss
## [1] "F2.0"
## 
## $labels
##                                      Married                       In a civil partnership 
##                                            1                                            2 
##            Separated (still legally married)     Separated (still in a civil partnership) 
##                                            3                                            4 
##                                     Divorced                                      Widowed 
##                                            5                                            6 
## Formerly in civil partnership, now dissolved  Formerly in civil partnership, partner died 
##                                            7                                            8 
## Never married and never in civil partnership                                      Refusal 
##                                            9                                           77 
##                                   Don't know                                    No answer 
##                                           88                                           99
```

``` r
table(d$maritala,useNA="always")
```

```
## 
##      1      2      3      4      5      6      7      8      9   <NA> 
##  51049   1620   1265     97   7694  10207    279     43  25194 443223
```

``` r
attributes(d$maritalb)
```

```
## $label
## [1] "Legal marital status, post coded"
## 
## $format.spss
## [1] "F2.0"
## 
## $labels
##                                                    Legally married 
##                                                                  1 
##                                In a legally registered civil union 
##                                                                  2 
##                                                  Legally separated 
##                                                                  3 
##                             Legally divorced/Civil union dissolved 
##                                                                  4 
##                                         Widowed/Civil partner died 
##                                                                  5 
## None of these (NEVER married or in legally registered civil union) 
##                                                                  6 
##                                                            Refusal 
##                                                                 77 
##                                                         Don't know 
##                                                                 88 
##                                                          No answer 
##                                                                 99
```

``` r
table(d$maritalb,useNA="always")
```

```
## 
##      1      2      3      4      5      6   <NA> 
## 160842   2861   3266  30068  31114  94412 218108
```

``` r
d$married<-case_when(d$marital==1 |
                    d$maritala==1 |
                    d$maritalb==1 ~1,
                    TRUE~0)
table(d$married,useNA="always")
```

```
## 
##      0      1   <NA> 
## 282344 258327      0
```

``` r
# mid-point centered marriage
d$married.c<-d$married-0.5
```

# Domicile


``` r
attributes(d$domicil)
```

```
## $label
## [1] "Domicile, respondent's description"
## 
## $format.spss
## [1] "F1.0"
## 
## $labels
##                       A big city Suburbs or outskirts of big city               Town or small city 
##                                1                                2                                3 
##                  Country village      Farm or home in countryside                          Refusal 
##                                4                                5                                7 
##                       Don't know                        No answer 
##                                8                                9
```

``` r
table(d$domicil,useNA="always")
```

```
## 
##      1      2      3      4      5   <NA> 
## 119151  60402 163463 165730  29809   2116
```

``` r
d$rural=case_when(
  d$domicil==4 | d$domicil==5~1,
  d$domicil<4~0,
  TRUE~NA_integer_
)
table(d$rural,useNA="always")
```

```
## 
##      0      1   <NA> 
## 343016 195539   2116
```

``` r
d$rural.c<-d$rural-0.5
table(d$rural.c,useNA="always")
```

```
## 
##   -0.5    0.5   <NA> 
## 343016 195539   2116
```

# Variable selection


``` r
fdat<-d %>%
  dplyr::select(
    essround,essround.c,idno,cntry,waves,pspwght,cntry_time,
    rlgdgr,rlgdgr.c,rlgdgr.cm,
    gndr,gndr.bin,gndr.c,
    agea,age.c,age_included,
    chldhhe,chldhm,childless1,childless2,childless3,
    eduyrs,eduyrs.c,
    married,married.c,
    rural,rural.c,
    same_gndr_partner,
    con,tra,ben,uni,sdi,sti,hed,ach,pow,sec,
    con.cz,tra.cz,ben.cz,uni.cz,sdi.cz,sti.cz,hed.cz,ach.cz,pow.cz,sec.cz,
    impdiff.R.cz,impenv.R.cz,impfree.R.cz,impfun.R.cz,imprich.R.cz,impsafe.R.cz,imptrad.R.cz,
    ipadvnt.R.cz,ipbhprp.R.cz,ipcrtiv.R.cz,ipeqopt.R.cz,ipfrule.R.cz,ipgdtim.R.cz,iphlppl.R.cz,
    iplylfr.R.cz,ipmodst.R.cz,iprspot.R.cz,ipshabt.R.cz,ipstrgv.R.cz,ipsuces.R.cz,ipudrst.R.cz
  )
```

# Data export


``` r
save(fdat,file="../data/fdat.rdata")
```

# Session information


``` r
s<-sessionInfo()
print(s,locale=F)
```

```
## R version 4.6.0 (2026-04-24 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 22631)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] dplyr_1.2.1    rio_1.3.0      knitr_1.51     rmarkdown_2.31
## 
## loaded via a namespace (and not attached):
##  [1] vctrs_0.7.3       cli_3.6.6         rlang_1.2.0       xfun_0.57         forcats_1.0.1    
##  [6] haven_2.5.5       generics_0.1.4    jsonlite_2.0.0    glue_1.8.1        htmltools_0.5.9  
## [11] sass_0.4.10       hms_1.1.4         evaluate_1.0.5    jquerylib_0.1.4   tibble_3.3.1     
## [16] tzdb_0.5.0        fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5   compiler_4.6.0   
## [21] pkgconfig_2.0.3   R.oo_1.27.1       R.utils_2.13.0    digest_0.6.39     R6_2.6.1         
## [26] utf8_1.2.6        readr_2.2.0       tidyselect_1.2.1  pillar_1.11.1     magrittr_2.0.5   
## [31] bslib_0.10.0      R.methodsS3_1.8.2 withr_3.0.2       tools_4.6.0       cachem_1.1.0
```
