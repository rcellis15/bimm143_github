# Class 12: HW Population Analysis Pt.2
Ryan Ellis (PID: A17673864)

- [Section 1. Proportion of G/G in a
  population](#section-1-proportion-of-gg-in-a-population)
- [Homework Problems](#homework-problems)

## Section 1. Proportion of G/G in a population

Downloaded a CSV file from
ensemble–\>https://www.ensembl.org/Homo_sapiens/Variation/Sample?db=core;r=17:39894595-39895595;v=rs8067378;vdb=variation;vf=959672880#\_

``` r
MXL<- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(MXL)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
round(table(MXL$Genotype..forward.strand.)/nrow(MXL)*100,2)
```


      A|A   A|G   G|A   G|G 
    34.38 32.81 18.75 14.06 

Lets look at the GBR dataset

``` r
GBR<-read.csv("373522-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

``` r
round(table(GBR$Genotype..forward.strand.)/nrow(GBR)*100,2)
```


      A|A   A|G   G|A   G|G 
    25.27 18.68 26.37 29.67 

This variant (GG) that is associated with childhood asthma is more
frequent in GBR populations than MXL populations.

Further dig.

## Homework Problems

> Q13: Read this file into R and determine the sample size for each
> genotype and their corresponding median expression levels for each of
> these genotypes.

``` r
pop<- read.table("rs8067378_ENSG00000172057.6.txt")

nrow(pop)
```

    [1] 462

``` r
table(pop$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
?median()
```

``` r
head(pop)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
pop|>
  arrange(geno)|>
  group_by(geno)|>
  summarize(MedianExp<-median(exp))
```

    # A tibble: 3 × 2
      geno  `MedianExp <- median(exp)`
      <chr>                      <dbl>
    1 A/A                         31.2
    2 A/G                         25.1
    3 G/G                         20.1

Sample size= 462 AA: 108 (median expression:31.2), AG:233(median
expression:25.1), GG:121(median expression: 20.1)

> Q14: Generate a boxplot with a box per genotype, what could you infer
> from the relative expression value between A/A and G/G displayed in
> this plot? Does the SNP effect the expression of ORMDL3?

``` r
library(ggplot2)
ggplot(pop)+
  aes(x=geno, y=exp, fill=geno)+
  geom_boxplot(notch=TRUE)
```

![](class12_files/figure-commonmark/unnamed-chunk-8-1.png)

You can infer using the above box plot, that genotyping has a
correlation with expression levels. With GG showing lower expression
levels than AA. Yes these SNP variation do indeed seem to play a role in
ORMDL3 expression, so depending on the variation this will either have
an increase or decrease in expression levels.
