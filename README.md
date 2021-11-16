
<!-- README.md is generated from README.Rmd. Please edit that file -->

# neopeptides

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The goal of **neopeptides** is to calculate and explore property indices
of peptides for cancer immunotherapy study.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ShixiangWang/neopeptides")
```

### Configuration

To use neopeptide, `blast` is required.

After installation of `blast`, you can set the path with

``` r
set_blast_path()
# e.g. 
# set_blast_path("/Users/wsx/anaconda3/bin/")
```

Then install required protein database with

``` r
install_database()
```

## Available features

  - `calc_geometric_mean()` - Calculate the Geometric Mean.
  - `calc_harmonic_mean()` - Calculate the Harmonic Mean.
  - `calc_iedb_score()` - Calculate IEDB Score for Peptides.
  - `calc_dissimilarity()` - Calculate Dissimilarity Value to Reference
    Proteome for Peptides.

## Examples

``` r
library(neopeptides)

calc_iedb_score(c("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP", "MTEYKLVVVG"))
#> => Running blastp for homology to IEDB antigens..
#> => Summing IEDB local alignments...
#> => Done.
#> => Removing temporary files...
#>                               peptide   iedb_score
#> 1: MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP 5.038247e-01
#> 2:                         MTEYKLVVVG 5.895595e-05
#>                                                                                    annotation
#> 1:                      18142|polyprotein precursor|NP_041724.2|West Nile virus|11082 DVGVSAL
#> 2: 32238|polyprotein [Hepatitis C virus subtype 1a]|ABV46251.2|Hepatitis C virus|11103 KLVVLG
calc_dissimilarity(c("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP", "MTEYKLVVVG"))
#> => Running blastp for homology to self antigens..
#> => Summing local alignments...
#> => Done.
#> => Removing temporary files...
#>                               peptide dissimilarity
#> 1: MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP             0
#> 2:                         MTEYKLVVVG             0
```

## Citation

Revisiting neoantigen depletion signal in the untreated cancer genome.
Shixiang Wang, Xuan Wang, Tao Wu, Zaoke He, Huimin Li, Xiaoqin Sun,
Xue-Song Liu bioRxiv 2020.05.11.089540; doi:
<https://doi.org/10.1101/2020.05.11.089540>

## Researches used this tool

- Wang X, Wang S, Han Y, Xu M, Li P, Ke M, Teng Z, Huang P, Diao Z, Yan Y, Meng Q, Kuang Y, Zheng W, Liu H, Liu X, Jia B. Association of CSMD1 with Tumor Mutation Burden and Other Clinical Outcomes in Gastric Cancer. Int J Gen Med. 2021;14:8293-8299
https://doi.org/10.2147/IJGM.S325910

## Reference

  - Richman LP, Vonderheide RH, Rech AJ. Neoantigen Dissimilarity to the
    Self-Proteome Predicts Immunogenicity and Response to Immune
    Checkpoint Blockade. Cell Systems. 2019 Oct;9(4):375-382.e4. DOI:
    10.1016/j.cels.2019.08.009.
  - R package antigen.garnish
    <https://github.com/immune-health/antigen.garnish>

## LICENSE

MIT © 2019-2020 Shixiang Wang, Xue-Song Liu

The software is made available for non commercial research purposes only
under the [MIT](LICENSE.md). However, notwithstanding any provision of
the MIT License, the software currently may not be used for commercial
purposes without explicit written permission after contacting Shixiang
Wang <wangshx@shanghaitech.edu.cn> or Xue-Song Liu
<liuxs@shanghaitech.edu.cn>.

-----

**[Cancer Biology Group](https://github.com/XSLiuLab) @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**
