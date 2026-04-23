
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gbspatial

<!-- badges: start -->

<!-- badges: end -->

The goal of gbspatial is to provide functions useful for spatial
transcriptomics data analysis.

## Installation

You can install the development version of gbspatial from
[GitHub](https://github.com/) with:

``` r
#devtools::install_github("gbarisano/gbspatial")
library(gbspatial)
```

## Examples

### PlotFOVSegmentation

``` r
PlotFOVSegmentation(seurat_obj = my_spatial_data,
                    fill = "updated_celltype",
                    img_dir= Sys.glob(file.path("data/DecodedFiles","*","*","CellStatsDir","CellComposite")),
                    fov_pos_file = Sys.glob(file.path("data/flatFiles", "*","*fov_positions_file.csv.gz")),
                    alpha=0.8)
#> Loading FOV positions for image 'TMA1' from: data/flatFiles/TMA1/TMA1_fov_positions_file.csv.gz
#> Loading FOV positions for image 'TMA2' from: data/flatFiles/TMA2/TMA2_fov_positions_file.csv.gz
#> Processing Image: TMA1
#>   -> Found 2 unique FOVs in image 'TMA1'.
#>   -> Rendering FOV: 1
#> Warning: Not validating FOV objects
#> Not validating FOV objects
#> Not validating FOV objects
#> Warning: Not validating Seurat objects
#> Not validating Seurat objects
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

    #> Saved: ./FOV_1_TMA1_segmentation.png
    #>   -> Rendering FOV: 2
    #> Warning: Not validating FOV objects
    #> Warning: Not validating FOV objects
    #> Not validating FOV objects
    #> Warning: Not validating Seurat objects
    #> Not validating Seurat objects
    #> Warning in PlotFOVSegmentation(seurat_obj = my_spatial_data, fill =
    #> "updated_celltype", : Could not find FOV 2 in the provided FOV position file
    #> for TMA1 - Guessing the coordinates, which may lead to segmentations being
    #> unaligned to the FOV image

<img src="man/figures/README-unnamed-chunk-3-2.png" alt="" width="100%" />

    #> Saved: ./FOV_2_TMA1_segmentation.png
    #> Processing Image: TMA2
    #>   -> Found 2 unique FOVs in image 'TMA2'.
    #>   -> Rendering FOV: 1
    #> Warning: Not validating Seurat objects
    #> Warning: Not validating FOV objects
    #> Not validating FOV objects
    #> Not validating FOV objects
    #> Warning: Not validating Seurat objects
    #> Warning in PlotFOVSegmentation(seurat_obj = my_spatial_data, fill =
    #> "updated_celltype", : Could not find FOV 1 in the provided FOV position file
    #> for TMA2 - Guessing the coordinates, which may lead to segmentations being
    #> unaligned to the FOV image

<img src="man/figures/README-unnamed-chunk-3-3.png" alt="" width="100%" />

    #> Saved: ./FOV_1_TMA2_segmentation.png
    #>   -> Rendering FOV: 2
    #> Warning: Not validating Seurat objects
    #> Warning: Not validating FOV objects
    #> Not validating FOV objects
    #> Not validating FOV objects
    #> Warning: Not validating Seurat objects

<img src="man/figures/README-unnamed-chunk-3-4.png" alt="" width="100%" />

    #> Saved: ./FOV_2_TMA2_segmentation.png
    #> Pipeline complete.
    knitr::include_graphics("FOV_1_TMA1_segmentation.png")

<img src="FOV_1_TMA1_segmentation.png" alt="" width="100%" />

``` r
knitr::include_graphics("FOV_2_TMA1_segmentation.png")
```

<img src="FOV_2_TMA1_segmentation.png" alt="" width="100%" />

``` r
knitr::include_graphics("FOV_1_TMA2_segmentation.png")
```

<img src="FOV_1_TMA2_segmentation.png" alt="" width="100%" />

``` r
knitr::include_graphics("FOV_2_TMA2_segmentation.png")
```

<img src="FOV_2_TMA2_segmentation.png" alt="" width="100%" />
