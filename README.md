# QuAPPro

### An R/Shiny web application for Quantification and Alignment of Polysome Profiles

The **QuAPPro** Shiny is an interactive web app for <ins>Qu</ins>antification and <ins>A</ins>lignment of <ins>P</ins>olysome <ins>Pro</ins>files as described by [Schiller et al. (BMC Bioinformatics 2026)](https://link.springer.com/article/10.1186/s12859-026-06379-2).
Polysome profiling is a key technique in the field of mRNA translation research. Ribosomal subunits are separated from monosomes and polysomes by ultracentrifugation on sucrose density gradients. During elution of the gradients, a UV absorbance profile is recorded, which can be used to study many different aspects of protein biosynthesis. In addition to UV absorbance, it is also possible to record fluorescence and to assess the association of fluorescently tagged proteins with ribosomes or polysomes. QuAPPro allows the import of many different text file formats, so that it can be used independently of the device or software that generated the profiles. 

<img src = "Figure2_Schiller_et_al.png" width="50%" height="50%">
Overview figure from [Schiller et al.](https://link.springer.com/article/10.1186/s12859-026-06379-2)

## Usage
QuAPPro be accessed [here](https://www.umm.uni-heidelberg.de/mi3/biochemie/tools/quappro) as a web app or you can use it locally by downloading and running the provided server version of the R code.
See the manual in the online app for details. 

## Example datasets
All profiles presented by [Schiller et al.](https://link.springer.com/article/10.1186/s12859-026-06379-2) can be found [here](https://github.com/johannaschott/QuAPPro/tree/main/profile_data).

## Remarks
This Shiny App was developed in R version 4.0.5. To run it locally, you need the following additional R packages:
- shiny
- shinythemes
- shinybusy
- colourpicker
- stringr
- colorspace
- shinyalert
- DT
- ggplot2
- ggbeeswarm
- dplyr
- cicerone

## Contact

Comments, suggestions and questions are very welcome.
- [Chiara Schiller](mailto:chiara.schiller@uni-heidelberg.de)
- [Johanna Schott](mailto:Johanna.Schott@medma.uni-heidelberg.de)

## Citation

When you use QuAPPro for a publication, please cite:

Schiller C, Lemmer M, Reitter S, Lehmann JA, Fenzl K and Schott J. QuAPPro: an R shiny app for quantification and alignment of polysome profiles. BMC Bioinformatics 27, 22 (2026). https://doi.org/10.1186/s12859-026-06379-2

