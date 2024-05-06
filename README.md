# QuAPPro

### An R/Shiny web application for Quantification and Alignment of Polysome Profiles

The **QuAPPro** Shiny is an interactive web app for <ins>Qu</ins>antification and <ins>A</ins>lignment of <ins>P</ins>olysome <ins>Pro</ins>files as described by [Schiller et al. (BioRxiv 2024)](https://www.biorxiv.org/content/10.1101/2024.05.02.592260v1).
Polysome profiling is a key technique in the field of mRNA translation research. Ribosomal subunits are separated from monosomes and polysomes by ultracentrifugation on sucrose density gradients. During elution of the gradients, a UV absorbance profile is recorded, which can be used to measure changes in global translation efficiency or ribosome biogenesis. In addition to UV absorbance, it is also possible to record fluorescence and to assess the association of fluorescently tagged proteins with ribosomes or polysomes. QuAPPro allows the import of many different text file formats, so that it can be used independently of the device or software that generated the profiles. 


![alt text](https://github.com/johannaschott/QuAPPro/blob/main/Figure1.png)



## Usage
QuAPPro be accessed [here](https://www.umm.uni-heidelberg.de/biochemie/shiny/) as a web app or you can use it locally by downloading and running the provided server version of the R code.
See the manual in the online app or [here](https://github.com/johannaschott/QuAPPro/blob/main/QuAPPro_v0-1-0/QuAPPro_v0-1-0_manual.Rmd) for details. 

## Example datasets
All profiles presented by [Schiller et al.](https://www.biorxiv.org/content/10.1101/2024.05.02.592260v1) can be found [here](https://github.com/johannaschott/QuAPPro/tree/main/profile_data).

## Remarks
This Shiny App was developed in R version 4.0.5. To run it locally, you need the following additional R packages:
- shiny
- shinythemes
- colourpicker
- stringr
- colorspace

## Contact

Comments, suggestions and questions are very welcome.
- [Chiara Schiller](mailto:chiara.schiller@uni-heidelberg.de)
- [Johanna Schott](mailto:Johanna.Schott@medma.uni-heidelberg.de)

