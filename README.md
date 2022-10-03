# QuAPPro

### An R/Shiny web application for Polysome Profile analysis


The **QuAPPro** Shiny web application allows interacive alignment, quantification and visualization of polysome profiles to study the tranlation status of cells. The app is hosted on a web server and is therefore accessible without locally running R. Profile analysis is performed by selecting a desired feature button and thereupon interacting with the visualized profile to individually select anchors, areas or peaks. For higher reproducibility, the user can apply several helper functions while interacting with the plotting area to precisely find i.e., local minima or maxima of the profiles. We additionally propose an alternative peak quantification approach to surface quantification of peak areas by enabling deconvolution of profiles into their underlying peaks. Generated ready-to-publish plots, quantification values or aligned and normalized profile values can be easily downloaded for further usage.




![alt text](https://github.com/johannaschott/QuAPPro/blob/main/Figure1.png)



## Usage

QuAPPro as online app in your browser can be accessed [here](https://www.umm.uni-heidelberg.de/biochemie/shiny/) or can be run locally by downloading and running the provided server version of the R code.
See the manual in the online app or [here](https://github.com/johannaschott/QuAPPro/blob/main/QuAPPro_manual_v0-0-5.Rmd) for details. 


## Example datasets
You can find example datasets [here](https://github.com/johannaschott/QuAPPro/tree/main/example_profiles).
The "Example_txt_format.txt" shows you an option to format your data in two columns (index and absorbance columns) if it is not coming as .pks file from the Peaktrak software of a Teledyne Isco fractionator or as .csv file from the gradient station from BioComp.


## Remarks
This Shiny App was developed in R version 4.0.5. 
QuAPPro is a work in progress web application as is this github repository.

## Contact

Comments, suggestions and questions are very welcome.
- [Chiara Schiller](mailto:chiara.schiller@stud.uni-heidelberg.de)
- [Johanna Schott](mailto:Johanna.Schott@medma.uni-heidelberg.de)

