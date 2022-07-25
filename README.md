# flashfm-ivis

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/461200048.svg)](https://zenodo.org/badge/latestdoi/461200048)
<!-- badges: end -->

[*flashfm-ivis*](http://shiny.mrc-bsu.cam.ac.uk/apps/flashfm-ivis/): interactive visualisation for fine-mapping of multiple quantitative traits.

Also available, [*finemap-ivis*](http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/) for interactive visualisation of single-trait fine-mapping results.

The **quick instruction README file** of this tool, its source code and pre-loaded dataset [can be downloaded publicly here](https://drive.google.com/drive/folders/167XScfXLekP5UmJcmM0_YK2bp0tCKcYW). Please Note: The latest updated *Version-2022-05-30* requires a few pre-loaded large datasets (i.e.: approx. 600MB), therefore we saved all source code and datasets in the [public online Google Drive](https://drive.google.com/drive/folders/167XScfXLekP5UmJcmM0_YK2bp0tCKcYW) for users to download easily. Users just need download the whole folders into the local machine and run the *app.R* (i.e.: Rshiny file) on the [RStudio](https://www.rstudio.com). 

Accompanying paper:
> *flashfm-ivis: interactive visualisation for fine-mapping of multiple quantitative traits* <br />
> Feng Zhou, Adam S Butterworth, Jennifer L Asimit <br />
> The paper and its supplementary material are published in *Bioinformatics*, https://doi.org/10.1093/bioinformatics/btac453

For more information, visit:
> flashfm-ivis: [http://shiny.mrc-bsu.cam.ac.uk/apps/flashfm-ivis/](http://shiny.mrc-bsu.cam.ac.uk/apps/flashfm-ivis/)   
> finemap-ivis: [http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/](http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/) 

## Overview
*flashfm-ivis* provides a suite of interactive visualisation plots to view potential causal genetic variants that underlie associations that are shared or distinct between multiple quantitative traits and compares results between single- and multi-trait fine-mapping. Unique features include network diagrams that show joint effects between variants for each trait and regional association plots that integrate fine-mapping results, all with user-controlled zoom features for an interactive exploration of potential causal variants across traits.

## Availability and Implementation
*flashfm-ivis* is an open-source software under the MIT license. It is available as an interactive web-based tool [http://shiny.mrc-bsu.cam.ac.uk/apps/flashfm-ivis/](http://shiny.mrc-bsu.cam.ac.uk/apps/flashfm-ivis/) and as an R package. Code and documentation are available at [https://github.com/fz-cambridge](https://github.com/fz-cambridge) and [https://github.com/jennasimit](https://github.com/jennasimit). Additional features can be downloaded as standalone R libraries to encourage reuse. 

## Installation
Users can directly use the online version of these tools without installing any packages or previous experience in programming languages. If users would like to run analyses in their local machines, then please download the codes from the [public online Google Drive](https://drive.google.com/drive/folders/167XScfXLekP5UmJcmM0_YK2bp0tCKcYW) and open the Rshiny app in the local [RStudio](https://www.rstudio.com).

## Authors
   - [Jennifer Asimit](https://www.mrc-bsu.cam.ac.uk/people/in-alphabetical-order/a-to-g/jennifer-asimit/) (University of Cambridge)
   - [Adam Butterworth](https://www.phpc.cam.ac.uk/people/ceu-group/ceu-senior-academic-staff/adam-butterworth/) (University of Cambridge)
   - [Feng Zhou](https://www.mrc-bsu.cam.ac.uk/people/in-alphabetical-order/t-to-z/feng-zhou/) (University of Cambridge)

## Acknowledgements
The authors thank Jana Soenksen and Inês Barroso, both from the Exeter Centre of Excellence for Diabetes Research (EXCEED), University of Exeter Medical School, Exeter, UK, for their valuable feedback and suggestions.

## Further references
> *FINEMAP: efficient variable selection using summary data from genome-wide association studies.* <br />
> Benner, C., Spencer, C. C., Havulinna, A. S., Salomaa, V., Ripatti, S., & Pirinen, M.(2016). <br />
> Bioinformatics, 32(10), 1493-1501. [https://doi.org/10.1093/bioinformatics/btw018](https://doi.org/10.1093/bioinformatics/btw018).

> *The flashfm approach for fine-mapping multiple quantitative traits.* <br />
> Hernandez, N., Soenksen, J., Newcombe, P., Sandhu, M., Barroso, I., Wallace, C., Asimit, J.L. (2021).<br />
> Nat Commun 12, 6147. [https://doi.org/10.1038/s41467-021-26364-y](https://doi.org/10.1038/s41467-021-26364-y).

> *Fine-mapping genetic associations.* <br />
> Hutchinson, A., Asimit, J., Wallace, C. (2020). <br />
> Human Molecular Genetics, Volume 29, Issue R1, Pages R81–R88. [https://doi.org/10.1093/hmg/ddaa148](https://doi.org/10.1093/hmg/ddaa148).

> *JAM: a scalable Bayesian framework for joint analysis of marginal SNP effects.* <br />
> Newcombe, P. J., Conti, D. V. & Richardson, S. (2016). <br />
> Genet. Epidemiol. 40, 188–201. [https://doi.org/10.1002/gepi.21953](https://doi.org/10.1002/gepi.21953).


