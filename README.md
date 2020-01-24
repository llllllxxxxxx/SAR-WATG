# Characterization of SAR Images with Weighted Amplitude Transition Graphs

#### Eduarda C. Chagas, Alejandro C. Frery, Heitor S. Ramos and Osvaldo A. Rosso 

</br>


### This repository contains all the data and code used to develop our research in the related paper submitted to LAGIRS 2020. 

---

#### Abstract

We propose a new technique for SAR image texture characterization based on ordinal pattern transition graphs.
	The proposal consists in
	(i) transforming a 2-D patch of data into a time series using a Hilbert Space Filling Curve,
	(ii) building an Ordinal Pattern Transition Graph with weighted edges;
	(iii) obtaining a probability distribution function from this graph;
	(iv) computing the Entropy and Statistical Complexity of this distribution.
	The weight of the edges is related to the absolute difference of observations.
	This modification takes into account the scattering properties of the target, and leads to a good characterization of several types of textures.
	Experiments with data from Munich urban areas, Guatemala forest regions, and Cape Canaveral ocean samples demonstrate the effectiveness of our technique, which achieves satisfactory levels of separability.

#### The repository is organized as follows:
- `/Code` - the scripts used to develop our research; 
- `/Common` - the BibTex files used in the reports developed; 
- `/Data` - the auxiliary data used during analysis; 
- `/Figures`- Illustrations used in final report; 
- `/Images`- Illustration of the results obtained throughout the research, alongside the methodology files corresponding to our *overview* figure; 
- `/Reports`- the reports developed during the study. 

#### In the code folder, we have the following scripts:
- `Generate_polsar_image.R`- Contains the functions of reading and analysis of SAR image.
- `Plot_SAR_TimeSeries.R`- It executes the plotting of SAR data as one-dimensional data
- `Plot_Transition_graph.R`- It executes the plotting of SAR image results using Hilbert curves and WATG.
- `Theory_Information.R`- Contains the implementation of the information theory descriptors used throughout the work.
- `Transition_graph.R`- It executes the analysis of SAR images using Hilbert curves and WATG.

#### Datasets

For this analysis, three SAR images with different regions were used, available at <a href="https://uavsar.jpl.nasa.gov/cgi-bin">jet propulsion laboratory</a>.:

- Parque Nacional Sierra del Lacandon, Guatemala (acquired April 10, 2015).

- Cape Canaveral Ocean Regions (acquired September 22, 2016).

- Urban area of the city of Munich, Germany (acquired June 5, 2015).

The images used in this experiment are results from the HHHH SAR band and each sample is represented by a 128 × 128 subimage.

A total of 160 samples were considered during the investigation, 40 samples from each category of regions: Guatemalan forest regions; oceanic regions of Cape Canaveral with behavior 1; Cape Canaveral Behavioral Ocean Regions 2 and Urban City of Munich.

#### Methodology

Our procedure consists of the following steps:
	1. Transforming a 2-D patch of data into a time series using a Hilbert Space Filling Curve,
	2. Building an Ordinal Pattern Transition Graph with weighted edges;
	3. Obtaining a probability distribution function from this graph;
	4. Computing the Entropy and Statistical Complexity of this distribution

![Methodology used in the characterization of SAR image textures](https://github.com/EduardaChagas/PolSAR-from-IT/blob/master/Figures/WATG.pdf)

### Software requirements

This code version is tested on the Linux operating system Ubuntu 18.10.

**Installing R 3.6.0 on Ubuntu 18.10**

```sh
$ sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/'
$ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
$ sudo apt update
$ sudo apt install r-base-dev
```

### Installation Guide

Prior to running the script experiments, we need to install the following required packages for R: 

``
install.packages(c('devtools', 'ggplot2', 'ggthemes', 'ggpubr', 'gtools', 'igraph', 'statcomp'))
``

The latest version of these packages was used by October 2019:

```
devtools       2.2.1   
ggplot2        3.2.1       
ggthemes       4.2.0      
ggpubr         0.2     
gtools         3.8.1      
igraph         1.2.4.1       
statcomp       0.0.1.1000      
```

### Running the scripts

Later, you can run the scripts using `Rscript` to call the experiment execution, this will automatically generate the data results and their plots. 

*This repository already has the actual results described in the article, if you are going to generate new ones, we recommend you to delete the original, as you can rollback further, if you need to*

```sh
$ Rscript 
```

This experiment is very computationally expensive and takes time to run, as we evaluate 160 samples of SAR images. 

Finally, if you have any questions or you want to report anaything, feel free to reach me at: eduarda.chagas@dcc.ufmg.br. 





