## Texture Analysis with Information Theory

</br>

#### Eduarda C. Chagas, Roger de A. Matos Júnior, Alejandro C. Frery, Heitor S. Ramos and Osvaldo A. Rosso 

</br>

---

### Abstract

We present a new approach for characterization and analysis textures from Information Theory.

The proposal is to receive a texture, to perform the linearization process by applying space-filling curves and through the Bandt & Pompe symbolization, to use the discriminatory power of Information Theory to perform the characterization.

Two analyses were used to validate the methodology.
In the first one, we evaluate how we can use causal descriptors of Information Theory to characterize Brodatz textures.
In the second, we use the plane complexity-entropy together with the alpha estimator in different regions extracted from SAR image textures.


### Datasets

#### - *Brodatz Textures*

For reference, we use an arbitrary subset of 45 Brodatz texture images, available at <a href="http://sipi.usc.edu/database/database.php?volume=textures">http://sipi.usc.edu</a>. 
Images have a dimension of 512 × 512 and in the studied subset 13 of these are equalized versions of histograms, ie versions with contrast alterations of other images present.

#### - *SAR images*

For this analysis, three SAR images with different regions were used:

- Parque Nacional Sierra del Lacandon, Guatemala (acquired April 10, 2015), available at <a href="https://uavsar.jpl.nasa.gov/cgi-bin/product.pl?jobName=Lacand_30202_15043_
006_150410_L090_CX_01#dados">jet propulsion laboratory</a>.

- Cape Canaveral Ocean Regions (acquired September 22, 2016).

- Urban area of the city of Munich, Germany (acquired June 5, 2015).

The images used in this experiment are results from the HHHH SAR band and like the Brodatz texture samples, each SAR sample is represented by a 128 × 128 subimage.

A total of 160 samples were considered during the investigation, 40 samples from each category of regions: Guatemalan forest regions; oceanic regions of Cape Canaveral with behavior 1; Cape Canaveral Behavioral Ocean Regions 2 and Urban City of Munich.

### Methodology

The proposed characterization algorithm for natural texture analysis consists of two modules, which are the linearization of the image intensity matrix data and the representation of the corresponding signals in the Entropy-Complexity Plane.
Linearization was performed with each of the three space filling curves studied (raster-1, raster-2 and hilbert) and the influence of their respective mappings on the Entropy-Complexity plane was analyzed.

The general implementation of texture image characterization can be described according to the following steps:

1. <b>Processing of input data</b>. Due to Hilbert curve constraints, images received as input must have power dimensions of 2, so this restriction must first be verified.
2. <b>Mapping Function Selection</b>. Only one scan method is used in each analysis. The step performs a two-dimensional data transformation to one-dimensional signals.
3. <b>Bandt & Pompe Symbolization</b>. We calculate the probability distribution of the data so that we can apply the quantifiers of Information Theory.
4. <b>Entropy-Complexity Plane</b>. Representing our main characterization tool, it is at this stage that it verifies the behavior of the mapping data dynamics and, consequently, their discriminating power.

![Methodology used in the characterization of SAR image textures](https://github.com/EduardaChagas/PolSAR-from-IT/blob/master/Images/Methodology/MethodologySAR.png)


