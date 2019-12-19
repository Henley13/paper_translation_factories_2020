# Method, plots and results from Chouaib et al. (2020)

This repository gathers the code used to explore and analyze a large part of the images from the following paper: 

>__Paper:__ A localization screen reveals translation factories and widespread co-translational protein targeting  
__Authors:__ Racha Chouaib<sup>1,2,3,11,+</sup>, Adham Safieddine<sup>1,2,3,+</sup>, Xavier Pichon<sup>1,2,+</sup>, Arthur Imbert<sup>4,5,6,+</sup>, Oh Sung Kwon<sup>7</sup>, Aubin Samacoits<sup>8,9</sup>, Abdel-Meneem Traboulsi<sup>1,2</sup>, Marie-Cécile Robert<sup>1,2</sup>, Nikolay Tsanov<sup>1,2</sup>, Emeline Coleno<sup>1,2</sup>, Ina Poser<sup>10</sup>, Christophe Zimmer<sup>8,9</sup>, Anthony Hyman<sup>10</sup>, Hervé Le Hir<sup>7</sup>, Kazem Zibara<sup>3,11</sup>, Marion Peter<sup>1,2</sup>, [Florian Mueller](mailto:muellerf.research@gmail.com)<sup>8,9,* </sup>, [Thomas Walter](mailto:thomas.walter@mines-paristech.fr)<sup>4,5,6,* </sup>, [Edouard Bertrand](mailto:edouard.bertrand@igmm.cnrs.fr)<sup>1,2,*</sup>
>
>><sup>1</sup>Institut de Génétique Moléculaire de Montpellier, University of Montpellier, CNRS, Montpellier, France  
<sup>2</sup>Equipe labélisée Ligue Nationale Contre le Cancer, University of Montpellier, CNRS, Montpellier, France  
<sup>3</sup>ER045, PRASE, DSST, Lebanese University, Beirut, Lebanon  
<sup>4</sup>MINES ParisTech, PSL-Research University, CBIO-Centre for Computational Biology, 77300 Fontainebleau, France  
<sup>5</sup>Institut Curie, 75248 Paris Cedex, France  
<sup>6</sup>INSERM, U900, 75248 Paris Cedex, France  
<sup>7</sup>Institut de biologie de l'Ecole normale supérieure (IBENS), Ecole normale supérieure, CNRS, INSERM, PSL Research University, 46 rue d'Ulm, 75005, Paris, France  
<sup>8</sup>Unité Imagerie et Modélisation, Institut Pasteur and CNRS UMR 3691, 28 rue du Docteur Roux, 75015 Paris; France  
<sup>9</sup>C3BI, USR 3756 IP CNRS – Paris, France  
<sup>10</sup>MPI-CBG, Pfotenhauer Str. 108, 01307 Dresden, Germany  
<sup>11</sup>Biology Department, Faculty of Sciences-I, Lebanese University, Beirut, Lebanon  
>
><sup>+</sup>Equal contributions  
<sup>*</sup>To whom correspondence should be addressed.

This paper provides cell-level qualitative and quantitative evidence about the non-random localization of several mRNAs within the cytoplasm. More specifically it emphasizes the **co-localization of several mRNAs with relative proteins** and the fact that **some mRNAs are translated in factories**.

## Data

Our entire dataset consists in 527 images. Three different channels were used:
- Dapi channel to label nuclei.
- FISH channel to spot mRNA molecules of a targeted gene.
- GFP channel to spot related proteins.

Each image is a 4D tensor (one channel dimension and three spatial dimensions). The following pipeline mainly exploits dapi and FISH channel, although some measures and observations required the GFP channel. Ultimately, we identified 9710 cells within these images for 32 different genes. 

| 2D projection of dapi channel | 2D projection of FISH channel |
| ------------- | ------------- |
| ![](images/dapi_2d_crop.png "Dapi channel") | ![](images/fish_2d_crop.png "FISH channel") |

If you have any question relative to the image acquisition, please contact [Edouard Bertrand](mailto:edouard.bertrand@igmm.cnrs.fr)

## Pipeline

The pipeline is made up of three different resources:
- **BigFISH**, a python library to manipulate FISH images, apply segmentation and detection algorithms, then compute spatial features at the cell-level. Except for nuclei segmentation and final results computation, the full pipeline is based on BigFISH. As the library is not public yet, the actual version used for this paper is directly integrated in this repository.
- [**NucleAIzer**](http://nucleaizer.org/), an online tool for nuclei segmentation. We actually scale it using a modified version of their open-sourced code.
- A more general environment with classic data science libraries to train classification models, perform statistical tests and plot results. 

![](images/pipeline.png "Pipeline")

### Projections



### Filtering



### Nuclei segmentation

![](images/nuc_segmentation.png "Nuclei segmentation")


#### NucleAIzer

#### Two-round segmentation

### Cell segmentation

![](images/cyt_segmentation.png "Cell segmentation")


### mRNAs detection

#### Spot detection

![](images/spot_detection.png "Spot detection")

#### Cluster decomposition

#### Foci detection

![](images/foci_detection.png "Foci detection")

### Postprocessing


### Cell extraction

![](images/plot_random.png "Cell coordinates")

### Hand-crafted features

![](images/plot_topography.png "mRNAs localize at different regions of the cell")

### Localization pattern classification

| Foci | Intranuclear | Nuclear | Perinuclear | Protrusion |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| ![](images/plot_foci.png "Cell with foci pattern") | ![](images/plot_intranuclear.png "Cell with intranuclear pattern") |  ![](images/plot_nuclear.png "Cell with nuclear pattern") | ![](images/plot_perinuclear.png "Cell with perinuclear pattern") | ![](images/plot_protrusion.png "Cell with protrusion pattern") |

### Visualization

![](images/tsne_annotation_legend.png "t-SNE of annotated cells")


![](images/heatmap.png "Proportion of cell classified with a localization pattern")

## Results



## References


