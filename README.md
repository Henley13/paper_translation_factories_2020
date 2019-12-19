
TODO

- [] add bigfish v0
- [] add nucleAIzer
- [] add final notebook
- [] add 2019_racha code and scripts
- [] complete README.md


# Method, plots and results from Chouaib et al. (2020)

This repository gathers the code used to explore and analyze a large part of the images from the following paper: 

__Paper:__ A localization screen reveals translation factories and widespread co-translational protein targeting  
__Authors:__ Racha Chouaib<sup>1,2,3,11,+</sup>, Adham Safieddine<sup>1,2,3,+</sup>, Xavier Pichon<sup>1,2,+</sup>, Arthur Imbert<sup>4,5,6,+</sup>, Oh Sung Kwon<sup>7</sup>, Aubin Samacoits<sup>8,9</sup>, Abdel-Meneem Traboulsi<sup>1,2</sup>, Marie-Cécile Robert<sup>1,2</sup>, Nikolay Tsanov<sup>1,2</sup>, Emeline Coleno<sup>1,2</sup>, Ina Poser<sup>10</sup>, Christophe Zimmer<sup>8,9</sup>, Anthony Hyman<sup>10</sup>, Hervé Le Hir<sup>7</sup>, Kazem Zibara<sup>3,11</sup>, Marion Peter<sup>1,2</sup>, [Florian Mueller](mailto:muellerf.research@gmail.com)<sup>8,9,* </sup>, [Thomas Walter](mailto:thomas.walter@mines-paristech.fr)<sup>4,5,6,* </sup>, [Edouard Bertrand](mailto:edouard.bertrand@igmm.cnrs.fr)<sup>1,2,*</sup>

><sup>1</sup>Institut de Génétique Moléculaire de Montpellier, University of Montpellier, CNRS, Montpellier, France  
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
| ![](images/dapi_2d_all.png "Dapi channel") | ![](images/fish_2d_all.png "FISH channel") |

If you have any question relative to the image acquisition, please contact [Edouard Bertrand](mailto:edouard.bertrand@igmm.cnrs.fr)

## Pipeline

The pipeline is made up of three different resources:
- **BigFISH**, a python library to manipulate FISH images, apply segmentation and detection algorithms, then compute spatial features at the cell-level. Except for nuclei segmentation and final results computation, the full pipeline is based on BigFISH. As the library is not public yet, the actual version used for this paper is directly integrated in this repository.
- [**NucleAIzer**](http://nucleaizer.org/), an online tool for nuclei segmentation. We actually scale it using a modified version of their open-sourced code.
- A more general environment with classic data science libraries to train classification models, perform statistical tests and plot results. 

![](images/pipeline.png "Pipeline")

If you have any question relative to the image analysis, please contact [Florian Mueller](mailto:muellerf.research@gmail.com) or [Thomas Walter](mailto:thomas.walter@mines-paristech.fr).

### 1. Projections

```python
import bigfish.stack as stack

# nuc : np.ndarray, np.uint
#   3-d dapi image with shape (z, y, x).

# focus projection and enhanced contrast for the dapi channel
nuc_focus = stack.focus_projection_fast(
    tensor=nuc,
    proportion=0.7,
    neighborhood_size=7)
nuc_focus = stack.rescale(
    tensor=nuc_focus, 
    channel_to_stretch=0)
```

![](images/dapi_2d_crop.png "zoomed in 2D focus projection of dapi channel")

```python
# cyt : np.ndarray, np.uint
#   3-d FISH image with shape (z, y, x).

# focus projection for FISH channel
cyt_focus = stack.focus_projection_fast(
    tensor=cyt,
    proportion=0.75,
    neighborhood_size=7)

# maximum intensity projection for FISH channel
cyt_in_focus = stack.in_focus_selection(
    image=cyt,
    proportion=0.80,
    neighborhood_size=30)
cyt_mip = stack.maximum_projection(cyt_in_focus)
```

![](images/fish_2d_crop.png "zoomed in 2D maximum intensity projection of FISH channel")


### 2. Filtering

```python
import bigfish.stack as stack
import bigfish.detection as detection

# cyt : np.ndarray, np.uint
#   3-d FISH image with shape (z, y, x).

# get sigma value for gaussian filters
sigma_z, sigma_yx = detection.get_sigma(
    resolution_z=300,
    resolution_yx=103,
    psf_z=350, 
    psf_yx=150)
sigma_log = (sigma_z, sigma_yx, sigma_yx)
sigma_background = (sigma_z*5, sigma_yx*5, sigma_yx*5)

# LoG filter
cyt_filtered_log = stack.log_filter(
    image=cyt, 
    sigma=sigma_log, 
    keep_dtype=True)

# large gaussian filter to estimate then remove noisy background
cyt_filtered_background = stack.remove_background_gaussian(
    image=cyt, 
    sigma=sigma_background)
```

### 3. Nuclei segmentation

![](images/nuc_segmentation.png "Nuclei segmentation")

#### NucleAIzer

#### Two-round segmentation

```python
import bigfish.segmentation as segmentation

# nuc_focus : np.ndarray, np.uint
#   2-d projection of dapi image with shape (y, x).
# nuc_mask : np.ndarray
#   Results of the nuclei segmentation, with shape (y, x).

# remove the nuclei we have already segmented in the image
unsegmented_nuclei = segmentation.remove_segmented_nuc(
    image=nuc_focus,
    mask=nuc_mask)
```

```python
import bigfish.segmentation as segmentation

# nuc_mask_1 : np.ndarray,
#   Results of the first nuclei segmentation, with shape (y, x).
# nuc_mask_2 : np.ndarray,
#   Results of the second nuclei segmentation, with shape (y, x).

# merge results of two successive segmentations
nuc_mask = segmentation.merge_labels(label_1=nuc_mask_1, label_2=nuc_mask_2)
nuc_mask = segmentation.dilate_erode_labels(label=nuc_mask)
```

### 4. Cell segmentation

```python
import bigfish.segmentation as segmentation

# cyt_focus : np.ndarray, np.uint 
#   2-d projection of FISH image with shape (y, x).
# threshold : int
#   Intensity pixel threshold to discriminate foreground (cells) from 
#   background
# nuc_mask : np.ndarray 
#   Result of the nuclei segmentation with shape (y, x).

# compute binary mask
mask = segmentation.build_cyt_binary_mask(
    image_projected=cyt_focus,
    threshold=threshold)
mask[nuc_mask > 0] = True

# compute relief
relief = segmentation.build_cyt_relief(
    image_projected=cyt_focus,
    nuc_labelled=nuc_mask,
    mask_cyt=mask,
    alpha=0.99)

# cell segmentation with watershed-based algorithm
cyt_mask = segmentation.cyt_watershed(
    relief=relief,
    nuc_labelled=nuc_mask,
    mask=mask,
    smooth=7)
```

![](images/cyt_segmentation.png "Cell segmentation")


### 5. mRNAs detection

#### Spot detection

```python
import bigfish.detection as detection

# cyt_filtered_log : np.ndarray, np.uint 
#   LoG filtered image with shape (z, y, x) or (y, x).
# threshold : int
#   Minimum threshold of a spot in the LoG filtered image to be kept.

# get sigma value for gaussian filters
sigma_z, sigma_yx = detection.get_sigma(
    resolution_z=300,
    resolution_yx=103,
    psf_z=350, 
    psf_yx=150)
sigma = (sigma_z, sigma_yx, sigma_yx)

# detect spots
mask_lm = detection.local_maximum_detection(
    image=cyt_filtered_log,
    minimum_distance=2)
spots, radius, _ = detection.spots_thresholding(
    image=cyt_filtered_log,
    sigma=sigma,
    mask_lm=mask_lm,
    threshold=threshold)
```

![](images/spot_detection.png "Spot detection")

#### Cluster decomposition

```python
import bigfish.detection as detection
import numpy as np

# cyt_filtered_background : np.ndarray
#   Image with shape (z, y, x) and filter with a large gaussian operator to 
#   estimate then remove background.
# spots : np.ndarray, np.int64
#   Coordinates of the detected spots with shape (nb_spots, 3).
# radius : Tuple[float]
#   Radius of the detected spots, one value for each dimension.

# cluster decomposition
(spots_out_cluster, spots_in_cluster, _, _) = detection.run_decomposition(
    image=cyt_no_background, 
    spots=spots, 
    radius=radius,
    min_area=2,
    resolution_z=300,
    resolution_yx=103,
    psf_z=350,
    psf_yx=150)
spots = np.concatenate((spots_out_cluster, spots_in_cluster[:, :3]), axis=0)
```

![](images/cluster_decomposition.png "Decomposition of dense mRNAs clusters")

#### Foci detection

```python
import bigfish.detection as detection

# spots : np.ndarray, np.int64
#   Coordinates of the detected spots with shape (nb_spots, 3).

# detect foci
clustered_spots = detection.cluster_spots(
    spots=spots,
    resolution_z=300,
    resolution_yx=103,
    radius=350,
    nb_min_spots=5)
foci = detection.extract_foci(clustered_spots=clustered_spots)
```

![](images/foci_detection.png "Foci detection")

### 6. Postprocessing and cell extraction

```python
import bigfish.stack as stack

# nuc_mask : np.ndarray
#   Results of the nuclei segmentation, with shape (y, x).
# clustered_spots : np.ndarray, np.int64
#   Coordinates of the detected spots with shape (nb_spots, 4). The last 
#   column is the cluster assigned to the spot. If no cluster was assigned,
#   value is -1.
# foci : np.ndarray, np.int64
#   Array with shape (nb_foci, 5). One coordinate per dimension for the 
#   foci centroid (zyx coordinates), the number of spots detected in the
#   foci and its index.

# binarize nuclei masks
binary_nuc_mask = nuc_mask > 0

# spots out of foci and inside foci
spots_out_foci = clustered_spots.copy()
spots_out_foci = spots_out_foci[spots_out_foci[:, 3] == -1, :]
spots_in_foci = clustered_spots.copy()
spots_in_foci = spots_in_foci[spots_in_foci[:, 3] != -1, :]

# remove foci inside nuclei
spots_in_foci_cleaned, foci_cleaned = stack.remove_transcription_site(
    mask_nuc=binary_nuc_mask,
    spots_in_foci=spots_in_foci,
    foci=foci)
```


```python
import bigfish.stack as stack

# nuc_mask : np.ndarray
#   Results of the nuclei segmentation, with shape (y, x).
# cyt_mask : np.ndarray
#   Results of the cell segmentation, with shape (y, x).
# spots_out_foci : np.ndarray, np.int64
#   Coordinate of the spots detected outside foci, with shape (nb_spots, 4). 
#   One coordinate per dimension (zyx coordinates) plus a default index 
#   (-1 for mRNAs spotted outside a foci).
# spots_in_foci_cleaned : np.ndarray, np.int64
#    Coordinate of the spots detected inside foci, with shape (nb_spots, 4). 
#    One coordinate per dimension (zyx coordinates) plus the index of the
#    foci. Spots from the transcription sites are removed.
# foci_cleaned : np.ndarray, np.int64
#   Array with shape (nb_foci, 5). One coordinate per dimension for the 
#   foci centroid (zyx coordinates), the number of spots detected in the
#   foci and its index. Transcription sites free.

# extract coordinates for each identified cell in the image
results = stack.extract_coordinates_image(
    cyt_labelled=cyt_mask,
    nuc_labelled=nuc_mask,
    spots_out=spots_out_foci,
    spots_in=spots_in_foci_cleaned,
    foci=foci_cleaned)

for i_cell, results_cell in enumerate(results):
    cyt_coord, nuc_coord, rna_coord, foci_coord, _ = results_cell
```

![](images/plot_random.png "Cell coordinates")

### 7. Hand-crafted features

```python
import bigfish.classification as classification

# cyt_coord : np.ndarray, np.int64
#   Coordinates of the cytoplasmic membrane (yx coordinates) with shape 
#   (nb_points, 2).
# nuc_coord : np.ndarray, np.int64
#   Coordinates of the nuclear membrane (yx coordinates) with shape 
#   (nb_points, 2).
# rna_coord : np.ndarray, np.int64
#   Coordinates of the mRNAs detected inside the cell, with shape (nb_spots, 4). 
#   One coordinate per dimension (zyx coordinates) plus the index of the
#   a potential foci (-1 if the mRNA is spotted outside a foci).

# get features names
features_name = classification.get_features_name()

# compute spatial features
features_cell = classification.get_features(
    cyt_coord=cyt_coord,
    nuc_coord=nuc_coord,
    rna_coord=rna_coord)
```

### 8. Localization patterns

| Foci | Intranuclear | Nuclear | Perinuclear | Protrusion |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| ![](images/plot_foci.png "Cell with foci pattern") | ![](images/plot_intranuclear.png "Cell with intranuclear pattern") |  ![](images/plot_nuclear.png "Cell with nuclear pattern") | ![](images/plot_perinuclear.png "Cell with perinuclear pattern") | ![](images/plot_protrusion.png "Cell with protrusion pattern") |

### 9. Visualization and results


## Results

![](images/tsne_annotation_legend.png "t-SNE of annotated cells")


![](images/heatmap.png "Proportion of cell classified with a localization pattern")


## References


