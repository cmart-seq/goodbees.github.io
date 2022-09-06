---
title: "gganimate flower"
date: 2022-09-05
draft: false
author: "Carly Martin"
categories: ["Art", "Data Viz"]
tags: ["digital-art", "animation"]
---

I made a version of this gif when I first considered creating a bioinformatics blog.

Two years later, I've finally got goodbees up and running, so to commemorate the occaison here is the code for this fun little guy. 

The animation simply bounces back and forth between a UMAP plot of a subset of root single-cell RNA sequencing [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4626007) from [Shahan et al. 2022](https://www.cell.com/developmental-cell/pdf/S1534-5807(22)00033-8.pdf) and the coordinates of a flower vector image. The number of points in the flower vector matches the number of cells in the RNA-sequencing dataset (more on that later).

With this outline, you should be able to re-create the gif yourself with the data in this github repo, or, make your own from scratch with a single cell dataset and vector image of your choice. 

The inputs are: 

1) a small scRNA-seq dataset (for example, 2796 root cell profiles)
2) a black and white vectorized image (I took a screenshot of a black and white flower drawing and opened it up in Adobe Illustrator, then performed "Rasterize as bitmap" (with transparent background) -> image trace -> save as .png

### Step 1: extract points from vector image

This portion of the tutorial is basically copy-and-pasted from Valentine Svensson's [blog post](https://www.nxn.se/valent/2020/4/7/converting-images-of-line-graphs-to-data), specifically the associated [jupyter notebook](https://github.com/vals/Blog/blob/master/200403-covid-mobility/Parse%20mobility%20screenshots.ipynb).

```python
import numpy as np
import matplotlib.pyplot as plt
import cv2
import csv
import os
```


```python
image = cv2.imread('./animate_umap_botanical/data/flower.png')
```


```python
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
```


```python
# Set threshold level
threshold_level = 50
```


```python
# Find coordinates of all pixels below threshold
coords = np.column_stack(np.where(gray < threshold_level))

print(coords)
```

    [[203 467]
     [203 468]
     [203 469]
     ...
     [678 646]
     [678 647]
     [678 648]]



```python
x = coords[:, 0]  # first column of coords

```


```python
y = coords[:, 1]  # first column of coords

```


```python
plt.scatter(x, y)
plt.show()
```


    
![Cluster identities](docs/img/main/desert_garden.jpg 'Cluster identities')

    



```python
#saving coords
np.savetxt("./animate_umap_botanical/data/flower_coors.csv", coords, delimiter=",")
``` 



