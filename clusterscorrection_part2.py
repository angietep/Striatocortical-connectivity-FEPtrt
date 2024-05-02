#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 15:27:41 2024

@author: angeles
"""

import numpy as np
import nibabel as nib
from scipy.ndimage import label

# Load the NIfTI image
nii_img = nib.load('your_nifti_image.nii.gz')
data = nii_img.get_fdata()

# Inverse p-values
inverse_p_values = 1 - data

# Threshold for p < 0.01
thresholded_data = inverse_p_values < 0.01

# Perform clustering using third-nearest neighbor criterion
clusters, num_clusters = label(thresholded_data, structure=np.ones((3, 3, 3)))

# Find clusters with 362 or more voxels
cluster_sizes = np.bincount(clusters.flatten())
large_clusters_indices = np.where(cluster_sizes >= 362)[0]

# Create a new image with only clusters > 362 and their corresponding inverse p-values
new_data = np.zeros_like(data)
for idx in large_clusters_indices:
    if idx == 0:
        continue  # Skip background cluster (label 0)
    cluster_mask = clusters == idx
    new_data[cluster_mask] = inverse_p_values[cluster_mask]

# Create a new NIfTI image with the modified data
new_nii_img = nib.Nifti1Image(new_data, affine=nii_img.affine, header=nii_img.header)

# Save the new NIfTI image
nib.save(new_nii_img, 'output_image_clusters_gt_362.nii.gz')
