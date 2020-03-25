# LRME-4DCT

## code

MATLAB code for "Lung Respiratory Motion Estimation Based on Kalman Filtering and 4D CT Image Registration"

## Dataset

DIR-lab dataset:https://www.dir-lab.com

POPI dataset: https://www.creatis.insa-lyon.fr/rio/popi-model?action=show&redirect=popi

## Implemented details:

- 4D CT data should be downloaded from the above datasets before implementation.
- In order to speed up the code running speed, the observation value can be pre calculated by the [isoPTV method](https://github.com/visva89/pTVreg) [1].
- For DIR-lab dataset, please start with the file /DIRlab_test/Registration_Kalman_DIR_lab.m
- For POPI dataset, please start with the file /POPI_test/Registration_Kalman_POPI.m

[1]Vishnevskiy V, Gass T, Szekely G, Tanner C, Goksel O. Isotropic total variation regularization of displacements in parametric image registration. IEEE TMI. 2017.
