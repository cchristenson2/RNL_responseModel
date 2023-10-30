# RNL_responseModel
### Companion code for publication in Brain Multiphysics
Citation: Chase Christenson, Chengyue Wu, David A. Hormuth, Shiliang Huang, Ande Bao, Andrew Brenner, Thomas E. Yankeelov,
Predicting the spatio-temporal response of recurrent glioblastoma treated with rhenium-186 labelled nanoliposomes,
Brain Multiphysics,
2023,
100084,
ISSN 2666-5220,
https://doi.org/10.1016/j.brain.2023.100084.

Code authors: Chase Christenson, Chengyuye Wu, David A. Hormuth II

The provided MATLAB functions can be used to predict glioblsatoma growth with either the patient-specific, or cohort based methods described in the manuscript. The clinical data is not publically available so replacement examples with similar formatting have been provided.

### Features
- Patient-specific calibration

  - Inputs MRI and segmentations with similar format to processed clinical data

  - Calibrates the provided model (M0 or M1 from manuscript)
$$\frac{\partial N}{\partial t}=\nabla \cdot \left(d\nabla N\right)+k_pN\left(1-\frac{N}{\theta}\right) \tag{M0}$$
$$\frac{\partial N}{\partial t}=\nabla \cdot \left(d\nabla N\right)+k_p(\textbf{x})N\left (1-\frac{N}{\theta}\right) \tag{M1}$$

  - Note: details for radiation coupled models (M2-M9) are available upon request

- Cohort calibration

  - Samples parameters from _in silico_ generated distribution to predict response using M0

### References
(1) Hormuth II DA, Eldridge SL, Weis JA, et al. Mechanically coupled reaction-diffusion model to predict Glioma growth: methodological details. Methods Mol Biol. 2018;1711:225–241.
### License
MIT
