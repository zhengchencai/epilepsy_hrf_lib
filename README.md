# The Hemodynamic Response Function Varies Across Anatomical Location and Pathology in the Epileptic Brain (Code)

Zhengchen Cai¹, Nicolás von Ellenrieder¹, Thaera Arafat¹, Hui Ming Kho², Gang Chen³, Andreas Koupparis⁴, Chifaou Abdallah¹, Roy Dudley⁵, Dang Khoa Nguyen⁶, Jeffery Hall¹, Francois Dubeau¹, Jean Gotman¹*, Boris Bernhardt¹*

¹ The Neuro (Montreal Neurological Institute-Hospital), McGill University, Montreal, Canada  
² Department of Neurosurgery, Osaka University Graduate School of Medicine, Suita, Japan  
³ Scientific and Statistical Computing Core, National Institute of Mental Health, Bethesda, USA  
⁴ The Cyprus Institute of Neurology and Genetics, Nicosia, Cyprus  
⁵ Montreal Children's Hospital, McGill University, Montreal, Canada  
⁶ Centre Hospitalier de l'Université de Montréal, Montréal, Canada  

\* These authors jointly supervised this work and share senior authorship

## Overview

This repository contains the analysis code for the manuscript.

1. **HRF deconvolution**: AFNI-based anatomical and functional preprocessing and HRF estimation using TENT basis functions
   - Location: `HRF_deconvolution/`
   - Image processing and HRF deconvolution method code adapted from Chen et al., NeuroImage 2023, [GitHub repository](https://github.com/afni/apaper_hrf_profiles)

2. **HRF library construction**: Building a library of representative HRF shapes from epileptic brain activity
   - Location: `HRF_lib/`
   - HRF time decomposition method code adapted from Kay et al., Nat Methods 2020, [GitHub repository](https://github.com/cvnlab/TDM) and Allen et al., Nat Neuroscience 2022, [GitHub repository](https://github.com/cvnlab/nsddatapaper)

3. **Bayesian modeling**: Modeling the effect of anatomical location and pathology on HRF features
   - Location: `Bayesian_model/`
   - Quarto-rendered HTML report reproducing manuscript figures including data preparation, hierarchical Bayesian modeling, marginal effect size visualization and model comparison and classification.

## Repository Structure

```
├── HRF_deconvolution/
│   ├── checkOblique.py           # check and correct oblique anatomical T1w images
│   ├── SSwarper.py               # skull-stripping and normalization using @SSwarper
│   ├── FSSUMA.py                 # FreeSurfer and SUMA processing scripts
│   ├── generate_afni_proc.py     # AFNI preprocessing pipelines
│   ├── generate_3dMSS.py         # 3dMSS HRF smoothing scripts
│   └── postAFNI.py               # prepare data for 3dMSS (called by generate_3dMSS.py)
│
├── HRF_lib/
│   ├── HRF_PCA.m                 # perform PCA on extracted HRFs
│   ├── HRF_PCA_plot.m            # visualization of PCA variance across parameters
│   ├── convert2hrflib.m          # build HRF library and extract features
│   ├── manifold.m                # manifold learning and density estimation
│   └── plot_manifold.R           # visualization of HRF manifold
│
├── Bayesian_model/
│   └── EpilepsyHRF_Cai_Supplementary_Code.html  # Quarto-rendered HTML report with code
│
├── hrflib18.csv                  # HRF library 
└── HRFbasis_500Hz.csv            # 3 PCA basis functions at 500Hz (same as EEG)
```

## License

This code is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).


## Contact

For questions or issues, please open an issue on GitHub or contact zhengchen.cai@mcgill.ca.

