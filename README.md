# Taper Curve Interpolation and Stem Volume Modelling

## Project Overview

This exercise uses stem analysis data from plantation trees to reconstruct taper curves and estimate stem volume under a hierarchical data structure. The workflow combines spline interpolation of diameter measurements along the stem, numerical integration for volume estimation, and mixed-effects modelling to account for correlations within compartments and plots. 

## Objectives

- Reconstruct tree taper curves from stem analysis measurements.
- Estimate stem volume from interpolated taper curves.
- Compare manual numerical integration with R-based integration.
- Fit and evaluate fixed-effects and mixed-effects stem volume models.
- Assess model fit and reliability using diagnostic plots and summary statistics. 

## Data

- `Stem_analysis_data.txt`
- Variables: `Compt`, `Plot`, `Tree`, `ht.st.m`, `d.st.cm`, `d13.cm`, `ht.m`, `dl.cm`, `rel.h`, `hl.m` 

## Methodology

### Stem taper reconstruction
The stem analysis measurements were ordered by compartment, plot, and tree. Monotone spline interpolation was then used to reconstruct taper curves from diameter measurements along the stem. The exercise code also illustrates why monotone interpolation is preferable to overly flexible cubic splines for taper data. 

### Volume estimation
For each tree, taper curves were integrated numerically to estimate stem volume in dm³. The script also evaluates volume by integrating a spline of cross-sectional area, allowing a comparison between two integration approaches. 

### Modelling
The derived tree-level volume values were then used in regression and mixed-effects modelling. The exercise diary shows that several candidate models were compared, and that the mixed-effects model `f2.lmem` produced the best overall performance based on AIC, BIC, bias, and RMSE. :contentReference[oaicite:13]{index=13}

## Main Outputs

- Tree-specific taper curves
- Interpolated stem profiles
- Tree-level stem volume estimates
- Volume-versus-diameter and volume-versus-height plots
- Model comparison tables
- Fit and reliability statistics 

## Key Results

The diary reports that the mixed-effects model `f2.lmem` achieved the best fit, with the lowest bias and RMSE among the tested models. This makes the exercise a strong example of volume modelling under hierarchical structure, which is directly relevant to forest inventory applications. :contentReference[oaicite:15]{index=15}

## Skills Demonstrated

- Stem analysis
- Taper-curve interpolation
- Numerical integration
- Tree volume estimation
- Mixed-effects modelling
- Hierarchical forestry data analysis
- Model comparison and validation
- R-based graphical diagnostics 
