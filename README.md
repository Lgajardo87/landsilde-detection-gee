# Landslide Detection in Google Earth Engine

This repository contains the Google Earth Engine (GEE) script developed for:

**Rainfall-Induced Landslide Detection in the Central–Southern Andes of Chile Using Integrated SAR, Optical, and Machine Learning Approaches in Google Earth Engine**

## Description

The script implements a multi-sensor methodological framework for rainfall-induced landslide detection using:

- Sentinel-1 SAR data (backscatter changes in VV and VH)
- Sentinel-2 optical imagery (NDVI changes)
- Shuttle Radar Topography Mission (SRTM) derived slope and elevation
- Random Forest supervised classification

## Usage

1. Open the script `code/landslide_detection_GEE_RF.js`
2. Copy the script into the Google Earth Engine Code Editor
3. Define the study area
4. Set pre- and post-event periods
5. Run

## Requirements

- Google Earth Engine account

## Methodological Background

This implementation is partially adapted from:

> Lindsay, E., et al. (2022). *Multi-Temporal Satellite Image Composites in Google Earth Engine for Improved Landslide Visibility: A Case Study of a Glacial Landscape.*  
> *Remote Sensing*, 14(10), 2301.  
> https://doi.org/10.3390/rs14102301

## License

This project is licensed under the **MIT License** – see the `LICENSE` file for details.

## Citation

If you use this code, please cite the associated article:

Gajardo, L., Jaque, E., Lillo, M., and Castro, F. (2025).  
*Rainfall-Induced Landslide Detection in the Central–Southern Andes of Chile Using Integrated SAR, Optical, and Machine Learning Approaches in Google Earth Engine.*  
*IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing*. DOI: [to be added upon publication].
