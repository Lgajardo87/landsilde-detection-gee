/*******************************************************
Title: Landslide Detection Using SAR, Optical, and Random Forest
Platform: Google Earth Engine
Author: Luis Gajardo et al.,

Description:
This script implements a multi-sensor methodology for
rainfall-induced landslide detection using Sentinel-1,
Sentinel-2, and Random Forest classification.

Note:
This implementation is partially adapted from the
methodological framework described in:

Lindsay, E., et al. (2022).
“Multi-Temporal Satellite Image Composites in Google Earth Engine
for Improved Landslide Visibility: A Case Study of a Glacial Landscape.”
Remote Sensing, 14(10), 2301.
https://doi.org/10.3390/rs14102301

The present script was independently developed and modified
to meet the objectives of this study.
*******************************************************/

// ===================
// 1. INITIAL DEFINITIONS
// ===================

// Define  pre- and post-event dates (SAR)


        // PRE-EVENT IMAGE
        // Dates over which to create a composite.
        var pre_start = ee.Date('2022-06-01');
        var pre_end = ee.Date('2023-06-20');
        
        // POST-EVENT IMAGE
        // Dates over which to create a median composite.
        var post_start = ee.Date('2023-06-25');
        var post_end = ee.Date('2024-02-02');
        
// Set the name of the google drive-folder where you want the image to be exported to.
        
        var my_google_drive_folder = "earthengine";


// Define  pre- and post-event dates (S2)

        // PRE-EVENT IMAGE
        var startDatePreEvent = '2023-01-01';
        var endDatePreEvent = '2023-04-01';

        // POST-EVENT IMAGE
        var startDatePostEvent = '2023-10-01';
        var endDatePostEvent = '2023-12-30';


// ===================
// 2. SENTINEL-2 PROCESSING
// ===================

// Apply cloud maskin using QA60

      function maskS2cloudsQA60(image) {
        var cloudBitMask = 1 << 10;  // nubes
        var cirrusBitMask = 1 << 11; // cirros

        var qa = image.select('QA60');
        var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                    qa.bitwiseAnd(cirrusBitMask).eq(0));
        return image.updateMask(mask).copyProperties(image, image.propertyNames());
      }


// Create median composites:

// Pre-event S2 image (filter + mask + median)


        var S2PreEvent = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
         .filterDate(startDatePreEvent, endDatePreEvent)
         .filterBounds(geometry)
         .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)) // Filtro por porcentaje de nubes
         .map(maskS2cloudsQA60)
         .median()
         .clip(geometry);


// Post-event S2 image (filter + mask + median)

        var S2PostEvent = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
         .filterDate(startDatePostEvent, endDatePostEvent)
         .filterBounds(geometry)
         .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)) // Filtro por porcentaje de nubes
         .map(maskS2cloudsQA60)
         .median()
         .clip(geometry);


// Calculate pre- and post-event NDVI

        var ndviPreEvent = S2PreEvent.normalizedDifference(['B8', 'B4']).rename('NDVI_Pre');
        var ndviPostEvent = S2PostEvent.normalizedDifference(['B8', 'B4']).rename('NDVI_Post');

// Compute NDVI difference 

        var ndviDiff = ndviPostEvent.subtract(ndviPreEvent).rename('NDVI_Diff');

// Compute NDVI percent change

        var ndviPctChange = ndviDiff.divide(ndviPreEvent.where(ndviPreEvent.eq(0), 0.0001))
                           .multiply(100)
                           .rename('NDVI_Pct_Change');

// Create NDVI change binary mask (threshold = 0.2)

        var changeThreshold = 0.2;
        var ndviChangeMask = ndviDiff.abs().gt(changeThreshold).rename('NDVI_Change_Mask');

// ===================
// 3. TOPOGRAPHIC VARIABLES
// ===================

// Set DEM 

        //Global: NASA SRTM Digital Elevation 30m
        var DEM = ee.Image('USGS/SRTMGL1_003')


var demaoi = DEM.clip(geometry);  // geometry = tu AOI

// Calculate terrain slope

var slope = ee.Terrain.slope(demaoi);


// ####################### Calculate Geometry #################################

var region = (geometry);


// ===================
// 4. SENTINEL-1 PROCESSING
// ===================

// Filter pre-event collections (VV & VH, IW mode) 

var pre_data = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(region)
  .filterMetadata('transmitterReceiverPolarisation','equals',["VV", "VH"])
  .filterMetadata('instrumentMode','equals','IW')
  .filterDate(pre_start, pre_end)

// Filter pre-event collections (VV & VH, IW mode)

var post_data = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(region)
  .filterMetadata('transmitterReceiverPolarisation','equals',["VV", "VH"])
  .filterMetadata('instrumentMode','equals','IW')
  .filterDate(post_start, post_end)

// Print orbit direction statistics

var pre_s = pre_data.aggregate_histogram('transmitterReceiverPolarisation')
print('size of pre collection', pre_data.size())

var post_s = post_data.aggregate_histogram('transmitterReceiverPolarisation')
print('size of post collection', post_data.size())
/**
*  Count Frequency of Ascending and Descending Orbital Passes
*/
var identify_pre  = ee.Dictionary(pre_data.aggregate_histogram('orbitProperties_pass'))
print('Number of images in each orbital direction (pre)', identify_pre )

var identify_post  = ee.Dictionary(post_data.aggregate_histogram('orbitProperties_pass'))
print('Number of images in each orbital direction (post)', identify_post )

  // Apply topographic correction :
/**
* Radiometric slope correction algorithm for topographic correction
* Author: Andreas Vollrath, described in https://doi.org/10.3390/rs12111867
*/
var slope_correction = function (collection,
                                options
                                ){

    // set defaults if undefined options
    options = options || {};
    var model = options.model || 'volume';
    var elevation = options.elevation || DEM
    // var elevation = options.elevation || ee.Image(DEM);  
    var buffer = options.buffer || 0;

    // we need a 90 degree in radians image for a couple of calculations
    var ninetyRad = ee.Image.constant(90).multiply(Math.PI/180);

    // Volumetric Model Hoekman 1990
    function _volume_model(theta_iRad, alpha_rRad){

      var nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan();
      var denominator = (ninetyRad.subtract(theta_iRad)).tan();
      return nominator.divide(denominator);
    }

    // surface model Ulander et al. 1996
    function _surface_model(theta_iRad, alpha_rRad, alpha_azRad){

      var nominator = (ninetyRad.subtract(theta_iRad)).cos();
      var denominator = alpha_azRad.cos()
        .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos());
      return nominator.divide(denominator);
    }

    // buffer function (thanks Noel)
    function _erode(img, distance) {

      var d = (img.not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()));

      return img.updateMask(d.gt(distance));
    }

    // calculate masks
    function _masking(alpha_rRad, theta_iRad, proj, buffer){

        // layover, where slope > radar viewing angle
        var layover = alpha_rRad.lt(theta_iRad).rename('layover');

        // shadow
        var shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow');

        // combine layover and shadow
        var mask = layover.and(shadow);

        // add buffer to final mask
        if (buffer > 0)
            mask = _erode(mask, buffer);

        return mask.rename('no_data_mask');
  }

    function _correct(image){

        // get image geometry and projection
        var geom = image.geometry();
        var proj = image.select(1).projection();

        // get look direction angle
        var heading = (ee.Terrain.aspect(
            image.select('angle')).reduceRegion(ee.Reducer.mean(), geom, 1000).get('aspect')
            );

        // Sigma0 to Power of input image
        var sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0));

        // Radar geometry
        var theta_iRad = image.select('angle').multiply(Math.PI/180).clip(geom);
        var phi_iRad = ee.Image.constant(heading).multiply(Math.PI/180);

        // Terrain geometry
        var alpha_sRad = ee.Terrain.slope(elevation).select('slope')
            .multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);
        var phi_sRad = ee.Terrain.aspect(elevation).select('aspect')
            .multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);

        // Model geometry

        //reduce to 3 angle
        var phi_rRad = phi_iRad.subtract(phi_sRad);

        // slope steepness in range
        var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();

        // slope steepness in azimuth
        var alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();

        // Gamma_nought
        var gamma0 = sigma0Pow .divide(theta_iRad.cos());

              // models
        if (model == 'volume')
          var corrModel = _volume_model(theta_iRad, alpha_rRad);

        if (model == 'surface')
          var corrModel = _surface_model(theta_iRad, alpha_rRad, alpha_azRad);

        if (model == 'direct')
          var corrModel = _direct_model(theta_iRad, alpha_rRad, alpha_azRad);

        // apply model to derive gamma0_flat
        var gamma0_flat = gamma0.divide(corrModel);

        // transform to dB-scale
        var gamma0_flatDB = (ee.Image.constant(10)
            .multiply(gamma0_flat.log10()).select(['VV', 'VH'])
            );

        // get Layover/Shadow mask
        var mask = _masking(alpha_rRad, theta_iRad, proj, buffer);

        // return gamma_flat plus mask
        return gamma0_flatDB.addBands(mask).copyProperties(image);


    }

    // run correction function and return corrected collection
    return collection.map(_correct);

};

// export function
exports.slope_correction = slope_correction;

// ####################### Create terrain corrected Pre- and Post- Event images #################################

var pre_corrected = slope_correction(pre_data)                          // Apply terrain correction
  .map(function(im) {return im.updateMask(im.select('no_data_mask'))})  // Apply no data mask
  .reduce(ee.Reducer.mean())                                            // Take the mean of the image collection, to make a single image
  .rename(pre_data.first().bandNames())
  .clip(region)
  
var post_corrected = slope_correction(post_data)
  .map(function(im) {return im.updateMask(im.select('no_data_mask'))}) // Apply no data mask
  .reduce(ee.Reducer.mean())
  .rename(post_data.first().bandNames())
  .clip(region)
  
var pre_VV = pre_corrected.select('VV').rename('preVV');
var pre_VH = pre_corrected.select('VH').rename('preVH');

var post_VV = post_corrected.select('VV').rename('postVV');
var post_VH = post_corrected.select('VH').rename('postVH');


// Calculate difference images (post - pre)


var diff_VV = (post_VV.subtract(pre_VV)).rename('diffVV');
var diff_VH = (post_VH.subtract(pre_VH)).rename('diffVH');

// Binarize differences (threshold = -0.8 dB)

var thresholdVV = -0.8;
var thresholdVH = -0.8;

var diffVV_binary = diff_VV.lte(thresholdVV).rename('diffVV_binary'); // lte: Less Than or Equal (menor o igual)
var diffVH_binary = diff_VH.lte(thresholdVH).rename('diffVH_binary');
//var diffVV_binary = diff_VV.abs().gt(thresholdVV).rename('diffVV_binary');
//var diffVH_binary = diff_VH.abs().gt(thresholdVH).rename('diffVH_binary');



// ===================
// 5. FEATURE STACK CREATION
// ===================

//var selectedBands = ['diff_VV', 'diff_VH','DEM', 'Slope','ndviDiff']; 
var selectedBands = ['diff_VV', 'diff_VH','DEM', 'Slope','ndviPctChange']; 

//  Create multi-band image
var allBands = diff_VV.rename('diff_VV')
  .addBands(diff_VH.rename('diff_VH'))
  .addBands(demaoi.rename('DEM'))
  .addBands(slope.rename('Slope'))
  .addBands(ndviDiff.rename('ndviDiff'))
  .addBands(ndviPctChange.rename('ndviPctChange'));

// Select relevant input features
var inputFeatures = allBands.select(selectedBands);

// ===================
// 6. TRAINING DATA GENERATION
// ===================

// Import landslide polygons

var trainingPolygons = ee.FeatureCollection(landslide);

//  Randomly split polygons: 80% train, 20% test

var trainingPolygonsRandom = trainingPolygons.randomColumn('random', 123);
var split = 0.8;
var trainingPolygonsSet = trainingPolygonsRandom.filter(ee.Filter.lt('random', split));
var testPolygonsSet = trainingPolygonsRandom.filter(ee.Filter.gte('random', split));

// Generate random points inside polygons (class = 1)

//Training
var numPointsPerPolygon = 50; 
var landslidePointsTrain = trainingPolygonsSet.map(function(f) {
  return ee.FeatureCollection.randomPoints({
    region: f.geometry(),
    points: numPointsPerPolygon,
    seed: 123
  }).map(function(pt) {
    return pt.set('class', 1);
  });
}).flatten();

// Test
var landslidePointsTest = testPolygonsSet.map(function(f) {
  return ee.FeatureCollection.randomPoints({
    region: f.geometry(),
    points: numPointsPerPolygon,
    seed: 456
  }).map(function(pt) {
    return pt.set('class', 1);
  });
}).flatten();

// Generate random points outside polygons (class = 0)

//Training
var randomPointsTrain = ee.FeatureCollection.randomPoints({
  region: geometry,
  points: 1000,
  seed: 789
}).filter(ee.Filter.notNull(['system:index']))
  .filterBounds(geometry)
  .filter(ee.Filter.intersects({
    leftField: '.geo',
    rightValue: trainingPolygonsSet.geometry(),
    maxError: 1
  }).not())
  .map(function(f) { return f.set('class', 0); });

// Test
var randomPointsTest = ee.FeatureCollection.randomPoints({
  region: geometry,
  points: 500,
  seed: 987
}).filter(ee.Filter.notNull(['system:index']))
  .filterBounds(geometry)
  .filter(ee.Filter.intersects({
    leftField: '.geo',
    rightValue: testPolygonsSet.geometry(),
    maxError: 1
  }).not())
  .map(function(f) { return f.set('class', 0); });

// Merge training and test points

var trainingData = landslidePointsTrain.merge(randomPointsTrain);
var testData = landslidePointsTest.merge(randomPointsTest);

// Sample bands at training and test locations

var trainingSample = inputFeatures.sampleRegions({
  collection: trainingData,
  properties: ['class'],
  scale: 30,
  geometries: true
});

//var testSample = inputFeatures.sampleRegions({
 // collection: testData,
 // properties: ['class'],
 // scale: 30,
  //geometries: true
//});


// ===================
// 7. RANDOM FOREST CLASSIFICATION
// ===================

//Train RF model with training samples

var classifier = ee.Classifier.smileRandomForest({
  numberOfTrees: 100,
  seed: 123
}).train({
  features: trainingSample,
  classProperty: 'class',
  inputProperties: inputFeatures.bandNames()
});

// Classify input image

var classified = inputFeatures.classify(classifier);

//Clean classification result

// Focal Majority (suavizado)
//var cleanedClassified = classified.focal_mode({
  //radius: 1,
  //units: 'pixels'
//});

// Remove small objects (< 30 connected pixels)
var connected = classified.connectedPixelCount(100, true);
var maskSmallObjects = connected.gte(30);
var cleanedClassified = classified.updateMask(maskSmallObjects);

// ===================
// 8. MODEL VALIDATION
// ===================

// Apply classifier to test samples

var predicted = classified.sampleRegions({
  collection: testData,
  properties: ['class'],
  scale: 30
});


// Check the number of samples
print('Training samples: ', trainingSample.size());
print('Class 1 quantity:', trainingSample.filter(ee.Filter.eq('class', 1)).size());
print('Class 0 quantity:', trainingSample.filter(ee.Filter.eq('class', 0)).size());

// Assess feature importance
var importance = ee.Dictionary(classifier.explain().get('importance'));
print('Assess feature importance:', importance);


// Compute confusion matrix and metrics

var confusionMatrix = predicted.errorMatrix('class', 'classification');

//var matrixArray = confusionMatrix.array();
//print('Confusion Matrix', matrixArray);


// Accuracy, Kappa

print('Confusion Matrix:', confusionMatrix);
print('Overall Accuracy:', confusionMatrix.accuracy());
print('Kappa Index:', confusionMatrix.kappa());



// Recall, Precision, Specificity, F1-score


var TN = confusionMatrix.array().get([0, 0]);
var FP = confusionMatrix.array().get([0, 1]);
var FN = confusionMatrix.array().get([1, 0]);
var TP = confusionMatrix.array().get([1, 1]);

// Sensitivity (Recall)
var recall = ee.Number(TP).divide(ee.Number(TP).add(FN));
print('Sensitivity (Recall) Class 1:', recall);

// Precision 
var precision = ee.Number(TP).divide(ee.Number(TP).add(FP));
print('Precision Class 1:', precision);

// Specificity for class 0
var specificity = ee.Number(TN).divide(ee.Number(TN).add(FP));
print('Specificity Class 0:', specificity);


// F1-Score
var f1Score = ee.Number(2).multiply(precision).multiply(recall)
  .divide(precision.add(recall));
print('F1-Score Class 1:', f1Score);


// ===================
// 9. VISUALIZATION
// ===================

Map.centerObject(trainingPolygons, 12);
Map.addLayer(classified, {min: 0, max: 1, palette: ['white', 'red']}, 'RF Classification');
Map.addLayer(cleanedClassified, {min: 0, max: 1, palette: ['white', 'red']}, 'Cleaned Classification');
Map.addLayer(ndviChangeMask.updateMask(ndviChangeMask), {palette: ['black', 'yellow']}, 'NDVI Change Mask (>0.2)');
Map.addLayer(diffVV_binary, {min:0, max:1, palette:['white', 'red']}, 'Binarized VV Change');
Map.addLayer(diffVH_binary, {min:0, max:1, palette:['white', 'blue']}, 'Binarizad VH Change');
Map.addLayer(trainingPolygonsSet, {color: 'orange'}, 'Training Polygons');



// ===================
// 10. EXPORT TO GOOGLE DRIVE
// ===================


Export.image.toDrive({
  image: cleanedClassified,
  description: 'cleanedClassified',
  folder: 'earthengine',
  fileNamePrefix: 'clasificacion_rf',
  region: geometry,
  scale: 10,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: diffVV_binary,
  description: 'diffVV_binary',
  scale: 10,
  folder: my_google_drive_folder,
  region: region
});

Export.image.toDrive({
  image: diffVH_binary,
  description: 'diffVH_binary',
  scale: 10,
  folder: my_google_drive_folder,
  region: region
  
});


Export.image.toDrive({
  image: ndviChangeMask,
  description: 'NDVI_ChangeMask',
  folder: 'EarthEngineExports',
  fileNamePrefix: 'NDVI_ChangeMask',
  region: geometry,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


Export.image.toDrive({
  image: S2PreEvent,
  description: 'S2_PreEvent',
  folder: 'EarthEngineExports',
  fileNamePrefix: 'S2_PreEvent',
  region: geometry,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: S2PostEvent,
  description: 'S2_PostEvent',
  folder: 'EarthEngineExports',
  fileNamePrefix: 'S2_PostEvent',
  region: geometry,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

