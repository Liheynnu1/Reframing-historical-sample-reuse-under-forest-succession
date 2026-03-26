// 所有分类特征的组合
var Feature_merge = ee.Image.cat([FH_feature, environment, NDVI_feature, Greenness, Texture_feature, phenology, sentinel2_spectrum1, sentinel2_glcm, senitnel1_phenology, sentinel_spectrum]);
print("Feature_merge", Feature_merge);
// 植物区系1中二级类特征优选
//  'NDVI_meanL', 'NDVI_speedyL', 'Texture_speedy2023L', 'FH2023', 'FH_speedyL', 'Texture_speedyL', 'FH_stdL',  'FH_medianL', 'Texture_meanL', 'Texture_ZFL', 'NDVI_ZFL', 
// 'Premean', 'Temmean', 'slope', 'aspect', 'PreMS', 'TemMS', 'elevation', 'SMCI'
// 'REPI_meand', 'REPI_maxd','NDI_stdDev',  'mRVI_cv', 'REPI_stdd','VH_stdDev', 'RVI_stdDev', 
// 'B5_contrast',
//  'VH', 'mRVI', 'greenness', 
// 'MCR',  'NDPI', 'greenness', 'elevation'
// 植物区系1中二级类最佳分类组合
var band1 = ['B8', 'B4', 'B2', 'B3', 'greenness', 'MCR', 'NDPI'];
// 植物区系1中三级类最佳分类组合
// var band1 = ['B8', 'B4', 'B2', 'B3', 'greenness', 'NDPI', 'FH2023'];
// 植物区系2中二级类特征优选
// 'NDVI_speedy2023L', 'FH_zfl', 'NDVI_stdL', 'Texture_zf', 'Texture_speed', 'Texture_speedyL', 'FH2023', 'FH_speedy2023L', 'NDVI_mean', 'FH_medain',
// 'REPI_maxQ3d', 'NDVI_shixuw', 'REPI_minQ1d', 'AVE_std','mRVI_sumwinter',
// 'B5_dvar', 'ratio', 'VDDPI_std',
// , 'PreMS', 'elevation',  'Premean',  'slope', 'SMCI'
// 植物区系2中二级类最佳分类组合
// var band1 = ['NDPI', 'IRECI', 'EVI', 'B2', 'SVVI', 'B5', 'NDRE', 'ratio', 'REPI_maxQ3d', 'NDVI_shixuw', 'REPI_minQ1d', 'mRVI_sumwinter', 'slope','elevation', 'NDVI_speedy2023L'];
// 植物区系3中二级类特征优选
// 'NDVI_cvL', 'FH_cvL', 'FH_medianL','NDVI_median','Texture_speedy2023L', 'FH_speedy2023L', 'NDVI_speedyL', 'FH2023'
// 'REPI_mediand', 'REPI_stdd', 'NDI_stdDev', 'mRVI_zf',
// 'elevation', 'slope', 'ratio', 'VDDPI', 'B5_dvar', 
// 植物区系3中二级类最佳分类组合
// var band1 = ['B8A', 'NDPI', 'B11', 'NDRE', 'B2', 'REPI', 'B5_dvar', 'Premean', 'Temmean'];
// 植物区系6中二级类特征优选
// 'FH2023', 'FH_mean','FH_speed2023','NDVI_zf','NDVImedian','NDVI_speed','Texture_speed2023',
// 'REPI_minQ1d',  'REPI_median', 'REPI_maxQ3d', 'REPI_stdd', 'NDVI_shixuw','VDDPI_stdDev',
// 'elevation', 'slope',
// 'DIF', 'mRVI','RVI','B5_dvar',
// 植物区系6中二级类最佳分类组合
// var band1 =  [ 'SVVI', 'B11',  'greenness', 'MCARI', 'B8', 'NDVIR1', 'B9',  'B4', 'REPI', 'DIF', 'mRVI', 'RVI', 'REPI_minQ1d', 'REPI_maxQ3d', 'elevation', 'slope', 'Premean', 'Temmean'];
// 构建各植物区系中用于植被分类的特征工程
var feature1 = Feature_merge.select(band1).updateMask(forest_region.eq(10));
print("feature1", feature1);
// 训练分类器
var samplesData = feature1.sampleRegions({
  collection: sample1,
  properties: ['firstclass'],
  tileScale: 4,
  scale: 30
});
var trainingPart = samplesData;
var testingPart = feature1.sampleRegions({
  collection: Forestsam_test,
  properties: ['firstclass'],
  tileScale: 4,
  scale: 30
});
// (3) Divide samplesdata into the training and testing parts.
// var randomSamplesData = samplesData.randomColumn('random');
// var splitValue = 0.8;  
// var trainingPart = randomSamplesData.filter(ee.Filter.lt('random', splitValue));
// var testingPart = randomSamplesData.filter(ee.Filter.gte('random', splitValue));

// 3. Classify using Planet-NICFI imagery and test and export the classificate maps.
// (1) Configurate the RF model
var classifier = ee.Classifier.smileRandomForest(50).train({
  features: trainingPart,
  classProperty: 'firstclass',
  inputProperties: band1
});

var classified19 = feature1.classify(classifier);
// 分类精度验证
var test = testingPart.classify(classifier);
// Print the confusion matrix.
var confusionMatrix = test.errorMatrix('firstclass', 'classification');
Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, {
      matrix: confusionMatrix.array(),
      fscore: confusionMatrix.fscore(),
      kappa: confusionMatrix.kappa(),
      accuracy: confusionMatrix.accuracy(),
      producersaccuracy: confusionMatrix.producersAccuracy(),
      consumersAccuracy: confusionMatrix.consumersAccuracy()
    }
  )]),
  description: "sampleall_zone1test",
  folder:"forest2similarity",
  fileFormat: "CSV"
});
// print("kappa", confusionMatrix.kappa());
// print("OA", confusionMatrix.accuracy());
var resultImg = classified19.clip(roi).toByte();
// Map.addLayer(forestsam_guding, {color: "ff0000"}, "forestsam_guding");
Map.addLayer(resultImg.updateMask(forest_region.eq(10)), {min:10, max: 40, palette: ['#00FF00','#FFFF00','#FF8C00','#FF0000']},'classified_result');
// 分类结果上传Asset
// Export.image.toAsset({
//   image: resultImg.updateMask(forest_region.eq(10)).clip(table5), 
//   description: 'Forestclassification2_Zone6', 
//   assetId: 'Forestclassification2_Zone6', 
//   // fileNamePrefix: 'lt-gee_disturbance_map', 
//   region: table5, 
//   scale: 30, 
//   crs: 'EPSG:4326', 
//   maxPixels: 1e13
// });
