// 植物区系的矢量边界和各区域的样本点数据
var sample4 = mixedF.merge(bamboo);
var sample4buff = sample4.map(function(f){
  return f.buffer(100,10);
});
sample = broadF.merge(confierF);
// 植物区系矢量转栅格
var Forestzone_img  =  ForestZone.reduceToImage({
    properties: ['shengtaizo'],
    reducer: ee.Reducer.mean()
    });
Forestzone_img = Forestzone_img.rename("forestzone");
// print("Forestzone_img",Forestzone_img);
// 计算各植物区系的样本集
var samplesall = Forestzone_img.sampleRegions({
  collection: sample,  
  scale: 30 ,
  tileScale: 4,
  geometries: true
});
var sample1 = samplesall.filterMetadata({name:"forestzone", operator:"equals", value: 6});
var conF = sample1.filterMetadata('firstclass','equals','10');
var broadF1 = sample1.filterMetadata('firstclass','equals','20');
print("conF",conF);
print("broadF1",broadF1);
sample1 = sample1.merge(mixedF).merge(bamboo);
print("sample1", sample1);
// print("geometry", geometry.merge(geometry2));
Map.addLayer(sample1, {color:"ff0000"}, "sample");
// // 植物区系的矢量边界
var forestzone1 = ForestZone.filterMetadata('shengtaizo','equals','1');
roi = forestzone1.merge(sample4buff);
var empty = ee.Image().toByte();
var outline = empty.paint({
 featureCollection:roi,  // 筛选的colletion
 color:0, //颜色透明
 width:2  //边界宽度
});
//绘制红色边界
Map.centerObject(forestzone1,8);
Map.addLayer(outline, {palette: "ff0000"}, "outline");
// 森林区域数据的获取
var forest_region = FNF2023.clip(roi);
Map.addLayer(forest_region,{},"FNF2023", false);
// 树高特征
var FH_feature = ee.Image.cat([FH2023.rename("FH2023"), FH_stdL, FH_ZFL, FH_cvL, FH_meanL, FH_medianL, FH_speedy2023L, FH_speedyL]);
FH_feature = FH_feature.updateMask(forest_region.eq(10));
// 成林速度特征
var NDVI_feature = ee.Image.cat([NDVI_ZFL, NDVI_cvL, NDVI_meanL, NDVI_medianL, NDVI_shixuw, NDVI_speedy2023L, NDVI_speedyL, NDVI_stdL, NDVI_zf]);
// NDVI_feature = NDVI_feature.clip(roi);
var Texture_feature = ee.Image.cat([Texture_ZFL, Texture_cvL, Texture_meanL, Texture_medianL, Texture_speedy2023L, Texture_speedyL, Texture_stdL]);
// 地理环境特征（地形、降水和温度）
var DEM = ee.Image('USGS/SRTMGL1_003');
var dem =DEM.clip(roi);
//slope
var slope = ee.Terrain.slope(dem);
var aspect = ee.Terrain.aspect(dem);
var environment = ee.Image.cat([dem, slope, aspect, PreMS.rename("PreMS"),Premean.rename("Premean"), TemMS.rename('TemMS'), Temmean.rename("Temmean"), SMCI.rename("SMCI")]);
environment = environment.updateMask(forest_region.eq(10));
// 物候特征：REPI
var phenology = ee.Image.cat([REPI_Q3Q1zfd, REPI_mind, REPI_cvd, REPI_maxQ3d, REPI_maxd, REPI_meand, REPI_mediand, REPI_minQ1d, REPI_stdd, REPI_zhenfud]); 
phenology = phenology.updateMask(forest_region.eq(10));
print("phenology", phenology);
// Map.centerObject(roi, 8);
// Sentinel-2数据去云
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
//=====================增加特征值(NDVI\REI\EVI\SAVI\REPI)===============
var addVIs2 = function(image){
        
        var NDVI = image.normalizedDifference(['B8','B4']);
        
        var EVI = image.expression(
        '2.5*((NIR-Red)/(NIR+6*Red-7.5*blue+10000))',
        {
            blue: image.select('B2'),    // 0.452-0.512, BLUE
            Red: image.select('B4'),    // 0.636-0.673 μm, RED
            NIR: image.select('B8'),    //  NIR
        });
        var REPI = image.expression(
        '705+3.5*(((B4+B7)/2-B6)/(B6-B5))',
        {
            B4: image.select('B4'),    // 0.636-0.673 μm, RED
            B5: image.select('B8'),    //  NIR
            B6: image.select('B6'), 
            B7: image.select('B7'), 
        });
        var SVVI = image.expression(
        "float(sqrt((B-((B+G+R+NIR+swir1+swir2)/5))**2+(G-((B+G+R+NIR+swir1+swir2)/5))**2+(R-((B+G+R+NIR+swir1+swir2)/5))**2+(NIR-((B+G+R+NIR+swir1+swir2)/5))**2+(swir1-((B+G+R+NIR+swir1+swir2)/5))**2+(swir2-((B+G+R+NIR+swir1+swir2)/5))**2)-sqrt((NIR-((NIR+swir1+swir2)/3))**2+(swir1-((NIR+swir1+swir2)/3))**2+(swir2-((NIR+swir1+swir2)/3))**2))",
        {
                        "B": image.select("B2"),
                        "G": image.select("B3"),
                        "R": image.select("B4"),
                        "NIR": image.select("B8"),
                        "swir1": image.select("B11"),
                        "swir2": image.select("B12")
                      });
        // NDRE
        var NDRE = image.normalizedDifference(['B8', 'B7']);
        // MCARI
        var MCARI = image.expression(
        "((B5-B4)-0.2*(B5-B3)*(B5/B4))",
        {
                        "B4": image.select("B4"),
                        "B5": image.select("B5"),
                        "B3": image.select("B3")
                      });
        // MCR
        var MCR = image.expression(
        "(B8-B5-1)/(sqrt(B8-B5)+1)",
        {
                        "B5": image.select("B5"),
                        "B8": image.select("B8")
                      });
        // NDPI
        var NDPI = image.expression(
        "(B8-(0.74*B4+0.26*B11))/(B8+(0.74*B4+0.26*B11))",
        {
                        "B4": image.select("B4"),
                        "B8": image.select("B8"),
                        "B11": image.select("B11")
                      });

        // REI
        var REI = image.normalizedDifference(['B6', 'B5']);
        // IRECI
        var IRECI = image.expression(
        "(B7-B4)/(B5/B6)",
        {
                        "B4": image.select("B4"),
                        "B5": image.select("B8"),
                        "B6": image.select("B11"),
                        "B7": image.select("B7")
                      });
        var NDVIR1 = image.normalizedDifference(['B8','B5']);
        //加入波段信息
        return image.addBands(NDVI.rename("NDVI"))
                  .addBands(SVVI.rename('SVVI'))
                  .addBands(EVI.rename('EVI'))
                  .addBands(NDRE.rename('NDRE'))
                  .addBands(REPI.rename('REPI'))
                  .addBands(REI.rename('REI'))
                  //.addBands(MSAVI.rename('MSAVI'))
                  .addBands(MCARI.rename('MCARI'))
                  .addBands(NDPI.rename('NDPI'))                 
                  .addBands(MCR.rename('MCR'))
                  .addBands(IRECI.rename('IRECI'))
                  .addBands(NDVIR1.rename('NDVIR1'));                     
    };
// var bands = ['B1', 'B2', 'B3','B4', 'B5', 'B6','B7', 'B8', 'B8A', 'B9','B10', 'B11', 'B12'];
var dataset = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
                  .filterDate('2022-10-01', '2023-10-01')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
                  .map(maskS2clouds)
                  .filterBounds(roi.geometry())
                  .map(function(image){return image.updateMask(forest_region.eq(10));})
                  .map(addVIs2);
// var datasetm = dataset.reduce(ee.Reducer.stdDev()).clip(roi);
// Sentinel-2数据纹理特征提取光谱特征
var sentinel2_spectrum = dataset.median();
print("sentinel2_spectrum", sentinel2_spectrum);
// print("sentinel2_spectrum",sentinel2_spectrum);
var bands = ['B1', 'B2', 'B3', 'B4','B5', 'B6', 'B7',"B8", 'B8A', 'B9', 'B10', 'B11', 'B12','SVVI', 'NDVI', 'EVI','NDRE',
            'REPI', 'REI', 'MCARI', 'NDPI', 'MCR', 'IRECI', 'NDVIR1'];
var sentinel2_spectrum1 = sentinel2_spectrum.select(bands);
print("sentinel2_spectrum1", sentinel2_spectrum1);
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
Map.addLayer(sentinel2_spectrum1, visualization, "Sentinel2");
// Sentinel-2数据纹理特征提取
var B5 = sentinel2_spectrum.select("B5");
B5 = B5.add(127.5).multiply(127.5).toUint16();
var sentinel2_glcm = B5.glcmTexture({size: 1,kernel:null});
//sentinel-2缨帽变换
//Sentinel-2数据k-t变换系数
var coefficients = ee.Array([  
  [0.0356, 0.0822, 0.1360, 0.2611, 0.2964, 0.3338, 0.3877, 0.3895, 0.0949, 0.3882, 0.1366, 0.4750],  
  [-0.0635, -0.1128, -0.1680, -0.3480, -0.3303, 0.0852, 0.3302, 0.3165, 0.0467, -0.4578,  -0.4064, 0.3625],  
  [0.0649, 0.1363, 0.2802, 0.3072, 0.5288, 0.1379, -0.0001, -0.0807, -0.0302,-0.4064, -0.5602,-0.1389]
]); 
var image= dataset.select(['B1','B2','B3','B4','B5','B6','B7','B8','B9','B11','B12','B8A']);
var arrayImage1D = image.mean().toArray();  
//print("arrayImage1D",arrayImage1D)
// Map.addLayer(arrayImage1D,{},"arrayImage1D");
 
 
var arrayImage2D = arrayImage1D.toArray(1);  
//print("arrayImage2D",arrayImage2D)
// Map.addLayer(arrayImage2D,{},"arrayImage2D");
//相乘变换  
var componentsImage = ee.Image(coefficients)   
                        .matrixMultiply(arrayImage2D)   
                        .arrayProject([0])   
                        .arrayFlatten([[  
                          'brightness', 'greenness', 'wetness' 
                        ]]);  
var Greenness =componentsImage.select(["greenness", 'wetness']);
//=====================增加sentinel-1特征值(NDVI\REI\EVI\SAVI\REPI)===============
var addVI = function(img){
        var ratio = img.select('VV').divide(img.select('VH'));
        var DIF = img.select('VV').subtract(img.select('VH'));
        var AVE =img.expression(
          '(VV+VH)/2',
          {
            VV: img.select('VV'),
            VH: img.select('VH'),
          });
        var NDI = img.normalizedDifference(['VV','VH']);
        
        var RVI = img.expression(
          'sqrt(VV/(VV+VH))*(VV/VH)',
           {
            VV: img.select('VV'),
            VH: img.select('VH'),
          });
       
        var mRVI = img.expression(
          'sqrt(VV/(VV+VH))*(4*VH/(VV+VH))',
           {
            VV: img.select('VV'),
            VH: img.select('VH'),
          });
        var VDDPI = img.expression(
          '(VV+VH)/VV',
           {
            VV: img.select('VV'),
            VH: img.select('VH'),
          });
        //加入波段信息
        return img.addBands(ratio.rename("ratio"))
                  .addBands(DIF.rename('DIF'))
                  .addBands(AVE.rename('AVE'))
                  .addBands(NDI.rename('NDI'))
                  .addBands(RVI.rename('RVI'))
                  .addBands(mRVI.rename('mRVI'))
                  .addBands(VDDPI.rename('VDDPI'));
    };
var s1_img = ee.ImageCollection('COPERNICUS/S1_GRD')
               .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
               .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
               .filter(ee.Filter.eq('instrumentMode', 'IW'))
               //降轨数据 'DESCENDING' 'ASCENDING'
              // .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
               .filterBounds(roi)
               .filterDate('2022-08-31', '2023-08-31')
               .map(function(image){return image.updateMask(forest_region.eq(10));})
               .map(addVI);
              // .median();
var rvi =s1_img.select('mRVI');
var summer = rvi.filterDate('2023-06-01', '2023-08-31');
summer = summer.max();
var winter = rvi.filterDate('2022-12-01', '2023-03-01');
winter = winter.max();
var VV_cha = summer.subtract(winter);
VV_cha = VV_cha.rename("mRVI_sumwinter");
var zf = (rvi.max()).subtract(rvi.min());
zf = zf.rename("mRVI_zf");
var std = rvi.reduce(ee.Reducer.stdDev());
std = std.rename("mRVI_std");
var cv = std.divide(rvi.mean());
cv = cv.rename("mRVI_cv");
var Sentinel1std = s1_img.reduce(ee.Reducer.stdDev());
// Map.addLayer(s1_img.select('VV').median(), {}, "sentinel1", false);
// Sentinel-1的物候特征
var senitnel1_phenology = ee.Image.cat([Sentinel1std, cv, std, zf, VV_cha]);
// Sentinel-1的极化特征
var sentinel_spectrum = s1_img.median();
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
