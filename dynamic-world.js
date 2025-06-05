// Load GAUL regions for New Zealand

var nzBounds = ee.Geometry.BBox(166.223, -47.502, 178.615, -34.32);

// Load Dynamic World label band, reduce with mode across date range
var dw = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
  // .filterDate('2024-11-01', '2025-04-01')
  .filterDate('2023-11-01', '2024-04-01')
  .filterBounds(nzBounds)
  .select('label')
  .reduce(ee.Reducer.mode())
  .clip(nzBounds)
  .toUint8();

// Add layer to map
Map.centerObject(nzBounds, 5);
Map.addLayer(dw, {
  min: 0,
  max: 8,
  palette: [
    '419bdf', '397d49', '88b053', '7a87c6', 'e49635',
    'dfc35a', 'c4281b', 'a59b8f', 'b39fe1'
  ]
}, 'NZ Land Cover');

Export.image.toCloudStorage({
  image: dw,
  // description: 'NZ_DynamicWorld_Mode_COG-2024-11-01-P5M',
  description: 'NZ_DynamicWorld_Mode_COG-2023-11-01-P5M',
  bucket: 'cost-effective-land-cover',
  fileNamePrefix: 'landcover/DW_mode_NZ',
  region: nzBounds,
  scale: 10, // 10 m resolution
  crs: 'EPSG:2193',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
});