// Load GAUL regions for New Zealand
var gaul = ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level1");
var nzRegions = gaul.filter(ee.Filter.eq('ADM0_NAME', 'New Zealand'));

// Select three regions by name
var regionNames = ['Auckland', 'Canterbury', 'Wanganui-manawatu'];  // Use exact names
var selectedRegions = nzRegions.filter(ee.Filter.inList('ADM1_NAME', regionNames));

// Load the Dynamic World collection once and reduce it
var dw = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
  .filterDate('2024-11-01', '2025-04-01')
  // .filter(ee.Filter.notNull(['label']))
  .select('label')
  .reduce(ee.Reducer.mode());

// Function to export and display clipped land cover per region
var processRegion = function(feature) {
  var name = String(feature.properties.ADM1_NAME);
  var geom = ee.Geometry(feature.geometry);
  var clipped = dw.clip(geom).toUint8();

  // Add layer to map
  Map.addLayer(clipped, {
    min: 0,
    max: 8,
    palette: [
      '419bdf', '397d49', '88b053', '7a87c6', 'e49635',
      'dfc35a', 'c4281b', 'a59b8f', 'b39fe1'
    ]
  }, name + ' Land Cover');

  // Export to Google Drive
  Export.image.toDrive({
    image: clipped,
    description: name.replace(/\s+/g, '_') + '_DynamicWorld',
    folder: 'GEE_Exports',
    // fileNamePrefix: ee.String(name).cat('_DynamicWorld').getInfo(),
    region: geom,
    scale: 10, // 10m resolution
    crs: 'EPSG:4326',
    maxPixels: 1e9,
    fileFormat: 'GeoTIFF',
    formatOptions: {
      cloudOptimized: true
    }
  });

  return feature;
};

selectedRegions.evaluate(function(fc) {
  fc.features.forEach(processRegion);
});

Map.clear();
Map.centerObject(selectedRegions, 5);
