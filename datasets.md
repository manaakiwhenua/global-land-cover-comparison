# Datasets

> Although Global land cover datasets offer a convenient and standardized means of analysis across large areas, their applicability and accuracy can vary significantly depending on geographical regions and resolution. Our evaluation of currently available global land cover datasets for mapping land cover types over New Zealand will consider the country’s unique geographical features, climate, and land use patterns. We would also assess spatial and temporal resolution and identify limitations for applicability in decision-making.

> The study will focus on three widely used global land cover datasets:

<table>
  <thead>
    <tr>
      <th>Dataset</th>
      <th>Source</th>
      <th>Resolution</th>
      <th>Recurrence (ISO 8601)</th>
      <th>License & Attribution</th>
      <th>Cost</th>
      <th>Validation</th>
      <th>Further Information</th>
      <th>Number of classes</th>
      <!-- <th>Classes & Description</th> -->
      <th>Format</th>
      <th>Projection</th>
      <th>Input Data</th>
      <th>Citation</th>
      <th>Download Link</th>
      <th>GEE Image Collection Name</th>
      <th>Attribution</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>World Cover Map</b></td>
      <td>ESA</td>
      <td>10m</td>
      <td>R2/2020/P1Y</td>
      <td>CC BY 4.0</td>
      <td>Free</td>
      <td>The ESA WorldCover product has been independently validated by Wageningen University (statistical accuracy) and IIASA (spatial accuracy).
        <ul>
            <li><a href="https://worldcover2020.esa.int/data/docs/WorldCover_PVR_V1.1.pdf">75% (2020)</a></li>
            <li><a href="https://worldcover2021.esa.int/data/docs/WorldCover_PVR_V2.0.pdf">76.7% (2021)</a></li>
        </ul>
      </td>
      <td><a href="https://github.com/ESA-WorldCover/esa-worldcover-datasets">ESA World Cover Datasets</a></td>
      <td>11</td>
      <!-- <td>1. Tree cover #006400<br>2. Shrubland #ffbb22<br>3. Grassland #ffff4c<br>4. Cropland #f096ff<br>5. Built-up #fa0000<br>6. Bare/sparse vegetation #b4b4b4<br>7. Snow and ice #f0f0f0<br>8. Permanent water bodies #0064c8<br>9. Herbaceous wetland #0096a0<br>10. Mangroves #00cf75<br>11. Moss and lichen #fae6a0</td> -->
      <td>COG</td>
      <td>EPSG:4326</td>
      <td>Sentinel-1, Sentinel-2</td>
      <!-- <td>Citation: European Space Agency. (2021).</td> -->
      <td>
        <ul>
            <li>2020 v100: Zanaga, D., Van De Kerchove, R., De Keersmaecker, W., Souverijns, N., Brockmann, C., Quast, R., Wevers, J., Grosu, A., Paccini, A., Vergnaud, S., Cartus, O., Santoro, M., Fritz, S., Georgieva, I., Lesiv, M., Carter, S., Herold, M., Li, Linlin, Tsendbazar, N.E., Ramoino, F., Arino, O., 2021. ESA WorldCover 10 m 2020 v100. <a href="https://doi.org/10.5281/zenodo.5571936">doi.org/10.5281/zenodo.5571936</a>
            </li>
            <li>2021 v200: Zanaga, D., Van De Kerchove, R., Daems, D., De Keersmaecker, W., Brockmann, C., Kirches, G., Wevers, J., Cartus, O., Santoro, M., Fritz, S., Lesiv, M., Herold, M., Tsendbazar, N.E., Xu, P., Ramoino, F., Arino, O., 2022. ESA WorldCover 10 m 2021 v200. <a href="https://doi.org/10.5281/zenodo.7254221">doi.org/10.5281/zenodo.7254221</a></li>
        </ul>
      </td>
      <td>
        <ul>
            <li><a href="https://worldcover2020.esa.int/download">World Cover 2020</a></li>
            <li><a href="https://worldcover2021.esa.int/download">World Cover 2021</a></li>
        </ul>
      </td>
      <td>
        <ul>
            <li>ESA/WorldCover/v100</li>
            <li>ESA/WorldCover/v200</li>
        </ul>
      </td>
      <td>© ESA WorldCover project / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium</td>
    </tr>
    <tr>
      <td><b>Dynamic World</b></td>
      <td>Google and the <a href="https://www.wri.org/">World Resources Institute</a></td>
      <td>10m</td>
      <td>R/2015/P1D</td>
      <td>CC BY 4.0</td>
      <td>Free</td>
      <td></td>
      <td><a href="https://dynamicworld.app/about/">Dynamic World</a></td>
      <td>9</td>
      <!-- <td>1. Water #419BDF<br>2. Trees #397D49<br>3. Grass #88B053<br>4. Flooded vegetation #7A87C6<br>5. Crops #E49635<br>6. Shrub & Scrub #DFC35A<br>7. Built Area #C4281B<br>8. Bare Ground #A59B8F<br>9. Snow & Ice #B39FE1</td> -->
      <td>COG</td>
      <td>EPSG:4326</td>
      <td>Sentinel-2 Top of Atmosphere</td>
      <td>Brown, C.F., Brumby, S.P., Guzder-Williams, B. et al. Dynamic World, Near real-time global 10 m land use land cover mapping. Sci Data 9, 251 (2022). <a href="https://doi.org/10.1038/s41597-022-01307-4">doi.org/10.1038/s41597-022-01307-4</a></td>
      <td><a href="https://dynamicworld.app/download">Dynamic World</a></td>
      <td>GOOGLE/DYNAMICWORLD/V1</td>
      <td>This dataset is produced for the Dynamic World Project by Google in partnership with National Geographic Society and the World Resources Institute.</td>
    </tr>
    <tr>
      <td><b>Global Land Cover Change Dataset (GLC FCS30D)</b></td>
      <td>Aerospace Information Research Institute, Chinese Academy of Sciences</td>
      <td>30m</td>
      <td>R4/1985/P5Y, R/2000/P1Y</td>
      <td>CC BY 4.0</td>
      <td>Free</td>
      <td>80.88% overall accuracy, validated over 84,000 global samples</td>
      <td><a href="https://data.casearth.cn/thematic/glc_fcs30">GLC FCS30</a></td>
      <td>35</td>
      <!-- <td>1. Forest #008000<br>2. Grassland #32CD32<br>3. Cropland #FFD700<br>4. Built-up #FF6347<br>5. Wetland #1E90FF<br>6. Bareland #D3D3D3<br>7. Snow/Ice #FFFFFF<br>8. Water #0000FF</td> -->
      <td>COG</td>
      <td>EPSG:4326</td>
      <td>Landsat</td>
      <td>Zhang, X., Zhao, T., Xu, H., Liu, W., Wang, J., Chen, X., and Liu, L.: GLC_FCS30D: the first global 30 m land-cover dynamics monitoring product with a fine classification system for the period from 1985 to 2022 generated using dense-time-series Landsat imagery and the continuous change-detection method, Earth Syst. Sci. Data, 16, 1353–1381, <a href="https://doi.org/10.5194/essd-16-1353-2024">doi.org/10.5194/essd-16-1353-2024</a>, 2024.</td>
      <td><a href="https://doi.org/10.5281/zenodo.8239305">GLC30</a></td>
      <td>
        <ul>
            <li>projects/sat-io/open-datasets/GLC-FCS30D/annual</li>
            <li>projects/sat-io/open-datasets/GLC-FCS30D/five-years-map</li>
        </ul>
      </td>
      <td></td>
    </tr>
  </tbody>
</table>

<<<<<<< Updated upstream

<!-- - World Cover Map
    - InputQuality: Three-band GeoTIFF providing three per-pixel quality indicators of the Sentinel-1 and Sentinel-2 input data
    - "Since the WorldCover maps for 2020 and 2021 were generated with different algorithm versions (v100 and v200, respectively), changes between the maps should be treated with caution, as they include both changes in real land cover and changes due to the used algorithms."
=======
- World Cover Map
    - Downloading
        - NB as it uses session tokens, must be downloaded interactively. Downloads a ZIP archive that must be extracted. NZ images are those that match the following regular expression: `ESA_WorldCover_10m_(?:2020|2021)_V(?:100|200)_60deg_macrotile_S90E120/ESA_WorldCover_10m_(?:2020|2021)_V(?:100|200)_S\d{2}E(?:165|168|171|174|177)_Map.tif`
    - Format: 3x3 degree tiles as Cloud Optimized GeoTIFFs (COGs)
    - Projection: EPSG:4326
    - Accuracy: 75% overall accuracy (2020)
    - Classes:
        1. Tree cover
        2. Shrubland
        3. Grassland
        4. Cropland
        5. Built-up
        6. Bare / sparse vegetation
        7. Snow and Ice
        8. Permanent water bodies
        9. Herbaceous Wetland
        10. Mangrove
        11. Moss and lichen
    - Download link:
        - 2020: https://worldcover2020.esa.int/download
        - 2021: https://worldcover2021.esa.int/downloader
    - Data layers:
        - Map: land cover map with 11 classes
        - InputQuality: Three-band GeoTIFF providing three per-pixel quality indicators of the Sentinel-1 and Sentinel-2 input data
    - Attribution: © ESA WorldCover project / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium
    - Citation: Zanaga, D., Van De Kerchove, R., Daems, D., De Keersmaecker, W., Brockmann, C., Kirches, G., Wevers, J., Cartus, O., Santoro, M., Fritz, S., Lesiv, M., Herold, M., Tsendbazar, N.E., Xu, P., Ramoino, F., Arino, O., 2022. ESA WorldCover 10 m 2021 v200. https://doi.org/10.5281/zenodo.7254221
    - Notes:
        - "Since the WorldCover maps for 2020 and 2021 were generated with different algorithm versions (v100 and v200, respectively), changes between the maps should be treated with caution, as they include both changes in real land cover and changes due to the used algorithms."
>>>>>>> Stashed changes

- Dynamic World:
    - Not an officially supported Google Product
    - Deep learning methodology
    - "Dynamic World is intended to be used as a data product for users to add custom rules with which to assign final class values, producing derivative land cover maps."
    - Updates every 2-5 days depending on location
    - Provides per-pixel probabilities across 9 land cover classes

- GLC30_FCS30D
    - "Developed using continuous change detection methods and leveraging the extensive Landsat imagery archives within the Google Earth Engine platform" -->

> Analysis will focus on different ecological zones of New Zealand, including urban areas, forested regions, coastal zones, agricultural land, and alpine zones, to ensure a diverse and representative evaluation of land cover types. We would acquire global land cover datasets for the most recent years available and compare them based on resolution, classification scheme, temporal resolution, and cost and availability. We will review source imagery and method of map production. We will compare datasets to available local datasets at similar time periods (LCDB (1ha) and woody layer (10m)), and where different we would cross-reference the datasets with Sentinel-2 satellite imagery and regional LiDAR, to validate overall classification accuracy. These analyses will allow us to summarise the strengths and weaknesses of each method in the New Zealand context, including the feasibility of implementation. 