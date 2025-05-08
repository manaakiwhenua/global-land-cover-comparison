# Datasets

> Although Global land cover datasets offer a convenient and standardized means of analysis across large areas, their applicability and accuracy can vary significantly depending on geographical regions and resolution. Our evaluation of currently available global land cover datasets for mapping land cover types over New Zealand will consider the country’s unique geographical features, climate, and land use patterns. We would also assess spatial and temporal resolution and identify limitations for applicability in decision-making.

> The study will focus on four widely used global land cover datasets:

| Dataset                         | Source                     | Resolution   | Recurrence (ISO 8601)          | License & Attribution                                   | Cost                             | Further Information                                                                                                 |
|---------------------------------|----------------------------|--------------|-------------------------|--------------------------------------------------------|----------------------------------|----------------------------------------------------------------------------------------------------------------------|
| **World Cover Map**             | ESA                        | 10m          | R2/2020/P1Y       | Creative Commons Attribution 4.0 International License. [More info](https://github.com/ESA-WorldCover/esa-worldcover-datasets) | Free of charge.                  | [World Cover 2021](https://worldcover2021.esa.int/download)                           |
| **Dynamic World**               | Google                     | 10m          | R/2015/P1D        | Apache License 2.0. [More info](https://github.com/google/dynamicworld/blob/master/LICENSE) | Freely available.                | [Dynamic World](https://dynamicworld.app/about/)                                   |
| **Global Land Cover 30m (GLC30)**| Chinese Academy of Sciences | 30m          | R/1985/P5Y       | License details not specified in available sources. | Available for download; cost details not specified. | [GLC30 Data](https://data.casearth.cn/thematic/glc_fcs30) |
| **PlanetScope Land Cover Dataset** | Planet Lab              | 20m          | R/2020/P1D        | Commercial licensing. [More info](https://www.planet.com/pricing/) | Pricing varies; e.g., $2.25 per km² for archive data with a minimum order of 250 km². | [PlanetScope Pricing](https://apollomapping.com/planetscope-satellite-imagery) |

- World Cover Map
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

- Dynamic World:
    - Notes:
        - Result of partnership between Google and the [World Resources Institute](https://www.wri.org/)
        - Not an officially supported Google Product
        - Deep learning
        - "Dynamic World is intended to be used as a data product for users to add custom rules with which to assign final class values, producing derivative land cover maps."
        - Based on [Sentinel-2 Top of Atmosphere](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2)
        - Updating every 2-5 days depending on location
        - Per-pixel probabilities across 9 land cover classes
        - Paper: https://doi.org/10.1038/s41597-022-01307-4
        - The Dynamic World dataset has been made available as an Earth Engine Image Collection under "GOOGLE/DYNAMICWORLD/V1". This is referenced in either the Earth Engine Python or JavaScript client library with: ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').

- GLC30
    - Citation: Chen, Jun et al.: 2015. Global land cover mapping at 30 m resolution: A POK-based operational approach. ISPRS Journal of Photogrammetry and Remote Sensing Volume 103, May 2015, Pages 7–27 http://dx.doi.org/10.1016/j.isprsjprs.2014.09.002

> Analysis will focus on different ecological zones of New Zealand, including urban areas, forested regions, coastal zones, agricultural land, and alpine zones, to ensure a diverse and representative evaluation of land cover types. We would acquire global land cover datasets for the most recent years available and compare them based on resolution, classification scheme, temporal resolution, and cost and availability. We will review source imagery and method of map production. We will compare datasets to available local datasets at similar time periods (LCDB (1ha) and woody layer (10m)), and where different we would cross-reference the datasets with Sentinel-2 satellite imagery and regional LiDAR, to validate overall classification accuracy. These analyses will allow us to summarise the strengths and weaknesses of each method in the New Zealand context, including the feasibility of implementation. 