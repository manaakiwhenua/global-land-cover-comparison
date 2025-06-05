import os
from dotenv import load_dotenv

# Load environment variables from the .env file
load_dotenv('secrets.env')

# Configuration
wfs_endpoint = config["lcdb"]["wfs_endpoint"]
layer = config["lcdb"]["layer"]
attribute_field = config["lcdb"]["attribute_field"]
nodata_value = config["lcdb"]["nodata_value"]

rule download_wfs:
    conda: "../envs/gdal.yaml"
    output:
        fgb=os.path.join(output_base, "lcdb", "lcdb.fgb")
    params:
        endpoint=wfs_endpoint.replace("{LRIS_KEY}", os.getenv('LRIS_KEY')),
        layer=layer
    shell: """
        ogr2ogr --config GDAL_HTTP_UNSAFESSL YES -progress -f FlatGeobuf -t_srs EPSG:2193 {output.fgb} WFS:\"{params.endpoint}\" {params.layer} -nln lcdb
        """

rule rasterize_fgb:
    conda: "../envs/gdal.yaml"
    input:
        fgb=os.path.join(output_base, "lcdb", "lcdb.fgb")
    output:
        raster=os.path.join(output_base, "lcdb", "lcdb_raster.tif")
    params:
        pixel_size=10,
        numeric_field=attribute_field,
        nodata_value=nodata_value
    shell: """
        gdal_rasterize -a {params.numeric_field} \
            -l lcdb \
            -of GTiff \
            -tr {params.pixel_size} {params.pixel_size} \
            -tap \
            -ot Byte \
            -co COMPRESS=DEFLATE \
            -co ZLEVEL=6 \
            -co TILED=YES \
            -co BLOCKXSIZE=256 \
            -co BLOCKYSIZE=256 \
            -a_srs EPSG:2193 \
            "{input.fgb}" \
            "{output.raster}"
        """

# # Rule: clip and reproject to existing grid
# rule clip_and_reproject:
#     conda: "envs/gdal.yaml"
#     input:
#         raster="data/{region}_raster.tif",
#         extent=grid_extent
#     output:
#         final_raster="output/{region}_clipped_reprojected.tif"
#     shell:
#         r"""
#         gdalwarp \
#             -t_srs EPSG:4326 \  # Or any desired output SRS
#             -cutline "{input.extent}" \
#             -crop_to_cutline \
#             -dstnodata {nodata_value} \
#             "{input.raster}" \
#             "{output.final_raster}"
#         """
