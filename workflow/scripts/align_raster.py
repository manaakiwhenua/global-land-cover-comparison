from osgeo import gdal

def align_raster_to_template(template_path, source_path, output_path):
    # Open template raster (reference)
    template_ds = gdal.Open(template_path)
    if template_ds is None:
        raise RuntimeError(f"Could not open template raster: {template_path}")

    # Open source raster to reproject
    source_ds = gdal.Open(source_path)
    if source_ds is None:
        raise RuntimeError(f"Could not open source raster: {source_path}")

    # Get template info
    proj = template_ds.GetProjection()
    gt = template_ds.GetGeoTransform()
    xsize = template_ds.RasterXSize
    ysize = template_ds.RasterYSize
    band_count = source_ds.RasterCount

    datatype = source_ds.GetRasterBand(1).DataType

    # Prepare creation options
    creation_options = [f"COMPRESS=DEFLATE", "TILED=YES", ]

    # Create output raster, with same size, geotransform and projection as template
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_path, xsize, ysize, band_count, datatype, creation_options)

    if out_ds is None:
        raise RuntimeError(f"Could not create output raster: {output_path}")

    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(proj)

    # Copy nodata from source if exists, for each band
    for i in range(band_count):
        src_band = source_ds.GetRasterBand(i+1)
        dst_band = out_ds.GetRasterBand(i+1)
        nodata = src_band.GetNoDataValue()
        if nodata is not None:
            dst_band.SetNoDataValue(nodata)

    # Reproject source to output raster grid
    gdal.ReprojectImage(source_ds, out_ds, source_ds.GetProjection(), proj, gdal.GRA_NearestNeighbour)

    # Flush and close datasets
    out_ds.FlushCache()
    del out_ds
    del source_ds
    del template_ds

if __name__ == "__main__":
    # snakemake.input and snakemake.output are lists
    template_path = snakemake.input.template
    source_path = snakemake.input.source
    output_path = snakemake.output[0]
    align_raster_to_template(template_path, source_path, output_path)