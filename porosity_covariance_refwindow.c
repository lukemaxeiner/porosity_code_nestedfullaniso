#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tiffio.h>

#include "porosity_covariance_refwindow.h"
#include "tiff.h"

static uint64_t porosity_field_index(uint32_t wx, uint32_t wy, uint32_t wz,
                                     uint32_t width, uint32_t height) {
    return (uint64_t)wz * (uint64_t)width * (uint64_t)height +
           (uint64_t)wy * (uint64_t)width + (uint64_t)wx;
}

static char *build_covariance_ref_window_output_path(const char *stack_path) {
    const char *suffix = "_porosity_covariance_ref_window.tif";
    const char *extension = strrchr(stack_path, '.');
    size_t base_length;
    char *output_path;

    if (stack_path == NULL) {
        return NULL;
    }

    base_length = extension != NULL ? (size_t)(extension - stack_path)
                                    : strlen(stack_path);
    output_path = (char *)malloc(base_length + strlen(suffix) + 1U);
    if (output_path == NULL) {
        return NULL;
    }

    memcpy(output_path, stack_path, base_length);
    memcpy(output_path + base_length, suffix, strlen(suffix) + 1U);
    return output_path;
}

static double compute_porosity_field_mean(const PorosityField *field) {
    uint64_t total_samples;
    uint64_t i;
    double sum = 0.0;

    total_samples = (uint64_t)field->width * (uint64_t)field->height *
                    (uint64_t)field->depth;
    if (total_samples == 0) {
        return 0.0;
    }

    for (i = 0; i < total_samples; ++i) {
        sum += field->values[i];
    }

    return sum / (double)total_samples;
}

int compute_reference_window_covariance_map(const PorosityField *field,
                                            uint32_t ref_x, uint32_t ref_y,
                                            uint32_t ref_z,
                                            PorosityCovarianceField *map) {
    uint32_t wz;
    uint32_t wy;
    uint32_t wx;
    double mean;
    double reference_value;

    if (field == NULL || map == NULL || field->values == NULL) {
        return EXIT_FAILURE;
    }

    if (ref_x >= field->width || ref_y >= field->height ||
        ref_z >= field->depth) {
        return EXIT_FAILURE;
    }

    memset(map, 0, sizeof(*map));

    mean = compute_porosity_field_mean(field);
    reference_value =
        field->values[porosity_field_index(ref_x, ref_y, ref_z, field->width,
                                           field->height)];

    map->width = field->width;
    map->height = field->height;
    map->depth = field->depth;
    map->ref_x = ref_x;
    map->ref_y = ref_y;
    map->ref_z = ref_z;
    map->reference_value = reference_value;
    map->field_mean = mean;
    map->values = (double *)malloc((size_t)map->width * (size_t)map->height *
                                   (size_t)map->depth * sizeof(double));
    if (map->values == NULL) {
        memset(map, 0, sizeof(*map));
        return EXIT_FAILURE;
    }

    for (wz = 0; wz < map->depth; ++wz) {
        for (wy = 0; wy < map->height; ++wy) {
            for (wx = 0; wx < map->width; ++wx) {
                uint64_t idx =
                    porosity_field_index(wx, wy, wz, map->width, map->height);
                double current_value = field->values[idx];

                map->values[idx] =
                    (reference_value - mean) * (current_value - mean);
            }
        }
    }

    return EXIT_SUCCESS;
}

void free_porosity_covariance_field(PorosityCovarianceField *map) {
    if (map == NULL) {
        return;
    }

    free(map->values);
    map->values = NULL;
    map->width = 0;
    map->height = 0;
    map->depth = 0;
    map->ref_x = 0;
    map->ref_y = 0;
    map->ref_z = 0;
    map->reference_value = 0.0;
    map->field_mean = 0.0;
}

int write_porosity_covariance_field_tiff(const char *path,
                                         const PorosityCovarianceField *map) {
    TIFF *tif;
    uint32_t wz;

    if (path == NULL || map == NULL || map->values == NULL) {
        return EXIT_FAILURE;
    }

    tif = TIFFOpen(path, "w8");
    if (tif == NULL) {
        return EXIT_FAILURE;
    }

    for (wz = 0; wz < map->depth; ++wz) {
        float *scanline;
        uint32_t wy;

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, map->width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, map->height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, map->height);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

        scanline = (float *)_TIFFmalloc((tmsize_t)map->width * sizeof(float));
        if (scanline == NULL) {
            TIFFClose(tif);
            return EXIT_FAILURE;
        }

        for (wy = 0; wy < map->height; ++wy) {
            uint32_t wx;

            for (wx = 0; wx < map->width; ++wx) {
                scanline[wx] = (float)map->values[porosity_field_index(
                    wx, wy, wz, map->width, map->height)];
            }

            if (TIFFWriteScanline(tif, scanline, wy, 0) < 0) {
                _TIFFfree(scanline);
                TIFFClose(tif);
                return EXIT_FAILURE;
            }
        }

        _TIFFfree(scanline);

        if (wz + 1U < map->depth) {
            TIFFWriteDirectory(tif);
        }
    }

    TIFFClose(tif);
    return EXIT_SUCCESS;
}

int process_porosity_covariance_ref_window(const char *stack_path,
                                           double voxel_resolution,
                                           double filter_resolution,
                                           int ref_x, int ref_y, int ref_z,
                                           FILE *out) {
    BinaryVolume volume;
    PorosityField field;
    PorosityCovarianceField map;
    char *output_path;

    if (stack_path == NULL || out == NULL) {
        return EXIT_FAILURE;
    }
    if (ref_x < 0 || ref_y < 0 || ref_z < 0) {
        fprintf(stderr,
                "REF_WINDOW_VALUES must be non-negative for %s\n",
                stack_path);
        return EXIT_FAILURE;
    }

    volume.width = 0;
    volume.height = 0;
    volume.depth = 0;
    volume.voxels = NULL;
    field.width = 0;
    field.height = 0;
    field.depth = 0;
    field.values = NULL;
    memset(&map, 0, sizeof(map));
    output_path = NULL;

    if (load_binary_volume_from_tiff(stack_path, &volume) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to load TIFF volume for %s\n", stack_path);
        return EXIT_FAILURE;
    }

    if (compute_porosity_field(volume.voxels, volume.width, volume.height,
                               volume.depth, filter_resolution,
                               voxel_resolution, &field) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to compute porosity field for %s\n",
                stack_path);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (compute_reference_window_covariance_map(
            &field, (uint32_t)ref_x, (uint32_t)ref_y, (uint32_t)ref_z,
            &map) != EXIT_SUCCESS) {
        fprintf(stderr,
                "Failed to compute reference-window covariance for %s\n",
                stack_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    output_path = build_covariance_ref_window_output_path(stack_path);
    if (output_path == NULL) {
        fprintf(stderr,
                "Failed to build reference-window covariance output path for %s\n",
                stack_path);
        free_porosity_covariance_field(&map);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (write_porosity_covariance_field_tiff(output_path, &map) !=
        EXIT_SUCCESS) {
        fprintf(stderr,
                "Failed to write reference-window covariance TIFF for %s\n",
                stack_path);
        free(output_path);
        free_porosity_covariance_field(&map);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "Reference-window covariance: %s\n"
                "  Output: %s\n"
                "  Field dimensions: %u x %u x %u\n"
                "  Reference window: %d %d %d\n\n",
                stack_path, output_path, map.width, map.height, map.depth,
                ref_x, ref_y, ref_z) < 0) {
        fprintf(stderr,
                "Failed to write reference-window covariance summary for %s\n",
                stack_path);
        free(output_path);
        free_porosity_covariance_field(&map);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    free(output_path);
    free_porosity_covariance_field(&map);
    free_porosity_field(&field);
    free_binary_volume(&volume);
    return EXIT_SUCCESS;
}
