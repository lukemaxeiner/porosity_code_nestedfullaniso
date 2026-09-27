#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tiffio.h>

#include "porosity_field.h"
#include "tiff.h"

#ifndef POROSITY_FIELD_STATIC
#define POROSITY_FIELD_STATIC static
#endif

POROSITY_FIELD_STATIC uint64_t volume_index(uint32_t x, uint32_t y, uint32_t z,
                                            uint32_t width,
                                            uint32_t height) {
    return (uint64_t)z * (uint64_t)width * (uint64_t)height +
           (uint64_t)y * (uint64_t)width + (uint64_t)x;
}

POROSITY_FIELD_STATIC uint64_t integral_index(uint32_t x, uint32_t y,
                                              uint32_t z,
                                              uint32_t width_plus_one,
                                              uint32_t height_plus_one) {
    return (uint64_t)z * (uint64_t)width_plus_one *
               (uint64_t)height_plus_one +
           (uint64_t)y * (uint64_t)width_plus_one + (uint64_t)x;
}

POROSITY_FIELD_STATIC void porosity_field_die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

POROSITY_FIELD_STATIC uint32_t filter_size_to_voxels(double filter_size,
                                                     double voxel_size) {
    uint32_t size_in_voxels;

    if (voxel_size <= 0.0) {
        porosity_field_die("voxel_size must be positive");
    }
    if (filter_size <= 0.0) {
        porosity_field_die("filter_size must be positive");
    }

    size_in_voxels = (uint32_t)llround(filter_size / voxel_size);
    if (size_in_voxels == 0) {
        size_in_voxels = 1;
    }
    if ((size_in_voxels % 2U) == 0U) {
        size_in_voxels += 1;
    }

    return size_in_voxels;
}

POROSITY_FIELD_STATIC uint64_t *build_void_integral_volume(
    const uint8_t *binary, uint32_t width, uint32_t height, uint32_t depth) {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t width_plus_one = width + 1;
    uint32_t height_plus_one = height + 1;
    uint32_t depth_plus_one = depth + 1;
    uint64_t *integral = (uint64_t *)calloc(
        (size_t)width_plus_one * (size_t)height_plus_one *
            (size_t)depth_plus_one,
        sizeof(uint64_t));

    if (integral == NULL) {
        porosity_field_die("Unable to allocate integral volume");
    }

    for (z = 1; z <= depth; ++z) {
        for (y = 1; y <= height; ++y) {
            for (x = 1; x <= width; ++x) {
                uint64_t idx = integral_index(x, y, z, width_plus_one,
                                              height_plus_one);
                uint64_t source_idx =
                    volume_index(x - 1, y - 1, z - 1, width, height);
                uint64_t voxel = binary[source_idx] == 0 ? 1 : 0;

                integral[idx] =
                    voxel +
                    integral[integral_index(x - 1, y, z, width_plus_one,
                                            height_plus_one)] +
                    integral[integral_index(x, y - 1, z, width_plus_one,
                                            height_plus_one)] +
                    integral[integral_index(x, y, z - 1, width_plus_one,
                                            height_plus_one)] -
                    integral[integral_index(x - 1, y - 1, z, width_plus_one,
                                            height_plus_one)] -
                    integral[integral_index(x - 1, y, z - 1, width_plus_one,
                                            height_plus_one)] -
                    integral[integral_index(x, y - 1, z - 1, width_plus_one,
                                            height_plus_one)] +
                    integral[integral_index(x - 1, y - 1, z - 1, width_plus_one,
                                            height_plus_one)];
            }
        }
    }

    return integral;
}

POROSITY_FIELD_STATIC uint64_t sum_voids_in_box(const uint64_t *integral,
                                                uint32_t width,
                                                uint32_t height,
                                                uint32_t depth, uint32_t x0,
                                                uint32_t y0, uint32_t z0,
                                                uint32_t x1, uint32_t y1,
                                                uint32_t z1) {
    uint32_t width_plus_one = width + 1;
    uint32_t height_plus_one = height + 1;

    (void)depth;

    return integral[integral_index(x1, y1, z1, width_plus_one,
                                   height_plus_one)] -
           integral[integral_index(x0, y1, z1, width_plus_one,
                                   height_plus_one)] -
           integral[integral_index(x1, y0, z1, width_plus_one,
                                   height_plus_one)] -
           integral[integral_index(x1, y1, z0, width_plus_one,
                                   height_plus_one)] +
           integral[integral_index(x0, y0, z1, width_plus_one,
                                   height_plus_one)] +
           integral[integral_index(x0, y1, z0, width_plus_one,
                                   height_plus_one)] +
           integral[integral_index(x1, y0, z0, width_plus_one,
                                   height_plus_one)] -
           integral[integral_index(x0, y0, z0, width_plus_one,
                                   height_plus_one)];
}

POROSITY_FIELD_STATIC char *build_porosity_field_output_path(
    const char *stack_path) {
    const char *suffix = "_porosity_field.tif";
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

int compute_porosity_field(const uint8_t *binary, uint32_t width,
                           uint32_t height, uint32_t depth,
                           double filter_size, double voxel_size,
                           PorosityField *field) {
    uint32_t wx;
    uint32_t wy;
    uint32_t wz;
    uint32_t window_size;
    uint32_t field_width;
    uint32_t field_height;
    uint32_t field_depth;
    uint64_t *integral;

    if (binary == NULL || field == NULL) {
        return EXIT_FAILURE;
    }

    window_size = filter_size_to_voxels(filter_size, voxel_size);
    integral = build_void_integral_volume(binary, width, height, depth);

    field_width = (width + window_size - 1U) / window_size;
    field_height = (height + window_size - 1U) / window_size;
    field_depth = (depth + window_size - 1U) / window_size;

    field->width = field_width;
    field->height = field_height;
    field->depth = field_depth;
    field->values = (double *)malloc((size_t)field_width *
                                     (size_t)field_height *
                                     (size_t)field_depth * sizeof(double));
    if (field->values == NULL) {
        free(integral);
        porosity_field_die("Unable to allocate porosity field");
    }

    for (wz = 0; wz < field_depth; ++wz) {
        uint32_t z0 = wz * window_size;
        uint32_t z1 = z0 + window_size < depth ? z0 + window_size : depth;

        for (wy = 0; wy < field_height; ++wy) {
            uint32_t y0 = wy * window_size;
            uint32_t y1 = y0 + window_size < height ? y0 + window_size : height;

            for (wx = 0; wx < field_width; ++wx) {
                uint32_t x0 = wx * window_size;
                uint32_t x1 = x0 + window_size < width ? x0 + window_size : width;
                uint64_t voids = sum_voids_in_box(integral, width, height, depth,
                                                  x0, y0, z0, x1, y1, z1);
                uint64_t voxel_count = (uint64_t)(x1 - x0) * (uint64_t)(y1 - y0) *
                                       (uint64_t)(z1 - z0);

                field->values[volume_index(wx, wy, wz, field_width,
                                           field_height)] =
                    voxel_count == 0 ? 0.0 : (double)voids / (double)voxel_count;
            }
        }
    }

    free(integral);
    return EXIT_SUCCESS;
}

void free_porosity_field(PorosityField *field) {
    if (field == NULL) {
        return;
    }

    free(field->values);
    field->values = NULL;
    field->width = 0;
    field->height = 0;
    field->depth = 0;
}

int write_porosity_field_tiff(const char *path, const PorosityField *field) {
    TIFF *tif;
    uint32_t wz;

    if (path == NULL || field == NULL || field->values == NULL) {
        return EXIT_FAILURE;
    }

    tif = TIFFOpen(path, "w");
    if (tif == NULL) {
        return EXIT_FAILURE;
    }

    for (wz = 0; wz < field->depth; ++wz) {
        float *scanline;
        uint32_t wy;

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, field->width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, field->height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, field->height);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

        scanline = (float *)_TIFFmalloc((tmsize_t)field->width * sizeof(float));
        if (scanline == NULL) {
            TIFFClose(tif);
            return EXIT_FAILURE;
        }

        for (wy = 0; wy < field->height; ++wy) {
            uint32_t wx;

            for (wx = 0; wx < field->width; ++wx) {
                scanline[wx] =
                    (float)field->values[volume_index(wx, wy, wz, field->width,
                                                      field->height)];
            }

            if (TIFFWriteScanline(tif, scanline, wy, 0) < 0) {
                _TIFFfree(scanline);
                TIFFClose(tif);
                return EXIT_FAILURE;
            }
        }

        _TIFFfree(scanline);

        if (wz + 1U < field->depth) {
            TIFFWriteDirectory(tif);
        }
    }

    TIFFClose(tif);
    return EXIT_SUCCESS;
}

int process_porosity_field(const char *stack_path, double voxel_resolution,
                           double filter_resolution, FILE *out) {
    BinaryVolume volume;
    PorosityField field;
    char *output_path;

    if (stack_path == NULL || out == NULL) {
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

    output_path = build_porosity_field_output_path(stack_path);
    if (output_path == NULL) {
        fprintf(stderr, "Failed to build porosity field output path for %s\n",
                stack_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (write_porosity_field_tiff(output_path, &field) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to write porosity field TIFF for %s\n",
                stack_path);
        free(output_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "Porosity field: %s\n"
                "  Output: %s\n"
                "  Dimensions: %u x %u x %u\n"
                "  Voxel resolution: %.6f\n"
                "  Filter resolution: %.6f\n\n",
                stack_path, output_path, field.width, field.height,
                field.depth, voxel_resolution, filter_resolution) < 0) {
        fprintf(stderr, "Failed to write porosity field summary for %s\n",
                stack_path);
        free(output_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    free(output_path);
    free_porosity_field(&field);
    free_binary_volume(&volume);
    return EXIT_SUCCESS;
}
