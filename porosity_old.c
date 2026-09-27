#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "porosity.h"

#ifndef POROSITY_STATIC
#define POROSITY_STATIC static
#endif

POROSITY_STATIC uint64_t volume_index(uint32_t x, uint32_t y, uint32_t z,
                                      uint32_t width, uint32_t height) {
    return (uint64_t)z * (uint64_t)width * (uint64_t)height +
           (uint64_t)y * (uint64_t)width + (uint64_t)x;
}

POROSITY_STATIC uint64_t integral_index(uint32_t x, uint32_t y, uint32_t z,
                                        uint32_t width_plus_one,
                                        uint32_t height_plus_one) {
    return (uint64_t)z * (uint64_t)width_plus_one * (uint64_t)height_plus_one +
           (uint64_t)y * (uint64_t)width_plus_one + (uint64_t)x;
}

POROSITY_STATIC void porosity_die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

POROSITY_STATIC uint32_t filter_size_to_voxels(double filter_size,
                                               double voxel_size) {
    uint32_t size_in_voxels;

    if (voxel_size <= 0.0) {
        porosity_die("voxel_size must be positive");
    }
    if (filter_size <= 0.0) {
        porosity_die("filter_size must be positive");
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

POROSITY_STATIC uint64_t *build_void_integral_volume(const uint8_t *binary,
                                                     uint32_t width,
                                                     uint32_t height,
                                                     uint32_t depth) {
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
        porosity_die("Unable to allocate integral volume");
    }

    for (z = 1; z <= depth; ++z) {
        for (y = 1; y <= height; ++y) {
            for (x = 1; x <= width; ++x) {
                uint64_t idx = integral_index(x, y, z, width_plus_one,
                                              height_plus_one);
                uint64_t source_idx =
                    volume_index(x - 1, y - 1, z - 1, width, height);
                uint64_t voxel = binary[source_idx] == 0 ? 0 : 1;

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

POROSITY_STATIC uint64_t sum_voids_in_box(const uint64_t *integral,
                                          uint32_t width, uint32_t height,
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

int compute_porosity_field(const uint8_t *binary, uint32_t width,
                           uint32_t height, uint32_t depth,
                           double filter_size, double voxel_size,
                           PorosityField *field) {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t window_size;
    uint32_t radius;
    uint64_t *integral;

    if (binary == NULL || field == NULL) {
        return EXIT_FAILURE;
    }

    window_size = filter_size_to_voxels(filter_size, voxel_size);
    radius = window_size / 2U;
    integral = build_void_integral_volume(binary, width, height, depth);

    field->width = width;
    field->height = height;
    field->depth = depth;
    field->values = (double *)malloc((size_t)width * (size_t)height *
                                     (size_t)depth * sizeof(double));
    if (field->values == NULL) {
        free(integral);
        porosity_die("Unable to allocate porosity field");
    }

    for (z = 0; z < depth; ++z) {
        for (y = 0; y < height; ++y) {
            for (x = 0; x < width; ++x) {
                uint32_t x0 = x > radius ? x - radius : 0;
                uint32_t y0 = y > radius ? y - radius : 0;
                uint32_t z0 = z > radius ? z - radius : 0;
                uint32_t x1 = x + radius + 1U < width ? x + radius + 1U : width;
                uint32_t y1 =
                    y + radius + 1U < height ? y + radius + 1U : height;
                uint32_t z1 = z + radius + 1U < depth ? z + radius + 1U : depth;
                uint64_t voids = sum_voids_in_box(integral, width, height, depth,
                                                  x0, y0, z0, x1, y1, z1);
                uint64_t voxel_count = (uint64_t)(x1 - x0) * (uint64_t)(y1 - y0) *
                                       (uint64_t)(z1 - z0);

                field->values[volume_index(x, y, z, width, height)] =
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

CovarianceSample compute_porosity_field_covariance(const PorosityField *field,
                                                   int lag_x, int lag_y,
                                                   int lag_z) {
    CovarianceSample sample;
    uint32_t x;
    uint32_t y;
    uint32_t z;
    double sum_a = 0.0;
    double sum_b = 0.0;
    double sum_ab = 0.0;

    memset(&sample, 0, sizeof(sample));
    sample.lag_x = lag_x;
    sample.lag_y = lag_y;
    sample.lag_z = lag_z;

    if (field == NULL || field->values == NULL) {
        return sample;
    }

    for (z = 0; z < field->depth; ++z) {
        int shifted_z = (int)z + lag_z;
        if (shifted_z < 0 || shifted_z >= (int)field->depth) {
            continue;
        }

        for (y = 0; y < field->height; ++y) {
            int shifted_y = (int)y + lag_y;
            if (shifted_y < 0 || shifted_y >= (int)field->height) {
                continue;
            }

            for (x = 0; x < field->width; ++x) {
                int shifted_x = (int)x + lag_x;
                if (shifted_x < 0 || shifted_x >= (int)field->width) {
                    continue;
                }

                double a = field->values[volume_index(x, y, z, field->width,
                                                      field->height)];
                double b = field->values[volume_index((uint32_t)shifted_x,
                                                      (uint32_t)shifted_y,
                                                      (uint32_t)shifted_z,
                                                      field->width,
                                                      field->height)];

                sum_a += a;
                sum_b += b;
                sum_ab += a * b;
                sample.samples += 1;
            }
        }
    }

    if (sample.samples > 0) {
        double mean_a = sum_a / (double)sample.samples;
        double mean_b = sum_b / (double)sample.samples;

        sample.mean = 0.5 * (mean_a + mean_b);
        sample.covariance =
            (sum_ab / (double)sample.samples) - (mean_a * mean_b);
    }

    return sample;
}

int compute_covariance_curve_along_axis(const PorosityField *field, char axis,
                                        uint32_t max_lag,
                                        CovarianceSample *samples_out) {
    uint32_t lag;

    if (field == NULL || samples_out == NULL) {
        return EXIT_FAILURE;
    }

    for (lag = 0; lag <= max_lag; ++lag) {
        if (axis == 'x') {
            samples_out[lag] =
                compute_porosity_field_covariance(field, (int)lag, 0, 0);
        } else if (axis == 'y') {
            samples_out[lag] =
                compute_porosity_field_covariance(field, 0, (int)lag, 0);
        } else if (axis == 'z') {
            samples_out[lag] =
                compute_porosity_field_covariance(field, 0, 0, (int)lag);
        } else {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
