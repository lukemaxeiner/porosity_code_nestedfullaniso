#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "porosity_covariance_lag.h"
#include "tiff.h"

static uint64_t porosity_field_index(uint32_t wx, uint32_t wy, uint32_t wz,
                                     uint32_t width, uint32_t height) {
    return (uint64_t)wz * (uint64_t)width * (uint64_t)height +
           (uint64_t)wy * (uint64_t)width + (uint64_t)wx;
}

static size_t porosity_covariance_grid_sample_count(const PorosityField *field) {
    if (field == NULL) {
        return 0U;
    }

    return (size_t)field->width * (size_t)field->height * (size_t)field->depth;
}

static int ensure_build_covariance_dir(const char *stack_path) {
    const char *slash;
    size_t dir_length;
    char *folder_path;
    int status;

    if (stack_path == NULL) {
        return EXIT_FAILURE;
    }

    slash = strrchr(stack_path, '/');
    dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    folder_path = (char *)malloc(dir_length + strlen("/BUILD_COVARIANCE") + 1U);
    if (folder_path == NULL) {
        return EXIT_FAILURE;
    }

    if (slash != NULL) {
        memcpy(folder_path, stack_path, dir_length);
    } else {
        folder_path[0] = '.';
    }
    memcpy(folder_path + dir_length, "/BUILD_COVARIANCE",
           strlen("/BUILD_COVARIANCE") + 1U);

    status = mkdir(folder_path, 0777);
    if (status != 0 && errno != EEXIST) {
        free(folder_path);
        return EXIT_FAILURE;
    }

    free(folder_path);
    return EXIT_SUCCESS;
}

static char *build_covariance_output_path(const char *stack_path,
                                          const char *suffix) {
    const char *extension = strrchr(stack_path, '.');
    const char *slash = strrchr(stack_path, '/');
    const char *file_name = slash != NULL ? slash + 1 : stack_path;
    size_t dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    size_t file_base_length =
        extension != NULL && extension > file_name
            ? (size_t)(extension - file_name)
            : strlen(file_name);
    char *output_path;
    size_t total_length;

    if (stack_path == NULL) {
        return NULL;
    }
    total_length = dir_length + strlen("/BUILD_COVARIANCE/") + file_base_length +
                   strlen(suffix) + 1U;
    output_path = (char *)malloc(total_length);
    if (output_path == NULL) {
        return NULL;
    }

    if (slash != NULL) {
        memcpy(output_path, stack_path, dir_length);
    } else {
        output_path[0] = '.';
    }
    memcpy(output_path + dir_length, "/BUILD_COVARIANCE/",
           strlen("/BUILD_COVARIANCE/"));
    memcpy(output_path + dir_length + strlen("/BUILD_COVARIANCE/"), file_name,
           file_base_length);
    memcpy(output_path + dir_length + strlen("/BUILD_COVARIANCE/") +
               file_base_length,
           suffix, strlen(suffix) + 1U);
    return output_path;
}

static char *build_covariance_lag_output_path(const char *stack_path,
                                              PorosityLagMode lag_mode,
                                              PorosityAxis axis,
                                              double filter_resolution) {
    const char *axis_name = "";
    char suffix[128];

    if (lag_mode == POROSITY_LAG_MODE_CURVE) {
        if (axis == POROSITY_AXIS_X) {
            axis_name = "x";
        } else if (axis == POROSITY_AXIS_Y) {
            axis_name = "y";
        } else {
            axis_name = "z";
        }

        if (snprintf(suffix, sizeof(suffix),
                     "_porosity_covariance_lag_axis_%s_filter_%g.txt",
                     axis_name, filter_resolution) >= (int)sizeof(suffix)) {
            return NULL;
        }
    } else {
        if (snprintf(suffix, sizeof(suffix), "%s",
                     "_porosity_covariance_lag.txt") >= (int)sizeof(suffix)) {
            return NULL;
        }
    }

    return build_covariance_output_path(stack_path, suffix);
}

static char *build_covariance_lag_paraview_path(const char *stack_path,
                                                double filter_resolution) {
    char suffix[128];

    if (snprintf(suffix, sizeof(suffix),
                 "_porosity_covariance_3dlag_filter_%g.vtk",
                 filter_resolution) >= (int)sizeof(suffix)) {
        return NULL;
    }

    return build_covariance_output_path(stack_path, suffix);
}

static double compute_porosity_field_mean(const PorosityField *field) {
    uint32_t wz;
    uint32_t wy;
    uint32_t wx;
    double sum = 0.0;
    uint64_t count = 0;

    if (field == NULL || field->values == NULL) {
        return 0.0;
    }

    for (wz = 0; wz < field->depth; ++wz) {
        for (wy = 0; wy < field->height; ++wy) {
            for (wx = 0; wx < field->width; ++wx) {
                sum += field->values[porosity_field_index(wx, wy, wz,
                                                          field->width,
                                                          field->height)];
                count += 1;
            }
        }
    }

    return count == 0 ? 0.0 : sum / (double)count;
}

static double compute_porosity_field_variance(const PorosityField *field) {
    uint32_t wz;
    uint32_t wy;
    uint32_t wx;
    uint64_t count = 0;
    double mean;
    double sum_sq = 0.0;

    if (field == NULL || field->values == NULL) {
        return 0.0;
    }

    mean = compute_porosity_field_mean(field);

    for (wz = 0; wz < field->depth; ++wz) {
        for (wy = 0; wy < field->height; ++wy) {
            for (wx = 0; wx < field->width; ++wx) {
                double phi = field->values[porosity_field_index(
                    wx, wy, wz, field->width, field->height)];
                double centered = phi - mean;

                sum_sq += centered * centered;
                count += 1;
            }
        }
    }

    return count == 0 ? 0.0 : sum_sq / (double)count;
}

int compute_porosity_covariance_sample(const PorosityField *field, int lag_x,
                                       int lag_y, int lag_z,
                                       PorosityCovarianceSample *sample) {
    uint32_t wz;
    uint32_t wy;
    uint32_t wx;
    double mu;
    double variance;
    double sum_centered_product = 0.0;

    if (field == NULL || sample == NULL || field->values == NULL) {
        return EXIT_FAILURE;
    }

    memset(sample, 0, sizeof(*sample));
    sample->lag_x = lag_x;
    sample->lag_y = lag_y;
    sample->lag_z = lag_z;
    mu = compute_porosity_field_mean(field);
    variance = compute_porosity_field_variance(field);

    for (wz = 0; wz < field->depth; ++wz) {
        int shifted_z = (int)wz + lag_z;

        if (shifted_z < 0 || shifted_z >= (int)field->depth) {
            continue;
        }

        for (wy = 0; wy < field->height; ++wy) {
            int shifted_y = (int)wy + lag_y;

            if (shifted_y < 0 || shifted_y >= (int)field->height) {
                continue;
            }

            for (wx = 0; wx < field->width; ++wx) {
                int shifted_x = (int)wx + lag_x;
                double phi_x;
                double phi_xr;

                if (shifted_x < 0 || shifted_x >= (int)field->width) {
                    continue;
                }

                phi_x = field->values[porosity_field_index(wx, wy, wz,
                                                           field->width,
                                                           field->height)];
                phi_xr = field->values[porosity_field_index(
                    (uint32_t)shifted_x, (uint32_t)shifted_y,
                    (uint32_t)shifted_z, field->width, field->height)];

                sum_centered_product += (phi_x - mu) * (phi_xr - mu);
                sample->samples += 1;
            }
        }
    }

    sample->mu = mu;

    if (sample->samples > 0) {
        sample->covariance = sum_centered_product / (double)sample->samples;
        if (variance > 0.0) {
            sample->covariance /= variance;
        } else {
            sample->covariance = 0.0;
        }
    }

    return EXIT_SUCCESS;
}

uint32_t porosity_covariance_axis_lag_count(const PorosityField *field,
                                            PorosityAxis axis) {
    if (field == NULL) {
        return 0;
    }

    if (axis == POROSITY_AXIS_X) {
        return field->width;
    }
    if (axis == POROSITY_AXIS_Y) {
        return field->height;
    }
    if (axis == POROSITY_AXIS_Z) {
        return field->depth;
    }

    return 0;
}

int compute_porosity_covariance_curve(const PorosityField *field,
                                      PorosityAxis axis,
                                      PorosityCovarianceSample *samples_out) {
    uint32_t lag;
    uint32_t lag_count;

    if (field == NULL || samples_out == NULL) {
        return EXIT_FAILURE;
    }

    lag_count = porosity_covariance_axis_lag_count(field, axis);
    if (lag_count == 0) {
        return EXIT_FAILURE;
    }

    for (lag = 0; lag < lag_count; ++lag) {
        int lag_x = 0;
        int lag_y = 0;
        int lag_z = 0;

        if (axis == POROSITY_AXIS_X) {
            lag_x = (int)lag;
        } else if (axis == POROSITY_AXIS_Y) {
            lag_y = (int)lag;
        } else {
            lag_z = (int)lag;
        }

        if (compute_porosity_covariance_sample(field, lag_x, lag_y, lag_z,
                                               &samples_out[lag]) !=
            EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static int compute_porosity_covariance_grid(
    const PorosityField *field, PorosityCovarianceSample *samples_out) {
    uint32_t lag_z;
    uint32_t lag_y;
    uint32_t lag_x;

    if (field == NULL || samples_out == NULL) {
        return EXIT_FAILURE;
    }

    for (lag_z = 0; lag_z < field->depth; ++lag_z) {
        for (lag_y = 0; lag_y < field->height; ++lag_y) {
            for (lag_x = 0; lag_x < field->width; ++lag_x) {
                uint64_t sample_index = porosity_field_index(
                    lag_x, lag_y, lag_z, field->width, field->height);

                if (compute_porosity_covariance_sample(
                        field, (int)lag_x, (int)lag_y, (int)lag_z,
                        &samples_out[sample_index]) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
            }
        }
    }

    return EXIT_SUCCESS;
}

int write_porosity_covariance_table(FILE *out,
                                    const PorosityCovarianceSample *samples,
                                    size_t sample_count,
                                    double filter_resolution) {
    size_t i;

    if (out == NULL || samples == NULL) {
        return EXIT_FAILURE;
    }

    if (fprintf(out, "lag_x lag_y lag_z samples mu covariance\n") < 0) {
        return EXIT_FAILURE;
    }

    for (i = 0; i < sample_count; ++i) {
        if (fprintf(out, "%.12f %.12f %.12f %llu %.12f %.12f\n",
                    (double)samples[i].lag_x * filter_resolution,
                    (double)samples[i].lag_y * filter_resolution,
                    (double)samples[i].lag_z * filter_resolution,
                    (unsigned long long)samples[i].samples, samples[i].mu,
                    samples[i].covariance) < 0) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static int write_porosity_covariance_vtk(FILE *out,
                                         const PorosityCovarianceSample *samples,
                                         uint32_t width, uint32_t height,
                                         uint32_t depth,
                                         double filter_resolution) {
    size_t point_count;
    size_t i;

    if (out == NULL || samples == NULL) {
        return EXIT_FAILURE;
    }

    point_count = (size_t)width * (size_t)height * (size_t)depth;

    if (fprintf(out,
                "# vtk DataFile Version 3.0\n"
                "Porosity covariance 3D lag grid\n"
                "ASCII\n"
                "DATASET STRUCTURED_POINTS\n"
                "DIMENSIONS %u %u %u\n"
                "ORIGIN 0 0 0\n"
                "SPACING %.12f %.12f %.12f\n"
                "POINT_DATA %zu\n"
                "SCALARS covariance double 1\n"
                "LOOKUP_TABLE default\n",
                width, height, depth, filter_resolution, filter_resolution,
                filter_resolution, point_count) < 0) {
        return EXIT_FAILURE;
    }

    for (i = 0; i < point_count; ++i) {
        if (fprintf(out, "%.12f\n", samples[i].covariance) < 0) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int process_porosity_covariance_lag(const char *stack_path,
                                    double voxel_resolution,
                                    double filter_resolution,
                                    PorosityLagMode lag_mode,
                                    PorosityAxis axis, FILE *out) {
    BinaryVolume volume;
    PorosityField field;
    PorosityCovarianceSample *samples = NULL;
    size_t sample_count = 0U;
    char *output_path;
    char *paraview_path;
    FILE *output_file;
    FILE *paraview_file;

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
    paraview_path = NULL;
    output_file = NULL;
    paraview_file = NULL;

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

    if (ensure_build_covariance_dir(stack_path) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to create BUILD_COVARIANCE for %s\n",
                stack_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    output_path = build_covariance_lag_output_path(stack_path, lag_mode, axis,
                                                   filter_resolution);
    if (output_path == NULL) {
        fprintf(stderr, "Failed to build lag covariance output path for %s\n",
                stack_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    output_file = fopen(output_path, "w");
    if (output_file == NULL) {
        fprintf(stderr, "Failed to open lag covariance output for %s\n",
                stack_path);
        free(output_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (lag_mode == POROSITY_LAG_MODE_CURVE) {
        sample_count = porosity_covariance_axis_lag_count(&field, axis);
        if (sample_count == 0) {
            fprintf(stderr, "Invalid covariance axis for %s\n", stack_path);
            fclose(output_file);
            free(output_path);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        samples = (PorosityCovarianceSample *)malloc((size_t)sample_count *
                                                     sizeof(*samples));
        if (samples == NULL) {
            fprintf(stderr, "Failed to allocate covariance curve for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        if (compute_porosity_covariance_curve(&field, axis, samples) !=
            EXIT_SUCCESS) {
            fprintf(stderr, "Failed to compute covariance curve for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        if (write_porosity_covariance_table(output_file, samples, sample_count,
                                            filter_resolution) !=
            EXIT_SUCCESS) {
            fprintf(stderr, "Failed to write covariance curve output for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }
    } else if (lag_mode == POROSITY_LAG_MODE_3DLAG) {
        sample_count = porosity_covariance_grid_sample_count(&field);
        if (sample_count == 0U) {
            fprintf(stderr, "Invalid covariance grid dimensions for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        paraview_path =
            build_covariance_lag_paraview_path(stack_path, filter_resolution);
        if (paraview_path == NULL) {
            fprintf(stderr, "Failed to build ParaView output path for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        samples = (PorosityCovarianceSample *)malloc(sample_count *
                                                     sizeof(*samples));
        if (samples == NULL) {
            fprintf(stderr, "Failed to allocate covariance grid for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free(paraview_path);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        if (compute_porosity_covariance_grid(&field, samples) != EXIT_SUCCESS) {
            fprintf(stderr, "Failed to compute lag covariance grid for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free(paraview_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        if (write_porosity_covariance_table(output_file, samples, sample_count,
                                            filter_resolution) !=
            EXIT_SUCCESS) {
            fprintf(stderr, "Failed to write lag covariance grid for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free(paraview_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        paraview_file = fopen(paraview_path, "w");
        if (paraview_file == NULL) {
            fprintf(stderr, "Failed to open ParaView output for %s\n",
                    stack_path);
            fclose(output_file);
            free(output_path);
            free(paraview_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        if (write_porosity_covariance_vtk(paraview_file, samples, field.width,
                                          field.height, field.depth,
                                          filter_resolution) != EXIT_SUCCESS) {
            fprintf(stderr, "Failed to write ParaView covariance grid for %s\n",
                    stack_path);
            fclose(paraview_file);
            fclose(output_file);
            free(output_path);
            free(paraview_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }

        fclose(paraview_file);
        paraview_file = NULL;
    } else {
        fprintf(stderr, "Invalid lag covariance mode for %s\n", stack_path);
        fclose(output_file);
        free(output_path);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    fclose(output_file);

    if (lag_mode == POROSITY_LAG_MODE_CURVE) {
        const char *axis_name = axis == POROSITY_AXIS_X
                                    ? "X"
                                    : (axis == POROSITY_AXIS_Y ? "Y" : "Z");

        if (fprintf(out,
                    "Lag covariance curve: %s\n"
                    "  Output: %s\n"
                    "  Field dimensions: %u x %u x %u\n"
                    "  Axis: %s\n"
                    "  Samples: %zu\n\n",
                    stack_path, output_path, field.width, field.height,
                    field.depth, axis_name, sample_count) < 0) {
            fprintf(stderr,
                    "Failed to write lag covariance curve summary for %s\n",
                    stack_path);
            free(output_path);
            free(samples);
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return EXIT_FAILURE;
        }
    } else if (fprintf(out,
                       "Lag covariance: %s\n"
                       "  Output: %s\n"
                       "  ParaView: %s\n"
                       "  Field dimensions: %u x %u x %u\n"
                       "  Lag grid dimensions: %u x %u x %u\n"
                       "  Lag spacing: %.12f\n"
                       "  Covariance samples: %llu\n\n",
                       stack_path, output_path,
                       paraview_path != NULL ? paraview_path : "(not written)",
                       field.width, field.height, field.depth, field.width,
                       field.height, field.depth,
                       filter_resolution, (unsigned long long)sample_count) < 0) {
        fprintf(stderr, "Failed to write lag covariance summary for %s\n",
                stack_path);
        free(output_path);
        free(paraview_path);
        free(samples);
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    free(output_path);
    free(paraview_path);
    free(samples);
    free_porosity_field(&field);
    free_binary_volume(&volume);
    return EXIT_SUCCESS;
}
