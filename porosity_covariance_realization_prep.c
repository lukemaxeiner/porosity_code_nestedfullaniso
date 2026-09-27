#include <errno.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "porosity_covariance_realization_prep.h"
#include "porosity_covariance_realization_generation.h"
#include "porosity_covariance_tail_noise.h"
#include "porosity_field.h"
#include "porosity_normal_score.h"
#include "tiff.h"

static size_t covariance_grid_index(uint32_t x, uint32_t y, uint32_t z,
                                    uint32_t width, uint32_t height) {
    return (size_t)z * (size_t)width * (size_t)height +
           (size_t)y * (size_t)width + (size_t)x;
}

static void free_covariance_grid_3d(CovarianceGrid3D *grid) {
    if (grid == NULL) {
        return;
    }

    free(grid->values);
    grid->values = NULL;
    grid->width = 0;
    grid->height = 0;
    grid->depth = 0;
    grid->mean = 0.0;
    grid->variance = 0.0;
}

static void free_symmetric_covariance_grid_3d(SymmetricCovarianceGrid3D *grid) {
    if (grid == NULL) {
        return;
    }

    free(grid->values);
    grid->values = NULL;
    grid->width = 0;
    grid->height = 0;
    grid->depth = 0;
    grid->origin_x = 0;
    grid->origin_y = 0;
    grid->origin_z = 0;
    grid->mean = 0.0;
    grid->variance = 0.0;
}

static void free_realization_prep_result(PorosityRealizationPrepResult *result) {
    if (result == NULL) {
        return;
    }

    free(result->samples);
    result->samples = NULL;
    result->lag_width = 0;
    result->lag_height = 0;
    result->lag_depth = 0;
    result->filter_resolution = 0.0;
    result->variance = 0.0;
    free_covariance_grid_3d(&result->covariance_grid);
    free_symmetric_covariance_grid_3d(&result->symmetric_grid);
    memset(&result->tail_noise, 0, sizeof(result->tail_noise));
    memset(&result->negative_covariance, 0,
           sizeof(result->negative_covariance));
    memset(&result->anisotropy, 0, sizeof(result->anisotropy));
    memset(&result->exponential_model, 0, sizeof(result->exponential_model));
    memset(&result->exponential_validation, 0,
           sizeof(result->exponential_validation));
    memset(&result->simulation_grid, 0, sizeof(result->simulation_grid));
    memset(&result->simulation_grid_compatibility, 0,
           sizeof(result->simulation_grid_compatibility));
    free_covariance_grid_3d(&result->evaluated_simulation_covariance_grid);
    memset(&result->evaluated_grid_validation, 0,
           sizeof(result->evaluated_grid_validation));
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

static char *build_build_covariance_output_path(const char *stack_path,
                                                const char *suffix) {
    const char *extension = strrchr(stack_path, '.');
    const char *slash = strrchr(stack_path, '/');
    const char *file_name = slash != NULL ? slash + 1 : stack_path;
    size_t dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    size_t file_base_length =
        extension != NULL && extension > file_name
            ? (size_t)(extension - file_name)
            : strlen(file_name);
    size_t total_length = dir_length + strlen("/BUILD_COVARIANCE/") +
                          file_base_length + strlen(suffix) + 1U;
    char *output_path = (char *)malloc(total_length);

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

static char *build_3dlag_table_path(const char *stack_path) {
    return build_build_covariance_output_path(stack_path,
                                              "_porosity_covariance_lag.txt");
}

static char *build_fitted_model_path(const char *stack_path) {
    return build_build_covariance_output_path(
        stack_path, "_fitted_covariance_model.txt");
}

//rebuilds the block porosity field from the tiff and returns its variance
//if normal_score_enabled, the field is rank-transformed first so the variance
//matches the field of the 3DLAG table
static double compute_porosity_field_variance_from_stack(const char *stack_path,
                                                         double voxel_resolution,
                                                         double filter_resolution,
                                                         int normal_score_enabled,
                                                         double *mean_out) {
    BinaryVolume volume;
    PorosityField field;
    uint64_t count = 0U;
    double mean = 0.0;
    double sum_sq = 0.0;

    volume.width = 0U;
    volume.height = 0U;
    volume.depth = 0U;
    volume.voxels = NULL;
    field.width = 0U;
    field.height = 0U;
    field.depth = 0U;
    field.values = NULL;

    if (load_binary_volume_from_tiff(stack_path, &volume) != EXIT_SUCCESS) {
        return -1.0;
    }
    if (compute_porosity_field(volume.voxels, volume.width, volume.height,
                               volume.depth, filter_resolution,
                               voxel_resolution, &field) != EXIT_SUCCESS) {
        free_binary_volume(&volume);
        return -1.0;
    }

    if (normal_score_enabled &&
        normal_score_transform_field(&field, NULL) != EXIT_SUCCESS){
            free_porosity_field(&field);
            free_binary_volume(&volume);
            return -1.0;
    }

    for (uint32_t z = 0; z < field.depth; ++z) {
        for (uint32_t y = 0; y < field.height; ++y) {
            for (uint32_t x = 0; x < field.width; ++x) {
                double value = field.values[covariance_grid_index(
                    x, y, z, field.width, field.height)];
                mean += value;
                count += 1U;
            }
        }
    }

    if (count > 0U) {
        mean /= (double)count;
    }

    if (mean_out != NULL) {
        *mean_out = mean;
    }

    for (uint32_t z = 0; z < field.depth; ++z) {
        for (uint32_t y = 0; y < field.height; ++y) {
            for (uint32_t x = 0; x < field.width; ++x) {
                double value = field.values[covariance_grid_index(
                    x, y, z, field.width, field.height)];
                double centered = value - mean;
                sum_sq += centered * centered;
            }
        }
    }

    free_porosity_field(&field);
    free_binary_volume(&volume);
    return count > 0U ? sum_sq / (double)count : -1.0;
}

static double covariance_radius_from_index(uint32_t x, uint32_t y, uint32_t z,
                                           double dx, double dy, double dz) {
    double rx = (double)x * dx;
    double ry = (double)y * dy;
    double rz = (double)z * dz;

    return sqrt(rx * rx + ry * ry + rz * rz);
}

static double interpolate_efold_distance(const CovarianceGrid3D *grid,
                                         double spacing_x, double spacing_y,
                                         double spacing_z, char axis,
                                         double target_value) {
    uint32_t limit = 0U;
    uint32_t i;
    double prev_distance = 0.0;
    double prev_value;

    if (grid == NULL || grid->values == NULL) {
        return 0.0;
    }

    prev_value = grid->values[covariance_grid_index(0U, 0U, 0U, grid->width,
                                                    grid->height)];

    if (axis == 'x') {
        limit = grid->width;
    } else if (axis == 'y') {
        limit = grid->height;
    } else if (axis == 'z') {
        limit = grid->depth;
    }

    for (i = 1U; i < limit; ++i) {
        uint32_t x = axis == 'x' ? i : 0U;
        uint32_t y = axis == 'y' ? i : 0U;
        uint32_t z = axis == 'z' ? i : 0U;
        double distance =
            covariance_radius_from_index(x, y, z, spacing_x, spacing_y, spacing_z);
        double value =
            grid->values[covariance_grid_index(x, y, z, grid->width,
                                               grid->height)];

        if (value <= target_value && prev_value > target_value) {
            double denom = prev_value - value;
            if (fabs(denom) < 1.0e-18) {
                return distance;
            }
            return prev_distance +
                   (distance - prev_distance) * (prev_value - target_value) /
                       denom;
        }

        prev_distance = distance;
        prev_value = value;
    }

    return prev_distance;
}

/*
 * Step 1: Reconstruct C(r) from the 3DLAG output samples.
 *
 * The current 3DLAG code stores a normalized covariance in sample.covariance,
 * i.e. covariance divided by the field variance. This function reconstructs the
 * dimensional covariance field:
 *
 *     C(r) = variance * sample_covariance(r)
 *
 * It expects the samples to be laid out exactly as the current 3DLAG writer
 * produces them: x-first, then y, then z.
 */
int reconstruct_covariance_grid_from_3dlag_samples(
    const PorosityCovarianceSample *samples, uint32_t width, uint32_t height,
    uint32_t depth, double variance, CovarianceGrid3D *grid_out) {
    size_t point_count;
    uint32_t lag_z;
    uint32_t lag_y;
    uint32_t lag_x;

    if (samples == NULL || grid_out == NULL) {
        return EXIT_FAILURE;
    }
    if (variance < 0.0) {
        return EXIT_FAILURE;
    }

    memset(grid_out, 0, sizeof(*grid_out));
    point_count = (size_t)width * (size_t)height * (size_t)depth;
    grid_out->values = (double *)malloc(point_count * sizeof(double));
    if (grid_out->values == NULL) {
        return EXIT_FAILURE;
    }

    grid_out->width = width;
    grid_out->height = height;
    grid_out->depth = depth;
    grid_out->variance = variance;
    grid_out->mean = point_count == 0U ? 0.0 : samples[0].mu;

    for (lag_z = 0; lag_z < depth; ++lag_z) {
        for (lag_y = 0; lag_y < height; ++lag_y) {
            for (lag_x = 0; lag_x < width; ++lag_x) {
                size_t index = covariance_grid_index(lag_x, lag_y, lag_z,
                                                     width, height);
                const PorosityCovarianceSample *sample = &samples[index];

                if (sample->lag_x != (int)lag_x || sample->lag_y != (int)lag_y ||
                    sample->lag_z != (int)lag_z) {
                    free_covariance_grid_3d(grid_out);
                    return EXIT_FAILURE;
                }

                /* sample->covariance already represents mean-centered,
                 * variance-normalized covariance, so reconstruction only
                 * needs C(r) = variance * sample->covariance.
                 * The mean is carried through grid_out for later realization
                 * generation, not for this reconstruction step. */
                grid_out->values[index] = variance * sample->covariance;
            }
        }
    }

    return EXIT_SUCCESS;
}

/*
 * Step 2: Symmetrize the covariance field onto a signed lag grid.
 *
 * IMPORTANT:
 * The available 3DLAG output only covers the nonnegative lag octant. Extending
 * that to all signed lags requires an additional symmetry assumption. This
 * function uses reflection symmetry across each coordinate axis:
 *
 *     C(rx, ry, rz) = C(|rx|, |ry|, |rz|)
 *
 * This is stronger than the general stationarity identity C(r) = C(-r), but it
 * is the most direct way to expand the measured octant into a full grid for
 * later spectral work.
 * Note from Poovathingal: Function written by Codex on 04/30/26 
 * and has not been verified by Poovathingal for accuracy.
 */
int symmetrize_covariance_grid_by_reflection(
    const CovarianceGrid3D *source, SymmetricCovarianceGrid3D *grid_out) {
    uint32_t z;
    uint32_t y;
    uint32_t x;
    uint32_t out_width;
    uint32_t out_height;
    uint32_t out_depth;
    size_t point_count;

    if (source == NULL || grid_out == NULL || source->values == NULL) {
        return EXIT_FAILURE;
    }

    out_width = source->width == 0U ? 0U : (2U * source->width - 1U);
    out_height = source->height == 0U ? 0U : (2U * source->height - 1U);
    out_depth = source->depth == 0U ? 0U : (2U * source->depth - 1U);
    point_count =
        (size_t)out_width * (size_t)out_height * (size_t)out_depth;

    memset(grid_out, 0, sizeof(*grid_out));
    grid_out->values = (double *)malloc(point_count * sizeof(double));
    if (grid_out->values == NULL) {
        return EXIT_FAILURE;
    }

    grid_out->width = out_width;
    grid_out->height = out_height;
    grid_out->depth = out_depth;
    grid_out->origin_x = (int)source->width - 1;
    grid_out->origin_y = (int)source->height - 1;
    grid_out->origin_z = (int)source->depth - 1;
    grid_out->mean = source->mean;
    grid_out->variance = source->variance;

    for (z = 0; z < out_depth; ++z) {
        int signed_z = (int)z - grid_out->origin_z;
        uint32_t reflected_z = (uint32_t)(signed_z < 0 ? -signed_z : signed_z);

        for (y = 0; y < out_height; ++y) {
            int signed_y = (int)y - grid_out->origin_y;
            uint32_t reflected_y =
                (uint32_t)(signed_y < 0 ? -signed_y : signed_y);

            for (x = 0; x < out_width; ++x) {
                int signed_x = (int)x - grid_out->origin_x;
                uint32_t reflected_x =
                    (uint32_t)(signed_x < 0 ? -signed_x : signed_x);
                size_t source_index = covariance_grid_index(
                    reflected_x, reflected_y, reflected_z, source->width,
                    source->height);
                size_t out_index = covariance_grid_index(x, y, z, out_width,
                                                        out_height);

                grid_out->values[out_index] = source->values[source_index];
            }
        }
    }

    return EXIT_SUCCESS;
}

int analyze_covariance_tail_noise(
    const PorosityCovarianceSample *samples, const CovarianceGrid3D *grid,
    double filter_resolution, CovarianceTailNoiseSummary *summary_out) {
    return check_covariance_tail_noise(samples, grid, filter_resolution,
                                       summary_out);
}

/*
 * Axis-aligned fully anisotropic exponential covariance model:
 *
 *     C(rx, ry, rz) = sill * exp(-sqrt(
 *         (rx / Lx)^2 +
 *         (ry / Ly)^2 +
 *         (rz / Lz)^2
 *     ))
 *
 * The function name is kept for now so the existing call sites still compile.
 */
static double evaluate_anisotropic_exponential_component(
    double rx, double ry, double rz,
    double lx, double ly, double lz) {
    double scaled_radius;

    if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0) {
        return 0.0;
    }

    scaled_radius = sqrt((rx * rx) / (lx * lx) +
                         (ry * ry) / (ly * ly) +
                         (rz * rz) / (lz * lz));

    return exp(-scaled_radius);
}

static double evaluate_exponential_covariance_model_at_lag(
    const ExponentialCovarianceModel *model,
    double rx, double ry, double rz) {
    double normalized = 0.0;

    if (model == NULL || model->sill <= 0.0) {
        return 0.0;
    }

    if (model->component_count > 0U) {
        uint32_t k;
        uint32_t count = model->component_count;

        if (count > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
            count = POROSITY_NESTED_EXPONENTIAL_COMPONENTS;
        }

        for (k = 0U; k < count; ++k) {
            double weight = model->component_weight[k];

            if (weight <= 0.0) {
                continue;
            }

            normalized +=
                weight *
                evaluate_anisotropic_exponential_component(
                    rx, ry, rz,
                    model->component_x_correlation_length[k],
                    model->component_y_correlation_length[k],
                    model->component_z_correlation_length[k]);
        }

        return model->sill * normalized;
    }

    /*
     * Backward-compatible fallback for old single-component model files.
     */
    return model->sill *
           evaluate_anisotropic_exponential_component(
               rx, ry, rz,
               model->x_correlation_length,
               model->y_correlation_length,
               model->z_correlation_length);
}

static int compute_covariance_fit_metrics(
    const CovarianceGrid3D *grid,
    const PorosityCovarianceSample *samples,
    double filter_resolution,
    double recommended_max_radius,
    ExponentialCovarianceModel *model) {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    double rss = 0.0;
    double weighted_rss = 0.0;
    double weight_sum = 0.0;
    double rss_x = 0.0;
    double rss_y = 0.0;
    double rss_z = 0.0;
    double rss_small = 0.0;
    double rss_recommended = 0.0;
    uint64_t eval_count = 0U;
    uint64_t count_x = 0U;
    uint64_t count_y = 0U;
    uint64_t count_z = 0U;
    uint64_t count_small = 0U;
    uint64_t count_recommended = 0U;

    if (grid == NULL || grid->values == NULL || model == NULL ||
        filter_resolution <= 0.0 || grid->width == 0U ||
        grid->height == 0U || grid->depth == 0U) {
        return EXIT_FAILURE;
    }

    for (z = 0U; z < grid->depth; ++z) {
        for (y = 0U; y < grid->height; ++y) {
            for (x = 0U; x < grid->width; ++x) {
                size_t index = covariance_grid_index(x, y, z, grid->width,
                                                     grid->height);
                double empirical = grid->values[index];
                double rx = (double)x * filter_resolution;
                double ry = (double)y * filter_resolution;
                double rz = (double)z * filter_resolution;
                double radius = sqrt(rx * rx + ry * ry + rz * rz);
                double fitted;
                double residual;
                double weight;

                if (!isfinite(empirical)) {
                    continue;
                }

                fitted = evaluate_exponential_covariance_model_at_lag(
                    model, rx, ry, rz);
                residual = empirical - fitted;

                rss += residual * residual;
                eval_count += 1U;

                weight = samples != NULL ? (double)samples[index].samples : 1.0;
                if (weight <= 0.0) {
                    weight = 1.0;
                }

                weighted_rss += weight * residual * residual;
                weight_sum += weight;

                if (x == 0U && y == 0U && z == 0U) {
                    double scale = fabs(empirical) > 1.0e-12
                                       ? fabs(empirical)
                                       : 1.0;
                    model->zero_lag_error = fabs(residual);
                    model->zero_lag_relative_error =
                        model->zero_lag_error / scale;
                }

                if (x > 0U && y == 0U && z == 0U) {
                    rss_x += residual * residual;
                    count_x += 1U;
                } else if (x == 0U && y > 0U && z == 0U) {
                    rss_y += residual * residual;
                    count_y += 1U;
                } else if (x == 0U && y == 0U && z > 0U) {
                    rss_z += residual * residual;
                    count_z += 1U;
                }

                if (radius > 0.0 && radius <= 2.0 * filter_resolution) {
                    rss_small += residual * residual;
                    count_small += 1U;
                }

                if (recommended_max_radius > 0.0 &&
                    radius <= recommended_max_radius) {
                    rss_recommended += residual * residual;
                    count_recommended += 1U;
                }
            }
        }
    }

    if (eval_count == 0U) {
        return EXIT_FAILURE;
    }

    model->rms_error = sqrt(rss / (double)eval_count);
    model->weighted_rms_error =
        weight_sum > 0.0 ? sqrt(weighted_rss / weight_sum) : 0.0;

    model->directional_rms_x =
        count_x > 0U ? sqrt(rss_x / (double)count_x) : 0.0;
    model->directional_rms_y =
        count_y > 0U ? sqrt(rss_y / (double)count_y) : 0.0;
    model->directional_rms_z =
        count_z > 0U ? sqrt(rss_z / (double)count_z) : 0.0;

    model->small_lag_rms_error =
        count_small > 0U ? sqrt(rss_small / (double)count_small) : 0.0;

    model->recommended_radius_rms_error =
        count_recommended > 0U
            ? sqrt(rss_recommended / (double)count_recommended)
            : 0.0;

    {
        const double parameter_count =
            model->model_type ==
                    POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL
                ? 8.0
                : 4.0;
        double mean_squared_error = rss / (double)eval_count;

        if (mean_squared_error <= 1.0e-300) {
            mean_squared_error = 1.0e-300;
        }

        model->aic =
            (double)eval_count * log(mean_squared_error) +
            2.0 * parameter_count;

        model->bic =
            (double)eval_count * log(mean_squared_error) +
            parameter_count * log((double)eval_count);
    }

    return EXIT_SUCCESS;
}

static double covariance_fit_priority_score(
    const ExponentialCovarianceModel *model) {
    double scale;
    double small;
    double weighted;
    double recommended;

    if (model == NULL || model->sill <= 0.0) {
        return INFINITY;
    }

    scale = model->sill;
    small = model->small_lag_rms_error / scale;
    weighted = model->weighted_rms_error / scale;
    recommended = model->recommended_radius_rms_error / scale;

    return 4.0 * small * small +
           weighted * weighted +
           0.5 * recommended * recommended;
}

int fit_nested_full_anisotropic_exponential_covariance_model(
    const CovarianceGrid3D *grid,
    const PorosityCovarianceSample *samples,
    double filter_resolution,
    double recommended_max_radius,
    ExponentialCovarianceModel *model_out) {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    double sill;
    double sum_x_d2 = 0.0;
    double sum_x_dlog = 0.0;
    double sum_y_d2 = 0.0;
    double sum_y_dlog = 0.0;
    double sum_z_d2 = 0.0;
    double sum_z_dlog = 0.0;
    double slope_x;
    double slope_y;
    double slope_z;
    uint64_t fit_count_x = 0U;
    uint64_t fit_count_y = 0U;
    uint64_t fit_count_z = 0U;

    if (grid == NULL || grid->values == NULL || model_out == NULL ||
        filter_resolution <= 0.0 || grid->width == 0U ||
        grid->height == 0U || grid->depth == 0U) {
        return EXIT_FAILURE;
    }

    memset(model_out, 0, sizeof(*model_out));

    sill = grid->values[covariance_grid_index(0U, 0U, 0U, grid->width,
                                              grid->height)];
    if (sill <= 0.0 || !isfinite(sill)) {
        return EXIT_FAILURE;
    }

    for (x = 1U; x < grid->width; ++x) {
        size_t index = covariance_grid_index(x, 0U, 0U, grid->width,
                                             grid->height);
        double covariance = grid->values[index];

        if (covariance > 0.0 && isfinite(covariance)) {
            double distance = (double)x * filter_resolution;
            double log_ratio = log(covariance / sill);
            double weight = samples != NULL ? (double)samples[index].samples : 1.0;

            if (weight <= 0.0) {
                weight = 1.0;
            }

            if (isfinite(log_ratio)) {
                sum_x_d2 += weight * distance * distance;
                sum_x_dlog += weight * distance * log_ratio;
                fit_count_x += 1U;
            }
        }
    }

    for (y = 1U; y < grid->height; ++y) {
        size_t index = covariance_grid_index(0U, y, 0U, grid->width,
                                             grid->height);
        double covariance = grid->values[index];

        if (covariance > 0.0 && isfinite(covariance)) {
            double distance = (double)y * filter_resolution;
            double log_ratio = log(covariance / sill);
            double weight = samples != NULL ? (double)samples[index].samples : 1.0;

            if (weight <= 0.0) {
                weight = 1.0;
            }

            if (isfinite(log_ratio)) {
                sum_y_d2 += weight * distance * distance;
                sum_y_dlog += weight * distance * log_ratio;
                fit_count_y += 1U;
            }
        }
    }

    for (z = 1U; z < grid->depth; ++z) {
        size_t index = covariance_grid_index(0U, 0U, z, grid->width,
                                             grid->height);
        double covariance = grid->values[index];

        if (covariance > 0.0 && isfinite(covariance)) {
            double distance = (double)z * filter_resolution;
            double log_ratio = log(covariance / sill);
            double weight = samples != NULL ? (double)samples[index].samples : 1.0;

            if (weight <= 0.0) {
                weight = 1.0;
            }

            if (isfinite(log_ratio)) {
                sum_z_d2 += weight * distance * distance;
                sum_z_dlog += weight * distance * log_ratio;
                fit_count_z += 1U;
            }
        }
    }

    if (fit_count_x == 0U || fit_count_y == 0U || fit_count_z == 0U ||
        sum_x_d2 <= 0.0 || sum_y_d2 <= 0.0 || sum_z_d2 <= 0.0) {
        return EXIT_FAILURE;
    }

    slope_x = sum_x_dlog / sum_x_d2;
    slope_y = sum_y_dlog / sum_y_d2;
    slope_z = sum_z_dlog / sum_z_d2;

    if (slope_x >= 0.0 || slope_y >= 0.0 || slope_z >= 0.0) {
        return EXIT_FAILURE;
    }

        {
        const double short_scales[] = {
            0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70
        };
        const double long_scales[] = {
            1.00, 1.25, 1.50, 2.00, 3.00, 4.00, 6.00
        };
        const double short_weights[] = {
            0.50, 0.60, 0.70, 0.80, 0.90, 0.95
        };
        const double base_lx = -1.0 / slope_x;
        const double base_ly = -1.0 / slope_y;
        const double base_lz = -1.0 / slope_z;
        ExponentialCovarianceModel best_model;
        double best_score = INFINITY;
        int have_best = 0;
        size_t is;
        size_t il;
        size_t iw;

        memset(&best_model, 0, sizeof(best_model));

        for (is = 0U; is < sizeof(short_scales) / sizeof(short_scales[0]); ++is) {
            for (il = 0U; il < sizeof(long_scales) / sizeof(long_scales[0]); ++il) {
                for (iw = 0U; iw < sizeof(short_weights) / sizeof(short_weights[0]); ++iw) {
                    ExponentialCovarianceModel candidate;
                    double score;

                    memset(&candidate, 0, sizeof(candidate));

                    candidate.model_type =
                        POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL;
                    candidate.sill = sill;
                    candidate.component_count = 2U;

                    candidate.component_weight[0] = short_weights[iw];
                    candidate.component_weight[1] = 1.0 - short_weights[iw];

                    candidate.component_x_correlation_length[0] =
                        base_lx * short_scales[is];
                    candidate.component_y_correlation_length[0] =
                        base_ly * short_scales[is];
                    candidate.component_z_correlation_length[0] =
                        base_lz * short_scales[is];

                    candidate.component_x_correlation_length[1] =
                        base_lx * long_scales[il];
                    candidate.component_y_correlation_length[1] =
                        base_ly * long_scales[il];
                    candidate.component_z_correlation_length[1] =
                        base_lz * long_scales[il];

                    candidate.x_correlation_length =
                        candidate.component_x_correlation_length[1];
                    candidate.y_correlation_length =
                        candidate.component_y_correlation_length[1];
                    candidate.z_correlation_length =
                        candidate.component_z_correlation_length[1];
                    candidate.correlation_length =
                        candidate.x_correlation_length;
                    candidate.in_plane_correlation_length =
                        candidate.x_correlation_length;
                    candidate.through_thickness_correlation_length =
                        candidate.z_correlation_length;

                    candidate.sample_count =
                        fit_count_x + fit_count_y + fit_count_z;

                    if (compute_covariance_fit_metrics(
                            grid, samples, filter_resolution,
                            recommended_max_radius,
                            &candidate) != EXIT_SUCCESS) {
                        continue;
                    }

                    score = covariance_fit_priority_score(&candidate);

                    if (!have_best || score < best_score) {
                        best_model = candidate;
                        best_score = score;
                        have_best = 1;
                    }
                }
            }
        }

        if (!have_best) {
            return EXIT_FAILURE;
        }

        *model_out = best_model;
    }

    return EXIT_SUCCESS;
}

int validate_exponential_covariance_model(
    const CovarianceGrid3D *grid, const ExponentialCovarianceModel *model,
    ExponentialCovarianceValidation *validation_out) {
    double zero_lag_covariance;
    double fitted_zero_lag;
    double scale;
    double error_for_ratio;
    uint32_t k;
    uint32_t count;
    double weight_sum = 0.0;

    if (grid == NULL || grid->values == NULL || model == NULL ||
        validation_out == NULL || grid->width == 0U || grid->height == 0U ||
        grid->depth == 0U) {
        return EXIT_FAILURE;
    }

    memset(validation_out, 0, sizeof(*validation_out));

    validation_out->sill_is_positive = model->sill > 0.0;
    validation_out->component_weights_are_valid = 1;
    validation_out->correlation_length_is_positive = 1;

    count = model->component_count;

    if (count > 0U) {
        if (count > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
            validation_out->component_weights_are_valid = 0;
            count = POROSITY_NESTED_EXPONENTIAL_COMPONENTS;
        }

        for (k = 0U; k < count; ++k) {
            double weight = model->component_weight[k];

            if (!isfinite(weight) || weight < -1.0e-12) {
                validation_out->component_weights_are_valid = 0;
            }

            if (weight > 0.0) {
                if (model->component_x_correlation_length[k] <= 0.0 ||
                    model->component_y_correlation_length[k] <= 0.0 ||
                    model->component_z_correlation_length[k] <= 0.0 ||
                    !isfinite(model->component_x_correlation_length[k]) ||
                    !isfinite(model->component_y_correlation_length[k]) ||
                    !isfinite(model->component_z_correlation_length[k])) {
                    validation_out->correlation_length_is_positive = 0;
                }
            }

            weight_sum += weight;
        }

        if (fabs(weight_sum - 1.0) > 1.0e-6) {
            validation_out->component_weights_are_valid = 0;
        }
    } else {
        validation_out->correlation_length_is_positive =
            model->x_correlation_length > 0.0 &&
            model->y_correlation_length > 0.0 &&
            model->z_correlation_length > 0.0;
    }

    validation_out->has_enough_samples =
        model->model_type ==
                POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL
            ? model->sample_count >= 8U
            : model->sample_count >= 3U;

    zero_lag_covariance = grid->values[covariance_grid_index(
        0U, 0U, 0U, grid->width, grid->height)];

    fitted_zero_lag =
        evaluate_exponential_covariance_model_at_lag(model, 0.0, 0.0, 0.0);

    scale = fabs(zero_lag_covariance) > 1.0e-12
                ? fabs(zero_lag_covariance)
                : 1.0;

    validation_out->zero_lag_relative_error =
        fabs(fitted_zero_lag - zero_lag_covariance) / scale;

    validation_out->zero_lag_matches_sill =
        validation_out->zero_lag_relative_error <= 0.10;

    error_for_ratio = model->weighted_rms_error > 0.0
                          ? model->weighted_rms_error
                          : model->rms_error;

    validation_out->rms_error_ratio = error_for_ratio / scale;

    validation_out->rms_error_is_reasonable =
        validation_out->rms_error_ratio <= 0.25;

    validation_out->is_valid =
        validation_out->sill_is_positive &&
        validation_out->correlation_length_is_positive &&
        validation_out->component_weights_are_valid &&
        validation_out->has_enough_samples &&
        validation_out->zero_lag_matches_sill &&
        validation_out->rms_error_is_reasonable;

    return EXIT_SUCCESS;
}

int analyze_negative_covariance_lobes(const CovarianceGrid3D *grid,
                                      double filter_resolution,
                                      NegativeCovarianceSummary *summary_out) {
    uint32_t z;
    uint32_t y;
    uint32_t x;

    if (grid == NULL || grid->values == NULL || summary_out == NULL ||
        filter_resolution <= 0.0) {
        return EXIT_FAILURE;
    }

    memset(summary_out, 0, sizeof(*summary_out));
    summary_out->min_covariance = 0.0;

    for (z = 0; z < grid->depth; ++z) {
        for (y = 0; y < grid->height; ++y) {
            for (x = 0; x < grid->width; ++x) {
                size_t index =
                    covariance_grid_index(x, y, z, grid->width, grid->height);
                double value = grid->values[index];

                if (x == 0U && y == 0U && z == 0U) {
                    continue;
                }

                if (value < 0.0) {
                    double radius =
                        covariance_radius_from_index(x, y, z,
                                                     filter_resolution,
                                                     filter_resolution,
                                                     filter_resolution);
                    summary_out->negative_count += 1U;
                    if (!summary_out->has_negative_lobes ||
                        value < summary_out->min_covariance) {
                        summary_out->min_covariance = value;
                        summary_out->most_negative_radius = radius;
                    }
                    summary_out->has_negative_lobes = 1;
                }
            }
        }
    }

    return EXIT_SUCCESS;
}

int analyze_covariance_anisotropy(const CovarianceGrid3D *grid,
                                  double filter_resolution,
                                  CovarianceAnisotropySummary *summary_out) {
    double target;
    double min_efold;
    double max_efold;

    if (grid == NULL || grid->values == NULL || summary_out == NULL ||
        filter_resolution <= 0.0) {
        return EXIT_FAILURE;
    }

    memset(summary_out, 0, sizeof(*summary_out));
    summary_out->zero_lag_covariance = grid->values[covariance_grid_index(
        0U, 0U, 0U, grid->width, grid->height)];
    target = summary_out->zero_lag_covariance / exp(1.0);

    summary_out->x_efold_distance = interpolate_efold_distance(
        grid, filter_resolution, filter_resolution, filter_resolution, 'x',
        target);
    summary_out->y_efold_distance = interpolate_efold_distance(
        grid, filter_resolution, filter_resolution, filter_resolution, 'y',
        target);
    summary_out->z_efold_distance = interpolate_efold_distance(
        grid, filter_resolution, filter_resolution, filter_resolution, 'z',
        target);

    min_efold = summary_out->x_efold_distance;
    max_efold = summary_out->x_efold_distance;
    if (summary_out->y_efold_distance < min_efold) {
        min_efold = summary_out->y_efold_distance;
    }
    if (summary_out->z_efold_distance < min_efold) {
        min_efold = summary_out->z_efold_distance;
    }
    if (summary_out->y_efold_distance > max_efold) {
        max_efold = summary_out->y_efold_distance;
    }
    if (summary_out->z_efold_distance > max_efold) {
        max_efold = summary_out->z_efold_distance;
    }

    summary_out->max_to_min_ratio =
        min_efold > 1.0e-12 ? max_efold / min_efold : 0.0;
    summary_out->is_anisotropic = summary_out->max_to_min_ratio > 1.25;

    return EXIT_SUCCESS;
}

int assess_realization_grid_compatibility(
    const PorosityRealizationGridSpec *simulation_grid,
    const ExponentialCovarianceModel *model,
    SimulationGridCompatibilitySummary *summary_out) {
    if (simulation_grid == NULL || model == NULL || summary_out == NULL ||
        model->x_correlation_length <= 0.0 ||
        model->y_correlation_length <= 0.0 ||
        model->z_correlation_length <= 0.0) {
        return EXIT_FAILURE;
    }

    memset(summary_out, 0, sizeof(*summary_out));

    summary_out->domain_x_over_correlation_length =
        simulation_grid->dim_x / model->x_correlation_length;
    summary_out->domain_y_over_correlation_length =
        simulation_grid->dim_y / model->y_correlation_length;
    summary_out->domain_z_over_correlation_length =
        simulation_grid->dim_z / model->z_correlation_length;

    summary_out->enough_domain_x =
        summary_out->domain_x_over_correlation_length >= 4.0;
    summary_out->enough_domain_y =
        summary_out->domain_y_over_correlation_length >= 4.0;
    summary_out->enough_domain_z =
        summary_out->domain_z_over_correlation_length >= 4.0;

    summary_out->overall_adequate =
        summary_out->enough_domain_x && summary_out->enough_domain_y &&
        summary_out->enough_domain_z;

    return EXIT_SUCCESS;
}

int evaluate_exponential_covariance_on_realization_grid(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid,
    CovarianceGrid3D *grid_out) {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    double dx;
    double dy;
    double dz;
    size_t point_count;

    if (model == NULL || simulation_grid == NULL || grid_out == NULL ||
        model->sill <= 0.0 ||
        model->x_correlation_length <= 0.0 ||
        model->y_correlation_length <= 0.0 ||
        model->z_correlation_length <= 0.0 ||
        simulation_grid->nx == 0U || simulation_grid->ny == 0U ||
        simulation_grid->nz == 0U) {
        return EXIT_FAILURE;
    }

    memset(grid_out, 0, sizeof(*grid_out));

    point_count = (size_t)simulation_grid->nx *
                  (size_t)simulation_grid->ny *
                  (size_t)simulation_grid->nz;

    grid_out->values = (double *)malloc(point_count * sizeof(double));
    if (grid_out->values == NULL) {
        return EXIT_FAILURE;
    }

    grid_out->width = simulation_grid->nx;
    grid_out->height = simulation_grid->ny;
    grid_out->depth = simulation_grid->nz;
    grid_out->mean = 0.0;
    grid_out->variance = model->sill;

    dx = simulation_grid->nx > 1U
             ? simulation_grid->dim_x / (double)(simulation_grid->nx - 1U)
             : simulation_grid->dim_x;
    dy = simulation_grid->ny > 1U
             ? simulation_grid->dim_y / (double)(simulation_grid->ny - 1U)
             : simulation_grid->dim_y;
    dz = simulation_grid->nz > 1U
             ? simulation_grid->dim_z / (double)(simulation_grid->nz - 1U)
             : simulation_grid->dim_z;

    for (z = 0U; z < simulation_grid->nz; ++z) {
        for (y = 0U; y < simulation_grid->ny; ++y) {
            for (x = 0U; x < simulation_grid->nx; ++x) {
                double rx = (double)x * dx;
                double ry = (double)y * dy;
                double rz = (double)z * dz;
                size_t index = covariance_grid_index(x, y, z,
                                                     simulation_grid->nx,
                                                     simulation_grid->ny);

                grid_out->values[index] =
                    evaluate_exponential_covariance_model_at_lag(
                        model, rx, ry, rz);
            }
        }
    }

    return EXIT_SUCCESS;
}

int validate_evaluated_covariance_grid(
    const CovarianceGrid3D *grid, const ExponentialCovarianceModel *model,
    EvaluatedCovarianceGridValidation *validation_out) {
    uint32_t z;
    uint32_t y;
    uint32_t x;
    double zero_lag;
    double scale;

    if (grid == NULL || grid->values == NULL || model == NULL ||
        validation_out == NULL || grid->width == 0U || grid->height == 0U ||
        grid->depth == 0U) {
        return EXIT_FAILURE;
    }

    memset(validation_out, 0, sizeof(*validation_out));
    validation_out->all_finite = 1;
    validation_out->all_nonnegative = 1;
    validation_out->min_value = grid->values[0];
    validation_out->max_value = grid->values[0];

    for (z = 0; z < grid->depth; ++z) {
        for (y = 0; y < grid->height; ++y) {
            for (x = 0; x < grid->width; ++x) {
                size_t index =
                    covariance_grid_index(x, y, z, grid->width, grid->height);
                double value = grid->values[index];

                if (!isfinite(value)) {
                    validation_out->all_finite = 0;
                }
                if (value < 0.0) {
                    validation_out->all_nonnegative = 0;
                }
                if (value < validation_out->min_value) {
                    validation_out->min_value = value;
                }
                if (value > validation_out->max_value) {
                    validation_out->max_value = value;
                }
            }
        }
    }

    zero_lag = grid->values[covariance_grid_index(0U, 0U, 0U, grid->width,
                                                  grid->height)];
    scale = fabs(model->sill) > 1.0e-12 ? fabs(model->sill) : 1.0;
    validation_out->zero_lag_relative_error =
        fabs(zero_lag - model->sill) / scale;
    validation_out->zero_lag_matches_sill =
        validation_out->zero_lag_relative_error <= 1.0e-8;
    validation_out->is_valid = validation_out->all_finite &&
                               validation_out->all_nonnegative &&
                               validation_out->zero_lag_matches_sill;

    return EXIT_SUCCESS;
}

static const char *covariance_model_type_label(
    const ExponentialCovarianceModel *model) {
    if (model != NULL &&
        model->model_type ==
            POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL) {
        return "nested fully anisotropic exponential";
    }

    return "fully anisotropic exponential";
}

static uint32_t covariance_model_component_count(
    const ExponentialCovarianceModel *model) {
    uint32_t count;

    if (model == NULL || model->component_count == 0U) {
        return 1U;
    }

    count = model->component_count;
    if (count > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
        count = POROSITY_NESTED_EXPONENTIAL_COMPONENTS;
    }

    return count;
}

static int write_covariance_model_components(
    FILE *out,
    const ExponentialCovarianceModel *model,
    const char *indent) {
    uint32_t count;
    uint32_t k;

    if (out == NULL || model == NULL || indent == NULL) {
        return EXIT_FAILURE;
    }

    count = covariance_model_component_count(model);

    for (k = 0U; k < count; ++k) {
        if (fprintf(out,
                    "%sComponent %u weight: %.12f\n"
                    "%sComponent %u X/Y/Z lengths: %.12f %.12f %.12f\n",
                    indent, k + 1U, model->component_weight[k],
                    indent, k + 1U,
                    model->component_x_correlation_length[k],
                    model->component_y_correlation_length[k],
                    model->component_z_correlation_length[k]) < 0) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int write_realization_prep_report(FILE *out,
                                  const PorosityRealizationPrepResult *result) {
    if (out == NULL || result == NULL) {
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "Realization Prep Report\n"
                "  Lag grid: %u x %u x %u\n"
                "  Filter resolution: %.12f\n"
                "  Variance: %.12f\n"
                "  Mean: %.12f\n"
                "  Tail noisy: %s\n"
                "  Tail noise fraction: %.12f\n"
                "  Recommended max radius: %.12f\n"
                "  Negative lobes: %s\n"
                "  Negative count: %llu\n"
                "  Min covariance: %.12f\n"
                "  X/Y/Z e-fold distances: %.12f %.12f %.12f\n"
                "  Anisotropy ratio: %.12f\n"
                "  Exponential model type: %s\n"
                "  Exponential sill: %.12f\n"
                "  Legacy X/Y/Z correlation lengths: %.12f %.12f %.12f\n"
                "  Exponential RMS error: %.12f\n"
                "  Exponential fit valid: %s\n"
                "  Simulation grid: %u x %u x %u\n"
                "  Simulation dimensions: %.12f %.12f %.12f\n"
                "  Domain/correlation-length ratios: %.12f %.12f %.12f\n"
                "  Simulation grid adequate: %s\n"
                "  Evaluated covariance grid valid: %s\n",
                result->lag_width, result->lag_height, result->lag_depth,
                result->filter_resolution, result->variance,
                result->covariance_grid.mean,
                result->tail_noise.tail_is_noisy ? "YES" : "NO",
                result->tail_noise.tail_noise_fraction,
                result->tail_noise.recommended_max_radius,
                result->negative_covariance.has_negative_lobes ? "YES" : "NO",
                (unsigned long long)result->negative_covariance.negative_count,
                result->negative_covariance.min_covariance,
                result->anisotropy.x_efold_distance,
                result->anisotropy.y_efold_distance,
                result->anisotropy.z_efold_distance,
                result->anisotropy.max_to_min_ratio,
                covariance_model_type_label(&result->exponential_model),
                result->exponential_model.sill,
                result->exponential_model.x_correlation_length,
                result->exponential_model.y_correlation_length,
                result->exponential_model.z_correlation_length,
                result->exponential_model.rms_error,
                result->exponential_validation.is_valid ? "YES" : "NO",
                result->simulation_grid.nx, result->simulation_grid.ny,
                result->simulation_grid.nz, result->simulation_grid.dim_x,
                result->simulation_grid.dim_y, result->simulation_grid.dim_z,
                result->simulation_grid_compatibility
                    .domain_x_over_correlation_length,
                result->simulation_grid_compatibility
                    .domain_y_over_correlation_length,
                result->simulation_grid_compatibility
                    .domain_z_over_correlation_length,
                result->simulation_grid_compatibility.overall_adequate ? "YES"
                                                                       : "NO",
                result->evaluated_grid_validation.is_valid ? "YES" : "NO") <
        0) {
        return EXIT_FAILURE;
    }
    if (write_covariance_model_components(
            out, &result->exponential_model, "  ") != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int write_directional_covariance_fit_table(
    FILE *out,
    const CovarianceGrid3D *grid,
    const PorosityCovarianceSample *samples,
    double filter_resolution,
    const ExponentialCovarianceModel *model) {
    uint32_t lag;

    if (out == NULL || grid == NULL || grid->values == NULL ||
        model == NULL || filter_resolution <= 0.0) {
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "axis lag distance samples empirical_covariance "
                "model_covariance residual\n") < 0) {
        return EXIT_FAILURE;
    }

    for (lag = 0U; lag < grid->width; ++lag) {
        size_t index = covariance_grid_index(lag, 0U, 0U, grid->width,
                                             grid->height);
        double distance = (double)lag * filter_resolution;
        double empirical = grid->values[index];
        double fitted = evaluate_exponential_covariance_model_at_lag(
            model, distance, 0.0, 0.0);
        double residual = empirical - fitted;
        uint64_t sample_count = samples != NULL ? samples[index].samples : 0U;

        if (fprintf(out, "x %u %.12f %llu %.12f %.12f %.12f\n",
                    lag, distance, (unsigned long long)sample_count,
                    empirical, fitted, residual) < 0) {
            return EXIT_FAILURE;
        }
    }

    for (lag = 0U; lag < grid->height; ++lag) {
        size_t index = covariance_grid_index(0U, lag, 0U, grid->width,
                                             grid->height);
        double distance = (double)lag * filter_resolution;
        double empirical = grid->values[index];
        double fitted = evaluate_exponential_covariance_model_at_lag(
            model, 0.0, distance, 0.0);
        double residual = empirical - fitted;
        uint64_t sample_count = samples != NULL ? samples[index].samples : 0U;

        if (fprintf(out, "y %u %.12f %llu %.12f %.12f %.12f\n",
                    lag, distance, (unsigned long long)sample_count,
                    empirical, fitted, residual) < 0) {
            return EXIT_FAILURE;
        }
    }

    for (lag = 0U; lag < grid->depth; ++lag) {
        size_t index = covariance_grid_index(0U, 0U, lag, grid->width,
                                             grid->height);
        double distance = (double)lag * filter_resolution;
        double empirical = grid->values[index];
        double fitted = evaluate_exponential_covariance_model_at_lag(
            model, 0.0, 0.0, distance);
        double residual = empirical - fitted;
        uint64_t sample_count = samples != NULL ? samples[index].samples : 0U;

        if (fprintf(out, "z %u %.12f %llu %.12f %.12f %.12f\n",
                    lag, distance, (unsigned long long)sample_count,
                    empirical, fitted, residual) < 0) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int write_fitted_covariance_model(FILE *out,
                                  const FittedCovarianceModelArtifact *model) {
    const char *model_type_name;
    uint32_t component_count;
    uint32_t k;

    if (out == NULL || model == NULL) {
        return EXIT_FAILURE;
    }

    model_type_name =
        model->exponential_model.model_type ==
                POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL
            ? "NESTED_FULLY_ANISOTROPIC_EXPONENTIAL"
            : "FULLY_ANISOTROPIC_EXPONENTIAL";

    component_count = model->exponential_model.component_count;
    if (component_count == 0U) {
        component_count = 1U;
    }
    if (component_count > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
        component_count = POROSITY_NESTED_EXPONENTIAL_COMPONENTS;
    }

    if (fprintf(out,
                "MODEL_TYPE = %s\n"
                "NORMAL_SCORE_TRANSFORM = %d\n"
                "MODELED_FIELD = %s\n"
                "MEAN_POROSITY = %.17g\n"
                "VARIANCE = %.17g\n"
                "FILTER_RESOLUTION = %.17g\n"
                "SILL = %.17g\n"
                "COMPONENT_COUNT = %u\n",
                model_type_name,
                model->normal_score_transform ? 1 : 0,
                model->normal_score_transform ? "NORMAL_SCORE" : "POROSITY",
                model->mean_porosity,
                model->variance,
                model->filter_resolution,
                model->exponential_model.sill,
                component_count) < 0) {
        return EXIT_FAILURE;
    }

    for (k = 0U; k < component_count; ++k) {
        if (fprintf(out,
                    "COMPONENT_%u_WEIGHT = %.17g\n"
                    "COMPONENT_%u_X_CORRELATION_LENGTH = %.17g\n"
                    "COMPONENT_%u_Y_CORRELATION_LENGTH = %.17g\n"
                    "COMPONENT_%u_Z_CORRELATION_LENGTH = %.17g\n",
                    k + 1U,
                    model->exponential_model.component_weight[k],
                    k + 1U,
                    model->exponential_model.component_x_correlation_length[k],
                    k + 1U,
                    model->exponential_model.component_y_correlation_length[k],
                    k + 1U,
                    model->exponential_model.component_z_correlation_length[k]) < 0) {
            return EXIT_FAILURE;
        }
    }

    if (fprintf(out,
                "CORRELATION_LENGTH = %.17g\n"
                "X_CORRELATION_LENGTH = %.17g\n"
                "Y_CORRELATION_LENGTH = %.17g\n"
                "Z_CORRELATION_LENGTH = %.17g\n"
                "IN_PLANE_CORRELATION_LENGTH = %.17g\n"
                "THROUGH_THICKNESS_CORRELATION_LENGTH = %.17g\n"
                "RMS_ERROR = %.17g\n"
                "WEIGHTED_RMS_ERROR = %.17g\n"
                "DIRECTIONAL_RMS_X = %.17g\n"
                "DIRECTIONAL_RMS_Y = %.17g\n"
                "DIRECTIONAL_RMS_Z = %.17g\n"
                "SMALL_LAG_RMS_ERROR = %.17g\n"
                "RECOMMENDED_RADIUS_RMS_ERROR = %.17g\n"
                "AIC = %.17g\n"
                "BIC = %.17g\n"
                "SAMPLE_COUNT = %llu\n"
                "FIT_VALID = %d\n"
                "ZERO_LAG_RELATIVE_ERROR = %.17g\n"
                "RMS_ERROR_RATIO = %.17g\n",
                model->exponential_model.correlation_length,
                model->exponential_model.x_correlation_length,
                model->exponential_model.y_correlation_length,
                model->exponential_model.z_correlation_length,
                model->exponential_model.in_plane_correlation_length,
                model->exponential_model.through_thickness_correlation_length,
                model->exponential_model.rms_error,
                model->exponential_model.weighted_rms_error,
                model->exponential_model.directional_rms_x,
                model->exponential_model.directional_rms_y,
                model->exponential_model.directional_rms_z,
                model->exponential_model.small_lag_rms_error,
                model->exponential_model.recommended_radius_rms_error,
                model->exponential_model.aic,
                model->exponential_model.bic,
                (unsigned long long)model->exponential_model.sample_count,
                model->exponential_validation.is_valid,
                model->exponential_validation.zero_lag_relative_error,
                model->exponential_validation.rms_error_ratio) < 0) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int read_fitted_covariance_model(const char *model_path,
                                 FittedCovarianceModelArtifact *model_out) {
    FILE *file;
    char line[512];
    int seen_mean = 0;
    int seen_variance = 0;
    int seen_filter_resolution = 0;
    int seen_sill = 0;
    int seen_correlation_length = 0;
    int seen_x_length = 0;
    int seen_y_length = 0;
    int seen_z_length = 0;
    int seen_in_plane_length = 0;
    int seen_through_thickness_length = 0;
    int seen_rms_error = 0;
    int seen_sample_count = 0;
    int seen_component_count = 0;
    int seen_component_weight[POROSITY_NESTED_EXPONENTIAL_COMPONENTS] = {0};
    int seen_component_x[POROSITY_NESTED_EXPONENTIAL_COMPONENTS] = {0};
    int seen_component_y[POROSITY_NESTED_EXPONENTIAL_COMPONENTS] = {0};
    int seen_component_z[POROSITY_NESTED_EXPONENTIAL_COMPONENTS] = {0};

    if (model_path == NULL || model_out == NULL) {
        return EXIT_FAILURE;
    }

    memset(model_out, 0, sizeof(*model_out));
    model_out->exponential_model.model_type =
        POROSITY_COVARIANCE_MODEL_FULLY_ANISOTROPIC_EXPONENTIAL;

    file = fopen(model_path, "r");
    if (file == NULL) {
        return EXIT_FAILURE;
    }

    while (fgets(line, sizeof(line), file) != NULL) {
        double value;
        unsigned long long ull_value;
        int int_value;
        unsigned int component_index;

        if (strstr(line, "MODEL_TYPE = NESTED_FULLY_ANISOTROPIC_EXPONENTIAL") != NULL) {
            model_out->exponential_model.model_type =
                POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL;
        } else if (strstr(line, "MODEL_TYPE = FULLY_ANISOTROPIC_EXPONENTIAL") != NULL) {
            model_out->exponential_model.model_type =
                POROSITY_COVARIANCE_MODEL_FULLY_ANISOTROPIC_EXPONENTIAL;
        } else if (sscanf(line, "NORMAL_SCORE_TRANSFORM = %d", &int_value) == 1){
            model_out->normal_score_transform = int_value != 0;
        } else if (sscanf(line, "MEAN_POROSITY = %lf", &value) == 1) {
            model_out->mean_porosity = value;
            seen_mean = 1;
        } else if (sscanf(line, "VARIANCE = %lf", &value) == 1) {
            model_out->variance = value;
            seen_variance = 1;
        } else if (sscanf(line, "FILTER_RESOLUTION = %lf", &value) == 1) {
            model_out->filter_resolution = value;
            seen_filter_resolution = 1;
        } else if (sscanf(line, "SILL = %lf", &value) == 1) {
            model_out->exponential_model.sill = value;
            seen_sill = 1;
        } else if (sscanf(line, "COMPONENT_COUNT = %d", &int_value) == 1) {
            if (int_value <= 0 ||
                int_value > (int)POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
                fclose(file);
                return EXIT_FAILURE;
            }
            model_out->exponential_model.component_count = (uint32_t)int_value;
            seen_component_count = 1;
        } else if (sscanf(line, "COMPONENT_%u_WEIGHT = %lf",
                          &component_index, &value) == 2) {
            if (component_index == 0U ||
                component_index > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
                fclose(file);
                return EXIT_FAILURE;
            }
            component_index -= 1U;
            model_out->exponential_model.component_weight[component_index] = value;
            seen_component_weight[component_index] = 1;
        } else if (sscanf(line, "COMPONENT_%u_X_CORRELATION_LENGTH = %lf",
                          &component_index, &value) == 2) {
            if (component_index == 0U ||
                component_index > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
                fclose(file);
                return EXIT_FAILURE;
            }
            component_index -= 1U;
            model_out->exponential_model
                .component_x_correlation_length[component_index] = value;
            seen_component_x[component_index] = 1;
        } else if (sscanf(line, "COMPONENT_%u_Y_CORRELATION_LENGTH = %lf",
                          &component_index, &value) == 2) {
            if (component_index == 0U ||
                component_index > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
                fclose(file);
                return EXIT_FAILURE;
            }
            component_index -= 1U;
            model_out->exponential_model
                .component_y_correlation_length[component_index] = value;
            seen_component_y[component_index] = 1;
        } else if (sscanf(line, "COMPONENT_%u_Z_CORRELATION_LENGTH = %lf",
                          &component_index, &value) == 2) {
            if (component_index == 0U ||
                component_index > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
                fclose(file);
                return EXIT_FAILURE;
            }
            component_index -= 1U;
            model_out->exponential_model
                .component_z_correlation_length[component_index] = value;
            seen_component_z[component_index] = 1;
        } else if (sscanf(line, "CORRELATION_LENGTH = %lf", &value) == 1) {
            model_out->exponential_model.correlation_length = value;
            seen_correlation_length = 1;
        } else if (sscanf(line, "X_CORRELATION_LENGTH = %lf", &value) == 1) {
            model_out->exponential_model.x_correlation_length = value;
            seen_x_length = 1;
        } else if (sscanf(line, "Y_CORRELATION_LENGTH = %lf", &value) == 1) {
            model_out->exponential_model.y_correlation_length = value;
            seen_y_length = 1;
        } else if (sscanf(line, "Z_CORRELATION_LENGTH = %lf", &value) == 1) {
            model_out->exponential_model.z_correlation_length = value;
            seen_z_length = 1;
        } else if (sscanf(line, "IN_PLANE_CORRELATION_LENGTH = %lf",
                          &value) == 1) {
            model_out->exponential_model.in_plane_correlation_length = value;
            seen_in_plane_length = 1;
        } else if (sscanf(line,
                          "THROUGH_THICKNESS_CORRELATION_LENGTH = %lf",
                          &value) == 1) {
            model_out->exponential_model
                .through_thickness_correlation_length = value;
            seen_through_thickness_length = 1;
        } else if (sscanf(line, "RMS_ERROR = %lf", &value) == 1) {
            model_out->exponential_model.rms_error = value;
            seen_rms_error = 1;
        } else if (sscanf(line, "WEIGHTED_RMS_ERROR = %lf", &value) == 1) {
            model_out->exponential_model.weighted_rms_error = value;
        } else if (sscanf(line, "DIRECTIONAL_RMS_X = %lf", &value) == 1) {
            model_out->exponential_model.directional_rms_x = value;
        } else if (sscanf(line, "DIRECTIONAL_RMS_Y = %lf", &value) == 1) {
            model_out->exponential_model.directional_rms_y = value;
        } else if (sscanf(line, "DIRECTIONAL_RMS_Z = %lf", &value) == 1) {
            model_out->exponential_model.directional_rms_z = value;
        } else if (sscanf(line, "SMALL_LAG_RMS_ERROR = %lf", &value) == 1) {
            model_out->exponential_model.small_lag_rms_error = value;
        } else if (sscanf(line, "RECOMMENDED_RADIUS_RMS_ERROR = %lf",
                          &value) == 1) {
            model_out->exponential_model.recommended_radius_rms_error = value;
        } else if (sscanf(line, "AIC = %lf", &value) == 1) {
            model_out->exponential_model.aic = value;
        } else if (sscanf(line, "BIC = %lf", &value) == 1) {
            model_out->exponential_model.bic = value;
        } else if (sscanf(line, "SAMPLE_COUNT = %llu", &ull_value) == 1) {
            model_out->exponential_model.sample_count = (uint64_t)ull_value;
            seen_sample_count = 1;
        } else if (sscanf(line, "FIT_VALID = %d", &int_value) == 1) {
            model_out->exponential_validation.is_valid = int_value;
        } else if (sscanf(line, "ZERO_LAG_RELATIVE_ERROR = %lf",
                          &value) == 1) {
            model_out->exponential_validation.zero_lag_relative_error = value;
        } else if (sscanf(line, "RMS_ERROR_RATIO = %lf", &value) == 1) {
            model_out->exponential_validation.rms_error_ratio = value;
        }
    }

    fclose(file);

    if (!seen_x_length && seen_in_plane_length) {
        model_out->exponential_model.x_correlation_length =
            model_out->exponential_model.in_plane_correlation_length;
        seen_x_length = 1;
    }

    if (!seen_y_length && seen_in_plane_length) {
        model_out->exponential_model.y_correlation_length =
            model_out->exponential_model.in_plane_correlation_length;
        seen_y_length = 1;
    }

    if (!seen_z_length && seen_through_thickness_length) {
        model_out->exponential_model.z_correlation_length =
            model_out->exponential_model.through_thickness_correlation_length;
        seen_z_length = 1;
    }

    if (!seen_x_length && seen_correlation_length) {
        model_out->exponential_model.x_correlation_length =
            model_out->exponential_model.correlation_length;
        seen_x_length = 1;
    }

    if (!seen_y_length && seen_correlation_length) {
        model_out->exponential_model.y_correlation_length =
            model_out->exponential_model.correlation_length;
        seen_y_length = 1;
    }

    if (!seen_z_length && seen_correlation_length) {
        model_out->exponential_model.z_correlation_length =
            model_out->exponential_model.correlation_length;
        seen_z_length = 1;
    }

    if (!seen_correlation_length && seen_x_length) {
        model_out->exponential_model.correlation_length =
            model_out->exponential_model.x_correlation_length;
        seen_correlation_length = 1;
    }

    if (!seen_in_plane_length && seen_x_length) {
        model_out->exponential_model.in_plane_correlation_length =
            model_out->exponential_model.x_correlation_length;
        seen_in_plane_length = 1;
    }

    if (!seen_through_thickness_length && seen_z_length) {
        model_out->exponential_model.through_thickness_correlation_length =
            model_out->exponential_model.z_correlation_length;
        seen_through_thickness_length = 1;
    }

    if (!seen_component_count) {
        model_out->exponential_model.component_count = 1U;
        model_out->exponential_model.component_weight[0] = 1.0;
        model_out->exponential_model.component_x_correlation_length[0] =
            model_out->exponential_model.x_correlation_length;
        model_out->exponential_model.component_y_correlation_length[0] =
            model_out->exponential_model.y_correlation_length;
        model_out->exponential_model.component_z_correlation_length[0] =
            model_out->exponential_model.z_correlation_length;
        seen_component_weight[0] = 1;
        seen_component_x[0] = seen_x_length;
        seen_component_y[0] = seen_y_length;
        seen_component_z[0] = seen_z_length;
    }

    for (uint32_t k = 0U; k < model_out->exponential_model.component_count; ++k) {
        if (!seen_component_weight[k] ||
            !seen_component_x[k] ||
            !seen_component_y[k] ||
            !seen_component_z[k]) {
            return EXIT_FAILURE;
        }
    }

    if (!seen_mean || !seen_variance || !seen_filter_resolution ||
        !seen_sill || !seen_x_length || !seen_y_length ||
        !seen_z_length || !seen_rms_error || !seen_sample_count) {
        return EXIT_FAILURE;
    }

    if (model_out->variance <= 0.0 ||
        model_out->filter_resolution <= 0.0 ||
        model_out->exponential_model.sill <= 0.0 ||
        model_out->exponential_model.x_correlation_length <= 0.0 ||
        model_out->exponential_model.y_correlation_length <= 0.0 ||
        model_out->exponential_model.z_correlation_length <= 0.0) {
        return EXIT_FAILURE;
    }

    for (uint32_t k = 0U; k < model_out->exponential_model.component_count; ++k) {
        if (model_out->exponential_model.component_weight[k] < 0.0 ||
            model_out->exponential_model.component_x_correlation_length[k] <= 0.0 ||
            model_out->exponential_model.component_y_correlation_length[k] <= 0.0 ||
            model_out->exponential_model.component_z_correlation_length[k] <= 0.0) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int read_3dlag_covariance_table(const char *table_path,
                                double filter_resolution,
                                PorosityCovarianceSample **samples_out,
                                uint32_t *width_out, uint32_t *height_out,
                                uint32_t *depth_out) {
    FILE *file;
    char line[8192];
    uint32_t max_x = 0U;
    uint32_t max_y = 0U;
    uint32_t max_z = 0U;
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    PorosityCovarianceSample *samples;

    if (table_path == NULL || samples_out == NULL || width_out == NULL ||
        height_out == NULL || depth_out == NULL || filter_resolution <= 0.0) {
        return EXIT_FAILURE;
    }

    *samples_out = NULL;
    *width_out = 0U;
    *height_out = 0U;
    *depth_out = 0U;

    file = fopen(table_path, "r");
    if (file == NULL) {
        return EXIT_FAILURE;
    }

    if (fgets(line, sizeof(line), file) == NULL) {
        fclose(file);
        return EXIT_FAILURE;
    }

    while (fgets(line, sizeof(line), file) != NULL) {
        double lag_x_phys;
        double lag_y_phys;
        double lag_z_phys;
        unsigned long long sample_count;
        double mu;
        double covariance;
        int parsed;
        long lag_x;
        long lag_y;
        long lag_z;

        parsed = sscanf(line, " %lf %lf %lf %llu %lf %lf", &lag_x_phys,
                        &lag_y_phys, &lag_z_phys, &sample_count, &mu,
                        &covariance);
        if (parsed != 6) {
            fclose(file);
            return EXIT_FAILURE;
        }

        (void)sample_count;
        (void)mu;
        (void)covariance;

        lag_x = lround(lag_x_phys / filter_resolution);
        lag_y = lround(lag_y_phys / filter_resolution);
        lag_z = lround(lag_z_phys / filter_resolution);
        if (lag_x < 0L || lag_y < 0L || lag_z < 0L) {
            fclose(file);
            return EXIT_FAILURE;
        }

        if ((uint32_t)lag_x > max_x) {
            max_x = (uint32_t)lag_x;
        }
        if ((uint32_t)lag_y > max_y) {
            max_y = (uint32_t)lag_y;
        }
        if ((uint32_t)lag_z > max_z) {
            max_z = (uint32_t)lag_z;
        }
    }

    width = max_x + 1U;
    height = max_y + 1U;
    depth = max_z + 1U;
    samples = (PorosityCovarianceSample *)malloc((size_t)width * (size_t)height *
                                                 (size_t)depth *
                                                 sizeof(*samples));
    if (samples == NULL) {
        fclose(file);
        return EXIT_FAILURE;
    }

    for (uint32_t z = 0; z < depth; ++z) {
        for (uint32_t y = 0; y < height; ++y) {
            for (uint32_t x = 0; x < width; ++x) {
                size_t index = covariance_grid_index(x, y, z, width, height);
                samples[index].lag_x = -1;
                samples[index].lag_y = -1;
                samples[index].lag_z = -1;
                samples[index].samples = 0U;
                samples[index].mu = 0.0;
                samples[index].covariance = 0.0;
            }
        }
    }

    rewind(file);
    if (fgets(line, sizeof(line), file) == NULL) {
        fclose(file);
        free(samples);
        return EXIT_FAILURE;
    }

    while (fgets(line, sizeof(line), file) != NULL) {
        double lag_x_phys;
        double lag_y_phys;
        double lag_z_phys;
        unsigned long long sample_count;
        double mu;
        double covariance;
        int parsed;
        uint32_t lag_x;
        uint32_t lag_y;
        uint32_t lag_z;
        size_t index;

        parsed = sscanf(line, " %lf %lf %lf %llu %lf %lf", &lag_x_phys,
                        &lag_y_phys, &lag_z_phys, &sample_count, &mu,
                        &covariance);
        if (parsed != 6) {
            fclose(file);
            free(samples);
            return EXIT_FAILURE;
        }

        lag_x = (uint32_t)lround(lag_x_phys / filter_resolution);
        lag_y = (uint32_t)lround(lag_y_phys / filter_resolution);
        lag_z = (uint32_t)lround(lag_z_phys / filter_resolution);
        index = covariance_grid_index(lag_x, lag_y, lag_z, width, height);
        if (samples[index].lag_x >= 0) {
            fclose(file);
            free(samples);
            return EXIT_FAILURE;
        }

        samples[index].lag_x = (int)lag_x;
        samples[index].lag_y = (int)lag_y;
        samples[index].lag_z = (int)lag_z;
        samples[index].samples = (uint64_t)sample_count;
        samples[index].mu = mu;
        samples[index].covariance = covariance;
    }

    fclose(file);

    for (uint32_t z = 0; z < depth; ++z) {
        for (uint32_t y = 0; y < height; ++y) {
            for (uint32_t x = 0; x < width; ++x) {
                size_t index = covariance_grid_index(x, y, z, width, height);
                if (samples[index].lag_x < 0 || samples[index].lag_y < 0 ||
                    samples[index].lag_z < 0) {
                    free(samples);
                    return EXIT_FAILURE;
                }
            }
        }
    }

    *samples_out = samples;
    *width_out = width;
    *height_out = height;
    *depth_out = depth;
    return EXIT_SUCCESS;
}

int prepare_realization_from_3dlag_table(
    const char *table_path, double filter_resolution, double variance,
    const PorosityRealizationGridSpec *simulation_grid,
    PorosityRealizationPrepResult *result_out) {
    if (table_path == NULL || filter_resolution <= 0.0 || variance < 0.0 ||
        simulation_grid == NULL || result_out == NULL) {
        return EXIT_FAILURE;
    }

    memset(result_out, 0, sizeof(*result_out));
    result_out->filter_resolution = filter_resolution;
    result_out->variance = variance;
    result_out->simulation_grid = *simulation_grid;

    if (read_3dlag_covariance_table(table_path, filter_resolution,
                                    &result_out->samples,
                                    &result_out->lag_width,
                                    &result_out->lag_height,
                                    &result_out->lag_depth) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (reconstruct_covariance_grid_from_3dlag_samples(
            result_out->samples, result_out->lag_width, result_out->lag_height,
            result_out->lag_depth, variance,
            &result_out->covariance_grid) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (symmetrize_covariance_grid_by_reflection(
            &result_out->covariance_grid,
            &result_out->symmetric_grid) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (analyze_covariance_tail_noise(result_out->samples,
                                      &result_out->covariance_grid,
                                      filter_resolution,
                                      &result_out->tail_noise) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (analyze_negative_covariance_lobes(
            &result_out->covariance_grid, filter_resolution,
            &result_out->negative_covariance) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (analyze_covariance_anisotropy(&result_out->covariance_grid,
                                      filter_resolution,
                                      &result_out->anisotropy) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (fit_nested_full_anisotropic_exponential_covariance_model(
            &result_out->covariance_grid,
            result_out->samples,
            filter_resolution,
            result_out->tail_noise.recommended_max_radius,
            &result_out->exponential_model) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (validate_exponential_covariance_model(
            &result_out->covariance_grid, &result_out->exponential_model,
            &result_out->exponential_validation) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (assess_realization_grid_compatibility(
            &result_out->simulation_grid, &result_out->exponential_model,
            &result_out->simulation_grid_compatibility) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (evaluate_exponential_covariance_on_realization_grid(
            &result_out->exponential_model, &result_out->simulation_grid,
            &result_out->evaluated_simulation_covariance_grid) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    if (validate_evaluated_covariance_grid(
            &result_out->evaluated_simulation_covariance_grid,
            &result_out->exponential_model,
            &result_out->evaluated_grid_validation) != EXIT_SUCCESS) {
        free_realization_prep_result(result_out);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static int write_exponential_fit_quality_report(
    FILE *out,
    const PorosityRealizationPrepResult *result) {
    if (out == NULL || result == NULL) {
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "Exponential Fit Quality\n"
                "  Model type: %s\n"
                "  Sill: %.12f\n"
                "  Legacy X correlation length: %.12f\n"
                "  Legacy Y correlation length: %.12f\n"
                "  Legacy Z correlation length: %.12f\n",
                covariance_model_type_label(&result->exponential_model),
                result->exponential_model.sill,
                result->exponential_model.x_correlation_length,
                result->exponential_model.y_correlation_length,
                result->exponential_model.z_correlation_length) < 0) {
        return EXIT_FAILURE;
    }

    if (write_covariance_model_components(
            out, &result->exponential_model, "  ") != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "  RMS error: %.12f\n"
                "  Weighted RMS error: %.12f\n"
                "  Directional RMS x/y/z: %.12f %.12f %.12f\n"
                "  Small-lag RMS error: %.12f\n"
                "  Recommended-radius RMS error: %.12f\n"
                "  AIC: %.12f\n"
                "  BIC: %.12f\n"
                "  Fit sample count: %llu\n"
                "  Fit valid: %s\n"
                "  Zero-lag relative error: %.12f\n"
                "  RMS error ratio: %.12f\n",
                result->exponential_model.rms_error,
                result->exponential_model.weighted_rms_error,
                result->exponential_model.directional_rms_x,
                result->exponential_model.directional_rms_y,
                result->exponential_model.directional_rms_z,
                result->exponential_model.small_lag_rms_error,
                result->exponential_model.recommended_radius_rms_error,
                result->exponential_model.aic,
                result->exponential_model.bic,
                (unsigned long long)result->exponential_model.sample_count,
                result->exponential_validation.is_valid ? "YES" : "NO",
                result->exponential_validation.zero_lag_relative_error,
                result->exponential_validation.rms_error_ratio) < 0) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int process_realization_prep_from_stack(
    const char *stack_path, double voxel_resolution, double filter_resolution,
    const PorosityRealizationGridSpec *simulation_grid,
    uint32_t realization_count, FILE *out) {
    PorosityRealizationPrepResult result;
    char *table_path = NULL;
    char *tail_report_path = NULL;
    char *fit_report_path = NULL;
    char *model_path = NULL;
    char *summary_report_path = NULL;
    FILE *tail_report = NULL;
    FILE *fit_report = NULL;
    FILE *model_file = NULL;
    FILE *summary_report = NULL;
    double variance;

    if (stack_path == NULL || simulation_grid == NULL || out == NULL ||
        realization_count == 0U) {
        return EXIT_FAILURE;
    }

    variance = compute_porosity_field_variance_from_stack(
        stack_path, voxel_resolution, filter_resolution, 0, NULL);
    if (variance < 0.0) {
        fprintf(stderr, "Failed to compute porosity-field variance for %s\n",
                stack_path);
        return EXIT_FAILURE;
    }

    if (ensure_build_covariance_dir(stack_path) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to create BUILD_COVARIANCE for %s\n",
                stack_path);
        return EXIT_FAILURE;
    }

    table_path = build_3dlag_table_path(stack_path);
    tail_report_path = build_build_covariance_output_path(
        stack_path, "_tail_noise_analysis.txt");
    fit_report_path = build_build_covariance_output_path(
        stack_path, "_exponential_fit_quality.txt");
    model_path = build_fitted_model_path(stack_path);
    summary_report_path = build_build_covariance_output_path(
        stack_path, "_covariance_build_summary.txt");
    if (table_path == NULL || tail_report_path == NULL || fit_report_path == NULL ||
        model_path == NULL || summary_report_path == NULL) {
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (prepare_realization_from_3dlag_table(table_path, filter_resolution,
                                             variance, simulation_grid,
                                             &result) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to prepare realization covariance for %s\n",
                stack_path);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    tail_report = fopen(tail_report_path, "w");
    fit_report = fopen(fit_report_path, "w");
    model_file = fopen(model_path, "w");
    summary_report = fopen(summary_report_path, "w");
    if (tail_report == NULL || fit_report == NULL || model_file == NULL ||
        summary_report == NULL) {
        if (tail_report != NULL) fclose(tail_report);
        if (fit_report != NULL) fclose(fit_report);
        if (model_file != NULL) fclose(model_file);
        if (summary_report != NULL) fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (fprintf(tail_report,
                "Tail Noise Analysis\n"
                "  Tail noisy: %s\n"
                "  Total tail bins: %u\n"
                "  Flagged tail bins: %u\n"
                "  Tail noise fraction: %.12f\n"
                "  Max relative fluctuation: %.12f\n"
                "  Min mean samples: %.12f\n"
                "  Recommended max radius: %.12f\n",
                result.tail_noise.tail_is_noisy ? "YES" : "NO",
                result.tail_noise.total_tail_bins,
                result.tail_noise.flagged_tail_bins,
                result.tail_noise.tail_noise_fraction,
                result.tail_noise.max_relative_fluctuation,
                result.tail_noise.min_mean_samples,
                result.tail_noise.recommended_max_radius) < 0) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (write_exponential_fit_quality_report(fit_report, &result) !=
        EXIT_SUCCESS) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    {
        FittedCovarianceModelArtifact model_artifact;

        memset(&model_artifact, 0, sizeof(model_artifact));
        model_artifact.mean_porosity = result.covariance_grid.mean;
        model_artifact.variance = result.variance;
        model_artifact.filter_resolution = result.filter_resolution;
        model_artifact.exponential_model = result.exponential_model;
        model_artifact.exponential_validation = result.exponential_validation;

        if (write_fitted_covariance_model(model_file, &model_artifact) !=
            EXIT_SUCCESS) {
            fclose(tail_report);
            fclose(fit_report);
            fclose(model_file);
            fclose(summary_report);
            free_realization_prep_result(&result);
            free(table_path);
            free(tail_report_path);
            free(fit_report_path);
            free(model_path);
            free(summary_report_path);
            return EXIT_FAILURE;
        }
    }

    if (write_realization_prep_report(summary_report, &result) != EXIT_SUCCESS) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    fclose(tail_report);
    fclose(fit_report);
    fclose(model_file);
    fclose(summary_report);

    if (generate_and_write_gaussian_porosity_realizations(
            stack_path, &result, realization_count, 12345U) != EXIT_SUCCESS) {
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "Covariance build: %s\n"
                "  3DLAG table: %s\n"
                "  Tail noise report: %s\n"
                "  Exponential fit report: %s\n"
                "  Fitted model: %s\n"
                "  Summary: %s\n"
                "  Realizations written: %u\n\n",
                stack_path, table_path, tail_report_path, fit_report_path,
                model_path, summary_report_path, realization_count) < 0) {
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    free_realization_prep_result(&result);
    free(table_path);
    free(tail_report_path);
    free(fit_report_path);
    free(model_path);
    free(summary_report_path);
    return EXIT_SUCCESS;
}

int process_realization_fit_from_stack(
    const char *stack_path, double voxel_resolution, double filter_resolution,
    const PorosityRealizationGridSpec *simulation_grid, 
    int normal_score_enabled, FILE *out) {
    double field_mean = 0.0;
    PorosityRealizationPrepResult result;
    char *table_path = NULL;
    char *tail_report_path = NULL;
    char *fit_report_path = NULL;
    char *model_path = NULL;
    char *summary_report_path = NULL;
    char *directional_fit_path = NULL;
    FILE *directional_fit_report = NULL;
    FILE *tail_report = NULL;
    FILE *fit_report = NULL;
    FILE *model_file = NULL;
    FILE *summary_report = NULL;
    double variance;

    if (stack_path == NULL || simulation_grid == NULL || out == NULL) {
        return EXIT_FAILURE;
    }

    variance = compute_porosity_field_variance_from_stack(
        stack_path, voxel_resolution, filter_resolution, normal_score_enabled,
        &field_mean);
    if (variance < 0.0) {
        fprintf(stderr, "Failed to compute porosity-field variance for %s\n",
                stack_path);
        return EXIT_FAILURE;
    }

    if (ensure_build_covariance_dir(stack_path) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to create BUILD_COVARIANCE for %s\n",
                stack_path);
        return EXIT_FAILURE;
    }

    table_path = build_3dlag_table_path(stack_path);
    tail_report_path = build_build_covariance_output_path(
        stack_path, "_tail_noise_analysis.txt");
    fit_report_path = build_build_covariance_output_path(
        stack_path, "_exponential_fit_quality.txt");
    model_path = build_fitted_model_path(stack_path);
    summary_report_path = build_build_covariance_output_path(
        stack_path, "_covariance_build_summary.txt");
    directional_fit_path = build_build_covariance_output_path(
        stack_path, "_directional_covariance_fit.txt");
    if (table_path == NULL || tail_report_path == NULL || fit_report_path == NULL ||
        model_path == NULL || summary_report_path == NULL ||
        directional_fit_path == NULL) {
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        free(directional_fit_path);
        return EXIT_FAILURE;
    }

    if (prepare_realization_from_3dlag_table(table_path, filter_resolution,
                                             variance, simulation_grid,
                                             &result) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to prepare realization covariance for %s\n",
                stack_path);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    //3D LAG table stores mean of input data set, it it doesn't match rebuild,
    //it came from different settings i.e theres an error
    if (fabs(result.covariance_grid.mean - field_mean) >
        1.0e-6 * (1.0 + fabs(field_mean))) {
        fprintf(stderr,
                "3DLAG table mean (%.9f) does not match the current field "
                "mean (%.9f) for %s.\n"
                "Rerun EMPIRICAL_COVARIANCE with the same "
                "NORMAL_SCORE_TRANSFORM and FILTER_RESOLUTION settings.\n",
                result.covariance_grid.mean, field_mean, stack_path);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        free(directional_fit_path);
        return EXIT_FAILURE;
    }

    tail_report = fopen(tail_report_path, "w");
    fit_report = fopen(fit_report_path, "w");
    model_file = fopen(model_path, "w");
    summary_report = fopen(summary_report_path, "w");
    directional_fit_report = fopen(directional_fit_path, "w");
    if (tail_report == NULL || fit_report == NULL || model_file == NULL ||
        summary_report == NULL || directional_fit_report == NULL) {
        if (tail_report != NULL) fclose(tail_report);
        if (fit_report != NULL) fclose(fit_report);
        if (model_file != NULL) fclose(model_file);
        if (summary_report != NULL) fclose(summary_report);
        if (directional_fit_report != NULL) fclose(directional_fit_report);
        free(directional_fit_path);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (fprintf(tail_report,
                "Tail Noise Analysis\n"
                "  Tail noisy: %s\n"
                "  Total tail bins: %u\n"
                "  Flagged tail bins: %u\n"
                "  Tail noise fraction: %.12f\n"
                "  Max relative fluctuation: %.12f\n"
                "  Min mean samples: %.12f\n"
                "  Recommended max radius: %.12f\n",
                result.tail_noise.tail_is_noisy ? "YES" : "NO",
                result.tail_noise.total_tail_bins,
                result.tail_noise.flagged_tail_bins,
                result.tail_noise.tail_noise_fraction,
                result.tail_noise.max_relative_fluctuation,
                result.tail_noise.min_mean_samples,
                result.tail_noise.recommended_max_radius) < 0) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (write_exponential_fit_quality_report(fit_report, &result) !=
        EXIT_SUCCESS) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    {
        FittedCovarianceModelArtifact model_artifact;

        memset(&model_artifact, 0, sizeof(model_artifact));
        model_artifact.mean_porosity = result.covariance_grid.mean;
        model_artifact.variance = result.variance;
        model_artifact.normal_score_transform = normal_score_enabled ? 1 : 0;
        model_artifact.filter_resolution = result.filter_resolution;
        model_artifact.exponential_model = result.exponential_model;
        model_artifact.exponential_validation = result.exponential_validation;

        if (write_fitted_covariance_model(model_file, &model_artifact) !=
            EXIT_SUCCESS) {
            fclose(tail_report);
            fclose(fit_report);
            fclose(model_file);
            fclose(summary_report);
            free_realization_prep_result(&result);
            free(table_path);
            free(tail_report_path);
            free(fit_report_path);
            free(model_path);
            free(summary_report_path);
            return EXIT_FAILURE;
        }
    }

    if (write_realization_prep_report(summary_report, &result) != EXIT_SUCCESS) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        return EXIT_FAILURE;
    }

    if (write_directional_covariance_fit_table(
            directional_fit_report,
            &result.covariance_grid,
            result.samples,
            result.filter_resolution,
            &result.exponential_model) != EXIT_SUCCESS) {
        fclose(tail_report);
        fclose(fit_report);
        fclose(model_file);
        fclose(summary_report);
        fclose(directional_fit_report);
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        free(directional_fit_path);
        return EXIT_FAILURE;
    }    
    fclose(tail_report);
    fclose(fit_report);
    fclose(model_file);
    fclose(summary_report);
    fclose(directional_fit_report);

    if (fprintf(out,
                "Covariance fit: %s\n"
                "  3DLAG table: %s\n"
                "  Tail noise report: %s\n"
                "  Exponential fit report: %s\n"
                "  Fitted model: %s\n"
                "  Summary: %s\n"
                "  Directional fit table: %s \n\n",
                stack_path, table_path, tail_report_path, fit_report_path,
                model_path, summary_report_path, directional_fit_path) < 0) {
        free_realization_prep_result(&result);
        free(table_path);
        free(tail_report_path);
        free(fit_report_path);
        free(model_path);
        free(summary_report_path);
        free(directional_fit_path);
        return EXIT_FAILURE;
    }

    free_realization_prep_result(&result);
    free(table_path);
    free(tail_report_path);
    free(fit_report_path);
    free(model_path);
    free(summary_report_path);
    free(directional_fit_path);
    return EXIT_SUCCESS;
}
