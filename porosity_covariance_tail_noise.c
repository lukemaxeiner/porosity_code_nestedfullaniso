#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "porosity_covariance_tail_noise.h"

#define COVARIANCE_TAIL_BIN_COUNT 12U
#define COVARIANCE_TAIL_START_FRACTION 0.75
#define COVARIANCE_TAIL_RELATIVE_FLUCTUATION_LIMIT 1.0
#define COVARIANCE_TAIL_MIN_MEAN_SAMPLES 32.0

typedef struct {
    double radius_min;
    double radius_max;
    double covariance_sum;
    double covariance_sq_sum;
    double sample_sum;
    uint32_t point_count;
    int flagged;
} CovarianceTailBinStats;

static size_t covariance_grid_index(uint32_t x, uint32_t y, uint32_t z,
                                    uint32_t width, uint32_t height) {
    return (size_t)z * (size_t)width * (size_t)height +
           (size_t)y * (size_t)width + (size_t)x;
}

static double covariance_lag_radius(uint32_t x, uint32_t y, uint32_t z,
                                    double spacing) {
    double dx = (double)x * spacing;
    double dy = (double)y * spacing;
    double dz = (double)z * spacing;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

int check_covariance_tail_noise(const PorosityCovarianceSample *samples,
                                const CovarianceGrid3D *grid,
                                double filter_resolution,
                                CovarianceTailNoiseSummary *summary_out) {
    CovarianceTailBinStats bins[COVARIANCE_TAIL_BIN_COUNT];
    uint32_t x;
    uint32_t y;
    uint32_t z;
    double max_radius;
    double tail_start;
    double bin_width;
    double recommended_max_radius;

    if (samples == NULL || grid == NULL || grid->values == NULL ||
        summary_out == NULL || filter_resolution <= 0.0) {
        return EXIT_FAILURE;
    }

    memset(summary_out, 0, sizeof(*summary_out));
    memset(bins, 0, sizeof(bins));

    max_radius = covariance_lag_radius(grid->width > 0U ? grid->width - 1U : 0U,
                                       grid->height > 0U ? grid->height - 1U : 0U,
                                       grid->depth > 0U ? grid->depth - 1U : 0U,
                                       filter_resolution);
    tail_start = COVARIANCE_TAIL_START_FRACTION * max_radius;
    recommended_max_radius = max_radius;
    bin_width = max_radius > tail_start
                    ? (max_radius - tail_start) / (double)COVARIANCE_TAIL_BIN_COUNT
                    : 1.0;

    for (x = 0; x < COVARIANCE_TAIL_BIN_COUNT; ++x) {
        bins[x].radius_min = tail_start + (double)x * bin_width;
        bins[x].radius_max = tail_start + (double)(x + 1U) * bin_width;
    }

    for (z = 0; z < grid->depth; ++z) {
        for (y = 0; y < grid->height; ++y) {
            for (x = 0; x < grid->width; ++x) {
                double radius =
                    covariance_lag_radius(x, y, z, filter_resolution);
                size_t index;
                uint32_t bin_index;

                if (radius < tail_start) {
                    continue;
                }

                if (bin_width <= 0.0) {
                    bin_index = 0U;
                } else {
                    double normalized = (radius - tail_start) / bin_width;
                    bin_index = (uint32_t)normalized;
                    if (bin_index >= COVARIANCE_TAIL_BIN_COUNT) {
                        bin_index = COVARIANCE_TAIL_BIN_COUNT - 1U;
                    }
                }

                index = covariance_grid_index(x, y, z, grid->width, grid->height);
                bins[bin_index].covariance_sum += grid->values[index];
                bins[bin_index].covariance_sq_sum +=
                    grid->values[index] * grid->values[index];
                bins[bin_index].sample_sum += (double)samples[index].samples;
                bins[bin_index].point_count += 1U;
            }
        }
    }

    for (x = 0; x < COVARIANCE_TAIL_BIN_COUNT; ++x) {
        if (bins[x].point_count > 0U) {
            double mean_covariance =
                bins[x].covariance_sum / (double)bins[x].point_count;
            double mean_covariance_sq =
                bins[x].covariance_sq_sum / (double)bins[x].point_count;
            double variance =
                mean_covariance_sq - mean_covariance * mean_covariance;
            double stddev = variance > 0.0 ? sqrt(variance) : 0.0;
            double relative_fluctuation =
                stddev / fmax(fabs(mean_covariance), 1.0e-12);
            double mean_samples = bins[x].sample_sum / (double)bins[x].point_count;

            summary_out->total_tail_bins += 1U;
            if (relative_fluctuation > summary_out->max_relative_fluctuation) {
                summary_out->max_relative_fluctuation = relative_fluctuation;
            }
            if (summary_out->min_mean_samples == 0.0 ||
                mean_samples < summary_out->min_mean_samples) {
                summary_out->min_mean_samples = mean_samples;
            }

            if (relative_fluctuation >
                    COVARIANCE_TAIL_RELATIVE_FLUCTUATION_LIMIT &&
                mean_samples < COVARIANCE_TAIL_MIN_MEAN_SAMPLES) {
                bins[x].flagged = 1;
                summary_out->flagged_tail_bins += 1U;
                if (bins[x].radius_min < recommended_max_radius) {
                    recommended_max_radius = bins[x].radius_min;
                }
            }
        }
    }

    if (summary_out->total_tail_bins > 0U) {
        summary_out->tail_noise_fraction =
            (double)summary_out->flagged_tail_bins /
            (double)summary_out->total_tail_bins;
    }
    summary_out->tail_is_noisy = summary_out->flagged_tail_bins > 0U;
    summary_out->recommended_max_radius = recommended_max_radius;

    return EXIT_SUCCESS;
}
