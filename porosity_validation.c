#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "porosity_covariance_lag.h"
#include "tiff.h"
#include "porosity_validation.h"

/* ======================================================================
 * Step 2a: one-point statistics
 * ====================================================================== */

static const double quantile_probabilities[VALIDATION_QUANTILE_COUNT] = {
    0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99
};

double validation_quantile_probability(int index) {
    if (index < 0 || index >= VALIDATION_QUANTILE_COUNT) {
        return NAN;
    }
    return quantile_probabilities[index];
}

/* qsort needs a function that says whether a < b (-1), a == b (0), a > b (1).
   (x > y) - (x < y) produces exactly that without subtraction overflow tricks. */
static int compare_doubles(const void *a, const void *b) {
    double x = *(const double *)a;
    double y = *(const double *)b;
    return (x > y) - (x < y);
}

/* Quantile of already-sorted data, "linear interpolation" definition.
 * This is the same rule numpy.quantile uses by default, so you can check
 * the C output against Python.
 *
 * Picture the n sorted values sitting at positions 0, 1, ..., n-1.
 * Probability p lands at position (n-1)*p. If that falls between two values,
 * blend them by how far along it is. */
static double sorted_quantile(const double *sorted, uint64_t n, double p) {
    double position = (double)(n - 1U) * p;
    uint64_t lower = (uint64_t)floor(position);
    double t = position - (double)lower;   /* 0 = exactly on lower, ->1 = near upper */

    if (lower + 1U >= n) {
        return sorted[n - 1U];
    }
    return (1.0 - t) * sorted[lower] + t * sorted[lower + 1U];
}

static int compute_one_point_statistics(const PorosityField *field,
                                        FieldValidationMetrics *out) {
    const uint64_t n = out->cell_count;
    double sum = 0.0;
    double sum_sq_centered = 0.0;
    double sum_cube_centered = 0.0;
    uint64_t empty_count = 0U;
    double *sorted;

    /* Pass 1: mean, and count empty cells while we're looping anyway. */
    for (uint64_t i = 0U; i < n; ++i) {
        double v = field->values[i];
        sum += v;
        if (v >= VALIDATION_EMPTY_THRESHOLD) {
            empty_count += 1U;
        }
    }
    out->mean = sum / (double)n;
    out->empty_fraction = (double)empty_count / (double)n;

    /* Pass 2: spread (variance) and lopsidedness (skewness) around the mean.
       variance = average of (v - mean)^2
       skewness = average of (v - mean)^3, divided by variance^(3/2) so it has no units */
    for (uint64_t i = 0U; i < n; ++i) {
        double d = field->values[i] - out->mean;
        sum_sq_centered += d * d;
        sum_cube_centered += d * d * d;
    }
    out->variance = sum_sq_centered / (double)n;
    out->skewness = out->variance > 0.0
                        ? (sum_cube_centered / (double)n) / pow(out->variance, 1.5)
                        : 0.0;

    /* Quantiles need the values in order. Sort a COPY -- the realization
       field is still going to be written to TIFF afterwards, and sorting the
       original would scramble its spatial layout. */
    sorted = (double *)malloc((size_t)n * sizeof(double));
    if (sorted == NULL) {
        fprintf(stderr, "Validation: unable to allocate %llu values for sorting\n",
                (unsigned long long)n);
        return EXIT_FAILURE;
    }
    memcpy(sorted, field->values, (size_t)n * sizeof(double));
    qsort(sorted, (size_t)n, sizeof(double), compare_doubles);

    out->min = sorted[0];
    out->max = sorted[n - 1U];
    for (int q = 0; q < VALIDATION_QUANTILE_COUNT; ++q) {
        out->quantiles[q] = sorted_quantile(sorted, n, quantile_probabilities[q]);
    }

    free(sorted);
    return EXIT_SUCCESS;
}

/* ======================================================================
 * Step 2b: two-point statistics (axis covariance)
 * ====================================================================== */

int validation_resolve_max_lag(const PorosityField *field, int requested_max_lag) {
    uint32_t shortest = field->width;
    int largest_possible;

    if (field->height < shortest) shortest = field->height;
    if (field->depth < shortest) shortest = field->depth;

    /* A lag of h cells needs at least one pair of cells h apart, so the
       biggest lag that exists on every axis is (shortest axis - 1). */
    largest_possible = (int)shortest - 1;
    if (largest_possible < 0) {
        largest_possible = 0;
    }

    if (requested_max_lag == 0) {
        /* Automatic: half the shortest axis. Past that, fewer than half the
           cells have a partner, and the estimate gets noisy. (Common
           geostatistics rule of thumb.) */
        return (int)(shortest / 2U);
    }
    if (requested_max_lag > largest_possible) {
        fprintf(stderr,
                "Warning: VALIDATION_MAX_LAG = %d is larger than this field allows; using %d\n",
                requested_max_lag, largest_possible);
        return largest_possible;
    }
    return requested_max_lag;
}

/* Fill cov[0..max_lag] for one axis.
 *
 * We reuse compute_porosity_covariance_sample() from the 3DLAG stage, so the
 * data and the realizations are measured with EXACTLY the estimator the
 * covariance fit was built on. That function returns a CORRELATION
 * (covariance divided by the field variance, "rho"), so we multiply the
 * variance back in to get covariance C(h) in porosity^2 units -- the same
 * units the fitted model's sill uses. */
static int compute_axis_covariance(const PorosityField *field, int step_x,
                                   int step_y, int step_z, double variance,
                                   int max_lag, double *cov) {
    for (int h = 0; h <= max_lag; ++h) {
        PorosityCovarianceSample sample;

        if (compute_porosity_covariance_sample(field, h * step_x, h * step_y,
                                               h * step_z, &sample) !=
            EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
        cov[h] = sample.covariance * variance;
    }
    return EXIT_SUCCESS;
}

/* ======================================================================
 * Step 2c: mean porosity of each z layer
 * ====================================================================== */

/* If the sample has a porosity gradient along z (more empty space at the top
 * than the bottom), the field is not stationary and that trend inflates the
 * long-range covariance the model is fit to. A realization from a stationary
 * model has no trend, so comparing these two profiles tells us whether
 * detrending the data is worth doing. */
static void compute_z_layer_means(const PorosityField *field, double *z_mean) {
    const uint64_t cells_per_layer = (uint64_t)field->width * field->height;

    for (uint32_t z = 0; z < field->depth; ++z) {
        /* Layers are contiguous in memory: index = z*width*height + y*width + x,
           so one whole layer starts at z * cells_per_layer. */
        const double *layer = field->values + (size_t)(z * cells_per_layer);
        double sum = 0.0;

        for (uint64_t i = 0U; i < cells_per_layer; ++i) {
            sum += layer[i];
        }
        z_mean[z] = sum / (double)cells_per_layer;
    }
}

/* ======================================================================
 * Step 2d: clusters of completely empty cells
 * ====================================================================== */

/* The data has isolated pores: blobs of touching cells that are all porosity
 * 1.0. A Gaussian random field reproduces the FRACTION of such cells (the
 * normal-score transform guarantees that), but not necessarily their SHAPE --
 * Gaussian upper-tail excursions tend to be small and scattered. Counting
 * clusters and measuring how big they get is how we check that.
 *
 * Algorithm: flood fill (breadth-first search). Walk every cell; when an
 * unvisited empty cell turns up, it is the seed of a new cluster, so push it
 * on a stack and repeatedly pop a cell and push its unvisited empty neighbors
 * until the stack empties. Everything reached that way is one cluster.
 *
 * "6-connected" means only the six face neighbors count as touching
 * (+-x, +-y, +-z). Cells meeting at an edge or corner are NOT one cluster.
 * That is the standard convention for pore networks: a corner touch is not a
 * path a gas molecule could travel through.
 *
 * The stack is explicit rather than using recursion: a single cluster can hold
 * thousands of cells, and that much recursion depth risks a stack overflow. */
static int compute_empty_clusters(const PorosityField *field,
                                 FieldValidationMetrics *out) {
    const int64_t nx = (int64_t)field->width;
    const int64_t ny = (int64_t)field->height;
    const int64_t nz = (int64_t)field->depth;
    const uint64_t n = out->cell_count;
    /* The 6 face neighbors as (dx, dy, dz) offsets. */
    static const int neighbor_offsets[6][3] = {
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}
    };
    uint8_t *visited;
    uint64_t *stack;
    uint64_t stack_size = 0U;

    visited = (uint8_t *)calloc((size_t)n, sizeof(uint8_t));
    /* Worst case every empty cell sits on the stack at once, so size it for all
       cells and we can never overflow it. */
    stack = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    if (visited == NULL || stack == NULL) {
        fprintf(stderr, "Validation: unable to allocate cluster workspace\n");
        free(visited);
        free(stack);
        return EXIT_FAILURE;
    }

    out->empty_cell_count = 0U;
    out->cluster_count = 0U;
    out->cluster_largest_size = 0U;

    for (int64_t z0 = 0; z0 < nz; ++z0) {
        for (int64_t y0 = 0; y0 < ny; ++y0) {
            for (int64_t x0 = 0; x0 < nx; ++x0) {
                uint64_t seed = (uint64_t)((z0 * ny + y0) * nx + x0);
                uint64_t cluster_size = 0U;

                if (visited[seed] ||
                    field->values[seed] < VALIDATION_EMPTY_THRESHOLD) {
                    continue;
                }

                /* New cluster: flood fill outward from this seed. */
                out->cluster_count += 1U;
                visited[seed] = 1U;
                stack[stack_size++] = seed;

                while (stack_size > 0U) {
                    uint64_t current = stack[--stack_size];
                    /* Recover (x, y, z) from the flat index. */
                    int64_t x = (int64_t)(current % (uint64_t)nx);
                    int64_t y = (int64_t)((current / (uint64_t)nx) % (uint64_t)ny);
                    int64_t z = (int64_t)(current / (uint64_t)(nx * ny));

                    cluster_size += 1U;

                    for (int k = 0; k < 6; ++k) {
                        int64_t xn = x + neighbor_offsets[k][0];
                        int64_t yn = y + neighbor_offsets[k][1];
                        int64_t zn = z + neighbor_offsets[k][2];
                        uint64_t neighbor;

                        /* Skip neighbors that fall outside the field. */
                        if (xn < 0 || xn >= nx || yn < 0 || yn >= ny ||
                            zn < 0 || zn >= nz) {
                            continue;
                        }
                        neighbor = (uint64_t)((zn * ny + yn) * nx + xn);
                        if (visited[neighbor] ||
                            field->values[neighbor] < VALIDATION_EMPTY_THRESHOLD) {
                            continue;
                        }
                        /* Mark on PUSH, not on pop: otherwise the same cell can
                           be pushed by several neighbors and counted twice. */
                        visited[neighbor] = 1U;
                        stack[stack_size++] = neighbor;
                    }
                }

                out->empty_cell_count += cluster_size;
                if (cluster_size > out->cluster_largest_size) {
                    out->cluster_largest_size = cluster_size;
                }
            }
        }
    }

    out->cluster_mean_size =
        out->cluster_count > 0U
            ? (double)out->empty_cell_count / (double)out->cluster_count
            : 0.0;

    free(visited);
    free(stack);
    return EXIT_SUCCESS;
}

/* ======================================================================
 * Public entry point
 * ====================================================================== */

void free_field_validation_metrics(FieldValidationMetrics *metrics) {
    if (metrics == NULL) {
        return;
    }
    free(metrics->cov_x);
    free(metrics->cov_y);
    free(metrics->cov_z);
    free(metrics->z_layer_mean);
    metrics->cov_x = NULL;
    metrics->cov_y = NULL;
    metrics->cov_z = NULL;
    metrics->z_layer_mean = NULL;
}

int compute_field_validation_metrics(const PorosityField *field,
                                     int requested_max_lag,
                                     FieldValidationMetrics *out) {
    size_t lag_array_length;

    if (field == NULL || field->values == NULL || out == NULL ||
        field->width == 0U || field->height == 0U || field->depth == 0U) {
        return EXIT_FAILURE;
    }

    /* Start from all zeros so every pointer is NULL. Then, if we fail halfway,
       free_field_validation_metrics() only frees what was really allocated. */
    memset(out, 0, sizeof(*out));
    out->nx = field->width;
    out->ny = field->height;
    out->nz = field->depth;
    out->cell_count = (uint64_t)field->width * field->height * field->depth;

    /* Step 2a */
    if (compute_one_point_statistics(field, out) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    /* Step 2b */
    out->max_lag = validation_resolve_max_lag(field, requested_max_lag);
    lag_array_length = (size_t)out->max_lag + 1U;
    out->cov_x = (double *)calloc(lag_array_length, sizeof(double));
    out->cov_y = (double *)calloc(lag_array_length, sizeof(double));
    out->cov_z = (double *)calloc(lag_array_length, sizeof(double));
    if (out->cov_x == NULL || out->cov_y == NULL || out->cov_z == NULL) {
        fprintf(stderr, "Validation: unable to allocate covariance arrays\n");
        free_field_validation_metrics(out);
        return EXIT_FAILURE;
    }

    /* (1,0,0) walks along x, (0,1,0) along y, (0,0,1) along z. */
    if (compute_axis_covariance(field, 1, 0, 0, out->variance, out->max_lag, out->cov_x) != EXIT_SUCCESS ||
        compute_axis_covariance(field, 0, 1, 0, out->variance, out->max_lag, out->cov_y) != EXIT_SUCCESS ||
        compute_axis_covariance(field, 0, 0, 1, out->variance, out->max_lag, out->cov_z) != EXIT_SUCCESS) {
        free_field_validation_metrics(out);
        return EXIT_FAILURE;
    }

    /* Step 2c */
    out->z_layer_mean = (double *)calloc((size_t)field->depth, sizeof(double));
    if (out->z_layer_mean == NULL) {
        fprintf(stderr, "Validation: unable to allocate z layer profile\n");
        free_field_validation_metrics(out);
        return EXIT_FAILURE;
    }
    compute_z_layer_means(field, out->z_layer_mean);

    /* Step 2d */
    if (compute_empty_clusters(field, out) != EXIT_SUCCESS) {
        free_field_validation_metrics(out);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void running_statistic_reset(RunningStatistic *stat) {
    if (stat == NULL) {
        return;
    }
    stat->count = 0U;
    stat->mean = 0.0;
    stat->m2 = 0.0;
}

void running_statistic_add(RunningStatistic *stat, double value) {
    double delta;

    if (stat == NULL) {
        return;
    }

    stat->count += 1U;
    delta = value - stat->mean;
    stat->mean += delta / (double)stat->count;
    stat->m2 += delta * (value - stat->mean);
}

double running_statistic_mean(const RunningStatistic *stat) {
    if (stat == NULL || stat->count == 0U) {
        return 0.0;
    }
    return stat->mean;
}

double running_statistic_variance(const RunningStatistic *stat) {
    if (stat == NULL || stat->count < 2U) {
        return 0.0;
    }
    return stat->m2 / (double)(stat->count - 1U);
}

double running_statistic_stddev(const RunningStatistic *stat) {
    return sqrt(running_statistic_variance(stat));
}

int ensemble_validation_accumulator_init(EnsembleValidationAccumulator *acc,
                                        int max_lag, uint32_t nz) {
    size_t cov_length;
    size_t i;

    if (acc == NULL) {
        return EXIT_FAILURE;
    }

    memset(acc, 0, sizeof(*acc));
    acc->max_lag = max_lag;
    acc->nz = nz;
    acc->realization_count = 0U;

    cov_length = (size_t)max_lag + 1U;
    acc->cov_x = (RunningStatistic *)calloc(cov_length, sizeof(RunningStatistic));
    acc->cov_y = (RunningStatistic *)calloc(cov_length, sizeof(RunningStatistic));
    acc->cov_z = (RunningStatistic *)calloc(cov_length, sizeof(RunningStatistic));
    acc->z_layer_mean = (RunningStatistic *)calloc((size_t)nz, sizeof(RunningStatistic));
    if (acc->cov_x == NULL || acc->cov_y == NULL || acc->cov_z == NULL ||
        (nz > 0U && acc->z_layer_mean == NULL)) {
        free(acc->cov_x);
        free(acc->cov_y);
        free(acc->cov_z);
        free(acc->z_layer_mean);
        acc->cov_x = NULL;
        acc->cov_y = NULL;
        acc->cov_z = NULL;
        acc->z_layer_mean = NULL;
        fprintf(stderr, "Validation: unable to allocate ensemble accumulator arrays\n");
        return EXIT_FAILURE;
    }

    for (i = 0U; i < VALIDATION_QUANTILE_COUNT; ++i) {
        running_statistic_reset(&acc->quantiles[i]);
    }
    running_statistic_reset(&acc->mean);
    running_statistic_reset(&acc->variance);
    running_statistic_reset(&acc->skewness);
    running_statistic_reset(&acc->min_value);
    running_statistic_reset(&acc->max_value);
    running_statistic_reset(&acc->empty_fraction);
    running_statistic_reset(&acc->cluster_count);
    running_statistic_reset(&acc->cluster_mean_size);
    running_statistic_reset(&acc->cluster_largest_size);
    if (acc->z_layer_mean != NULL) {
        for (i = 0U; i < (size_t)nz; ++i) {
            running_statistic_reset(&acc->z_layer_mean[i]);
        }
    }
    for (i = 0U; i < cov_length; ++i) {
        running_statistic_reset(&acc->cov_x[i]);
        running_statistic_reset(&acc->cov_y[i]);
        running_statistic_reset(&acc->cov_z[i]);
    }
    return EXIT_SUCCESS;
}

int ensemble_validation_accumulator_add(EnsembleValidationAccumulator *acc,
                                       const FieldValidationMetrics *metrics) {
    size_t i;

    if (acc == NULL || metrics == NULL) {
        return EXIT_FAILURE;
    }

    if (acc->max_lag != metrics->max_lag) {
        fprintf(stderr,
                "Validation: ensemble max_lag mismatch (%d vs %d); cannot average realizations\n",
                acc->max_lag, metrics->max_lag);
        return EXIT_FAILURE;
    }
    if (acc->nz != metrics->nz) {
        fprintf(stderr,
                "Validation: ensemble z-depth mismatch (%u vs %u); cannot average realizations\n",
                acc->nz, metrics->nz);
        return EXIT_FAILURE;
    }

    running_statistic_add(&acc->mean, metrics->mean);
    running_statistic_add(&acc->variance, metrics->variance);
    running_statistic_add(&acc->skewness, metrics->skewness);
    running_statistic_add(&acc->min_value, metrics->min);
    running_statistic_add(&acc->max_value, metrics->max);
    running_statistic_add(&acc->empty_fraction, metrics->empty_fraction);
    running_statistic_add(&acc->cluster_count, (double)metrics->cluster_count);
    running_statistic_add(&acc->cluster_mean_size, metrics->cluster_mean_size);
    running_statistic_add(&acc->cluster_largest_size,
                          (double)metrics->cluster_largest_size);

    for (i = 0U; i < VALIDATION_QUANTILE_COUNT; ++i) {
        running_statistic_add(&acc->quantiles[i], metrics->quantiles[i]);
    }

    for (i = 0U; i < (size_t)acc->max_lag + 1U; ++i) {
        running_statistic_add(&acc->cov_x[i], metrics->cov_x[i]);
        running_statistic_add(&acc->cov_y[i], metrics->cov_y[i]);
        running_statistic_add(&acc->cov_z[i], metrics->cov_z[i]);
    }

    if (acc->z_layer_mean != NULL && metrics->z_layer_mean != NULL) {
        for (i = 0U; i < (size_t)acc->nz; ++i) {
            running_statistic_add(&acc->z_layer_mean[i], metrics->z_layer_mean[i]);
        }
    }

    acc->realization_count += 1U;
    return EXIT_SUCCESS;
}

void free_ensemble_validation_accumulator(EnsembleValidationAccumulator *acc) {
    if (acc == NULL) {
        return;
    }
    free(acc->cov_x);
    free(acc->cov_y);
    free(acc->cov_z);
    free(acc->z_layer_mean);
    acc->cov_x = NULL;
    acc->cov_y = NULL;
    acc->cov_z = NULL;
    acc->z_layer_mean = NULL;
    acc->realization_count = 0U;
    acc->max_lag = 0;
    acc->nz = 0U;
}

int compute_data_field_validation_metrics(const char *stack_path,
                                         double voxel_resolution,
                                         double filter_resolution,
                                         int requested_max_lag,
                                         FieldValidationMetrics *metrics_out,
                                         double *block_spacing_out) {
    BinaryVolume volume = {0U, 0U, 0U, NULL};
    PorosityField field = {0U, 0U, 0U, NULL};
    uint32_t block_voxels;
    double block_spacing = 0.0;
    int status = EXIT_FAILURE;

    if (stack_path == NULL || metrics_out == NULL || voxel_resolution <= 0.0 ||
        filter_resolution <= 0.0) {
        return EXIT_FAILURE;
    }

    if (load_binary_volume_from_tiff(stack_path, &volume) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    block_voxels = (uint32_t)llround(filter_resolution / voxel_resolution);
    if (block_voxels == 0U) {
        block_voxels = 1U;
    }
    if ((block_voxels % 2U) == 0U) {
        block_voxels += 1U;
    }
    block_spacing = (double)block_voxels * voxel_resolution;

    if (compute_porosity_field(volume.voxels, volume.width, volume.height,
                               volume.depth, filter_resolution, voxel_resolution,
                               &field) != EXIT_SUCCESS) {
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (compute_field_validation_metrics(&field, requested_max_lag, metrics_out) !=
        EXIT_SUCCESS) {
        free_porosity_field(&field);
        free_binary_volume(&volume);
        return EXIT_FAILURE;
    }

    if (block_spacing_out != NULL) {
        *block_spacing_out = block_spacing;
    }

    free_porosity_field(&field);
    free_binary_volume(&volume);
    status = EXIT_SUCCESS;
    return status;
}

void warn_if_validation_grids_differ(const FieldValidationMetrics *data_metrics,
                                     double data_block_spacing,
                                     const PorosityRealizationGridSpec *simulation_grid,
                                     FILE *out) {
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double spacing_x;
    double spacing_y;
    double spacing_z;
    int mismatch = 0;

    if (data_metrics == NULL || simulation_grid == NULL || out == NULL) {
        return;
    }

    spacing_x = simulation_grid->nx > 1U
                    ? simulation_grid->dim_x / (double)(simulation_grid->nx - 1U)
                    : simulation_grid->dim_x;
    spacing_y = simulation_grid->ny > 1U
                    ? simulation_grid->dim_y / (double)(simulation_grid->ny - 1U)
                    : simulation_grid->dim_y;
    spacing_z = simulation_grid->nz > 1U
                    ? simulation_grid->dim_z / (double)(simulation_grid->nz - 1U)
                    : simulation_grid->dim_z;

    dx = fabs(spacing_x - data_block_spacing) / data_block_spacing;
    dy = fabs(spacing_y - data_block_spacing) / data_block_spacing;
    dz = fabs(spacing_z - data_block_spacing) / data_block_spacing;
    if (dx > 0.01 || dy > 0.01 || dz > 0.01 ||
        data_metrics->nx != simulation_grid->nx ||
        data_metrics->ny != simulation_grid->ny ||
        data_metrics->nz != simulation_grid->nz) {
        mismatch = 1;
    }

    if (mismatch) {
        fprintf(out,
                "Warning: realization grid support differs from data support; variance and covariance comparisons are only meaningful on matching support. "
                "data_field=(%u,%u,%u) block_spacing=%.12g, simulation_grid=(%u,%u,%u) spacing=(%.12g,%.12g,%.12g)\n",
                data_metrics->nx, data_metrics->ny, data_metrics->nz,
                data_block_spacing, simulation_grid->nx, simulation_grid->ny,
                simulation_grid->nz, spacing_x, spacing_y, spacing_z);
    }
}

static double resolve_stat_value(const RunningStatistic *stat) {
    if (stat == NULL) {
        return 0.0;
    }
    return running_statistic_mean(stat);
}

static double resolve_stat_stddev(const RunningStatistic *stat) {
    if (stat == NULL) {
        return 0.0;
    }
    return running_statistic_stddev(stat);
}

static double compute_z_score(double data_value, const RunningStatistic *ens_stat) {
    double mean = resolve_stat_value(ens_stat);
    double stddev = resolve_stat_stddev(ens_stat);

    if (ens_stat == NULL || ens_stat->count < 2U || stddev == 0.0) {
        return NAN;
    }
    return (data_value - mean) / stddev;
}

static double compute_efold_lag(const double *cov_curve, int max_lag, double c0) {
    double threshold;
    int h;

    if (cov_curve == NULL || max_lag <= 0 || c0 <= 0.0) {
        return NAN;
    }

    threshold = c0 / exp(1.0);
    for (h = 1; h <= max_lag; ++h) {
        double value = cov_curve[h];
        if (value <= threshold) {
            if (h == 1) {
                return 1.0;
            }
            if (cov_curve[h - 1] <= threshold) {
                return (double)h;
            }
            return (double)(h - 1) +
                   ((threshold - cov_curve[h - 1]) /
                    (cov_curve[h] - cov_curve[h - 1]));
        }
    }
    return NAN;
}

int write_realization_validation_report(FILE *out,
                                       const char *stack_path,
                                       const FieldValidationMetrics *data_metrics,
                                       double data_block_spacing,
                                       const PorosityRealizationGridSpec *simulation_grid,
                                       const EnsembleValidationAccumulator *acc,
                                       int include_z_profile) {
    size_t i;
    double cov_rms_x = 0.0;
    double cov_rms_y = 0.0;
    double cov_rms_z = 0.0;
    double quantile_rms = 0.0;
    double quantile_max_abs = 0.0;
    double empty_fraction_error = 0.0;
    double largest_cluster_ratio = NAN;
    double efold_data_x = NAN;
    double efold_data_y = NAN;
    double efold_data_z = NAN;
    double efold_ens_x = NAN;
    double efold_ens_y = NAN;
    double efold_ens_z = NAN;
    double sim_spacing_x = 0.0;
    double sim_spacing_y = 0.0;
    double sim_spacing_z = 0.0;

    if (out == NULL || stack_path == NULL || data_metrics == NULL || acc == NULL) {
        return EXIT_FAILURE;
    }
    if (simulation_grid != NULL) {
        sim_spacing_x = simulation_grid->nx > 1U
                            ? simulation_grid->dim_x / (double)(simulation_grid->nx - 1U)
                            : simulation_grid->dim_x;
        sim_spacing_y = simulation_grid->ny > 1U
                            ? simulation_grid->dim_y / (double)(simulation_grid->ny - 1U)
                            : simulation_grid->dim_y;
        sim_spacing_z = simulation_grid->nz > 1U
                            ? simulation_grid->dim_z / (double)(simulation_grid->nz - 1U)
                            : simulation_grid->dim_z;
    }

    if (data_metrics->max_lag > 0) {
        for (int h = 1; h <= data_metrics->max_lag; ++h) {
            double dx = data_metrics->cov_x[h] - resolve_stat_value(&acc->cov_x[h]);
            double dy = data_metrics->cov_y[h] - resolve_stat_value(&acc->cov_y[h]);
            double dz = data_metrics->cov_z[h] - resolve_stat_value(&acc->cov_z[h]);
            cov_rms_x += dx * dx;
            cov_rms_y += dy * dy;
            cov_rms_z += dz * dz;
        }
        cov_rms_x = sqrt(cov_rms_x / (double)data_metrics->max_lag) /
                   data_metrics->variance;
        cov_rms_y = sqrt(cov_rms_y / (double)data_metrics->max_lag) /
                   data_metrics->variance;
        cov_rms_z = sqrt(cov_rms_z / (double)data_metrics->max_lag) /
                   data_metrics->variance;
    }

    for (i = 0U; i < VALIDATION_QUANTILE_COUNT; ++i) {
        double diff = data_metrics->quantiles[i] - resolve_stat_value(&acc->quantiles[i]);
        quantile_rms += diff * diff;
        if (fabs(diff) > quantile_max_abs) {
            quantile_max_abs = fabs(diff);
        }
    }
    quantile_rms = sqrt(quantile_rms / (double)VALIDATION_QUANTILE_COUNT);

    empty_fraction_error = fabs(data_metrics->empty_fraction -
                                resolve_stat_value(&acc->empty_fraction));
    if (data_metrics->cluster_largest_size != 0U) {
        largest_cluster_ratio = resolve_stat_value(&acc->cluster_largest_size) /
                                (double)data_metrics->cluster_largest_size;
    }

    efold_data_x = compute_efold_lag(data_metrics->cov_x, data_metrics->max_lag,
                                     data_metrics->cov_x[0]);
    efold_data_y = compute_efold_lag(data_metrics->cov_y, data_metrics->max_lag,
                                     data_metrics->cov_y[0]);
    efold_data_z = compute_efold_lag(data_metrics->cov_z, data_metrics->max_lag,
                                     data_metrics->cov_z[0]);
    efold_ens_x = compute_efold_lag(acc->cov_x != NULL ? &acc->cov_x[0].mean : NULL,
                                   data_metrics->max_lag,
                                   data_metrics->cov_x[0]);
    efold_ens_y = compute_efold_lag(acc->cov_y != NULL ? &acc->cov_y[0].mean : NULL,
                                   data_metrics->max_lag,
                                   data_metrics->cov_y[0]);
    efold_ens_z = compute_efold_lag(acc->cov_z != NULL ? &acc->cov_z[0].mean : NULL,
                                   data_metrics->max_lag,
                                   data_metrics->cov_z[0]);

    fprintf(out,
            "# Realization validation\n"
            "# stack_path %s\n"
            "# realization_count %llu\n"
            "# max_lag %d\n"
            "# data_field_dimensions %u %u %u\n"
            "# data_block_spacing %.12g\n"
            "# simulation_grid_dimensions %u %u %u\n"
            "# simulation_grid_spacing %.12g %.12g %.12g\n"
            "# z_score_definition z = (data - ensemble_mean) / ensemble_stddev; n/a when stddev=0 or count<2\n",
            stack_path, (unsigned long long)acc->realization_count, data_metrics->max_lag,
            data_metrics->nx, data_metrics->ny, data_metrics->nz, data_block_spacing,
            simulation_grid->nx, simulation_grid->ny, simulation_grid->nz,
            sim_spacing_x, sim_spacing_y, sim_spacing_z);

    fprintf(out, "# One-point statistics\n# metric data ensemble_mean ensemble_stddev z_score\n");
    fprintf(out, "mean %.6g %.6g %.6g %s\n",
            data_metrics->mean,
            resolve_stat_value(&acc->mean),
            resolve_stat_stddev(&acc->mean),
            isfinite(compute_z_score(data_metrics->mean, &acc->mean)) ? "" : "n/a");
    fprintf(out, "variance %.6g %.6g %.6g %s\n",
            data_metrics->variance,
            resolve_stat_value(&acc->variance),
            resolve_stat_stddev(&acc->variance),
            isfinite(compute_z_score(data_metrics->variance, &acc->variance)) ? "" : "n/a");
    fprintf(out, "skewness %.6g %.6g %.6g %s\n",
            data_metrics->skewness,
            resolve_stat_value(&acc->skewness),
            resolve_stat_stddev(&acc->skewness),
            isfinite(compute_z_score(data_metrics->skewness, &acc->skewness)) ? "" : "n/a");
    fprintf(out, "min %.6g %.6g %.6g %s\n",
            data_metrics->min,
            resolve_stat_value(&acc->min_value),
            resolve_stat_stddev(&acc->min_value),
            isfinite(compute_z_score(data_metrics->min, &acc->min_value)) ? "" : "n/a");
    fprintf(out, "max %.6g %.6g %.6g %s\n",
            data_metrics->max,
            resolve_stat_value(&acc->max_value),
            resolve_stat_stddev(&acc->max_value),
            isfinite(compute_z_score(data_metrics->max, &acc->max_value)) ? "" : "n/a");
    fprintf(out, "empty_fraction %.6g %.6g %.6g %s\n",
            data_metrics->empty_fraction,
            resolve_stat_value(&acc->empty_fraction),
            resolve_stat_stddev(&acc->empty_fraction),
            isfinite(compute_z_score(data_metrics->empty_fraction, &acc->empty_fraction)) ? "" : "n/a");
    for (i = 0U; i < VALIDATION_QUANTILE_COUNT; ++i) {
        double z = compute_z_score(data_metrics->quantiles[i], &acc->quantiles[i]);
        fprintf(out, "%s %.6g %.6g %.6g %s\n",
                i == 0 ? "q01" : i == 1 ? "q05" : i == 2 ? "q25" : i == 3 ? "q50" : i == 4 ? "q75" : i == 5 ? "q95" : "q99",
                data_metrics->quantiles[i],
                resolve_stat_value(&acc->quantiles[i]),
                resolve_stat_stddev(&acc->quantiles[i]),
                isfinite(z) ? "" : "n/a");
    }
    fprintf(out, "cluster_count %.6g %.6g %.6g %s\n",
            (double)data_metrics->cluster_count,
            resolve_stat_value(&acc->cluster_count),
            resolve_stat_stddev(&acc->cluster_count),
            isfinite(compute_z_score((double)data_metrics->cluster_count, &acc->cluster_count)) ? "" : "n/a");
    fprintf(out, "cluster_mean_size %.6g %.6g %.6g %s\n",
            data_metrics->cluster_mean_size,
            resolve_stat_value(&acc->cluster_mean_size),
            resolve_stat_stddev(&acc->cluster_mean_size),
            isfinite(compute_z_score(data_metrics->cluster_mean_size, &acc->cluster_mean_size)) ? "" : "n/a");
    fprintf(out, "cluster_largest_size %.6g %.6g %.6g %s\n",
            (double)data_metrics->cluster_largest_size,
            resolve_stat_value(&acc->cluster_largest_size),
            resolve_stat_stddev(&acc->cluster_largest_size),
            isfinite(compute_z_score((double)data_metrics->cluster_largest_size, &acc->cluster_largest_size)) ? "" : "n/a");

    fprintf(out, "# Axis covariance\n# lag C_data_x C_ens_mean_x C_ens_sd_x C_data_y C_ens_mean_y C_ens_sd_y C_data_z C_ens_mean_z C_ens_sd_z\n");
    for (int h = 0; h <= data_metrics->max_lag; ++h) {
        fprintf(out, "%d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
                h,
                data_metrics->cov_x[h],
                resolve_stat_value(&acc->cov_x[h]),
                resolve_stat_stddev(&acc->cov_x[h]),
                data_metrics->cov_y[h],
                resolve_stat_value(&acc->cov_y[h]),
                resolve_stat_stddev(&acc->cov_y[h]),
                data_metrics->cov_z[h],
                resolve_stat_value(&acc->cov_z[h]),
                resolve_stat_stddev(&acc->cov_z[h]));
    }

    if (include_z_profile) {
        fprintf(out, "# Z layer profile\n# layer mean_data mean_ens_mean mean_ens_sd\n");
        for (i = 0U; i < (size_t)data_metrics->nz; ++i) {
            fprintf(out, "%zu %.9g %.9g %.9g\n",
                    i, data_metrics->z_layer_mean[i],
                    resolve_stat_value(&acc->z_layer_mean[i]),
                    resolve_stat_stddev(&acc->z_layer_mean[i]));
        }
    }

    fprintf(out, "# Summary scores\nscore value\n");
    fprintf(out, "cov_rms_x %.6g\n", cov_rms_x);
    fprintf(out, "cov_rms_y %.6g\n", cov_rms_y);
    fprintf(out, "cov_rms_z %.6g\n", cov_rms_z);
    fprintf(out, "quantile_rms %.6g\n", quantile_rms);
    fprintf(out, "quantile_max_abs_error %.6g\n", quantile_max_abs);
    fprintf(out, "empty_fraction_error %.6g\n", empty_fraction_error);
    fprintf(out, "largest_cluster_ratio %.6g\n", isfinite(largest_cluster_ratio) ? largest_cluster_ratio : NAN);
    fprintf(out, "efold_data_x %.6g\n", efold_data_x);
    fprintf(out, "efold_data_y %.6g\n", efold_data_y);
    fprintf(out, "efold_data_z %.6g\n", efold_data_z);
    fprintf(out, "efold_ens_x %.6g\n", efold_ens_x);
    fprintf(out, "efold_ens_y %.6g\n", efold_ens_y);
    fprintf(out, "efold_ens_z %.6g\n", efold_ens_z);
    if (simulation_grid != NULL) {
        fprintf(out, "efold_data_x_phys %.6g\n", isfinite(efold_data_x) ? efold_data_x * data_block_spacing : NAN);
        fprintf(out, "efold_data_y_phys %.6g\n", isfinite(efold_data_y) ? efold_data_y * data_block_spacing : NAN);
        fprintf(out, "efold_data_z_phys %.6g\n", isfinite(efold_data_z) ? efold_data_z * data_block_spacing : NAN);
    }

    return EXIT_SUCCESS;
}
