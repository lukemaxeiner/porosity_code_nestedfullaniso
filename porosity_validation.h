#ifndef POROSITY_VALIDATION_H
#define POROSITY_VALIDATION_H

/*
 * porosity_validation: measure a porosity field so the DATA and each generated
 * REALIZATION can be compared number-for-number.
 *
 * The rule for this module: every metric is computed by the same function for
 * both kinds of field. If the data and a realization disagree, it must be
 * because the fields differ, never because they were measured differently.
 */

#include <stdint.h>
#include <stdio.h>

#include "porosity_covariance_realization_prep.h"
#include "porosity_field.h"

/* Probabilities for the quantiles we report: 1%, 5%, 25%, 50%, 75%, 95%, 99%. */
#define VALIDATION_QUANTILE_COUNT 7

/* A cell counts as "completely empty" (porosity 1.0) if it is within this
   distance of 1. Exact == on doubles is unreliable after arithmetic. */
#define VALIDATION_EMPTY_THRESHOLD (1.0 - 1e-9)

typedef struct {
    /* Size of the field that was measured. */
    uint32_t nx;
    uint32_t ny;
    uint32_t nz;
    uint64_t cell_count;

    /* --- Step 2a: one-point statistics (ignore where cells are) --- */
    double mean;
    double variance;   /* population variance: divide by N, same as the 3DLAG code */
    double skewness;   /* 0 = symmetric histogram; negative = long tail toward low porosity */
    double min;
    double max;
    double quantiles[VALIDATION_QUANTILE_COUNT]; /* same order as the probabilities */
    double empty_fraction;                       /* fraction of cells >= VALIDATION_EMPTY_THRESHOLD */

    /* --- Step 2b: two-point statistics (how neighbors relate) --- */
    int max_lag;       /* covariance computed for lags 0..max_lag cells */
    double *cov_x;     /* length max_lag + 1; cov_x[h] = C(h cells along x) */
    double *cov_y;
    double *cov_z;

    /* --- Step 2c: mean porosity of each z layer (length nz) --- */
    double *z_layer_mean;

    /* --- Step 2d: 6-connected clusters of completely empty cells --- */
    uint64_t empty_cell_count;     /* cells >= VALIDATION_EMPTY_THRESHOLD */
    uint64_t cluster_count;        /* how many separate empty clusters */
    double cluster_mean_size;      /* empty_cell_count / cluster_count */
    uint64_t cluster_largest_size; /* cells in the biggest cluster */
} FieldValidationMetrics;

/* Welford's online algorithm: numerically stable running mean and variance. */
typedef struct {
    uint64_t count;
    double mean;
    double m2;
} RunningStatistic;

void running_statistic_reset(RunningStatistic *stat);
void running_statistic_add(RunningStatistic *stat, double value);
double running_statistic_mean(const RunningStatistic *stat);
double running_statistic_variance(const RunningStatistic *stat);
double running_statistic_stddev(const RunningStatistic *stat);

typedef struct {
    uint64_t realization_count;
    int max_lag;
    uint32_t nz;

    RunningStatistic mean;
    RunningStatistic variance;
    RunningStatistic skewness;
    RunningStatistic min_value;
    RunningStatistic max_value;
    RunningStatistic empty_fraction;
    RunningStatistic quantiles[VALIDATION_QUANTILE_COUNT];
    RunningStatistic cluster_count;
    RunningStatistic cluster_mean_size;
    RunningStatistic cluster_largest_size;

    RunningStatistic *cov_x;
    RunningStatistic *cov_y;
    RunningStatistic *cov_z;
    RunningStatistic *z_layer_mean;
} EnsembleValidationAccumulator;

typedef struct {
    int enabled;
    int requested_max_lag;
    double voxel_resolution;
    double filter_resolution;
} PorosityValidationRequest;

/* The probability that goes with quantiles[i] (e.g. index 3 -> 0.50). */
double validation_quantile_probability(int index);

/* Turn the VALIDATION_MAX_LAG setting into an actual lag count for this field.
   0 means automatic (half the shortest axis). Too-large values are clamped. */
int validation_resolve_max_lag(const PorosityField *field, int requested_max_lag);

/* Fill *out with every metric for *field. Allocates the covariance arrays, so
   call free_field_validation_metrics() when done. Returns EXIT_SUCCESS/FAILURE. */
int compute_field_validation_metrics(const PorosityField *field,
                                     int requested_max_lag,
                                     FieldValidationMetrics *out);

/* Release the arrays inside *metrics (safe to call on a zeroed struct). */
void free_field_validation_metrics(FieldValidationMetrics *metrics);

int ensemble_validation_accumulator_init(EnsembleValidationAccumulator *acc,
                                        int max_lag, uint32_t nz);
int ensemble_validation_accumulator_add(EnsembleValidationAccumulator *acc,
                                       const FieldValidationMetrics *metrics);
void free_ensemble_validation_accumulator(EnsembleValidationAccumulator *acc);

int compute_data_field_validation_metrics(const char *stack_path,
                                         double voxel_resolution,
                                         double filter_resolution,
                                         int requested_max_lag,
                                         FieldValidationMetrics *metrics_out,
                                         double *block_spacing_out);
void warn_if_validation_grids_differ(const FieldValidationMetrics *data_metrics,
                                     double data_block_spacing,
                                     const PorosityRealizationGridSpec *simulation_grid,
                                     FILE *out);
int write_realization_validation_report(FILE *out,
                                       const char *stack_path,
                                       const FieldValidationMetrics *data_metrics,
                                       double data_block_spacing,
                                       const PorosityRealizationGridSpec *simulation_grid,
                                       const EnsembleValidationAccumulator *acc,
                                       int include_z_profile);

#endif
