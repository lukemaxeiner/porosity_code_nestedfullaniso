#ifndef POROSITY_COVARIANCE_REALIZATION_PREP_H
#define POROSITY_COVARIANCE_REALIZATION_PREP_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "porosity_covariance_lag.h"

typedef struct {
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    double mean;
    double variance;
    double *values;
} CovarianceGrid3D;

typedef struct {
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    int origin_x;
    int origin_y;
    int origin_z;
    double mean;
    double variance;
    double *values;
} SymmetricCovarianceGrid3D;

typedef struct {
    int tail_is_noisy;
    uint32_t total_tail_bins;
    uint32_t flagged_tail_bins;
    double tail_noise_fraction;
    double max_relative_fluctuation;
    double min_mean_samples;
    double recommended_max_radius;
} CovarianceTailNoiseSummary;

#define POROSITY_NESTED_EXPONENTIAL_COMPONENTS 2U

typedef enum {
    POROSITY_COVARIANCE_MODEL_FULLY_ANISOTROPIC_EXPONENTIAL = 1,
    POROSITY_COVARIANCE_MODEL_NESTED_FULLY_ANISOTROPIC_EXPONENTIAL = 2
} PorosityCovarianceModelType;

typedef struct {
    PorosityCovarianceModelType model_type;
    double sill;

    /*
     * Nested fully anisotropic exponential model:
     *
     * C(rx,ry,rz) = sill * sum_k weight[k] * exp(-sqrt(
     *     (rx / Lx[k])^2 +
     *     (ry / Ly[k])^2 +
     *     (rz / Lz[k])^2
     * ))
     *
     * For the nested model, component_count should be 2 and the weights
     * should be nonnegative and sum to 1.
     */
    uint32_t component_count;
    double component_weight[POROSITY_NESTED_EXPONENTIAL_COMPONENTS];
    double component_x_correlation_length[POROSITY_NESTED_EXPONENTIAL_COMPONENTS];
    double component_y_correlation_length[POROSITY_NESTED_EXPONENTIAL_COMPONENTS];
    double component_z_correlation_length[POROSITY_NESTED_EXPONENTIAL_COMPONENTS];

    /*
     * Legacy compatibility:
     * keep these fields so older output and generation paths can still compile.
     * For the nested model, set these to the long/dominant component lengths
     * or the max component length in each direction.
     */
    double correlation_length;
    double x_correlation_length;
    double y_correlation_length;
    double z_correlation_length;

    /*
     * Legacy compatibility with the transverse-isotropic code.
     */
    double in_plane_correlation_length;
    double through_thickness_correlation_length;

    double rms_error;
    double weighted_rms_error;
    double directional_rms_x;
    double directional_rms_y;
    double directional_rms_z;
    double small_lag_rms_error;
    double recommended_radius_rms_error;
    double aic;
    double bic;
    uint64_t sample_count;
    double zero_lag_error;
    double zero_lag_relative_error;
} ExponentialCovarianceModel;

typedef struct {
    int is_valid;
    int sill_is_positive;
    int correlation_length_is_positive;
    int component_weights_are_valid;
    int has_enough_samples;
    int zero_lag_matches_sill;
    int rms_error_is_reasonable;
    double zero_lag_relative_error;
    double rms_error_ratio;
} ExponentialCovarianceValidation;

typedef struct {
    int has_negative_lobes;
    uint64_t negative_count;
    double min_covariance;
    double most_negative_radius;
} NegativeCovarianceSummary;

typedef struct {
    double zero_lag_covariance;
    double x_efold_distance;
    double y_efold_distance;
    double z_efold_distance;
    double max_to_min_ratio;
    int is_anisotropic;
} CovarianceAnisotropySummary;

typedef struct {
    uint32_t nx;
    uint32_t ny;
    uint32_t nz;
    double dim_x;
    double dim_y;
    double dim_z;
} PorosityRealizationGridSpec;

typedef struct {
    int is_valid;
    int all_finite;
    int zero_lag_matches_sill;
    int all_nonnegative;
    double zero_lag_relative_error;
    double min_value;
    double max_value;
} EvaluatedCovarianceGridValidation;

typedef struct {
    double domain_x_over_correlation_length;
    double domain_y_over_correlation_length;
    double domain_z_over_correlation_length;
    int enough_domain_x;
    int enough_domain_y;
    int enough_domain_z;
    int overall_adequate;
} SimulationGridCompatibilitySummary;

typedef struct {
    PorosityCovarianceSample *samples;
    uint32_t lag_width;
    uint32_t lag_height;
    uint32_t lag_depth;
    double filter_resolution;
    double variance;
    CovarianceGrid3D covariance_grid;
    SymmetricCovarianceGrid3D symmetric_grid;
    CovarianceTailNoiseSummary tail_noise;
    NegativeCovarianceSummary negative_covariance;
    CovarianceAnisotropySummary anisotropy;
    ExponentialCovarianceModel exponential_model;
    ExponentialCovarianceValidation exponential_validation;
    PorosityRealizationGridSpec simulation_grid;
    SimulationGridCompatibilitySummary simulation_grid_compatibility;
    CovarianceGrid3D evaluated_simulation_covariance_grid;
    EvaluatedCovarianceGridValidation evaluated_grid_validation;
} PorosityRealizationPrepResult;

typedef struct {
    double mean_porosity;
    double variance;
    double filter_resolution;
    ExponentialCovarianceModel exponential_model;
    ExponentialCovarianceValidation exponential_validation;
} FittedCovarianceModelArtifact;

int reconstruct_covariance_grid_from_3dlag_samples(
    const PorosityCovarianceSample *samples, uint32_t width, uint32_t height,
    uint32_t depth, double variance, CovarianceGrid3D *grid_out);

int symmetrize_covariance_grid_by_reflection(
    const CovarianceGrid3D *source, SymmetricCovarianceGrid3D *grid_out);

int analyze_covariance_tail_noise(
    const PorosityCovarianceSample *samples, const CovarianceGrid3D *grid,
    double filter_resolution, CovarianceTailNoiseSummary *summary_out);

int fit_nested_full_anisotropic_exponential_covariance_model(
const CovarianceGrid3D *grid,
const PorosityCovarianceSample *samples,
double filter_resolution,
double recommended_max_radius,
ExponentialCovarianceModel *model_out);

int validate_exponential_covariance_model(
    const CovarianceGrid3D *grid, const ExponentialCovarianceModel *model,
    ExponentialCovarianceValidation *validation_out);

int analyze_negative_covariance_lobes(const CovarianceGrid3D *grid,
                                      double filter_resolution,
                                      NegativeCovarianceSummary *summary_out);

int analyze_covariance_anisotropy(const CovarianceGrid3D *grid,
                                  double filter_resolution,
                                  CovarianceAnisotropySummary *summary_out);

int assess_realization_grid_compatibility(
    const PorosityRealizationGridSpec *simulation_grid,
    const ExponentialCovarianceModel *model,
    SimulationGridCompatibilitySummary *summary_out);

int evaluate_exponential_covariance_on_realization_grid(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid,
    CovarianceGrid3D *grid_out);

int validate_evaluated_covariance_grid(
    const CovarianceGrid3D *grid, const ExponentialCovarianceModel *model,
    EvaluatedCovarianceGridValidation *validation_out);

int write_realization_prep_report(FILE *out,
                                  const PorosityRealizationPrepResult *result);

int read_3dlag_covariance_table(const char *table_path,
                                double filter_resolution,
                                PorosityCovarianceSample **samples_out,
                                uint32_t *width_out, uint32_t *height_out,
                                uint32_t *depth_out);

int prepare_realization_from_3dlag_table(
    const char *table_path, double filter_resolution, double variance,
    const PorosityRealizationGridSpec *simulation_grid,
    PorosityRealizationPrepResult *result_out);

int process_realization_prep_from_stack(
    const char *stack_path, double voxel_resolution, double filter_resolution,
    const PorosityRealizationGridSpec *simulation_grid,
    uint32_t realization_count, FILE *out);

int write_fitted_covariance_model(FILE *out,
                                  const FittedCovarianceModelArtifact *model);

int read_fitted_covariance_model(const char *model_path,
                                 FittedCovarianceModelArtifact *model_out);

int process_realization_fit_from_stack(
    const char *stack_path, double voxel_resolution, double filter_resolution,
    const PorosityRealizationGridSpec *simulation_grid, FILE *out);

#endif
