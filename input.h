#ifndef INPUT_H
#define INPUT_H

#include "porosity_covariance_lag.h"

typedef struct {
    char *stack_path;
    double voxel_resolution;
} InputOptions;

typedef struct {
    int stack_count;
    InputOptions *stacks;
} InputBatch;

typedef enum {
    POROSITY_REALIZATION_METHOD_CHOLESKY = 0,
    POROSITY_REALIZATION_METHOD_FFT = 1
} PorosityRealizationMethod;

typedef struct {
    int global_porosity_enabled;
    int porosity_field_enabled;
    int run_full_pipeline;
    int empirical_covariance_enabled;
    int fit_covariance_model_enabled;
    int generate_porosity_realizations_enabled;
    int normal_score_transform_enabled; //normal score transform addition
    int validate_realization_enabled; //validate realizations yes/no
    int validation_max_lag; //validation max lag = N; N = 0 -> Auto
    PorosityRealizationMethod porosity_realization_method;
    double filter_resolution;
    PorosityLagMode porosity_lag_mode;
    PorosityAxis porosity_axis;
    int porosity_realization_count;
    int realization_grid_nx;
    int realization_grid_ny;
    int realization_grid_nz;
    double realization_dim_x;
    double realization_dim_y;
    double realization_dim_z;
    int porosity_ref_window_enabled;
    int ref_window_x;
    int ref_window_y;
    int ref_window_z;
    InputBatch input_data;
} RunConfig;

int parse_input(int argc, char **argv, RunConfig *config);
void free_input(InputBatch *batch);

#endif
