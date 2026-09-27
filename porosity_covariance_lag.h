#ifndef POROSITY_COVARIANCE_LAG_H
#define POROSITY_COVARIANCE_LAG_H

#include <stdint.h>
#include <stdio.h>

#include "porosity_field.h"

typedef enum {
    POROSITY_AXIS_X = 0,
    POROSITY_AXIS_Y = 1,
    POROSITY_AXIS_Z = 2
} PorosityAxis;

typedef enum {
    POROSITY_LAG_MODE_NO = 0,
    POROSITY_LAG_MODE_3DLAG = 1,
    POROSITY_LAG_MODE_CURVE = 2
} PorosityLagMode;

typedef struct {
    int lag_x;
    int lag_y;
    int lag_z;
    uint64_t samples;
    double mu;
    double covariance;
} PorosityCovarianceSample;

int compute_porosity_covariance_sample(const PorosityField *field, int lag_x,
                                       int lag_y, int lag_z,
                                       PorosityCovarianceSample *sample);
uint32_t porosity_covariance_axis_lag_count(const PorosityField *field,
                                            PorosityAxis axis);
int compute_porosity_covariance_curve(const PorosityField *field,
                                      PorosityAxis axis,
                                      PorosityCovarianceSample *samples_out);
int write_porosity_covariance_table(FILE *out,
                                    const PorosityCovarianceSample *samples,
                                    size_t sample_count,
                                    double filter_resolution);
int process_porosity_covariance_lag(const char *stack_path,
                                    double voxel_resolution,
                                    double filter_resolution,
                                    PorosityLagMode lag_mode,
                                    PorosityAxis axis, 
                                    int normal_score_enabled, FILE *out);

#endif
