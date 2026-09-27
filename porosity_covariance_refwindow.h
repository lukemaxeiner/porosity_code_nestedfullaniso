#ifndef POROSITY_COVARIANCE_REFWINDOW_H
#define POROSITY_COVARIANCE_REFWINDOW_H

#include <stdint.h>
#include <stdio.h>

#include "porosity_field.h"

typedef struct {
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    double *values;
    uint32_t ref_x;
    uint32_t ref_y;
    uint32_t ref_z;
    double reference_value;
    double field_mean;
} PorosityCovarianceField;

int compute_reference_window_covariance_map(const PorosityField *field,
                                            uint32_t ref_x, uint32_t ref_y,
                                            uint32_t ref_z,
                                            PorosityCovarianceField *map);
void free_porosity_covariance_field(PorosityCovarianceField *map);
int write_porosity_covariance_field_tiff(const char *path,
                                         const PorosityCovarianceField *map);
int process_porosity_covariance_ref_window(const char *stack_path,
                                           double voxel_resolution,
                                           double filter_resolution,
                                           int ref_x, int ref_y, int ref_z,
                                           FILE *out);

#endif
