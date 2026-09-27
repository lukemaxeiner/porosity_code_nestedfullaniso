#ifndef POROSITY_FIELD_H
#define POROSITY_FIELD_H

#include <stdint.h>
#include <stdio.h>

typedef struct {
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    double *values;
} PorosityField;

int compute_porosity_field(const uint8_t *binary, uint32_t width,
                           uint32_t height, uint32_t depth,
                           double filter_size, double voxel_size,
                           PorosityField *field);
void free_porosity_field(PorosityField *field);
int write_porosity_field_tiff(const char *path, const PorosityField *field);
int process_porosity_field(const char *stack_path, double voxel_resolution,
                           double filter_resolution, FILE *out);

#endif
