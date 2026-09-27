#ifndef TIFF_STACK_H
#define TIFF_STACK_H

#include <stdint.h>
#include <stdio.h>

typedef struct {
    uint64_t total_voxels;
    uint64_t void_voxels;
    uint32_t width;
    uint32_t height;
    uint16_t slices;
    int has_palette;
} StackStats;

typedef struct {
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    uint8_t *voxels;
} BinaryVolume;

StackStats analyze_stack(const char *path);
double compute_porosity(StackStats stats);
int write_stack_summary(FILE *out, const char *path, StackStats stats);
int load_binary_volume_from_tiff(const char *path, BinaryVolume *volume);
void free_binary_volume(BinaryVolume *volume);

#endif
