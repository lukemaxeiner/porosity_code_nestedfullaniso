#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tiffio.h>

#include "tiff.h"

static void tiff_die(const char *message, const char *path) {
    if (path != NULL) {
        fprintf(stderr, "%s: %s\n", message, path);
    } else {
        fprintf(stderr, "%s\n", message);
    }
    exit(EXIT_FAILURE);
}

static uint8_t palette_index_to_binary(uint8_t index, const uint16_t *red,
                                       const uint16_t *green,
                                       const uint16_t *blue) {
    if (red == NULL || green == NULL || blue == NULL) {
        return index == 0 ? 0 : 1;
    }

    return (red[index] == 0 && green[index] == 0 && blue[index] == 0) ? 0 : 1;
}

static uint64_t binary_volume_index(uint32_t x, uint32_t y, uint32_t z,
                                    uint32_t width, uint32_t height) {
    return (uint64_t)z * (uint64_t)width * (uint64_t)height +
           (uint64_t)y * (uint64_t)width + (uint64_t)x;
}

static int read_stack_shape(TIFF *tif, uint32_t *width, uint32_t *height,
                            uint32_t *depth) {
    uint32_t current_width = 0;
    uint32_t current_height = 0;
    uint16_t bits_per_sample = 0;
    uint16_t samples_per_pixel = 0;
    uint32_t page_count = 0;

    if (tif == NULL || width == NULL || height == NULL || depth == NULL) {
        return EXIT_FAILURE;
    }

    do {
        if (!TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &current_width) ||
            !TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &current_height) ||
            !TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample)) {
            return EXIT_FAILURE;
        }

        TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

        if (bits_per_sample != 8 || samples_per_pixel != 1) {
            return EXIT_FAILURE;
        }

        if (page_count == 0) {
            *width = current_width;
            *height = current_height;
        } else if (*width != current_width || *height != current_height) {
            return EXIT_FAILURE;
        }

        page_count += 1;
    } while (TIFFReadDirectory(tif));

    *depth = page_count;
    return page_count > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

static void accumulate_porosity(TIFF *tif, StackStats *stats) {
    do {
        uint32_t width = 0;
        uint32_t height = 0;
        uint16_t bits_per_sample = 0;
        uint16_t samples_per_pixel = 0;
        uint16_t photometric = 0;
        uint16_t *red = NULL;
        uint16_t *green = NULL;
        uint16_t *blue = NULL;

        if (!TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) ||
            !TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height) ||
            !TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample)) {
            tiff_die("Unable to read required TIFF metadata", NULL);
        }

        TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
        TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &photometric);

        if (bits_per_sample != 8 || samples_per_pixel != 1) {
            tiff_die("This program currently expects 8-bit single-channel TIFF pages",
                     NULL);
        }

        if (stats->slices == 0) {
            stats->width = width;
            stats->height = height;
        } else if (stats->width != width || stats->height != height) {
            tiff_die("All slices in a stack must have matching dimensions", NULL);
        }

        if (photometric == PHOTOMETRIC_PALETTE) {
            if (!TIFFGetField(tif, TIFFTAG_COLORMAP, &red, &green, &blue)) {
                tiff_die("Palette TIFF is missing its color map", NULL);
            }
            stats->has_palette = 1;
        }

        tsize_t scanline_size = TIFFScanlineSize(tif);
        uint8_t *scanline = (uint8_t *)_TIFFmalloc((tmsize_t)scanline_size);
        if (scanline == NULL) {
            tiff_die("Unable to allocate memory for scanline", NULL);
        }

        for (uint32_t row = 0; row < height; ++row) {
            if (TIFFReadScanline(tif, scanline, row, 0) < 0) {
                _TIFFfree(scanline);
                tiff_die("Failed while reading scanline", NULL);
            }

            for (uint32_t col = 0; col < width; ++col) {
                uint8_t is_void =
                    palette_index_to_binary(scanline[col], red, green, blue);
                stats->total_voxels += 1;
                if (is_void == 0) { /*Edited by Savio Poovathingal manually*/ 
                    stats->void_voxels += 1;
                }
            }
        }

        _TIFFfree(scanline);
        stats->slices += 1;
    } while (TIFFReadDirectory(tif));
}

StackStats analyze_stack(const char *path) {
    StackStats stats;

    memset(&stats, 0, sizeof(stats));

    TIFF *tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        tiff_die("Unable to open TIFF stack", path);
    }

    accumulate_porosity(tif, &stats);
    TIFFClose(tif);
    return stats;
}

double compute_porosity(StackStats stats) {
    if (stats.total_voxels == 0) {
        return 0.0;
    }

    return (double)stats.void_voxels / (double)stats.total_voxels;
}

int write_stack_summary(FILE *out, const char *path, StackStats stats) {
    double porosity = compute_porosity(stats);

    if (fprintf(out, "Stack: %s\n", path) < 0 ||
        fprintf(out, "  Dimensions: %u x %u x %u\n", stats.width, stats.height,
                stats.slices) < 0 ||
        fprintf(out, "  Total voxels: %llu\n",
                (unsigned long long)stats.total_voxels) < 0 ||
        fprintf(out, "  Void voxels: %llu\n",
                (unsigned long long)stats.void_voxels) < 0 ||
        fprintf(out, "  Porosity: %.6f (%.2f%%)\n\n", porosity,
                porosity * 100.0) < 0) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int load_binary_volume_from_tiff(const char *path, BinaryVolume *volume) {
    TIFF *tif = NULL;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t depth = 0;
    uint32_t z = 0;

    if (path == NULL || volume == NULL) {
        return EXIT_FAILURE;
    }

    memset(volume, 0, sizeof(*volume));

    tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        return EXIT_FAILURE;
    }

    if (read_stack_shape(tif, &width, &height, &depth) != EXIT_SUCCESS) {
        TIFFClose(tif);
        return EXIT_FAILURE;
    }

    TIFFClose(tif);
    tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        return EXIT_FAILURE;
    }

    volume->width = width;
    volume->height = height;
    volume->depth = depth;
    volume->voxels =
        (uint8_t *)malloc((size_t)width * (size_t)height * (size_t)depth);
    if (volume->voxels == NULL) {
        TIFFClose(tif);
        memset(volume, 0, sizeof(*volume));
        return EXIT_FAILURE;
    }

    do {
        uint16_t photometric = 0;
        uint16_t *red = NULL;
        uint16_t *green = NULL;
        uint16_t *blue = NULL;
        tsize_t scanline_size = TIFFScanlineSize(tif);
        uint8_t *scanline = (uint8_t *)_TIFFmalloc((tmsize_t)scanline_size);

        if (scanline == NULL) {
            TIFFClose(tif);
            free_binary_volume(volume);
            return EXIT_FAILURE;
        }

        TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &photometric);
        if (photometric == PHOTOMETRIC_PALETTE) {
            if (!TIFFGetField(tif, TIFFTAG_COLORMAP, &red, &green, &blue)) {
                _TIFFfree(scanline);
                TIFFClose(tif);
                free_binary_volume(volume);
                return EXIT_FAILURE;
            }
        }

        for (uint32_t row = 0; row < height; ++row) {
            if (TIFFReadScanline(tif, scanline, row, 0) < 0) {
                _TIFFfree(scanline);
                TIFFClose(tif);
                free_binary_volume(volume);
                return EXIT_FAILURE;
            }

            for (uint32_t col = 0; col < width; ++col) {
                volume->voxels[binary_volume_index(col, row, z, width, height)] =
                    palette_index_to_binary(scanline[col], red, green, blue);
            }
        }

        _TIFFfree(scanline);
        z += 1;
    } while (TIFFReadDirectory(tif));

    TIFFClose(tif);
    return EXIT_SUCCESS;
}

void free_binary_volume(BinaryVolume *volume) {
    if (volume == NULL) {
        return;
    }

    free(volume->voxels);
    volume->voxels = NULL;
    volume->width = 0;
    volume->height = 0;
    volume->depth = 0;
}
