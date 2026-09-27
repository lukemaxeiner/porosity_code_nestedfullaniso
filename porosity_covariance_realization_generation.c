#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "porosity_covariance_realization_generation.h"
#include "porosity_normal_score.h"
#define MAX_CHOLESKY_POINTS 4096U
#define MAX_FFT_EMBEDDED_POINTS 16777216U

typedef struct {
    double real;
    double imag;
} ComplexValue;

typedef struct {
    double mean;
    double variance;
    double min_value;
    double max_value;
    uint64_t count_below_zero;
    uint64_t count_above_one;
    uint64_t total_count;
    double fraction_outside_unit_interval;
} PorosityGenerationDiagnostics;

static int compute_grid_spacing(const PorosityRealizationGridSpec *simulation_grid,
                                double *dx, double *dy, double *dz);

static void free_porosity_field_storage(PorosityField *field) {
    if (field == NULL) {
        return;
    }

    free(field->values);
    field->values = NULL;
    field->width = 0U;
    field->height = 0U;
    field->depth = 0U;
}

static int ensure_generated_porosity_fields_dir(const char *stack_path) {
    const char *slash;
    size_t dir_length;
    char *folder_path;
    int status;

    if (stack_path == NULL) {
        return EXIT_FAILURE;
    }

    slash = strrchr(stack_path, '/');
    dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    folder_path =
        (char *)malloc(dir_length + strlen("/GENERATED_POROSITY_FIELDS") + 1U);
    if (folder_path == NULL) {
        return EXIT_FAILURE;
    }

    if (slash != NULL) {
        memcpy(folder_path, stack_path, dir_length);
    } else {
        folder_path[0] = '.';
    }
    memcpy(folder_path + dir_length, "/GENERATED_POROSITY_FIELDS",
           strlen("/GENERATED_POROSITY_FIELDS") + 1U);

    status = mkdir(folder_path, 0777);
    if (status != 0 && errno != EEXIST) {
        free(folder_path);
        return EXIT_FAILURE;
    }

    free(folder_path);
    return EXIT_SUCCESS;
}

static char *build_realization_output_path(const char *stack_path,
                                           uint32_t realization_index,
                                           const char *output_extension) {
    const char *input_extension = strrchr(stack_path, '.');
    const char *slash = strrchr(stack_path, '/');
    const char *file_name = slash != NULL ? slash + 1 : stack_path;
    size_t dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    size_t file_base_length =
        input_extension != NULL && input_extension > file_name
            ? (size_t)(input_extension - file_name)
            : strlen(file_name);
    char suffix[128];
    size_t total_length;
    char *output_path;

    if (stack_path == NULL || output_extension == NULL) {
        return NULL;
    }

    if (snprintf(suffix, sizeof(suffix), "_gaussian_realization_%03u.%s",
                 realization_index + 1U, output_extension) >=
        (int)sizeof(suffix)) {
        return NULL;
    }

    total_length = dir_length + strlen("/GENERATED_POROSITY_FIELDS/") +
                   file_base_length + strlen(suffix) + 1U;
    output_path = (char *)malloc(total_length);
    if (output_path == NULL) {
        return NULL;
    }

    if (slash != NULL) {
        memcpy(output_path, stack_path, dir_length);
    } else {
        output_path[0] = '.';
    }
    memcpy(output_path + dir_length, "/GENERATED_POROSITY_FIELDS/",
           strlen("/GENERATED_POROSITY_FIELDS/"));
    memcpy(output_path + dir_length + strlen("/GENERATED_POROSITY_FIELDS/"),
           file_name, file_base_length);
    memcpy(output_path + dir_length +
               strlen("/GENERATED_POROSITY_FIELDS/") + file_base_length,
           suffix, strlen(suffix) + 1U);

    return output_path;
}

static char *build_generation_diagnostics_output_path(const char *stack_path) {
    const char *input_extension = strrchr(stack_path, '.');
    const char *slash = strrchr(stack_path, '/');
    const char *file_name = slash != NULL ? slash + 1 : stack_path;
    size_t dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    size_t file_base_length =
        input_extension != NULL && input_extension > file_name
            ? (size_t)(input_extension - file_name)
            : strlen(file_name);
    const char *suffix = "_generation_diagnostics.txt";
    size_t total_length;
    char *output_path;

    if (stack_path == NULL) {
        return NULL;
    }

    total_length = dir_length + strlen("/GENERATED_POROSITY_FIELDS/") +
                   file_base_length + strlen(suffix) + 1U;

    output_path = (char *)malloc(total_length);
    if (output_path == NULL) {
        return NULL;
    }

    if (slash != NULL) {
        memcpy(output_path, stack_path, dir_length);
    } else {
        output_path[0] = '.';
    }

    memcpy(output_path + dir_length, "/GENERATED_POROSITY_FIELDS/",
           strlen("/GENERATED_POROSITY_FIELDS/"));

    memcpy(output_path + dir_length + strlen("/GENERATED_POROSITY_FIELDS/"),
           file_name, file_base_length);

    memcpy(output_path + dir_length +
               strlen("/GENERATED_POROSITY_FIELDS/") + file_base_length,
           suffix, strlen(suffix) + 1U);

    return output_path;
}

static int write_porosity_field_vtk(FILE *out, const PorosityField *field,
                                    const PorosityRealizationGridSpec *grid) {
    size_t point_count;
    size_t i;
    double dx;
    double dy;
    double dz;

    if (out == NULL || field == NULL || grid == NULL || field->values == NULL) {
        return EXIT_FAILURE;
    }
    if (field->width != grid->nx || field->height != grid->ny ||
        field->depth != grid->nz) {
        return EXIT_FAILURE;
    }
    if (compute_grid_spacing(grid, &dx, &dy, &dz) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    point_count = (size_t)field->width * (size_t)field->height *
                  (size_t)field->depth;
    if (fprintf(out,
                "# vtk DataFile Version 3.0\n"
                "Gaussian porosity realization\n"
                "ASCII\n"
                "DATASET STRUCTURED_POINTS\n"
                "DIMENSIONS %u %u %u\n"
                "ORIGIN 0 0 0\n"
                "SPACING %.17g %.17g %.17g\n"
                "POINT_DATA %zu\n"
                "SCALARS porosity double 1\n"
                "LOOKUP_TABLE default\n",
                field->width, field->height, field->depth, dx, dy, dz,
                point_count) < 0) {
        return EXIT_FAILURE;
    }

    for (i = 0; i < point_count; ++i) {
        if (fprintf(out, "%.17g\n", field->values[i]) < 0) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static int compute_porosity_generation_diagnostics(
    const PorosityField *field,
    PorosityGenerationDiagnostics *diagnostics_out) {
    size_t point_count;
    size_t i;
    double sum = 0.0;
    double sum_sq = 0.0;

    if (field == NULL || field->values == NULL || diagnostics_out == NULL) {
        return EXIT_FAILURE;
    }

    point_count = (size_t)field->width * (size_t)field->height *
                  (size_t)field->depth;

    if (point_count == 0U) {
        return EXIT_FAILURE;
    }

    memset(diagnostics_out, 0, sizeof(*diagnostics_out));

    diagnostics_out->min_value = field->values[0];
    diagnostics_out->max_value = field->values[0];
    diagnostics_out->total_count = (uint64_t)point_count;

    for (i = 0U; i < point_count; ++i) {
        double value = field->values[i];

        sum += value;
        sum_sq += value * value;

        if (value < diagnostics_out->min_value) {
            diagnostics_out->min_value = value;
        }
        if (value > diagnostics_out->max_value) {
            diagnostics_out->max_value = value;
        }

        if (value < 0.0) {
            diagnostics_out->count_below_zero += 1U;
        }
        if (value > 1.0) {
            diagnostics_out->count_above_one += 1U;
        }
    }

    diagnostics_out->mean = sum / (double)point_count;
    diagnostics_out->variance =
        sum_sq / (double)point_count -
        diagnostics_out->mean * diagnostics_out->mean;

    if (diagnostics_out->variance < 0.0 &&
        diagnostics_out->variance > -1.0e-15) {
        diagnostics_out->variance = 0.0;
    }

    diagnostics_out->fraction_outside_unit_interval =
        (double)(diagnostics_out->count_below_zero +
                 diagnostics_out->count_above_one) /
        (double)diagnostics_out->total_count;

    return EXIT_SUCCESS;
}

static int write_generation_diagnostics_header(FILE *out) {
    if (out == NULL) {
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "realization mean variance min_value max_value "
                "count_below_zero count_above_one total_count "
                "fraction_outside_0_1\n") < 0) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static int write_generation_diagnostics_row(
    FILE *out,
    uint32_t realization_index,
    const PorosityGenerationDiagnostics *diagnostics) {
    if (out == NULL || diagnostics == NULL) {
        return EXIT_FAILURE;
    }

    if (fprintf(out,
                "%u %.12f %.12f %.12f %.12f %llu %llu %llu %.12f\n",
                realization_index + 1U,
                diagnostics->mean,
                diagnostics->variance,
                diagnostics->min_value,
                diagnostics->max_value,
                (unsigned long long)diagnostics->count_below_zero,
                (unsigned long long)diagnostics->count_above_one,
                (unsigned long long)diagnostics->total_count,
                diagnostics->fraction_outside_unit_interval) < 0) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static uint32_t next_random_u32(unsigned int *seed) {
    uint32_t state = (uint32_t)*seed;

    state = state * UINT32_C(1664525) + UINT32_C(1013904223);
    *seed = (unsigned int)state;

    return state;
}

static double uniform_open01(unsigned int *seed) {
    return ((double)next_random_u32(seed) + 1.0) /
           ((double)UINT32_MAX + 2.0);
}

static double sample_standard_normal(unsigned int *seed) {
    const double two_pi = 6.283185307179586476925286766559;
    double u1 = uniform_open01(seed);
    double u2 = uniform_open01(seed);

    return sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
}

static int compute_grid_spacing(const PorosityRealizationGridSpec *simulation_grid,
                                double *dx, double *dy, double *dz) {
    if (simulation_grid == NULL || dx == NULL || dy == NULL || dz == NULL) {
        return EXIT_FAILURE;
    }

    *dx = simulation_grid->nx > 1U
              ? simulation_grid->dim_x / (double)(simulation_grid->nx - 1U)
              : simulation_grid->dim_x;
    *dy = simulation_grid->ny > 1U
              ? simulation_grid->dim_y / (double)(simulation_grid->ny - 1U)
              : simulation_grid->dim_y;
    *dz = simulation_grid->nz > 1U
              ? simulation_grid->dim_z / (double)(simulation_grid->nz - 1U)
              : simulation_grid->dim_z;
    return EXIT_SUCCESS;
}

static void linear_index_to_grid(size_t linear_index, uint32_t width,
                                 uint32_t height, uint32_t *x, uint32_t *y,
                                 uint32_t *z) {
    size_t plane = (size_t)width * (size_t)height;
    *z = (uint32_t)(linear_index / plane);
    linear_index %= plane;
    *y = (uint32_t)(linear_index / (size_t)width);
    *x = (uint32_t)(linear_index % (size_t)width);
}

static double evaluate_anisotropic_exponential_component(
    double rx, double ry, double rz,
    double lx, double ly, double lz) {
    double scaled_radius;

    if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0) {
        return 0.0;
    }

    scaled_radius = sqrt((rx * rx) / (lx * lx) +
                         (ry * ry) / (ly * ly) +
                         (rz * rz) / (lz * lz));

    return exp(-scaled_radius);
}

static int validate_generation_covariance_model(
    const ExponentialCovarianceModel *model) {
    uint32_t count;
    uint32_t k;
    double weight_sum = 0.0;

    if (model == NULL || model->sill <= 0.0) {
        return EXIT_FAILURE;
    }

    count = model->component_count;

    if (count == 0U) {
        return model->x_correlation_length > 0.0 &&
                       model->y_correlation_length > 0.0 &&
                       model->z_correlation_length > 0.0
                   ? EXIT_SUCCESS
                   : EXIT_FAILURE;
    }

    if (count > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
        return EXIT_FAILURE;
    }

    for (k = 0U; k < count; ++k) {
        double weight = model->component_weight[k];

        if (!isfinite(weight) || weight < 0.0) {
            return EXIT_FAILURE;
        }

        if (weight > 0.0) {
            if (model->component_x_correlation_length[k] <= 0.0 ||
                model->component_y_correlation_length[k] <= 0.0 ||
                model->component_z_correlation_length[k] <= 0.0 ||
                !isfinite(model->component_x_correlation_length[k]) ||
                !isfinite(model->component_y_correlation_length[k]) ||
                !isfinite(model->component_z_correlation_length[k])) {
                return EXIT_FAILURE;
            }
        }

        weight_sum += weight;
    }

    return fabs(weight_sum - 1.0) <= 1.0e-6 ? EXIT_SUCCESS : EXIT_FAILURE;
}

static double evaluate_generation_covariance_model_at_lag(
    const ExponentialCovarianceModel *model,
    double rx, double ry, double rz) {
    double normalized = 0.0;
    uint32_t count;
    uint32_t k;

    if (model == NULL || model->sill <= 0.0) {
        return 0.0;
    }

    count = model->component_count;

    if (count == 0U) {
        return model->sill *
               evaluate_anisotropic_exponential_component(
                   rx, ry, rz,
                   model->x_correlation_length,
                   model->y_correlation_length,
                   model->z_correlation_length);
    }

    if (count > POROSITY_NESTED_EXPONENTIAL_COMPONENTS) {
        count = POROSITY_NESTED_EXPONENTIAL_COMPONENTS;
    }

    for (k = 0U; k < count; ++k) {
        double weight = model->component_weight[k];

        if (weight <= 0.0) {
            continue;
        }

        normalized +=
            weight *
            evaluate_anisotropic_exponential_component(
                rx, ry, rz,
                model->component_x_correlation_length[k],
                model->component_y_correlation_length[k],
                model->component_z_correlation_length[k]);
    }

    return model->sill * normalized;
}

static int checked_multiply_size(size_t a, size_t b, size_t *out) {
    if (out == NULL) {
        return EXIT_FAILURE;
    }

    if (a != 0U && b > SIZE_MAX / a) {
        return EXIT_FAILURE;
    }

    *out = a * b;
    return EXIT_SUCCESS;
}

static int is_power_of_two_size(size_t value) {
    return value != 0U && (value & (value - 1U)) == 0U;
}

static size_t next_power_of_two_size(size_t value) {
    size_t result = 1U;

    while (result < value) {
        if (result > SIZE_MAX / 2U) {
            return 0U;
        }
        result *= 2U;
    }

    return result;
}

static size_t fft_embedding_size_for_grid_count(uint32_t count) {
    size_t needed;

    if (count <= 1U) {
        return 1U;
    }

    needed = 2U * (size_t)(count - 1U);
    return next_power_of_two_size(needed);
}

static size_t embedded_index(size_t x, size_t y, size_t z,
                             size_t nx, size_t ny) {
    return x + nx * (y + ny * z);
}

static void swap_complex_values(ComplexValue *a, ComplexValue *b) {
    ComplexValue tmp = *a;
    *a = *b;
    *b = tmp;
}

static int fft_1d_in_place(ComplexValue *data, size_t n, int inverse) {
    const double two_pi = 6.283185307179586476925286766559;
    size_t i;
    size_t j = 0U;
    size_t length;

    if (data == NULL || !is_power_of_two_size(n)) {
        return EXIT_FAILURE;
    }

    for (i = 1U; i < n; ++i) {
        size_t bit = n >> 1U;

        while ((j & bit) != 0U) {
            j ^= bit;
            bit >>= 1U;
        }
        j ^= bit;

        if (i < j) {
            swap_complex_values(&data[i], &data[j]);
        }
    }

    for (length = 2U; length <= n; length <<= 1U) {
        double angle = (inverse ? two_pi : -two_pi) / (double)length;
        ComplexValue w_length;
        size_t half_length = length >> 1U;
        size_t start;

        w_length.real = cos(angle);
        w_length.imag = sin(angle);

        for (start = 0U; start < n; start += length) {
            ComplexValue w;
            size_t offset;

            w.real = 1.0;
            w.imag = 0.0;

            for (offset = 0U; offset < half_length; ++offset) {
                ComplexValue u = data[start + offset];
                ComplexValue v;
                ComplexValue next_w;

                v.real = data[start + offset + half_length].real * w.real -
                         data[start + offset + half_length].imag * w.imag;
                v.imag = data[start + offset + half_length].real * w.imag +
                         data[start + offset + half_length].imag * w.real;

                data[start + offset].real = u.real + v.real;
                data[start + offset].imag = u.imag + v.imag;
                data[start + offset + half_length].real = u.real - v.real;
                data[start + offset + half_length].imag = u.imag - v.imag;

                next_w.real = w.real * w_length.real - w.imag * w_length.imag;
                next_w.imag = w.real * w_length.imag + w.imag * w_length.real;
                w = next_w;
            }
        }

        if (length == n) {
            break;
        }
    }

    if (inverse) {
        double scale = 1.0 / (double)n;
        for (i = 0U; i < n; ++i) {
            data[i].real *= scale;
            data[i].imag *= scale;
        }
    }

    return EXIT_SUCCESS;
}

static int fft_3d_in_place(ComplexValue *data,
                           size_t nx, size_t ny, size_t nz,
                           int inverse) {
    size_t x;
    size_t y;
    size_t z;
    size_t max_dim = nx;
    ComplexValue *line = NULL;

    if (data == NULL || !is_power_of_two_size(nx) ||
        !is_power_of_two_size(ny) || !is_power_of_two_size(nz)) {
        return EXIT_FAILURE;
    }

    if (ny > max_dim) {
        max_dim = ny;
    }
    if (nz > max_dim) {
        max_dim = nz;
    }

    line = (ComplexValue *)malloc(max_dim * sizeof(ComplexValue));
    if (line == NULL) {
        return EXIT_FAILURE;
    }

    for (z = 0U; z < nz; ++z) {
        for (y = 0U; y < ny; ++y) {
            if (fft_1d_in_place(&data[embedded_index(0U, y, z, nx, ny)],
                                nx, inverse) != EXIT_SUCCESS) {
                free(line);
                return EXIT_FAILURE;
            }
        }
    }

    for (z = 0U; z < nz; ++z) {
        for (x = 0U; x < nx; ++x) {
            for (y = 0U; y < ny; ++y) {
                line[y] = data[embedded_index(x, y, z, nx, ny)];
            }
            if (fft_1d_in_place(line, ny, inverse) != EXIT_SUCCESS) {
                free(line);
                return EXIT_FAILURE;
            }
            for (y = 0U; y < ny; ++y) {
                data[embedded_index(x, y, z, nx, ny)] = line[y];
            }
        }
    }

    for (y = 0U; y < ny; ++y) {
        for (x = 0U; x < nx; ++x) {
            for (z = 0U; z < nz; ++z) {
                line[z] = data[embedded_index(x, y, z, nx, ny)];
            }
            if (fft_1d_in_place(line, nz, inverse) != EXIT_SUCCESS) {
                free(line);
                return EXIT_FAILURE;
            }
            for (z = 0U; z < nz; ++z) {
                data[embedded_index(x, y, z, nx, ny)] = line[z];
            }
        }
    }

    free(line);
    return EXIT_SUCCESS;
}

static int fill_circulant_covariance_embedding(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid,
    ComplexValue *embedding,
    size_t embed_nx,
    size_t embed_ny,
    size_t embed_nz) {
    double dx;
    double dy;
    double dz;
    size_t x;
    size_t y;
    size_t z;

    if (model == NULL || simulation_grid == NULL || embedding == NULL) {
        return EXIT_FAILURE;
    }

    if (compute_grid_spacing(simulation_grid, &dx, &dy, &dz) !=
        EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    for (z = 0U; z < embed_nz; ++z) {
        size_t z_lag_index = z <= embed_nz - z ? z : embed_nz - z;
        double rz = (double)z_lag_index * dz;

        for (y = 0U; y < embed_ny; ++y) {
            size_t y_lag_index = y <= embed_ny - y ? y : embed_ny - y;
            double ry = (double)y_lag_index * dy;

            for (x = 0U; x < embed_nx; ++x) {
                size_t x_lag_index = x <= embed_nx - x ? x : embed_nx - x;
                double rx = (double)x_lag_index * dx;
                size_t index = embedded_index(x, y, z, embed_nx, embed_ny);

                embedding[index].real =
                    evaluate_generation_covariance_model_at_lag(
                        model, rx, ry, rz);
                embedding[index].imag = 0.0;
            }
        }
    }

    return EXIT_SUCCESS;
}


static double covariance_between_points(
    size_t i, size_t j, uint32_t width, uint32_t height, double dx, double dy,
    double dz, const ExponentialCovarianceModel *model) {
    uint32_t xi;
    uint32_t yi;
    uint32_t zi;
    uint32_t xj;
    uint32_t yj;
    uint32_t zj;
    double rx;
    double ry;
    double rz;

    linear_index_to_grid(i, width, height, &xi, &yi, &zi);
    linear_index_to_grid(j, width, height, &xj, &yj, &zj);

    rx = (double)((int)xj - (int)xi) * dx;
    ry = (double)((int)yj - (int)yi) * dy;
    rz = (double)((int)zj - (int)zi) * dz;

    return evaluate_generation_covariance_model_at_lag(model, rx, ry, rz);
}

static int fill_covariance_matrix(const ExponentialCovarianceModel *model,
                                  const PorosityRealizationGridSpec *simulation_grid,
                                  double *matrix, size_t point_count) {
    size_t i;
    size_t j;
    double dx;
    double dy;
    double dz;

    if (compute_grid_spacing(simulation_grid, &dx, &dy, &dz) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    for (i = 0; i < point_count; ++i) {
        for (j = i; j < point_count; ++j) {
            double value = covariance_between_points(
                i, j, simulation_grid->nx, simulation_grid->ny, dx, dy, dz,
                model);
            matrix[i * point_count + j] = value;
            matrix[j * point_count + i] = value;
        }
    }

    return EXIT_SUCCESS;
}

static int cholesky_decompose_in_place(double *matrix, size_t n) {
    size_t i;
    size_t j;
    size_t k;

    for (i = 0; i < n; ++i) {
        for (j = 0; j <= i; ++j) {
            double sum = matrix[i * n + j];

            for (k = 0; k < j; ++k) {
                sum -= matrix[i * n + k] * matrix[j * n + k];
            }

            if (i == j) {
                if (sum <= 0.0) {
                    return EXIT_FAILURE;
                }
                matrix[i * n + j] = sqrt(sum);
            } else {
                matrix[i * n + j] = sum / matrix[j * n + j];
            }
        }

        for (j = i + 1; j < n; ++j) {
            matrix[i * n + j] = 0.0;
        }
    }

    return EXIT_SUCCESS;
}

static void multiply_lower_by_vector(const double *lower, const double *vector,
                                     size_t n, double *result) {
    size_t i;
    size_t j;

    for (i = 0; i < n; ++i) {
        double sum = 0.0;
        for (j = 0; j <= i; ++j) {
            sum += lower[i * n + j] * vector[j];
        }
        result[i] = sum;
    }
}

int generate_gaussian_porosity_realization_from_model(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid, double mean_porosity,
    unsigned int seed, PorosityField *field_out) {
    size_t point_count;
    double *covariance_matrix = NULL;
    double *standard_normal = NULL;
    double *correlated = NULL;
    unsigned int rng_state = seed;
    size_t i;

    if (model == NULL || simulation_grid == NULL || field_out == NULL) {
        return EXIT_FAILURE;
    }
    if (validate_generation_covariance_model(model) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }
    
    if (simulation_grid->nx == 0U || simulation_grid->ny == 0U ||
        simulation_grid->nz == 0U) {
        return EXIT_FAILURE;
    }

    point_count = (size_t)simulation_grid->nx * (size_t)simulation_grid->ny *
                  (size_t)simulation_grid->nz;
    if (point_count == 0U || point_count > MAX_CHOLESKY_POINTS) {
        return EXIT_FAILURE;
    }

    memset(field_out, 0, sizeof(*field_out));
    field_out->values = (double *)malloc(point_count * sizeof(double));
    covariance_matrix = (double *)malloc(point_count * point_count *
                                         sizeof(double));
    standard_normal = (double *)malloc(point_count * sizeof(double));
    correlated = (double *)malloc(point_count * sizeof(double));
    if (field_out->values == NULL || covariance_matrix == NULL ||
        standard_normal == NULL || correlated == NULL) {
        free_porosity_field_storage(field_out);
        free(covariance_matrix);
        free(standard_normal);
        free(correlated);
        return EXIT_FAILURE;
    }

    if (fill_covariance_matrix(model, simulation_grid, covariance_matrix,
                               point_count) != EXIT_SUCCESS) {
        free_porosity_field_storage(field_out);
        free(covariance_matrix);
        free(standard_normal);
        free(correlated);
        return EXIT_FAILURE;
    }

    if (cholesky_decompose_in_place(covariance_matrix, point_count) !=
        EXIT_SUCCESS) {
        free_porosity_field_storage(field_out);
        free(covariance_matrix);
        free(standard_normal);
        free(correlated);
        return EXIT_FAILURE;
    }

    for (i = 0; i < point_count; ++i) {
        standard_normal[i] = sample_standard_normal(&rng_state);
    }

    multiply_lower_by_vector(covariance_matrix, standard_normal, point_count,
                             correlated);

    for (i = 0; i < point_count; ++i) {
        field_out->values[i] = mean_porosity + correlated[i];
    }

    field_out->width = simulation_grid->nx;
    field_out->height = simulation_grid->ny;
    field_out->depth = simulation_grid->nz;

    free(covariance_matrix);
    free(standard_normal);
    free(correlated);
    return EXIT_SUCCESS;
}

int generate_gaussian_porosity_realization_from_model_fft(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid,
    double mean_porosity,
    unsigned int seed,
    PorosityField *field_out) {
    size_t point_xy;
    size_t point_count;
    size_t embed_nx;
    size_t embed_ny;
    size_t embed_nz;
    size_t embed_xy;
    size_t embedded_count;
    ComplexValue *covariance_spectrum = NULL;
    ComplexValue *white_noise = NULL;
    unsigned int rng_state = seed;
    double max_eigenvalue = 0.0;
    double min_eigenvalue = 0.0;
    double negative_tolerance;
    size_t negative_count = 0U;
    size_t significant_negative_count = 0U;
    size_t i;
    size_t x;
    size_t y;
    size_t z;

    if (model == NULL || simulation_grid == NULL || field_out == NULL) {
        return EXIT_FAILURE;
    }
    if (validate_generation_covariance_model(model) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }
    if (simulation_grid->nx == 0U || simulation_grid->ny == 0U ||
        simulation_grid->nz == 0U) {
        return EXIT_FAILURE;
    }

    if (checked_multiply_size((size_t)simulation_grid->nx,
                              (size_t)simulation_grid->ny,
                              &point_xy) != EXIT_SUCCESS ||
        checked_multiply_size(point_xy,
                              (size_t)simulation_grid->nz,
                              &point_count) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    embed_nx = fft_embedding_size_for_grid_count(simulation_grid->nx);
    embed_ny = fft_embedding_size_for_grid_count(simulation_grid->ny);
    embed_nz = fft_embedding_size_for_grid_count(simulation_grid->nz);
    if (embed_nx == 0U || embed_ny == 0U || embed_nz == 0U) {
        return EXIT_FAILURE;
    }

    if (checked_multiply_size(embed_nx, embed_ny, &embed_xy) !=
            EXIT_SUCCESS ||
        checked_multiply_size(embed_xy, embed_nz, &embedded_count) !=
            EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    if (embedded_count == 0U || embedded_count > MAX_FFT_EMBEDDED_POINTS) {
        fprintf(stderr,
                "FFT embedded grid is too large: %zu points. Limit is %u.\n",
                embedded_count, MAX_FFT_EMBEDDED_POINTS);
        return EXIT_FAILURE;
    }

    memset(field_out, 0, sizeof(*field_out));
    field_out->values = (double *)malloc(point_count * sizeof(double));
    covariance_spectrum =
        (ComplexValue *)calloc(embedded_count, sizeof(ComplexValue));
    white_noise =
        (ComplexValue *)calloc(embedded_count, sizeof(ComplexValue));

    if (field_out->values == NULL || covariance_spectrum == NULL ||
        white_noise == NULL) {
        free_porosity_field_storage(field_out);
        free(covariance_spectrum);
        free(white_noise);
        return EXIT_FAILURE;
    }

    if (fill_circulant_covariance_embedding(
            model, simulation_grid, covariance_spectrum,
            embed_nx, embed_ny, embed_nz) != EXIT_SUCCESS) {
        free_porosity_field_storage(field_out);
        free(covariance_spectrum);
        free(white_noise);
        return EXIT_FAILURE;
    }

    if (fft_3d_in_place(covariance_spectrum,
                        embed_nx, embed_ny, embed_nz,
                        0) != EXIT_SUCCESS) {
        free_porosity_field_storage(field_out);
        free(covariance_spectrum);
        free(white_noise);
        return EXIT_FAILURE;
    }

    min_eigenvalue = covariance_spectrum[0].real;
    for (i = 0U; i < embedded_count; ++i) {
        double lambda = covariance_spectrum[i].real;

        if (lambda > max_eigenvalue) {
            max_eigenvalue = lambda;
        }
        if (lambda < min_eigenvalue) {
            min_eigenvalue = lambda;
        }
    }

    if (max_eigenvalue <= 0.0 || !isfinite(max_eigenvalue)) {
        free_porosity_field_storage(field_out);
        free(covariance_spectrum);
        free(white_noise);
        return EXIT_FAILURE;
    }

    negative_tolerance = 1.0e-8 * max_eigenvalue;
    for (i = 0U; i < embedded_count; ++i) {
        double lambda = covariance_spectrum[i].real;

        if (lambda < 0.0) {
            negative_count += 1U;
            if (lambda < -negative_tolerance) {
                significant_negative_count += 1U;
            }
            lambda = 0.0;
        }

        covariance_spectrum[i].real = lambda;
        covariance_spectrum[i].imag = 0.0;
    }

    if (significant_negative_count > 0U) {
        fprintf(stderr,
                "Warning: FFT circulant embedding clipped %zu negative "
                "eigenvalues (%zu below tolerance). Minimum eigenvalue was %.17g.\n",
                negative_count, significant_negative_count, min_eigenvalue);
    }

    for (i = 0U; i < embedded_count; ++i) {
        white_noise[i].real = sample_standard_normal(&rng_state);
        white_noise[i].imag = 0.0;
    }

    if (fft_3d_in_place(white_noise,
                        embed_nx, embed_ny, embed_nz,
                        0) != EXIT_SUCCESS) {
        free_porosity_field_storage(field_out);
        free(covariance_spectrum);
        free(white_noise);
        return EXIT_FAILURE;
    }

    for (i = 0U; i < embedded_count; ++i) {
        double scale = sqrt(covariance_spectrum[i].real);

        white_noise[i].real *= scale;
        white_noise[i].imag *= scale;
    }

    if (fft_3d_in_place(white_noise,
                        embed_nx, embed_ny, embed_nz,
                        1) != EXIT_SUCCESS) {
        free_porosity_field_storage(field_out);
        free(covariance_spectrum);
        free(white_noise);
        return EXIT_FAILURE;
    }

    for (z = 0U; z < (size_t)simulation_grid->nz; ++z) {
        for (y = 0U; y < (size_t)simulation_grid->ny; ++y) {
            for (x = 0U; x < (size_t)simulation_grid->nx; ++x) {
                size_t output_index =
                    x + (size_t)simulation_grid->nx *
                            (y + (size_t)simulation_grid->ny * z);
                size_t source_index =
                    embedded_index(x, y, z, embed_nx, embed_ny);

                field_out->values[output_index] =
                    mean_porosity + white_noise[source_index].real;
            }
        }
    }

    field_out->width = simulation_grid->nx;
    field_out->height = simulation_grid->ny;
    field_out->depth = simulation_grid->nz;

    free(covariance_spectrum);
    free(white_noise);
    return EXIT_SUCCESS;
}


int generate_gaussian_porosity_realization_from_prep(
    const PorosityRealizationPrepResult *prep_result, unsigned int seed,
    PorosityField *field_out) {
    if (prep_result == NULL || field_out == NULL) {
        return EXIT_FAILURE;
    }

    return generate_gaussian_porosity_realization_from_model(
        &prep_result->exponential_model, &prep_result->simulation_grid,
        prep_result->covariance_grid.mean, seed, field_out);
}

int generate_and_write_gaussian_porosity_realizations(
    const char *stack_path, const PorosityRealizationPrepResult *prep_result,
    uint32_t realization_count, unsigned int base_seed) {
    uint32_t i;
    char *diagnostics_path = NULL;
    FILE *diagnostics_file = NULL;

    if (stack_path == NULL || prep_result == NULL || realization_count == 0U) {
        return EXIT_FAILURE;
    }

    if (ensure_generated_porosity_fields_dir(stack_path) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    diagnostics_path = build_generation_diagnostics_output_path(stack_path);
    if (diagnostics_path == NULL) {
        return EXIT_FAILURE;
    }

    diagnostics_file = fopen(diagnostics_path, "w");
    if (diagnostics_file == NULL) {
        free(diagnostics_path);
        return EXIT_FAILURE;
    }

    if (write_generation_diagnostics_header(diagnostics_file) != EXIT_SUCCESS) {
        fclose(diagnostics_file);
        free(diagnostics_path);
        return EXIT_FAILURE;
    }

    for (i = 0U; i < realization_count; ++i) {
        PorosityField field;
        char *tiff_output_path;
        char *vtk_output_path;
        FILE *vtk_file;

        memset(&field, 0, sizeof(field));
        tiff_output_path = build_realization_output_path(stack_path, i, "tif");
        vtk_output_path = build_realization_output_path(stack_path, i, "vtk");
        if (tiff_output_path == NULL || vtk_output_path == NULL) {
            free(tiff_output_path);
            free(vtk_output_path);
            fclose(diagnostics_file);
            free(diagnostics_path); 
            return EXIT_FAILURE;
        }

        if (generate_gaussian_porosity_realization_from_prep(
                prep_result, base_seed + i, &field) != EXIT_SUCCESS) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        {
            PorosityGenerationDiagnostics diagnostics;

            if (compute_porosity_generation_diagnostics(
                    &field, &diagnostics) != EXIT_SUCCESS ||
                write_generation_diagnostics_row(
                    diagnostics_file, i, &diagnostics) != EXIT_SUCCESS) {
                free(tiff_output_path);
                free(vtk_output_path);
                free_porosity_field_storage(&field);
                fclose(diagnostics_file);
                free(diagnostics_path);
                return EXIT_FAILURE;
            }
        }

        if (write_porosity_field_tiff(tiff_output_path, &field) != EXIT_SUCCESS) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path); 
            return EXIT_FAILURE;
        }

        vtk_file = fopen(vtk_output_path, "w");
        if (vtk_file == NULL) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path); 
            return EXIT_FAILURE;
        }

        if (write_porosity_field_vtk(vtk_file, &field,
                                     &prep_result->simulation_grid) !=
            EXIT_SUCCESS) {
            fclose(vtk_file);
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path); 
            return EXIT_FAILURE;
        }

        fclose(vtk_file);
        free(tiff_output_path);
        free(vtk_output_path);
        free_porosity_field_storage(&field);
        fclose(diagnostics_file);
        free(diagnostics_path); 
    }

    return EXIT_SUCCESS;
}

int generate_and_write_gaussian_porosity_realizations_from_model(
    const char *stack_path, const FittedCovarianceModelArtifact *model,
    const PorosityRealizationGridSpec *simulation_grid,
    uint32_t realization_count, unsigned int base_seed) {
    uint32_t i;
    char *diagnostics_path = NULL;
    FILE *diagnostics_file = NULL;

    if (stack_path == NULL || model == NULL || simulation_grid == NULL ||
        realization_count == 0U) {
        return EXIT_FAILURE;
    }

    if (ensure_generated_porosity_fields_dir(stack_path) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    diagnostics_path = build_generation_diagnostics_output_path(stack_path);
    if (diagnostics_path == NULL) {
        return EXIT_FAILURE;
    }

    diagnostics_file = fopen(diagnostics_path, "w");
    if (diagnostics_file == NULL) {
        free(diagnostics_path);
        return EXIT_FAILURE;
    }

    if (write_generation_diagnostics_header(diagnostics_file) != EXIT_SUCCESS) {
        fclose(diagnostics_file);
        free(diagnostics_path);
        return EXIT_FAILURE;
    }

    for (i = 0U; i < realization_count; ++i) {
        PorosityField field;
        char *tiff_output_path;
        char *vtk_output_path;
        FILE *vtk_file;

        memset(&field, 0, sizeof(field));
        tiff_output_path = build_realization_output_path(stack_path, i, "tif");
        vtk_output_path = build_realization_output_path(stack_path, i, "vtk");
        if (tiff_output_path == NULL || vtk_output_path == NULL) {
            free(tiff_output_path);
            free(vtk_output_path);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        if (generate_gaussian_porosity_realization_from_model(
                &model->exponential_model, simulation_grid,
                model->mean_porosity, base_seed + i, &field) != EXIT_SUCCESS) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path); 
            return EXIT_FAILURE;
        }

        {
            PorosityGenerationDiagnostics diagnostics;

            if (compute_porosity_generation_diagnostics(
                    &field, &diagnostics) != EXIT_SUCCESS ||
                write_generation_diagnostics_row(
                    diagnostics_file, i, &diagnostics) != EXIT_SUCCESS) {
                free(tiff_output_path);
                free(vtk_output_path);
                free_porosity_field_storage(&field);
                fclose(diagnostics_file);
                free(diagnostics_path);
                return EXIT_FAILURE;
            }
        }

        if (write_porosity_field_tiff(tiff_output_path, &field) != EXIT_SUCCESS) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        vtk_file = fopen(vtk_output_path, "w");
        if (vtk_file == NULL) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        if (write_porosity_field_vtk(vtk_file, &field, simulation_grid) !=
            EXIT_SUCCESS) {
            fclose(vtk_file);
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        fclose(vtk_file);
        free(tiff_output_path);
        free(vtk_output_path);
        free_porosity_field_storage(&field);
        fclose(diagnostics_file);
        free(diagnostics_path);
    }

    return EXIT_SUCCESS;
}
/*
 * Generates, (optionally) back-transforms, and writes realizations.
 * gen_model/gen_mean describe the Gaussian field actually sampled.
 * score_table == NULL : the sampled field already is porosity.
 * score_table != NULL : the sampled field is normal scores; each value is
 *                       mapped back to porosity before diagnostics/output.
 */
static int write_gaussian_realizations_impl(
    const char *stack_path,
    const ExponentialCovarianceModel *gen_model,
    double gen_mean,
    const NormalScoreTable *score_table,
    const PorosityRealizationGridSpec *simulation_grid,
    PorosityRealizationMethod realization_method,
    uint32_t realization_count,
    unsigned int base_seed
) {
    uint32_t i;
    char *diagnostics_path = NULL;
    FILE *diagnostics_file = NULL;

    if (stack_path == NULL || gen_model == NULL || simulation_grid == NULL ||
        realization_count == 0U) {
        return EXIT_FAILURE;
    }

    if (ensure_generated_porosity_fields_dir(stack_path) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    diagnostics_path = build_generation_diagnostics_output_path(stack_path);
    if (diagnostics_path == NULL) {
        return EXIT_FAILURE;
    }

    diagnostics_file = fopen(diagnostics_path, "w");
    if (diagnostics_file == NULL) {
        free(diagnostics_path);
        return EXIT_FAILURE;
    }

    if (write_generation_diagnostics_header(diagnostics_file) != EXIT_SUCCESS) {
        fclose(diagnostics_file);
        free(diagnostics_path);
        return EXIT_FAILURE;
    }

    for (i = 0U; i < realization_count; ++i) {
        PorosityField field;
        char *tiff_output_path;
        char *vtk_output_path;
        FILE *vtk_file;

        memset(&field, 0, sizeof(field));
        tiff_output_path = build_realization_output_path(stack_path, i, "tif");
        vtk_output_path = build_realization_output_path(stack_path, i, "vtk");
        if (tiff_output_path == NULL || vtk_output_path == NULL) {
            free(tiff_output_path);
            free(vtk_output_path);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        if (realization_method == POROSITY_REALIZATION_METHOD_CHOLESKY) {
            if (generate_gaussian_porosity_realization_from_model(
                    gen_model, simulation_grid,
                    gen_mean, base_seed + i, &field) !=
                EXIT_SUCCESS) {
                free(tiff_output_path);
                free(vtk_output_path);
                free_porosity_field_storage(&field);
                fclose(diagnostics_file);
                free(diagnostics_path);
                return EXIT_FAILURE;
            }
        } else if (realization_method == POROSITY_REALIZATION_METHOD_FFT) {
            if (generate_gaussian_porosity_realization_from_model_fft(
                    gen_model, simulation_grid,
                    gen_mean, base_seed + i, &field) !=
                EXIT_SUCCESS) {
                free(tiff_output_path);
                free(vtk_output_path);
                free_porosity_field_storage(&field);
                fclose(diagnostics_file);
                free(diagnostics_path);
                return EXIT_FAILURE;
            }
        } else {
            fprintf(stderr, "Unknown porosity realization method.\n");
            free(tiff_output_path);
            free(vtk_output_path);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
            
        }
        //Map every cell back to porosity through the data's sorted values so
        //units are consistent
        if (score_table != NULL &&
            normal_score_back_transform_field(&field, score_table) !=
                EXIT_SUCCESS) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }


        {
            PorosityGenerationDiagnostics diagnostics;

            if (compute_porosity_generation_diagnostics(
                    &field, &diagnostics) != EXIT_SUCCESS ||
                write_generation_diagnostics_row(
                    diagnostics_file, i, &diagnostics) != EXIT_SUCCESS) {
                free(tiff_output_path);
                free(vtk_output_path);
                free_porosity_field_storage(&field);
                fclose(diagnostics_file);
                free(diagnostics_path);
                return EXIT_FAILURE;
            }
        }

        if (write_porosity_field_tiff(tiff_output_path, &field) !=
            EXIT_SUCCESS) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        vtk_file = fopen(vtk_output_path, "w");
        if (vtk_file == NULL) {
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        if (write_porosity_field_vtk(vtk_file, &field, simulation_grid) !=
            EXIT_SUCCESS) {
            fclose(vtk_file);
            free(tiff_output_path);
            free(vtk_output_path);
            free_porosity_field_storage(&field);
            fclose(diagnostics_file);
            free(diagnostics_path);
            return EXIT_FAILURE;
        }

        fclose(vtk_file);
        free(tiff_output_path);
        free(vtk_output_path);
        free_porosity_field_storage(&field);
    }

    fclose(diagnostics_file);
    free(diagnostics_path);
    return EXIT_SUCCESS;
}

int generate_and_write_gaussian_porosity_realizations_from_model_with_method(
                        const char *stack_path,
                        const FittedCovarianceModelArtifact *model,
                        const PorosityRealizationGridSpec *simulation_grid,
                        PorosityRealizationMethod realization_method,
                        uint32_t realization_count,
                        unsigned int base_seed
                    ) {
        NormalScoreTable score_table = {0U, NULL};
    ExponentialCovarianceModel generation_model;
    char *table_path;
    int status;

    if (stack_path == NULL || model == NULL) {
        return EXIT_FAILURE;
    }

    if (!model->normal_score_transform) {
        /* Original behavior: sample porosity directly. */
        return write_gaussian_realizations_impl(
            stack_path, &model->exponential_model, model->mean_porosity, NULL,
            simulation_grid, realization_method, realization_count, base_seed);
    }

    /*
     * The model was fit to normal scores. The back-transform
     *     phi = F_data^-1( Phi(y) )
     * assumes y is STANDARD normal: mean 0, variance 1. The fitted sill is the
     * empirical score variance, which sits slightly below 1 when porosity
     * values are tied (tied blocks share one score), so keep only the SHAPE of
     * the fit, the correlation C(r)/C(0), and sample with mean 0, sill 1.
     */
    generation_model = model->exponential_model;
    generation_model.sill = 1.0;

    table_path = build_normal_score_table_path(stack_path);
    if (table_path == NULL ||
        read_normal_score_table(table_path, &score_table) != EXIT_SUCCESS) {
        fprintf(stderr,
                "Model was fit with NORMAL_SCORE_TRANSFORM but the lookup table "
                "could not be read: %s\n",
                table_path != NULL ? table_path : "(path error)");
        free(table_path);
        return EXIT_FAILURE;
    }
    free(table_path);

    status = write_gaussian_realizations_impl(
        stack_path, &generation_model, 0.0, &score_table, simulation_grid,
        realization_method, realization_count, base_seed);

    free_normal_score_table(&score_table);
    return status;
}
