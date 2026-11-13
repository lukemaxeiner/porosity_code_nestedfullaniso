#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "input.h"
#include "porosity_covariance_lag.h"
#include "porosity_covariance_realization_generation.h"
#include "porosity_covariance_realization_prep.h"
#include "porosity_covariance_refwindow.h"
#include "porosity_field.h"
#include "porosity_global.h"

static char *build_fitted_model_path(const char *stack_path) {
    const char *extension = strrchr(stack_path, '.');
    const char *slash = strrchr(stack_path, '/');
    const char *file_name = slash != NULL ? slash + 1 : stack_path;
    size_t dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    size_t file_base_length =
        extension != NULL && extension > file_name
            ? (size_t)(extension - file_name)
            : strlen(file_name);
    const char *suffix = "_fitted_covariance_model.txt";
    size_t total_length = dir_length + strlen("/BUILD_COVARIANCE/") +
                          file_base_length + strlen(suffix) + 1U;
    char *output_path = (char *)malloc(total_length);

    if (output_path == NULL) {
        return NULL;
    }

    if (slash != NULL) {
        memcpy(output_path, stack_path, dir_length);
    } else {
        output_path[0] = '.';
    }
    memcpy(output_path + dir_length, "/BUILD_COVARIANCE/",
           strlen("/BUILD_COVARIANCE/"));
    memcpy(output_path + dir_length + strlen("/BUILD_COVARIANCE/"), file_name,
           file_base_length);
    memcpy(output_path + dir_length + strlen("/BUILD_COVARIANCE/") +
               file_base_length,
           suffix, strlen(suffix) + 1U);

    return output_path;
}

int main(int argc, char **argv) {
    RunConfig config;
    double porosity_sum = 0.0;

    if (parse_input(argc, argv, &config) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    for (int i = 0; i < config.input_data.stack_count; ++i) {
        InputOptions *stack = &config.input_data.stacks[i];
        PorosityRealizationGridSpec realization_grid;

        realization_grid.nx = (uint32_t)config.realization_grid_nx;
        realization_grid.ny = (uint32_t)config.realization_grid_ny;
        realization_grid.nz = (uint32_t)config.realization_grid_nz;
        realization_grid.dim_x = config.realization_dim_x;
        realization_grid.dim_y = config.realization_dim_y;
        realization_grid.dim_z = config.realization_dim_z;

        if (config.global_porosity_enabled &&
            process_global_porosity(stack->stack_path, &porosity_sum) !=
                EXIT_SUCCESS) {
            free_input(&config.input_data);
            return EXIT_FAILURE;
        }

        if (config.porosity_field_enabled &&
            process_porosity_field(stack->stack_path, stack->voxel_resolution,
                                   config.filter_resolution, stdout) !=
                EXIT_SUCCESS) {
            free_input(&config.input_data);
            return EXIT_FAILURE;
        }

        if (config.porosity_lag_mode != POROSITY_LAG_MODE_NO &&
            !config.empirical_covariance_enabled &&
            process_porosity_covariance_lag(
                stack->stack_path, stack->voxel_resolution,
                config.filter_resolution, config.porosity_lag_mode,
                config.porosity_axis, stdout) != EXIT_SUCCESS) {
            free_input(&config.input_data);
            return EXIT_FAILURE;
        }

        if (config.porosity_ref_window_enabled &&
            process_porosity_covariance_ref_window(
                stack->stack_path, stack->voxel_resolution,
                config.filter_resolution, config.ref_window_x,
                config.ref_window_y, config.ref_window_z, stdout) !=
                EXIT_SUCCESS) {
            free_input(&config.input_data);
            return EXIT_FAILURE;
        }

        if (config.empirical_covariance_enabled &&
            process_porosity_covariance_lag(
                stack->stack_path, stack->voxel_resolution,
                config.filter_resolution, POROSITY_LAG_MODE_3DLAG,
                config.porosity_axis, stdout) != EXIT_SUCCESS) {
            free_input(&config.input_data);
            return EXIT_FAILURE;
        }

        if (config.fit_covariance_model_enabled &&
            process_realization_fit_from_stack(
                stack->stack_path, stack->voxel_resolution,
                config.filter_resolution, &realization_grid, stdout) !=
                EXIT_SUCCESS) {
            free_input(&config.input_data);
            return EXIT_FAILURE;
        }

        if (config.generate_porosity_realizations_enabled) {
            FittedCovarianceModelArtifact model_artifact;
            char *model_path = build_fitted_model_path(stack->stack_path);

            if (model_path == NULL) {
                free_input(&config.input_data);
                return EXIT_FAILURE;
            }
            if (read_fitted_covariance_model(model_path, &model_artifact) !=
                EXIT_SUCCESS) {
                fprintf(stderr, "Failed to read fitted covariance model for %s\n",
                        stack->stack_path);
                free(model_path);
                free_input(&config.input_data);
                return EXIT_FAILURE;
            }

            if (generate_and_write_gaussian_porosity_realizations_from_model_with_method(
                    stack->stack_path,
                    &model_artifact,
                    &realization_grid,
                    config.porosity_realization_method,
                    (uint32_t)config.porosity_realization_count,
                    12345U) != EXIT_SUCCESS) {
                free(model_path);
                free_input(&config.input_data);
                return EXIT_FAILURE;
            }

            if (fprintf(stdout,
                        "Porosity realizations: %s\n"
                        "  Fitted model: %s\n"
                        "  Realizations written: %d\n\n",
                        stack->stack_path, model_path,
                        config.porosity_realization_count) < 0) {
                free(model_path);
                free_input(&config.input_data);
                return EXIT_FAILURE;
            }

            free(model_path);
        }
    }

    if (config.global_porosity_enabled && config.input_data.stack_count > 1) {
        double average_porosity =
            porosity_sum / (double)config.input_data.stack_count;
        printf("Average porosity across %d stacks: %.6f (%.2f%%)\n",
               config.input_data.stack_count, average_porosity,
               average_porosity * 100.0);
    }

    free_input(&config.input_data);
    return EXIT_SUCCESS;
}
