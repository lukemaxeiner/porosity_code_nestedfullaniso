#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "input.h"

static void release_batch(InputBatch *batch) {
    if (batch == NULL || batch->stacks == NULL) {
        return;
    }

    for (int i = 0; i < batch->stack_count; ++i) {
        free(batch->stacks[i].stack_path);
    }

    free(batch->stacks);
    batch->stacks = NULL;
    batch->stack_count = 0;
}

void free_input(InputBatch *batch) { release_batch(batch); }

static char *trim_whitespace(char *text) {
    char *end;

    if (text == NULL) {
        return NULL;
    }

    while (*text != '\0' && isspace((unsigned char)*text)) {
        ++text;
    }

    if (*text == '\0') {
        return text;
    }

    end = text + strlen(text) - 1;
    while (end > text && isspace((unsigned char)*end)) {
        *end = '\0';
        --end;
    }

    return text;
}

static int equals_ignore_case(const char *left, const char *right) {
    while (*left != '\0' && *right != '\0') {
        if (tolower((unsigned char)*left) != tolower((unsigned char)*right)) {
            return 0;
        }
        ++left;
        ++right;
    }

    return *left == '\0' && *right == '\0';
}

static int is_comment_or_blank(const char *line) {
    while (*line != '\0' && isspace((unsigned char)*line)) {
        ++line;
    }

    return *line == '\0' || *line == '#';
}

static char *duplicate_string(const char *text) {
    size_t length;
    char *copy;

    if (text == NULL) {
        return NULL;
    }

    length = strlen(text);
    copy = (char *)malloc(length + 1U);
    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy, text, length + 1U);
    return copy;
}

static int append_input_entry(InputBatch *batch, int *capacity,
                              const char *path, double voxel_resolution) {
    if (batch->stack_count == *capacity) {
        int new_capacity = *capacity == 0 ? 4 : (*capacity * 2);
        InputOptions *resized =
            (InputOptions *)realloc(batch->stacks,
                                    (size_t)new_capacity *
                                        sizeof(InputOptions));
        if (resized == NULL) {
            fprintf(stderr, "Unable to allocate input records\n");
            return EXIT_FAILURE;
        }
        batch->stacks = resized;
        *capacity = new_capacity;
    }

    batch->stacks[batch->stack_count].stack_path = duplicate_string(path);
    if (batch->stacks[batch->stack_count].stack_path == NULL) {
        fprintf(stderr, "Unable to allocate stack path\n");
        return EXIT_FAILURE;
    }
    batch->stacks[batch->stack_count].voxel_resolution = voxel_resolution;

    batch->stack_count += 1;
    return EXIT_SUCCESS;
}

static int parse_yes_no(const char *value, int *enabled) {
    if (value == NULL || enabled == NULL) {
        return EXIT_FAILURE;
    }

    if (equals_ignore_case(value, "YES")) {
        *enabled = 1;
        return EXIT_SUCCESS;
    }
    if (equals_ignore_case(value, "NO")) {
        *enabled = 0;
        return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
}

static int parse_porosity_realization_method(
    const char *value,
    PorosityRealizationMethod *method
) {
    if (value == NULL || method == NULL) {
        return EXIT_FAILURE;
    }

    if (equals_ignore_case(value, "CHOLESKY")) {
        *method = POROSITY_REALIZATION_METHOD_CHOLESKY;
        return EXIT_SUCCESS;
    }

    if (equals_ignore_case(value, "FFT") ||
        equals_ignore_case(value, "CIRCULANT") ||
        equals_ignore_case(value, "CIRCULANT_EMBEDDING")) {
        *method = POROSITY_REALIZATION_METHOD_FFT;
        return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
}

static int parse_porosity_lag_mode(const char *value,
                                   PorosityLagMode *mode) {
    if (value == NULL || mode == NULL) {
        return EXIT_FAILURE;
    }

    if (equals_ignore_case(value, "NO")) {
        *mode = POROSITY_LAG_MODE_NO;
        return EXIT_SUCCESS;
    }
    if (equals_ignore_case(value, "3DLAG")) {
        *mode = POROSITY_LAG_MODE_3DLAG;
        return EXIT_SUCCESS;
    }
    if (equals_ignore_case(value, "CURVE")) {
        *mode = POROSITY_LAG_MODE_CURVE;
        return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
}

static int parse_porosity_axis(const char *value, PorosityAxis *axis) {
    if (value == NULL || axis == NULL) {
        return EXIT_FAILURE;
    }

    if (equals_ignore_case(value, "X")) {
        *axis = POROSITY_AXIS_X;
        return EXIT_SUCCESS;
    }
    if (equals_ignore_case(value, "Y")) {
        *axis = POROSITY_AXIS_Y;
        return EXIT_SUCCESS;
    }
    if (equals_ignore_case(value, "Z")) {
        *axis = POROSITY_AXIS_Z;
        return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
}

static int parse_three_ints(const char *value, int *x, int *y, int *z) {
    int parsed;

    if (value == NULL || x == NULL || y == NULL || z == NULL) {
        return EXIT_FAILURE;
    }

    parsed = sscanf(value, " %d %d %d", x, y, z);
    if (parsed != 3) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static int parse_three_doubles(const char *value, double *x, double *y,
                               double *z) {
    int parsed;

    if (value == NULL || x == NULL || y == NULL || z == NULL) {
        return EXIT_FAILURE;
    }

    parsed = sscanf(value, " %lf %lf %lf", x, y, z);
    if (parsed != 3) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static int parse_input_data_value(const char *value, InputBatch *batch,
                                  int *capacity) {
    const char *cursor = value;

    if (value == NULL || batch == NULL || capacity == NULL) {
        return EXIT_FAILURE;
    }

    while (*cursor != '\0') {
        char entry_buffer[4096];
        char stack_path_buffer[4096];
        double voxel_resolution;
        size_t entry_length;
        const char *entry_end;
        int parsed_count;

        while (*cursor != '\0' && isspace((unsigned char)*cursor)) {
            ++cursor;
        }

        if (*cursor == '\0') {
            break;
        }

        if (*cursor != '"') {
            fprintf(stderr,
                    "Invalid INPUT_DATA format. Expected quoted entries separated by ';'\n");
            return EXIT_FAILURE;
        }
        ++cursor;

        entry_end = strchr(cursor, '"');
        if (entry_end == NULL) {
            fprintf(stderr, "Invalid INPUT_DATA format. Missing closing quote\n");
            return EXIT_FAILURE;
        }

        entry_length = (size_t)(entry_end - cursor);
        if (entry_length >= sizeof(entry_buffer)) {
            fprintf(stderr, "INPUT_DATA entry is too long\n");
            return EXIT_FAILURE;
        }

        memcpy(entry_buffer, cursor, entry_length);
        entry_buffer[entry_length] = '\0';

        parsed_count = sscanf(entry_buffer, " %4095s %lf", stack_path_buffer,
                              &voxel_resolution);
        if (parsed_count != 2) {
            fprintf(stderr,
                    "Invalid INPUT_DATA entry. Expected: \"<stack.tif> <voxel_resolution>\"\n");
            return EXIT_FAILURE;
        }

        if (voxel_resolution <= 0.0) {
            fprintf(stderr,
                    "INPUT_DATA voxel_resolution must be positive\n");
            return EXIT_FAILURE;
        }

        if (append_input_entry(batch, capacity, stack_path_buffer,
                               voxel_resolution) !=
            EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }

        cursor = entry_end + 1;
        while (*cursor != '\0' && isspace((unsigned char)*cursor)) {
            ++cursor;
        }

        if (*cursor == ';') {
            ++cursor;
            continue;
        }

        if (*cursor != '\0') {
            fprintf(stderr,
                    "Invalid INPUT_DATA format. Entries must end with ';'\n");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int parse_input(int argc, char **argv, RunConfig *config) {
    FILE *input_file;
    char line[8192];
    int capacity = 0;
    int input_data_seen = 0;
    int filter_resolution_seen = 0;
    int run_full_pipeline_seen = 0;
    int empirical_covariance_seen = 0;
    int fit_covariance_seen = 0;
    int generate_realizations_seen = 0;
    int porosity_axis_seen = 0;
    int realization_count_seen = 0;
    int realization_grid_res_seen = 0;
    int realization_grid_dim_seen = 0;
    int ref_window_values_seen = 0;
    int validation_max_lag_seen = 0; //realization validation

    if (config == NULL) {
        return EXIT_FAILURE;
    }

    memset(config, 0, sizeof(*config));
    config->porosity_lag_mode = POROSITY_LAG_MODE_NO;
    config->porosity_axis = POROSITY_AXIS_X;
    config->porosity_realization_method = POROSITY_REALIZATION_METHOD_CHOLESKY;

    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s <input_config.txt>\n"
                "Expected entries:\n"
                "  RUN_FULL_PIPELINE = YES|NO\n"
                "  EMPIRICAL_COVARIANCE = YES|NO\n"
                "  FIT_COVARIANCE_MODEL = YES|NO\n"
                "  GENERATE_POROSITY_REALIZATIONS = YES|NO\n"
                "  POROSITY_REALIZATION_METHOD = CHOLESKY|FFT\n"
                "  NORMAL_SCORE_TRANSFORM = YES|NO\n"
                "  VALIDATE REALIZATIONS = YES|NO\n"
                "  VALIDATION_MAX_LAG = <N> (0 = automatic)\n"
                "  GLOBAL_POROSITY = YES|NO\n"
                "  POROSITY_FIELD = YES|NO\n"
                "  FILTER_RESOLUTION = <value>\n"
                "  POROSITY_COVARIANCE_LAG = NO|3DLAG|CURVE\n"
                "  POROSITY_AXIS = X|Y|Z\n"
                "  POROSITY_REALIZATION_COUNT = <N>\n"
                "  POROSITY_REALIZATION_GRID_RES = <Nx> <Ny> <Nz>\n"
                "  POROSITY_REALIZATION_GRID_DIM = <x> <y> <z>\n"
                "  POROSITY_COVARIANCE_REF_WINDOW = YES|NO\n"
                "  REF_WINDOW_VALUES = <x> <y> <z>\n"
                "  INPUT_DATA = \"<stack.tif> <voxel_resolution>\" ;\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    input_file = fopen(argv[1], "r");
    if (input_file == NULL) {
        fprintf(stderr, "Unable to open input file: %s\n", argv[1]);
        return EXIT_FAILURE;
    }

    while (fgets(line, sizeof(line), input_file) != NULL) {
        char *key;
        char *value;
        char *equals_sign;

        key = trim_whitespace(line);
        if (is_comment_or_blank(key)) {
            continue;
        }

        equals_sign = strchr(key, '=');
        if (equals_sign == NULL) {
            fprintf(stderr, "Invalid input line. Expected KEY = VALUE\n");
            fclose(input_file);
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }

        *equals_sign = '\0';
        value = trim_whitespace(equals_sign + 1);
        key = trim_whitespace(key);

        if (equals_ignore_case(key, "RUN_FULL_PIPELINE")) {
            if (parse_yes_no(value, &config->run_full_pipeline) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid RUN_FULL_PIPELINE value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            run_full_pipeline_seen = 1;
        } else if (equals_ignore_case(key, "EMPIRICAL_COVARIANCE")) {
            if (parse_yes_no(value, &config->empirical_covariance_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid EMPIRICAL_COVARIANCE value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            empirical_covariance_seen = 1;
        } else if (equals_ignore_case(key, "FIT_COVARIANCE_MODEL")) {
            if (parse_yes_no(value, &config->fit_covariance_model_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid FIT_COVARIANCE_MODEL value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            fit_covariance_seen = 1;
        } else if (equals_ignore_case(key, "GENERATE_POROSITY_REALIZATIONS")) {
            if (parse_yes_no(value,
                             &config->generate_porosity_realizations_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid GENERATE_POROSITY_REALIZATIONS value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            generate_realizations_seen = 1;

        } else if (equals_ignore_case(key, "POROSITY_REALIZATION_METHOD")) {
            if (parse_porosity_realization_method(
                    value,
                    &config->porosity_realization_method
                ) != EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_REALIZATION_METHOD value. Expected CHOLESKY or FFT\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "NORMAL_SCORE_TRANSFORM")){
            if(parse_yes_no(value, &config->normal_score_transform_enabled) != EXIT_SUCCESS){
                fprintf(stderr,
                    "Invalid NORMAL_SCORE_TRANSFORM value. Expected YES or NO\n");
                    fclose(input_file);
                    release_batch(&config->input_data);
                    return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "VALIDATE_REALIZATIONS")) {
            if (parse_yes_no(value, &config->validate_realization_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid VALIDATE_REALIZATIONS value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "VALIDATION_MAX_LAG")) {
            if (sscanf(value, " %d", &config->validation_max_lag) != 1) {
                fprintf(stderr,
                        "Invalid VALIDATION_MAX_LAG. Expected: <N>\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            validation_max_lag_seen = 1;
        } else if (equals_ignore_case(key, "GLOBAL_POROSITY")) {
            if (parse_yes_no(value, &config->global_porosity_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid GLOBAL_POROSITY value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "POROSITY_FIELD")) {
            if (parse_yes_no(value, &config->porosity_field_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_FIELD value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "FILTER_RESOLUTION")) {
            if (sscanf(value, " %lf", &config->filter_resolution) != 1) {
                fprintf(stderr,
                        "Invalid FILTER_RESOLUTION. Expected: <value>\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            filter_resolution_seen = 1;
        } else if (equals_ignore_case(key, "POROSITY_COVARIANCE_LAG")) {
            if (parse_porosity_lag_mode(value, &config->porosity_lag_mode) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_COVARIANCE_LAG value. Expected NO, 3DLAG, or CURVE\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "POROSITY_AXIS")) {
            if (parse_porosity_axis(value, &config->porosity_axis) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_AXIS value. Expected X, Y, or Z\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            porosity_axis_seen = 1;
        } else if (equals_ignore_case(key, "POROSITY_REALIZATION_COUNT")) {
            if (sscanf(value, " %d", &config->porosity_realization_count) != 1) {
                fprintf(stderr,
                        "Invalid POROSITY_REALIZATION_COUNT. Expected: <N>\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            realization_count_seen = 1;
        } else if (equals_ignore_case(key, "POROSITY_REALIZATION_GRID_RES")) {
            if (parse_three_ints(value, &config->realization_grid_nx,
                                 &config->realization_grid_ny,
                                 &config->realization_grid_nz) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_REALIZATION_GRID_RES. Expected: <Nx> <Ny> <Nz>\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            realization_grid_res_seen = 1;
        } else if (equals_ignore_case(key, "POROSITY_REALIZATION_GRID_DIM")) {
            if (parse_three_doubles(value, &config->realization_dim_x,
                                    &config->realization_dim_y,
                                    &config->realization_dim_z) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_REALIZATION_GRID_DIM. Expected: <x> <y> <z>\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            realization_grid_dim_seen = 1;
        } else if (equals_ignore_case(key, "POROSITY_COVARIANCE_REF_WINDOW")) {
            if (parse_yes_no(value, &config->porosity_ref_window_enabled) !=
                EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid POROSITY_COVARIANCE_REF_WINDOW value. Expected YES or NO\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else if (equals_ignore_case(key, "REF_WINDOW_VALUES")) {
            if (parse_three_ints(value, &config->ref_window_x,
                                 &config->ref_window_y,
                                 &config->ref_window_z) != EXIT_SUCCESS) {
                fprintf(stderr,
                        "Invalid REF_WINDOW_VALUES. Expected: <x> <y> <z>\n");
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
            ref_window_values_seen = 1;
        } else if (equals_ignore_case(key, "INPUT_DATA")) {
            input_data_seen = 1;
            if (parse_input_data_value(value, &config->input_data, &capacity) !=
                EXIT_SUCCESS) {
                fclose(input_file);
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        } else {
            fprintf(stderr, "Unknown input key: %s\n", key);
            fclose(input_file);
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
    }

    fclose(input_file);

    if (!input_data_seen || config->input_data.stack_count == 0) {
        fprintf(stderr, "INPUT_DATA is required and must contain at least one entry\n");
        release_batch(&config->input_data);
        return EXIT_FAILURE;
    }

    if (!run_full_pipeline_seen) {
        config->run_full_pipeline = 0;
    }
    if (!empirical_covariance_seen) {
        config->empirical_covariance_enabled = 0;
    }
    if (!fit_covariance_seen) {
        config->fit_covariance_model_enabled = 0;
    }
    if (!generate_realizations_seen) {
        config->generate_porosity_realizations_enabled = 0;
    }

    if (config->run_full_pipeline) {
        config->empirical_covariance_enabled = 1;
        config->fit_covariance_model_enabled = 1;
        config->generate_porosity_realizations_enabled = 1;
    }

    if (!filter_resolution_seen) {
        fprintf(stderr, "FILTER_RESOLUTION is required\n");
        release_batch(&config->input_data);
        return EXIT_FAILURE;
    }

    if (config->filter_resolution <= 0.0) {
        fprintf(stderr, "FILTER_RESOLUTION must be positive\n");
        release_batch(&config->input_data);
        return EXIT_FAILURE;
    }

    if (config->generate_porosity_realizations_enabled ||
        config->fit_covariance_model_enabled) {
        if (!realization_count_seen) {
            if (config->generate_porosity_realizations_enabled) {
                fprintf(stderr,
                        "POROSITY_REALIZATION_COUNT is required when GENERATE_POROSITY_REALIZATIONS = YES\n");
                release_batch(&config->input_data);
                return EXIT_FAILURE;
            }
        }
        if (!realization_grid_res_seen) {
            fprintf(stderr,
                    "POROSITY_REALIZATION_GRID_RES is required when FIT_COVARIANCE_MODEL = YES or GENERATE_POROSITY_REALIZATIONS = YES\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
        if (!realization_grid_dim_seen) {
            fprintf(stderr,
                    "POROSITY_REALIZATION_GRID_DIM is required when FIT_COVARIANCE_MODEL = YES or GENERATE_POROSITY_REALIZATIONS = YES\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
        if (config->realization_grid_nx <= 0 || config->realization_grid_ny <= 0 ||
            config->realization_grid_nz <= 0) {
            fprintf(stderr,
                    "POROSITY_REALIZATION_GRID_RES values must all be positive\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
        if (config->generate_porosity_realizations_enabled &&
            config->porosity_realization_count <= 0) {
            fprintf(stderr,
                    "POROSITY_REALIZATION_COUNT must be positive\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
        if (config->realization_dim_x <= 0.0 || config->realization_dim_y <= 0.0 ||
            config->realization_dim_z <= 0.0) {
            fprintf(stderr,
                    "POROSITY_REALIZATION_GRID_DIM values must all be positive\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
    }

    if (config->empirical_covariance_enabled ||
        config->fit_covariance_model_enabled ||
        config->generate_porosity_realizations_enabled) {
        if (config->porosity_lag_mode == POROSITY_LAG_MODE_NO) {
            config->porosity_lag_mode = POROSITY_LAG_MODE_3DLAG;
        }
    }

    if ((config->empirical_covariance_enabled || config->run_full_pipeline) &&
        config->porosity_lag_mode != POROSITY_LAG_MODE_3DLAG) {
        fprintf(stderr,
                "EMPIRICAL_COVARIANCE requires POROSITY_COVARIANCE_LAG = 3DLAG\n");
        release_batch(&config->input_data);
        return EXIT_FAILURE;
    }

    if (config->porosity_lag_mode == POROSITY_LAG_MODE_CURVE &&
        !porosity_axis_seen) {
        fprintf(stderr,
                "POROSITY_AXIS is required when POROSITY_COVARIANCE_LAG = CURVE\n");
        release_batch(&config->input_data);
        return EXIT_FAILURE;
    }

    if (config->validate_realization_enabled) {
        if (!config->generate_porosity_realizations_enabled) {
            fprintf(stderr,
                    "VALIDATE_REALIZATIONS = YES requires GENERATE_POROSITY_REALIZATIONS = YES\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
        if (config->porosity_realization_count < 2) {
            fprintf(stderr,
                    "VALIDATE_REALIZATIONS = YES requires POROSITY_REALIZATION_COUNT >= 2 "
                    "(10 or more recommended)\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
        if (config->validation_max_lag < 0) {
            fprintf(stderr, "VALIDATION_MAX_LAG must be 0 (automatic) or positive\n");
            release_batch(&config->input_data);
            return EXIT_FAILURE;
        }
    } else if (validation_max_lag_seen) {
        fprintf(stderr,
                "Warning: VALIDATION_MAX_LAG is ignored because VALIDATE_REALIZATIONS = NO\n");
    }

    if (config->porosity_ref_window_enabled && !ref_window_values_seen) {
        fprintf(stderr,
                "REF_WINDOW_VALUES is required when POROSITY_COVARIANCE_REF_WINDOW = YES\n");
        release_batch(&config->input_data);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
