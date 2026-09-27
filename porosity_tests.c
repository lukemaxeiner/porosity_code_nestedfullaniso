#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#define POROSITY_STATIC
#define main porosity_entry_main
#include "porosity.c"
#undef main

typedef int (*test_fn)(void);

static const char *test_root_dir =
    "/Users/saviopoovathingal/Desktop/FiberFormTrials/test_artifacts";
static char test_run_dir[PATH_MAX];
static unsigned int test_file_counter = 0;

static int expect_int(const char *label, long long actual, long long expected) {
    if (actual != expected) {
        fprintf(stderr, "%s: expected %lld, got %lld\n", label, expected, actual);
        return 0;
    }
    return 1;
}

static int expect_contains(const char *label, const char *haystack,
                           const char *needle) {
    if (strstr(haystack, needle) == NULL) {
        fprintf(stderr, "%s: missing substring '%s'\n", label, needle);
        return 0;
    }
    return 1;
}

static void ensure_tmp_dir(void) {
    static int ready = 0;
    if (!ready) {
        int pid = (int)getpid();
        if (mkdir(test_root_dir, 0700) != 0 && errno != EEXIST) {
            fprintf(stderr, "mkdir failed for %s: %s\n", test_root_dir,
                    strerror(errno));
            exit(EXIT_FAILURE);
        }
        snprintf(test_run_dir, sizeof(test_run_dir), "%s/run_%d", test_root_dir,
                 pid);
        if (mkdir(test_run_dir, 0700) != 0 && errno != EEXIST) {
            fprintf(stderr, "mkdir failed for %s: %s\n", test_run_dir,
                    strerror(errno));
            exit(EXIT_FAILURE);
        }
        ready = 1;
    }
}

static void build_path(char *buffer, size_t size, const char *filename) {
    ensure_tmp_dir();
    snprintf(buffer, size, "%s/%s", test_run_dir, filename);
}

static void build_unique_path(char *buffer, size_t size, const char *prefix,
                              const char *extension) {
    char filename[128];

    snprintf(filename, sizeof(filename), "%s_%u%s", prefix, test_file_counter++,
             extension);
    build_path(buffer, size, filename);
}

static void create_grayscale_stack(const char *path, uint32_t width,
                                   uint32_t height, uint16_t slices,
                                   const uint8_t *voxels) {
    TIFF *tif = TIFFOpen(path, "w");
    if (tif == NULL) {
        fprintf(stderr, "Failed to create TIFF: %s\n", path);
        exit(EXIT_FAILURE);
    }

    for (uint16_t slice = 0; slice < slices; ++slice) {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

        const uint8_t *slice_data = voxels + ((size_t)slice * width * height);
        for (uint32_t row = 0; row < height; ++row) {
            if (TIFFWriteScanline(tif, (void *)(slice_data + row * width), row,
                                  0) < 0) {
                TIFFClose(tif);
                fprintf(stderr, "Failed to write grayscale scanline\n");
                exit(EXIT_FAILURE);
            }
        }

        if (slice + 1 < slices) {
            TIFFWriteDirectory(tif);
        }
    }

    TIFFClose(tif);
}

static void create_palette_stack(const char *path, uint32_t width,
                                 uint32_t height, uint16_t slices,
                                 const uint8_t *voxels) {
    TIFF *tif = TIFFOpen(path, "w");
    uint16_t red[256] = {0};
    uint16_t green[256] = {0};
    uint16_t blue[256] = {0};

    if (tif == NULL) {
        fprintf(stderr, "Failed to create palette TIFF: %s\n", path);
        exit(EXIT_FAILURE);
    }

    red[1] = 65535;
    green[1] = 65535;
    blue[1] = 65535;

    for (uint16_t slice = 0; slice < slices; ++slice) {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_COLORMAP, red, green, blue);

        const uint8_t *slice_data = voxels + ((size_t)slice * width * height);
        for (uint32_t row = 0; row < height; ++row) {
            if (TIFFWriteScanline(tif, (void *)(slice_data + row * width), row,
                                  0) < 0) {
                TIFFClose(tif);
                fprintf(stderr, "Failed to write palette scanline\n");
                exit(EXIT_FAILURE);
            }
        }

        if (slice + 1 < slices) {
            TIFFWriteDirectory(tif);
        }
    }

    TIFFClose(tif);
}

static int read_stream(FILE *stream, char *buffer, size_t size) {
    long length;
    if (fflush(stream) != 0) {
        return 0;
    }
    if (fseek(stream, 0, SEEK_END) != 0) {
        return 0;
    }
    length = ftell(stream);
    if (length < 0) {
        return 0;
    }
    if (fseek(stream, 0, SEEK_SET) != 0) {
        return 0;
    }
    if ((size_t)length >= size) {
        length = (long)size - 1;
    }
    if (fread(buffer, 1, (size_t)length, stream) != (size_t)length) {
        return 0;
    }
    buffer[length] = '\0';
    return 1;
}

static FILE *open_capture_file(char *path_buffer, size_t path_size,
                               const char *prefix) {
    build_unique_path(path_buffer, path_size, prefix, ".txt");
    return fopen(path_buffer, "w+");
}

static int capture_stdio_call(int target_fd, int (*fn)(int, char **), int argc,
                              char **argv, char *buffer, size_t size,
                              int *result) {
    int saved_fd = dup(target_fd);
    char capture_path[PATH_MAX];
    FILE *tmp = open_capture_file(capture_path, sizeof(capture_path), "capture");
    if (saved_fd < 0 || tmp == NULL) {
        return 0;
    }

    fflush(stdout);
    fflush(stderr);
    if (dup2(fileno(tmp), target_fd) < 0) {
        close(saved_fd);
        fclose(tmp);
        return 0;
    }

    *result = fn(argc, argv);
    fflush(stdout);
    fflush(stderr);

    if (dup2(saved_fd, target_fd) < 0) {
        close(saved_fd);
        fclose(tmp);
        return 0;
    }

    close(saved_fd);
    if (!read_stream(tmp, buffer, size)) {
        fclose(tmp);
        return 0;
    }

    fclose(tmp);
    return 1;
}

static int capture_die(const char *message, const char *path, char *buffer,
                       size_t size, int *exit_code) {
    int pipefd[2];
    pid_t pid;

    if (pipe(pipefd) != 0) {
        return 0;
    }

    pid = fork();
    if (pid < 0) {
        close(pipefd[0]);
        close(pipefd[1]);
        return 0;
    }

    if (pid == 0) {
        close(pipefd[0]);
        dup2(pipefd[1], STDERR_FILENO);
        close(pipefd[1]);
        die(message, path);
        _exit(EXIT_SUCCESS);
    }

    close(pipefd[1]);
    ssize_t count = read(pipefd[0], buffer, size - 1);
    int status = 0;
    close(pipefd[0]);
    waitpid(pid, &status, 0);

    if (count < 0) {
        return 0;
    }

    buffer[count] = '\0';
    if (WIFEXITED(status)) {
        *exit_code = WEXITSTATUS(status);
        return 1;
    }

    return 0;
}

static int test_die_emits_message_and_failure(void) {
    char stderr_buffer[256];
    int exit_code = 0;

    if (!capture_die("boom", "fixture.tif", stderr_buffer, sizeof(stderr_buffer),
                     &exit_code)) {
        return 0;
    }

    return expect_int("die exit code", exit_code, EXIT_FAILURE) &&
           expect_contains("die stderr", stderr_buffer, "boom: fixture.tif");
}

static int test_palette_index_to_binary_without_palette(void) {
    return expect_int("zero maps to solid",
                      palette_index_to_binary(0, NULL, NULL, NULL), 0) &&
           expect_int("non-zero maps to void",
                      palette_index_to_binary(9, NULL, NULL, NULL), 1);
}

static int test_palette_index_to_binary_with_palette(void) {
    uint16_t red[256] = {0};
    uint16_t green[256] = {0};
    uint16_t blue[256] = {0};

    red[7] = 65535;
    green[7] = 65535;
    blue[7] = 65535;

    return expect_int("black palette entry is solid",
                      palette_index_to_binary(0, red, green, blue), 0) &&
           expect_int("white palette entry is void",
                      palette_index_to_binary(7, red, green, blue), 1);
}

static int test_accumulate_porosity_grayscale_stack(void) {
    char path[PATH_MAX];
    TIFF *tif;
    StackStats stats;
    const uint8_t voxels[] = {
        0, 1, 0, 1,
        1, 1, 0, 0,
    };

    build_path(path, sizeof(path), "grayscale_stack.tif");
    create_grayscale_stack(path, 2, 2, 2, voxels);

    memset(&stats, 0, sizeof(stats));
    tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        return 0;
    }

    accumulate_porosity(tif, &stats);
    TIFFClose(tif);

    return expect_int("grayscale total voxels", (long long)stats.total_voxels,
                      8) &&
           expect_int("grayscale void voxels", (long long)stats.void_voxels, 4) &&
           expect_int("grayscale width", stats.width, 2) &&
           expect_int("grayscale height", stats.height, 2) &&
           expect_int("grayscale slices", stats.slices, 2) &&
           expect_int("grayscale has_palette", stats.has_palette, 0);
}

static int test_accumulate_porosity_palette_stack(void) {
    char path[PATH_MAX];
    TIFF *tif;
    StackStats stats;
    const uint8_t voxels[] = {
        0, 1,
        1, 0,
    };

    build_path(path, sizeof(path), "palette_stack.tif");
    create_palette_stack(path, 2, 1, 2, voxels);

    memset(&stats, 0, sizeof(stats));
    tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        return 0;
    }

    accumulate_porosity(tif, &stats);
    TIFFClose(tif);

    return expect_int("palette total voxels", (long long)stats.total_voxels, 4) &&
           expect_int("palette void voxels", (long long)stats.void_voxels, 2) &&
           expect_int("palette has_palette", stats.has_palette, 1);
}

static int test_analyze_stack_returns_expected_counts(void) {
    char path[PATH_MAX];
    const uint8_t voxels[] = {
        0, 0, 1,
        1, 0, 1,
    };
    StackStats stats;

    build_path(path, sizeof(path), "analyze_stack.tif");
    create_grayscale_stack(path, 3, 1, 2, voxels);
    stats = analyze_stack(path);

    return expect_int("analyze total voxels", (long long)stats.total_voxels, 6) &&
           expect_int("analyze void voxels", (long long)stats.void_voxels, 3) &&
           expect_int("analyze width", stats.width, 3) &&
           expect_int("analyze height", stats.height, 1) &&
           expect_int("analyze slices", stats.slices, 2);
}

static int test_compute_porosity_fraction(void) {
    StackStats stats;

    memset(&stats, 0, sizeof(stats));
    stats.total_voxels = 8;
    stats.void_voxels = 2;

    return expect_int("compute_porosity scaled value",
                      (long long)(compute_porosity(stats) * 1000000.0), 250000);
}

static int test_write_stack_summary_emits_expected_text(void) {
    char capture_path[PATH_MAX];
    FILE *tmp = open_capture_file(capture_path, sizeof(capture_path), "summary");
    char output[512];
    StackStats stats;

    if (tmp == NULL) {
        return 0;
    }

    memset(&stats, 0, sizeof(stats));
    stats.total_voxels = 8;
    stats.void_voxels = 2;
    stats.width = 2;
    stats.height = 2;
    stats.slices = 2;

    if (write_stack_summary(tmp, "fixture.tif", stats) != EXIT_SUCCESS) {
        fclose(tmp);
        return 0;
    }

    if (!read_stream(tmp, output, sizeof(output))) {
        fclose(tmp);
        return 0;
    }

    fclose(tmp);

    return expect_contains("summary path", output, "Stack: fixture.tif") &&
           expect_contains("summary dimensions", output, "Dimensions: 2 x 2 x 2") &&
           expect_contains("summary total", output, "Total voxels: 8") &&
           expect_contains("summary void", output, "Void voxels: 2") &&
           expect_contains("summary porosity", output, "Porosity: 0.250000 (25.00%)");
}

static int test_porosity_main_without_args_fails(void) {
    char output[512];
    char *argv[] = {"porosity"};
    int result = 0;

    if (!capture_stdio_call(STDERR_FILENO, porosity_main, 1, argv, output,
                            sizeof(output), &result)) {
        return 0;
    }

    return expect_int("porosity_main no-arg status", result, EXIT_FAILURE) &&
           expect_contains("porosity_main usage", output, "Usage: porosity") &&
           expect_contains("porosity_main binary wording", output,
                           "non-zero voxels in a binary stack");
}

static int test_porosity_main_with_inputs_reports_average(void) {
    char path_a[PATH_MAX];
    char path_b[PATH_MAX];
    char output[2048];
    const uint8_t voxels_a[] = {0, 1, 1, 0};
    const uint8_t voxels_b[] = {0, 0, 1, 0};
    char *argv[] = {"porosity", path_a, path_b};
    int result = 0;

    build_path(path_a, sizeof(path_a), "main_a.tif");
    build_path(path_b, sizeof(path_b), "main_b.tif");
    create_grayscale_stack(path_a, 2, 2, 1, voxels_a);
    create_grayscale_stack(path_b, 2, 2, 1, voxels_b);

    if (!capture_stdio_call(STDOUT_FILENO, porosity_main, 3, argv, output,
                            sizeof(output), &result)) {
        return 0;
    }

    return expect_int("porosity_main status", result, EXIT_SUCCESS) &&
           expect_contains("first stack reported", output, path_a) &&
           expect_contains("second stack reported", output, path_b) &&
           expect_contains("average reported", output,
                           "Average porosity across 2 stacks: 0.375000 (37.50%)");
}

static int test_entry_main_wrapper_executes(void) {
    char path[PATH_MAX];
    char output[1024];
    const uint8_t voxels[] = {0, 1, 1, 1};
    char *argv[] = {"porosity", path};
    int result = 0;

    build_path(path, sizeof(path), "entry_main.tif");
    create_grayscale_stack(path, 2, 2, 1, voxels);

    if (!capture_stdio_call(STDOUT_FILENO, porosity_entry_main, 2, argv, output,
                            sizeof(output), &result)) {
        return 0;
    }

    return expect_int("entry main status", result, EXIT_SUCCESS) &&
           expect_contains("entry main summary", output, "Porosity: 0.750000 (75.00%)");
}

int main(void) {
    const struct {
        const char *name;
        test_fn fn;
    } tests[] = {
        {"die emits message and failure", test_die_emits_message_and_failure},
        {"palette_index_to_binary without palette",
         test_palette_index_to_binary_without_palette},
        {"palette_index_to_binary with palette",
         test_palette_index_to_binary_with_palette},
        {"accumulate_porosity grayscale stack",
         test_accumulate_porosity_grayscale_stack},
        {"accumulate_porosity palette stack",
         test_accumulate_porosity_palette_stack},
        {"analyze_stack expected counts", test_analyze_stack_returns_expected_counts},
        {"compute_porosity fraction", test_compute_porosity_fraction},
        {"write_stack_summary output", test_write_stack_summary_emits_expected_text},
        {"porosity_main without args", test_porosity_main_without_args_fails},
        {"porosity_main with inputs", test_porosity_main_with_inputs_reports_average},
        {"entry main wrapper", test_entry_main_wrapper_executes},
    };
    size_t total = sizeof(tests) / sizeof(tests[0]);
    size_t passed = 0;

    for (size_t i = 0; i < total; ++i) {
        int ok = tests[i].fn();
        printf("[%s] %s\n", ok ? "PASS" : "FAIL", tests[i].name);
        if (ok) {
            passed += 1;
        }
    }

    printf("%zu/%zu tests passed\n", passed, total);
    return passed == total ? EXIT_SUCCESS : EXIT_FAILURE;
}
