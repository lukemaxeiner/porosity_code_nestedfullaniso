#include <stdio.h>
#include <stdlib.h>

#include "porosity_global.h"
#include "tiff.h"

int process_global_porosity(const char *stack_path, double *porosity_sum) {
    StackStats stats;
    double porosity;

    if (stack_path == NULL || porosity_sum == NULL) {
        return EXIT_FAILURE;
    }

    stats = analyze_stack(stack_path);
    porosity = compute_porosity(stats);

    if (write_stack_summary(stdout, stack_path, stats) != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to write summary for %s\n", stack_path);
        return EXIT_FAILURE;
    }

    *porosity_sum += porosity;
    return EXIT_SUCCESS;
}
