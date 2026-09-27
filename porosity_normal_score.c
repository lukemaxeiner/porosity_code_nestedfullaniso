#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "porosity_normal_score.h"

//Standard normal cdf and its inverse

//CDF of standard normal distribtion using erfc as its convenient in C
//Phi(x) = 0.5 * erfc(-x/sqrt(2))

double normal_cdf(double y) {return 0.5 * erfc(-y / sqrt(2.0)); }

//Peter Acklam's rational approximation of inverse normal CDF
//accuracy order of ~10^-9
//range is split in 3: A central region where (p-0.5) works, and tails of sqrt(-2log(p))

double normal_inverse_cdf(double p){
    //tabulated formula values
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                               -2.759285104469687e+02, 1.383577518672690e+02,
                               -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                               -1.556989798598866e+02, 6.680131188771972e+01,
                               -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                               -2.400758277161838e+00, -2.549732539343734e+00,
                               4.374664141464968e+00,  2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                               2.445134137142996e+00, 3.754408661907416e+00};

    const double p_low = 0.02425; //corresponds to bottom 2.4% of curve, formula value
    const double p_high = 1.0 - p_low;

    double q;
    double r;

    if (!(p>0.0)){
        return -INFINITY;
    }
    if (!(p<1.0)){
        return INFINITY;
    }

    if (p < p_low){
        //lower tail
        q = sqrt(-2.0 * log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / 
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p > p_high){
        //higher tail
        q = sqrt(-2.0 * log(1.0 - p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / 
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    //central region
    q = p - 0.5;
    r = q * q;

    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}


//forward transform. porosity -> normal scores

//necessary struct for data sorting
typedef struct{
    double value;
    size_t index;
} ValueWithIndex;

//aray sorting
static int compare_value_with_index(const void *left, const void *right){
    double a = ((const ValueWithIndex *)left)->value;
    double b = ((const ValueWithIndex *)right)->value;

    return (a > b) - (a < b);
}

int normal_score_transform_field(PorosityField *field, NormalScoreTable *table_out) {
    size_t n;
    size_t i;
    ValueWithIndex *order;

    if (field == NULL || field->values == NULL){
        return EXIT_FAILURE;
    }

    n = (size_t)field->width * (size_t)field->height * (size_t)field->depth;
    if (n == 0U){
        return EXIT_FAILURE;
    }

    //Sort the values, remembering where each came from
    order = (ValueWithIndex *)malloc(n * sizeof(*order));
    if (order == NULL){
        return EXIT_FAILURE;
    }
    for (i = 0U; i < n; ++i){
        if (!isfinite(field->values[i])){
            free(order);
            return EXIT_FAILURE;
        }
        order[i].value = field->values[i];
        order[i].index = i;
    }
    qsort(order, n, sizeof(*order), compare_value_with_index);

    //keep sorted porosity values for back-transform table
    if (table_out != NULL){
        table_out->count = n;
        table_out->sorted_values = (double *)malloc(n*sizeof(double));
        if(table_out->sorted_values == NULL){
            table_out->count = 0U;
            free(order);
            return EXIT_FAILURE;
        }
        for (i = 0U; i < n; ++i){
            table_out->sorted_values[i] = order[i].value;
        }
    }

    //Rank -> percentile -> gaussian score
    //percentile is p = (k + 0.5) / n

    i = 0U;
    while (i < n){
        size_t j = i;
        double p;
        double score;
        size_t k;

        while (j + 1U < n && order[j + 1U].value == order[i].value){
            ++j; //positions i and j hold same value
        }

        p = (0.5 * (double)(i+j) + 0.5) / (double)n;
        score = normal_inverse_cdf(p);

        for (k = i; k <= j; ++k){
            field->values[order[k].index] = score;
        }

        i = j + 1U;
    }

    free(order);
    return EXIT_SUCCESS;
}

//back transform. gaussian values -> porosity

double normal_score_back_transform_value(double y, const NormalScoreTable *table){
    const size_t n = table->count;
    double position;
    size_t lower;
    double t;

    //invert plotting position used in forward step
    position = normal_cdf(y) * (double)n - 0.5;

    //beyond smallest/largest observed block: clmap to data range
    if (position <= 0.0){
        return table->sorted_values[0];
    }
    if (position >= (double)(n - 1U)){
        return table->sorted_values[n - 1U];
    }

    lower = (size_t)position;
    t = position - (double)lower;
    return (1.0 - t) * table->sorted_values[lower] + t * table->sorted_values[lower + 1U];
}

int normal_score_back_transform_field(PorosityField *field, const NormalScoreTable *table){
    size_t n;

    if (field == NULL || field->values == NULL || table == NULL || table->sorted_values == NULL ||
        table->count == 0U){
            return EXIT_FAILURE;
        }
    
    n = (size_t)field->width * (size_t)field->height * (size_t)field->depth;
    for (size_t i = 0U; i < n; ++i){
        field->values[i] = normal_score_back_transform_value(field->values[i], table);
    }

    return EXIT_SUCCESS;
}

//Lookup table file I/O
//builds a file in the build_covariance subfolder like other outputs

char *build_normal_score_table_path(const char *stack_path){
    const char *suffix = "_normal_score_table.txt";
    const char *folder = "/BUILD_COVARIANCE";
    const char *slash;
    const char *file_name;
    const char *extension;
    size_t dir_length;
    size_t base_length;
    char *path;

    if (stack_path == NULL){
        return NULL;
    }

    slash = strrchr(stack_path, '/');
    file_name = slash != NULL ? slash + 1 : stack_path;
    extension = strrchr(file_name, '.');
    dir_length = slash != NULL ? (size_t)(slash - stack_path) : 1U;
    base_length = extension != NULL ? (size_t)(extension - file_name) : strlen(file_name);

    path = (char *)malloc(dir_length + strlen(folder) + base_length + strlen(suffix) + 1U);

    if (path == NULL){
        return NULL;
    }

    if (slash != NULL){
        memcpy(path, stack_path, dir_length);
    }
    else{
        path[0] = '.';
    }
    memcpy(path + dir_length, folder, strlen(folder));
    memcpy(path + dir_length + strlen(folder), file_name, base_length);
    memcpy(path + dir_length + strlen(folder) + base_length, suffix, strlen(suffix) + 1U);

    return path;
}

int write_normal_score_table(const char *path, const NormalScoreTable *table) {
    FILE *file;
 
    if (path == NULL || table == NULL || table->sorted_values == NULL) {
        return EXIT_FAILURE;
    }
 
    file = fopen(path, "w");
    if (file == NULL) {
        return EXIT_FAILURE;
    }
 
    if (fprintf(file, "NORMAL_SCORE_TABLE_COUNT = %zu\n", table->count) < 0) {
        fclose(file);
        return EXIT_FAILURE;
    }
    for (size_t i = 0U; i < table->count; ++i) {
        if (fprintf(file, "%.17g\n", table->sorted_values[i]) < 0) {
            fclose(file);
            return EXIT_FAILURE;
        }
    }
 
    return fclose(file) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
 
int read_normal_score_table(const char *path, NormalScoreTable *table_out) {
    FILE *file;
    size_t count;
 
    if (path == NULL || table_out == NULL) {
        return EXIT_FAILURE;
    }
    memset(table_out, 0, sizeof(*table_out));
 
    file = fopen(path, "r");
    if (file == NULL) {
        return EXIT_FAILURE;
    }
 
    if (fscanf(file, " NORMAL_SCORE_TABLE_COUNT = %zu", &count) != 1 ||
        count == 0U) {
        fclose(file);
        return EXIT_FAILURE;
    }
 
    table_out->sorted_values = (double *)malloc(count * sizeof(double));
    if (table_out->sorted_values == NULL) {
        fclose(file);
        return EXIT_FAILURE;
    }
 
    for (size_t i = 0U; i < count; ++i) {
        if (fscanf(file, " %lf", &table_out->sorted_values[i]) != 1 ||
            (i > 0U &&
             table_out->sorted_values[i] < table_out->sorted_values[i - 1U])) {
            /* missing value, or not ascending: the table is corrupt */
            fclose(file);
            free_normal_score_table(table_out);
            return EXIT_FAILURE;
        }
    }
 
    fclose(file);
    table_out->count = count;
    return EXIT_SUCCESS;
}
 
void free_normal_score_table(NormalScoreTable *table) {
    if (table == NULL) {
        return;
    }
    free(table->sorted_values);
    table->sorted_values = NULL;
    table->count = 0U;
}


