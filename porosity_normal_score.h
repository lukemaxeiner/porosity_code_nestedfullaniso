#ifndef POROSITY_NORMAL_SCORE_H
#define POROSITY_NORMAL_SCORE_H

#include <stddef.h>
#include <stdio.h>

#include "porosity_field.h"

typedef struct{
    size_t count; //number of porosity blocks in source field
    double *sorted_values; //porosity values
} NormalScoreTable;

//standard normal CDF
double normal_cdf(double y);

//inverse standard normal CDF
double normal_inverse_cdf(double p);

//replace every value in field with its normal score
int normal_score_transform_field(PorosityField *field, NormalScoreTable *table_out);

//map one gaussian value back to porosity through lookup table
double normal_score_back_transform_value(double y, const NormalScoreTable *table);

//make every value in field back to porosity in place
int normal_score_back_transform_field(PorosityField *field, const NormalScoreTable *table);

// I/O for file output
char *build_normal_score_table_path(const char *stack_path);

int write_normal_score_table(const char *path, const NormalScoreTable *table);
int read_normal_score_table(const char *path, NormalScoreTable *table_out);
void free_normal_score_table(NormalScoreTable *table);


#endif