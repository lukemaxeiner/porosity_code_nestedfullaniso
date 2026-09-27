#ifndef POROSITY_COVARIANCE_REALIZATION_GENERATION_H
#define POROSITY_COVARIANCE_REALIZATION_GENERATION_H
#include "input.h"
#include "porosity_covariance_realization_prep.h"
#include "porosity_field.h"
#include "porosity_validation.h"

int generate_gaussian_porosity_realization_from_model(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid, double mean_porosity,
    unsigned int seed, PorosityField *field_out);

int generate_gaussian_porosity_realization_from_model_fft(
    const ExponentialCovarianceModel *model,
    const PorosityRealizationGridSpec *simulation_grid,
    double mean_porosity,
    unsigned int seed,
    PorosityField *field_out);
    
int generate_gaussian_porosity_realization_from_prep(
    const PorosityRealizationPrepResult *prep_result, unsigned int seed,
    PorosityField *field_out);

int generate_and_write_gaussian_porosity_realizations(
    const char *stack_path, const PorosityRealizationPrepResult *prep_result,
    uint32_t realization_count, unsigned int base_seed);

int generate_and_write_gaussian_porosity_realizations_from_model(
    const char *stack_path, const FittedCovarianceModelArtifact *model,
    const PorosityRealizationGridSpec *simulation_grid,
    uint32_t realization_count, unsigned int base_seed);

int generate_and_write_gaussian_porosity_realizations_from_model_with_method(
    const char *stack_path,
    const FittedCovarianceModelArtifact *model,
    const PorosityRealizationGridSpec *simulation_grid,
    PorosityRealizationMethod realization_method,
    uint32_t realization_count,
    unsigned int base_seed,
    const PorosityValidationRequest *validation);

#endif