#ifndef POROSITY_COVARIANCE_TAIL_NOISE_H
#define POROSITY_COVARIANCE_TAIL_NOISE_H

#include "porosity_covariance_realization_prep.h"

int check_covariance_tail_noise(const PorosityCovarianceSample *samples,
                                const CovarianceGrid3D *grid,
                                double filter_resolution,
                                CovarianceTailNoiseSummary *summary_out);

#endif
