"""Epidemic modeling tutorials by Reza Sameni, Emory University."""

from .models import seirp, seirp_saturated_resource, si_controlled, si_alpha_controlled
from .growth import (
    rt_exp_fit_gen_ratios,
    rt_exp_fit_log_lin_reg,
    rt_exp_fit_nonlin_ls,
    exp_model,
)
from .kalman import (
    generic_extended_kalman_filter,
    si_alpha_model_ekf,
    si_alpha_model_ekf_opt_controlled,
    si_alpha_model_backward_ekf,
    si_alpha_model_backward_ekf_opt_controlled,
    new_case_ekf_estimator_with_optimal_npi,
    rt_exp_fit_ekf,
)
from .data import read_covid19_data, read_oxford_data, read_geo_table, prepare_cases

from .npi import (
    npi_cost,
    fit_npi_model,
    forecast_npi,
    optimal_npi,
    train_npi_prescriptor,
    prescribe_npi,
    train_predict_prescribe_npi,
    forecast_quality_assessment,
    nonnegative_least_squares,
    default_si_params,
    NPI_COLUMNS,
    NPI_MAXES,
)
from .layers import exp_layer, my_tanh_layer, torch_exp_layer, torch_my_tanh_layer
from .spatial import diffusion_2d, population_motion_2d

__version__ = "2.0.0"
