#!/usr/bin/env python
# coding: utf-8

"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2025-02-15 20:55:07
@desc        Main RIBOSS module
"""




import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import chi2_contingency
from types import SimpleNamespace
from statsmodels.stats.multitest import multipletests




def calculate_expected(observed):
    row_sums = observed.sum(axis=1)
    col_sums = observed.sum(axis=0)
    total = observed.sum()
    expected = np.outer(row_sums, col_sums) / total
    return expected


def bootstrap_chi_square(observed_data, num_simulations=10000):
    """
    Bootstrap Chi-squared statistic.
    """
    
    if isinstance(observed_data, pd.DataFrame):
        observed_data = observed_data.values

    expected_data = calculate_expected(observed_data)
    n_rows, n_cols = observed_data.shape
    bootstrap_stats = []

    for _ in range(num_simulations):
        bootstrap_sample = np.zeros_like(observed_data)
        for i in range(n_rows):
            row_total = observed_data[i, :].sum()
            row_probs = expected_data[i, :] / expected_data[i, :].sum()
            bootstrap_sample[i, :] = np.random.multinomial(row_total, row_probs)

        bootstrap_stats.append(chi2_contingency(bootstrap_sample)[0])  # Chi-squared statistic

    return np.array(bootstrap_stats)



def chi_square_posthoc(obs, num_simulations=1000, subset_indices=([0, 1], [0, 0])):
    """
    Chi-square test with bootstrap and post-hoc analysis (optional subset).
    """
    
    observed = obs + (obs == 0) * 1e-10
    chi2, pval, _, expected = chi2_contingency(observed)
    bootstrap_stats = bootstrap_chi_square(observed, num_simulations)
    bootstrap_pval = np.mean(bootstrap_stats >= chi2)

    # Post-hoc analysis (adjusted standardized residuals)
    row_totals = observed.sum(axis=1, keepdims=True)
    col_totals = observed.sum(axis=0, keepdims=True)
    grand_total = observed.sum()

    adjusted_residuals = (observed - expected) / np.sqrt(expected * (1 - row_totals / grand_total) * (1 - col_totals / grand_total))

    if subset_indices is not None:  # Handle subset
        adjusted_residuals_subset = adjusted_residuals[subset_indices[0], subset_indices[1]]
        posthoc_pvals = 2 * stats.norm.sf(np.abs(adjusted_residuals_subset))
        n_comparisons = len(adjusted_residuals_subset)  # Correct number of comparisons
    else:  # Full table
        posthoc_pvals = 2 * stats.norm.sf(np.abs(adjusted_residuals))
        n_comparisons = observed.size

    fdr_pvals = multipletests(posthoc_pvals.flatten(), method='fdr_bh')[1].reshape(posthoc_pvals.shape)

    ChiSquareResult = SimpleNamespace(
        statistic=chi2,
        pvalue=pval,
        bootstrap_pvalue=bootstrap_pval,
        adjusted_residuals=adjusted_residuals,
        posthoc_pvalues=fdr_pvals,
    )

    return ChiSquareResult




def calculate_g_statistic(observed, expected):
    """
    Calculates the G-statistic (log-likelihood ratio).
    """
    
    observed = observed.astype(float)
    expected = expected.astype(float)
    with np.errstate(divide='ignore', invalid='ignore'): # Suppress warnings locally
        g_values = 2 * observed * np.log(observed / expected)
    g_values = np.nan_to_num(g_values) # Convert nan to 0
    return np.sum(g_values)



def bootstrap_g_statistic(observed_data, num_simulations=10000):
    """
    Bootstrap G-statistic.
    """
    
    if isinstance(observed_data, pd.DataFrame):
        observed_data = observed_data.values

    expected_data = calculate_expected(observed_data)
    n_rows, n_cols = observed_data.shape
    bootstrap_stats = []

    for _ in range(num_simulations):
        bootstrap_sample = np.zeros_like(observed_data)
        for i in range(n_rows):
            row_total = observed_data[i, :].sum()
            row_probs = expected_data[i, :] / expected_data[i, :].sum()
            bootstrap_sample[i, :] = np.random.multinomial(row_total, row_probs)

        bootstrap_stats.append(calculate_g_statistic(bootstrap_sample, expected_data)) #Calculate G-statistic

    return np.array(bootstrap_stats)


def g_test_posthoc(obs, num_simulations=1000, subset_indices=([0, 1], [0, 0])):
    """
    G-test with bootstrap and post-hoc analysis (optional subset).
    """

    observed = obs + (obs == 0) * 1e-10
    g = calculate_g_statistic(observed, calculate_expected(observed))
    pval = chi2_contingency(observed, correction=False, lambda_="log-likelihood")[1]
    bootstrap_stats = bootstrap_g_statistic(observed, num_simulations)
    bootstrap_pval = np.mean(bootstrap_stats >= g)

    # Post-hoc analysis (adjusted standardized residuals)
    expected = calculate_expected(observed)
    row_totals = observed.sum(axis=1, keepdims=True)
    col_totals = observed.sum(axis=0, keepdims=True)
    grand_total = observed.sum()

    adjusted_residuals = (observed - expected) / np.sqrt(expected * (1 - row_totals / grand_total) * (1 - col_totals / grand_total))

    if subset_indices is not None:  # Handle subset
        adjusted_residuals_subset = adjusted_residuals[subset_indices[0], subset_indices[1]]
        posthoc_pvals = 2 * stats.norm.sf(np.abs(adjusted_residuals_subset))
        n_comparisons = len(adjusted_residuals_subset)  # Correct number of comparisons
    else:  # Full table
        posthoc_pvals = 2 * stats.norm.sf(np.abs(adjusted_residuals))
        n_comparisons = observed.size

    fdr_pvals = multipletests(posthoc_pvals.flatten(), method='fdr_bh')[1].reshape(posthoc_pvals.shape)

    ChiSquareResult = SimpleNamespace(
        statistic=g,
        pvalue=pval,
        bootstrap_pvalue=bootstrap_pval,
        adjusted_residuals=adjusted_residuals, 
        posthoc_pvalues=fdr_pvals,
    )

    return ChiSquareResult

