#!/usr/bin/env python
# coding: utf-8


"""
functions for bootstrap sampling for chi-square 
https://github.com/vishanth10/resamp

MIT License

Copyright (c) 2023 vishanth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import logging
import numpy as np
import pandas as pd
import scipy.stats as stats
from types import SimpleNamespace



def calculate_chi_square(observed, expected):
    chi_square = ((observed - expected) ** 2 / expected).sum().sum()
    return chi_square
    
    
def calculate_chi_abs(observed, expected=None):
    try:
        # Check if the input data is DataFrame and convert to numpy arrays if necessary
        if isinstance(observed, pd.DataFrame):
            observed = observed.values
        if expected is None:
            expected = calculate_expected(observed)
        elif isinstance(expected, pd.DataFrame):
            expected = expected.values

        # Calculate the chi absolute statistic
        chi_abs = (np.abs(observed - expected) / expected).sum().sum()
        return chi_abs
    except Exception as e:
        logging.error("Error calculating p-value: ", exc_info=True)
        return None


def chi_abs_stat(observed_data, expected_data=None):
    """
    Compute chi-abs statistic

    Input:
        observed data (array or data frame): Observed counts
        expected data (array or data frame, optional): Expected counts. If not provided,
        the function will calculate an expected table that assumes equal probabilities across columns.

    Output:
        chi-abs statistic
    """

    try:
        # If expected data is not provided, calculate it
        if expected_data is None:
            expected_data = calculate_expected(observed_data)

        # Ensure dimensions match
        if observed_data.shape != expected_data.shape:
            raise ValueError("Dimensions of observed and expected data do not match")

        chi_abs = calculate_chi_abs(observed_data, expected_data)
        return chi_abs
    except Exception as e:
        logging.error("Error calculating chi absolute: ", exc_info=True)
        return None




"""
Calculate residuals and standardised residuals
https://stackoverflow.com/a/20457483
https://creativecommons.org/licenses/by-sa/4.0/
"""

def residuals(observed, expected):
    return (observed - expected) / np.sqrt(expected)


def stdres(observed, expected):
    n = observed.sum()
    rsum, csum = stats.contingency.margins(observed)
    # With integers, the calculation
    #     csum * rsum * (n - rsum) * (n - csum)
    # might overflow, so convert rsum and csum to floating point.
    rsum = rsum.astype(np.float64)
    csum = csum.astype(np.float64)
    v = csum * rsum * (n - rsum) * (n - csum) / n**3
    return (observed - expected) / np.sqrt(v)




def bootstrap_chi_abs(observed_data, num_simulations=10000, with_replacement=True):
    
    """
    Generates a bootstrap distribution of the chi absolute statistic for an n*n contingency table.

    Parameters:
        observed_data (np.array or pd.DataFrame): n*n contingency table with observed frequencies.
        num_simulations (int): Number of bootstrap samples to generate.
        with_replacement (bool): Indicates whether sampling should be with replacement.
    """
    
    if isinstance(observed_data, pd.DataFrame):
        observed_data = observed_data.values

    total_rows, total_columns = observed_data.shape
    expected_data = calculate_expected(observed_data)

    results = np.zeros(num_simulations)
    total_counts_per_column = observed_data.sum(axis=0)

    # Create a pooled data array combining all categories across rows and columns
    pooled_data = np.concatenate([np.repeat(row, sum(observed_data[row, :])) for row in range(total_rows)])

    np.random.seed(2020)
    for i in range(num_simulations):
        sim = np.zeros_like(observed_data)

        for col in range(total_columns):
            column_sample = np.random.choice(pooled_data, total_counts_per_column[col], replace=with_replacement)

            for row in range(total_rows):
                # Count occurrences of each category in the column sample
                sim[row, col] = np.sum(column_sample == row)
        #print(sim)
        # Calculate the chi absolute statistic for the simulated data
        chi_abs = (np.abs(sim - expected_data) / expected_data).sum().sum()
        results[i] = chi_abs

    return results


def calculate_expected(observed):
    row_sums = observed.sum(axis=1)
    col_sums = observed.sum(axis=0)
    total = observed.sum().sum()
    expected = np.outer(row_sums, col_sums) / total
    return expected


def calculate_p_value_bootstrap(observed_data, simulated_data, two_tailed=False):
    """
    Calculates the p-value for the chi absolute statistic using bootstrap methods, 
    determining first if the observed statistic lies on the left or right side of the distribution's mean.

    Parameters:
        observed_data (float): The observed chi absolute statistic.
        simulated_data (np.array): The array of chi absolute statistics from bootstrap samples.
        two_tailed (bool): If True, perform a two-tailed test. Defaults to False (one-tailed test).

    Returns:
        float: The p-value.
    """
    try:
        # Determine the side of the distribution where the observed data lies
        mean_simulated_data = np.mean(simulated_data)
        is_right_side = observed_data > mean_simulated_data
        
        if two_tailed:
            if is_right_side:
                # For a two-tailed test, consider both tails of the distribution (right side logic)
                tail_proportion = np.mean(simulated_data >= observed_data)
            else:
                # For a two-tailed test, consider both tails of the distribution (left side logic)
                tail_proportion = np.mean(simulated_data <= observed_data)
            p_value = tail_proportion
        else:
            if is_right_side:
                # For a one-tailed test, only consider the tail of interest (right side logic)
                p_value = np.mean(simulated_data >= observed_data)
            else:
                # For a one-tailed test, only consider the tail of interest (left side logic)
                p_value = np.mean(simulated_data <= observed_data)
        
        return p_value
    except Exception as e:
        logging.error("Error in calculating p-value: ", exc_info=True)
        return None




def chi_square_posthoc(observed, num_simulations=1000):
    """
    Perform chi-square tests and bootstrap sampling using resamp functions above.
    Perform chi-square post-hoc tests according to Beasley and Schumacker (1995).
    
    Author: CS Lim
    
    Input:
        * observed: observed values as contingency table (required)
        * num_simulations: number of simulations (default=1000)

    Output:
        * ChiSquareResult: chi-square statistic, p-value, bootstrap p-value, 
            post-hoc statistic, posthoc p-value
    """
    
    expected = calculate_expected(observed)
    chi_square = calculate_chi_square(observed, expected)
    chi_abs = chi_abs_stat(observed)
    bootstrap_results = bootstrap_chi_abs(observed, num_simulations=num_simulations)
    
    # post-hoc test for 2x3 contingency table
    # convert adjusted standardised residuals (z-scores) to chi-square statisticss
    posthoc = stdres(observed, expected)**2
    
    # convert chi-square statistics to p-values, take the first 3 values (the remaining 3 values are identical)
    # https://stackoverflow.com/a/34965620
    posthoc_pval = [stats.chi2.sf(i, 1) for i in posthoc][0]
    
    # chi-square test p-value
    chi, pval, _, _ = stats.chi2_contingency(observed)
    # bootstrap p-value
    bootstrap_pval = calculate_p_value_bootstrap(chi_abs,bootstrap_results)
    
    ChiSquareResult = SimpleNamespace(statistic=chi, pvalue=pval, bootstrap_pvalue=bootstrap_pval, 
                                      posthoc_statistic=posthoc[0], posthoc_pvalue=posthoc_pval)

    return ChiSquareResult