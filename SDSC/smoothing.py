import numpy as np
import statsmodels.api as sm

# import pandas as pd
# import numpy.matlib as mp
# import scipy as sp
# import statsmodels


def two_step_STL(x, start_p=None, H=7, robust=True):
    # x should contain daily cases
    # start_p -for partial smoothing starting from the observation with the index start_p
    # H horizon
    # two step procedure: 1) outlier detection with robust STL, which are given by residual+seasonality of STL on non-log daily
    #                    2) non-robust STL
    z = x.copy()
    if start_p is None:
        start_p = np.where(np.cumsum(x) == 0)[0]
        if len(start_p) == 0:
            start_p = 0
        elif start_p[0] + H > len(x):
            return z
        else:
            start_p = start_p[0]

    if robust is True:
        stl_log_daily = sm.tsa.seasonal.STL(
            z[start_p:], robust=True, seasonal=H, period=H, trend=2 * H + 1
        ).fit()
        trend_no_outl = stl_log_daily.trend
        trend_no_outl = np.where(trend_no_outl > 0, trend_no_outl, 0)

    else:
        trend_no_outl = z[start_p:]

    stl_no_outliers = sm.tsa.seasonal.STL(
        trend_no_outl, seasonal=H, period=H, trend=2 * H + 1
    ).fit()
    # second more smooth non robust application
    smoothed = np.where(stl_no_outliers.trend > 0, stl_no_outliers.trend, 0)

    z[start_p:] = smoothed
    return z


def piecewise_STL(x, Ws=3, H=7):
    # x contains cumulative cases
    # piecewise STL: first the history is divided into the overlapping subhistories, each of 2*len_piece size
    # next STL is applied to the i-th one and i-1-the, in the overlapping region the smoothing result is the convex combination of two
    # smoothings with the sigmoid weights

    z = np.diff(x.copy(), axis=0)  # will be changed during the rescaling
    z0 = z.copy()  # daily cases
    len_piece = 3 * H

    # imputations switch
    imputations, thresh = True, 0.2
    if imputations:
        a = np.where(np.diff(x) > 0)[0]
        for i, j in enumerate(a):
            if (i >= 1) & (j > H):
                temp = a[i] - a[i - 1]
                if (temp >= 2) & (np.mean(z0[j - H - 1 : j - 1]) > -np.log(thresh)):
                    z[a[i - 1] + 1 : a[i]] = round(z0[j] / temp)
                    z[j] = z0[j] - round(z0[j] / temp) * (temp - 1)

    # for the small number of observations use the binomial filter
    # (maybe to change also to the number of non-zeros in the last n weeks)
    int_lims = np.unique([0] + list(np.sort(np.arange(len(z), -1, -len_piece))))
    ind = np.arange(len_piece)
    weights = 1.0 / (1.0 + np.exp(10.549 * (ind / len_piece - 0.518)))
    if len(int_lims) > 3:
        smoothed_z = z.copy()
        # last sub-interval
        next_sub_interval = z[int_lims[-3] :].copy()
        smoothed_2 = two_step_STL(next_sub_interval, 0, robust=True)
        smoothed_z[int_lims[-2] :] = smoothed_2[len_piece:].copy()
        # the (-excess) of observations, which due to outlier and seasonality detection
        outl_last_int = np.sum(z0[-(len_piece + 1) :]) - np.sum(smoothed_2[len_piece:])

        # if the excess is positive, scale backwards
        if (
            (outl_last_int > 0)
            & (np.sum(z[int_lims[-2] :]) > 0)
            & (np.count_nonzero(z[: int_lims[-2]]) > H)
        ):

            scaler = np.sum(z0[int_lims[-2] :]) / np.sum(smoothed_z[int_lims[-2] :])
            smoothed_z[int_lims[-2] :] = smoothed_z[int_lims[-2] :] * scaler
            smoothed_2 = smoothed_2 * scaler
        # if excess is negative, scale forward

        if (outl_last_int < 0) & (np.sum(smoothed_z[int_lims[-2] :]) > 0):

            scaler = np.sum(z0[int_lims[-2] :]) / np.sum(smoothed_z[int_lims[-2] :])
            smoothed_z[int_lims[-2] :] = smoothed_z[int_lims[-2] :] * scaler
            smoothed_2 = smoothed_2 * scaler

        # repreat backwards in subintervals
        for i in range(len(int_lims) - 4, -1, -1):
            first_sub_interval = z[int_lims[i] : int_lims[i + 2]].copy()
            smoothed_1 = two_step_STL(first_sub_interval, 0, robust=False)
            smoothed_z[int_lims[i + 1] : int_lims[i + 2]] = smoothed_2[:len_piece] * (
                1 - weights
            ) + smoothed_1[-len_piece:] * (weights)
            smoothed_2 = smoothed_1.copy()
            outl_last_int = np.sum(z0[int_lims[i + 1] :]) - np.sum(
                smoothed_z[int_lims[i + 1] :]
            )

            if (outl_last_int > 0) & (np.sum(z[: int_lims[i + 1]]) > 0):
                z[: int_lims[i + 1]] = (
                    z[: int_lims[i + 1]]
                    * (outl_last_int + np.sum(z0[: int_lims[i + 1]]))
                    / (np.sum(z[: int_lims[i + 1]]))
                )

            if (outl_last_int < 0) & (np.sum(smoothed_z[int_lims[i + 1] :]) > 0):
                smoothed_z[int_lims[i + 1] :] = (
                    smoothed_z[int_lims[i + 1] :]
                    * np.sum(z0[int_lims[i + 1] :])
                    / (np.sum(smoothed_z[int_lims[i + 1] :]))
                )

        smoothed_z[: int_lims[1]] = smoothed_2[: int_lims[1]]
    else:
        smoothed_z = two_step_STL(z, 0)

    cumsum_sm = np.cumsum(list([x[0]]) + list(smoothed_z))

    # final scaling
    if cumsum_sm[-1] > 0:
        cumsum_sm = cumsum_sm * x[-1] / cumsum_sm[-1]

    return cumsum_sm
