import os
import sys
import numpy as np
from numpy.testing import assert_allclose
import pandas as pd

sys.path.append(
    os.path.abspath(
        "d:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code/Rt_estimate_reconstruction"
    )
)

input_folder = "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code/Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/"
data_source = "_rtlive"


def tests_positive_infected():
    # Check that the time series for
    # infected individuals is always weakly positive
    df = pd.read_csv("{}dataset{}.csv".format(input_folder, data_source))
    for days_infectious in range(5, 10 + 1):
        if days_infectious == 5:
            mask = df["infected_{}".format(days_infectious)] < 0.0
        else:
            mask_temp = df["infected_{}".format(days_infectious)] < 0.0
            mask = mask | mask_temp
    test = df.loc[
        mask,
    ].shape[0]
    print(
        df.loc[
            mask,
        ]
    )
    assert_allclose(0, test)


def test_new_cases():
    # Check that new cases are always
    # weakly positive
    df = pd.read_csv("{}dataset{}.csv".format(input_folder, data_source))
    mask = df["new_cases"] < 0
    test = df.loc[
        mask,
    ].shape[0]
    print(
        df.loc[
            mask,
        ]
    )
    assert_allclose(0, test)


def tests_growth_rates():
    # Check that the growth rate of infected
    # individuals is equal to (-gamma) when
    # new cases are zero
    df = pd.read_csv("{}dataset{}.csv".format(input_folder, data_source))
    for days_infectious in range(5, 10 + 1):
        mask = (df["new_cases"] == 0.0) & ~np.isnan(
            (df["gr_infected_{}".format(days_infectious)])
        )
        assert_allclose(
            df.loc[mask, "gr_infected_{}".format(days_infectious)],
            -1 / float(days_infectious),
        )


def tests_growth_rates_2():
    # Check that the growth rate of infected
    # individuals is bounded below by gamma
    # when new cases are positive
    df = pd.read_csv("{}dataset{}.csv".format(input_folder, data_source))
    for days_infectious in range(5, 10 + 1):
        mask = (df["new_cases"] > 0.0) & df[
            "gr_infected_{}".format(days_infectious)
        ] <= -1 / float(days_infectious)
        assert_allclose(
            df.loc[
                mask,
            ].shape[0],
            0,
        )


def tests_growth_rates_zero_new_cases():
    # Additional check for growth rates
    # with zero new cases
    df = pd.read_csv("{}dataset{}.csv".format(input_folder, data_source))
    # mask = df["new_cases"] == 0.0
    for days_infectious in range(5, 10 + 1):
        df["temp"] = (
            df["infected_{}".format(days_infectious)]
            / df.groupby("Country/Region").shift()[
                "infected_{}".format(days_infectious)
            ]
            - 1
        )
        df_temp = df.loc[:, ["temp", "gr_infected_{}".format(days_infectious)]].dropna()
        assert_allclose(
            df_temp["temp"], df_temp["gr_infected_{}".format(days_infectious)]
        )


print("1")
tests_positive_infected()
print("2")
test_new_cases()
print("3")
tests_growth_rates()
print("4")
tests_growth_rates_2()
print("5")
tests_growth_rates_zero_new_cases()
