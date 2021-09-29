import os
import sys
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.dirname(__file__) + "/" + "../.."))


def construct_dataset(file_name, var_name):
    """Convenience function for constructing
    a clean Pandas dataframe from the CSV
    files provided by JH CSSE on their Github
    repo

    Args:
        file_name (str): File name / URL of CSV file
        var_name (name): Variable name

    Returns:
        df: Dataframe
    """
    df = pd.read_csv(file_name)
    del df["Lat"], df["Long"]

    # Melt to long format
    df = pd.melt(
        df,
        id_vars=["Province/State", "Country/Region"],
        value_vars=list(df.columns.values[2:]),
    )
    df.rename(columns={"variable": "Date", "value": var_name}, inplace=True)

    # For some countries, data are reported
    # by regions / states; aggregate to country level
    return df.groupby(["Country/Region", "Date"]).sum().reset_index()


##############
# Parameters #
##############

# methods = ["", "_ETH", "_ilmenau", "_RKI", "_sdsc", "_epiforecasts"]
methods = ["", "_RKI"]
# input_folder = "./relevant_scripts_adjusted/"
input_folder = "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code/Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/"
output_folder = input_folder
min_cases = 100
days_infectious_list = [5, 6, 7, 8, 9, 10]  # Values of (1 / gamma) used in constructing
# time series of infected individuals
# end_date             = '2020-05-06'        # End of sample
restrict_end_sample = True

for method in methods:
    print(
        "Construct dataset with data from ",
        method[1:] if method != "" else "original method",
    )

    #####################
    # Construct dataset #
    #####################

    # Read in data on total cases
    df = construct_dataset(
        file_name="{}time_series_covid19_confirmed_global{}.csv".format(
            input_folder, method
        ),
        var_name="total_cases",
    )

    # Only consider days after a minimum
    # number of total cases has been reached
    mask = df["total_cases"] >= min_cases
    df = df.loc[
        mask,
    ]

    # Clean up the dataframe
    df["Date"] = pd.to_datetime(df["Date"])
    df.reset_index(inplace=True)
    del df["index"]

    # Sort by date
    df.sort_values(by=["Country/Region", "Date"], ascending=True, inplace=True)

    # Construct derived flow variables (only new cases, removed recoveries / deaths)
    for var_name in ["cases"]:
        df["new_{}".format(var_name)] = (
            df["total_{}".format(var_name)]
            - df.groupby("Country/Region").shift()["total_{}".format(var_name)]
        )

    # Construct number of infected
    for days_infectious in days_infectious_list:
        df["infected_{}".format(days_infectious)] = np.nan
        for country in df["Country/Region"].unique():
            mask = df["Country/Region"] == country
            df_country = (
                df.loc[
                    mask,
                ]
                .copy()
                .reset_index()
            )
            T = df_country.shape[0]

            # Initialize number of infected
            infected = np.zeros(T) * np.nan
            infected[0] = df_country["total_cases"][0]

            # Main loop
            for tt in range(1, T):
                gamma = 1 / float(days_infectious)

                # Calculate number of infected recursively;
                # In the JH CSSE dataset, there are some
                # data problems whereby new cases are occasionally
                # reported to be negative; in these case, take zero
                # when constructing time series for # of invected,
                # and then change values to NaN's later on
                infected[tt] = (1 - gamma) * infected[tt - 1] + np.maximum(
                    df_country["new_cases"][tt], 0.0
                )
            df.loc[mask, "infected_{}".format(days_infectious)] = infected

    # In the original JH CSSE dataset, there are
    # some inconsistencies in the data
    # Replace with NaN's in these cases
    mask = df["new_cases"] < 0
    df.loc[mask, "new_cases"] = np.nan
    print(
        "     Inconsistent observations in new_cases in JH CSSE dataset: {:}".format(
            mask.sum()
        )
    )
    for days_infectious in days_infectious_list:
        df.loc[mask, "infected_{}".format(days_infectious)] = np.nan

    # Calculate growth rate of infected
    for days_infectious in days_infectious_list:
        df["gr_infected_{}".format(days_infectious)] = (
            df["infected_{}".format(days_infectious)]
            / df.groupby("Country/Region").shift(1)[
                "infected_{}".format(days_infectious)
            ]
        ) - 1
        mask = (
            df.groupby("Country/Region").shift(1)["infected_{}".format(days_infectious)]
            == 0.0
        )
        df.loc[mask, "gr_infected_{}".format(days_infectious)] = np.nan

    # Deal with potential consecutive zeros in the number of infected
    for days_infectious in days_infectious_list:
        mask = (df["infected_{}".format(days_infectious)] == 0.0) & (
            df.groupby("Country/Region").shift(1)["infected_{}".format(days_infectious)]
            == 0.0
        )
        df.loc[mask, "gr_infected_{}".format(days_infectious)] = -(1 / days_infectious)
        if mask.sum() > 0:
            print(
                "     Number of observations with zero infected (with {} infectious days) over two consecutive days: {:}".format(
                    days_infectious, mask.sum()
                )
            )

    # Set to NaN observations with very small
    # number of cases but very high growth rates
    # to avoid these observations acting as
    # large outliers
    for days_infectious in days_infectious_list:
        gamma = 1 / float(days_infectious)
        mask = (df["new_cases"] <= 25) & (
            df["gr_infected_{}".format(days_infectious)] >= gamma * (5 - 1)
        )  # Implicit upper bound on R
        df.loc[
            mask,
            [
                "infected_{}".format(days_infectious),
                "gr_infected_{}".format(days_infectious),
            ],
        ] = np.nan

    # Set to NaN observations implausibly
    # high growth rates that are likely
    # due to data issues
    for days_infectious in days_infectious_list:
        gamma = 1 / float(days_infectious)
        mask = df["gr_infected_{}".format(days_infectious)] >= gamma * (
            10 - 1
        )  # Implicit upper bound on R
        df.loc[
            mask,
            [
                "infected_{}".format(days_infectious),
                "gr_infected_{}".format(days_infectious),
            ],
        ] = np.nan

    # Remove initial NaN values for growth rates
    for country in df["Country/Region"].unique():
        mask = df["Country/Region"] == country
        T = df.loc[
            mask,
        ].shape[0]
        df.loc[mask, "days_since_min_cases"] = range(T)
    mask = df["days_since_min_cases"] >= 1
    df = df.loc[
        mask,
    ]
    del df["days_since_min_cases"]

    # Save final dataset
    df.to_csv("{}/dataset{}.csv".format(output_folder, method), index=False)
