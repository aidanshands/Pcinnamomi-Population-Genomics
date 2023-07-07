#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
from itertools import combinations
Data = pd.read_csv("Perc_Diseased_Average5WPI_for_AUDPC_for_Py.csv")
def calculate_random_pair_means(df, seed=None):
    np.random.seed(seed)
    grouped = df.groupby('Isolate')
    result_df = pd.DataFrame(columns=['Isolate', '5WPI_Random_Pair_Mean'])
    for treatment, group in grouped:
        if len(group) >= 2:
            pairs = list(combinations(group['WPI5'], 2))
            selected_pairs = np.random.choice(len(pairs), size=3, replace=False)
            pair_means = [np.mean(pair) for pair in np.array(pairs)[selected_pairs]]
            pair_means_df = pd.DataFrame({'Isolate': treatment, '5WPI_Random_Pair_Mean': pair_means})
            result_df = pd.concat([result_df, pair_means_df], ignore_index=True)
    return result_df

result = calculate_random_pair_means(Data, seed=42)
result.to_csv("5WPI_Random_Pair_Means_for_AUDPC.csv", index = None)
