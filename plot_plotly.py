#!/usr/bin/env python3

import plotly as py
import plotly.graph_objs as go
import numpy as np
import pandas
from scipy.stats import norm
import plotly.tools as tls


def get_csv_col(filepath, col_id):
    return np.array(pandas.read_csv(filepath, usecols=[col_id], delimiter=',')).transpose()[0]


def get_data_from_files(files, col_id):
    if files:
        data = []
        for file in files:
            data.append(get_csv_col(file, col_id))
        if len(data) == 1:
            return data[0]
        else:
            x = np.concatenate((data[0], data[1]))
            for d in data[2:]:
                x = np.concatenate((x, d))
            return x
    else:
        return None


def histogram(files, col_id, title):
    x = get_data_from_files(files, col_id)
    data = [go.Histogram(x=x)]
    layout = go.Layout(title=title)
    fig = go.Figure(data=data, layout=layout)
    py.offline.plot(fig, filename='histogram.html')


def histogram2(files):
    data_sets = []

    data_pp = get_data_from_files(files, 'distPP')
    x = get_data_from_files(files, 'distCP')
    # ['dist_P-P', 'dist_C-P']
    print(len(x))

    # for col_id in col_ids:
    #     print(col_id)
    #     data_sets.append(get_data_from_files(files, col_id))
    # print(len(data_sets))
    # col_nr = len(col_ids)
    # print(col_nr)

    # ## optional 2nd and 3rd column concatenate
    # data_sets = [data_sets[0], np.concatenate(data_sets[1], data_sets[2])]
    # col_nr -= 1
    #
    # fig = tls.make_subplots(rows=col_nr, cols=1, shared_xaxes=True)
    # print(fig.to_string())


histogram2(['Kyu_dists.csv', 'Ding_dists.csv'])
