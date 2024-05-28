        if Z_low == Z_high:
            x = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zindex['AGB+SNcc']['low']]].columns.to_numpy().astype(np.float64)
            y = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zindex['AGB+SNcc']['low']]].loc[elem].to_numpy().astype(np.float64)
            interp = interpolate.interp1d(x, y, kind='linear', fill_value='extrapolate')
        else:
            x_low = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zindex['AGB+SNcc']['low']]].columns.to_numpy().astype(np.float64)
            x_high = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zindex['AGB+SNcc']['high']]].columns.to_numpy().astype(np.float64)
            y_low = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zindex['AGB+SNcc']['low']]].loc[elem].to_numpy().astype(np.float64)
            y_high = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zindex['AGB+SNcc']['high']]].loc[elem].to_numpy().astype(np.float64)