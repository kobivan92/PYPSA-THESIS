# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import matplotlib.pyplot as plt
import pypsa
import seaborn as sns
import plotly.graph_objects as go

from scripts._helpers import configure_logging, set_scenario_config

sns.set_theme("paper", style="whitegrid")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_elec_statistics",
            opts="Ept-12h",
            clusters="37",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    n.loads.carrier = "load"
    n.carriers.loc["load", ["nice_name", "color"]] = "Load", "darkred"
    colors = n.carriers.set_index("nice_name").color.where(
        lambda s: s != "", "lightgrey"
    )

    def rename_index(ds):
        specific = ds.index.map(lambda x: f"{x[1]}\n({x[0]})")
        generic = ds.index.get_level_values("carrier")
        duplicated = generic.duplicated(keep=False)
        index = specific.where(duplicated, generic)
        return ds.set_axis(index)

    def plot_static_per_carrier(ds, ax, drop_zero=True):
        if drop_zero:
            ds = ds[ds != 0]
        ds = ds.dropna()
        c = colors[ds.index.get_level_values("carrier")]
        ds = ds.pipe(rename_index)
        label = f"{ds.attrs['name']} [{ds.attrs['unit']}]"
        ds.plot.barh(color=c.values, xlabel=label, ax=ax)
        ax.grid(axis="y")

    def save_plotly_bar(ds, filename, colors, title=None, orientation='h'):
        if ds.empty:
            with open(filename, 'w') as f:
                f.write("<html><body><h2>No data to plot</h2></body></html>")
            return
        if hasattr(ds, 'attrs') and 'name' in ds.attrs and 'unit' in ds.attrs:
            y_label = f"{ds.attrs['name']} [{ds.attrs['unit']}]"
        else:
            y_label = title or ''
        if orientation == 'h':
            fig = go.Figure(go.Bar(
                y=ds.index.astype(str),
                x=ds.values,
                orientation='h',
                marker_color=[colors.get(i[1] if isinstance(i, tuple) else i, 'lightgrey') for i in ds.index],
            ))
            fig.update_layout(yaxis_title='', xaxis_title=y_label, title=title)
        else:
            fig = go.Figure(go.Bar(
                x=ds.index.astype(str),
                y=ds.values,
                marker_color=[colors.get(i[1] if isinstance(i, tuple) else i, 'lightgrey') for i in ds.index],
            ))
            fig.update_layout(xaxis_title='', yaxis_title=y_label, title=title)
        fig.write_html(filename)

    fig, ax = plt.subplots()
    ds = n.statistics.capacity_factor().dropna()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.capacity_factor_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.capacity_factor_bar.replace('.pdf', '.html'), colors, title='Capacity Factor')

    fig, ax = plt.subplots()
    ds = n.statistics.installed_capacity().dropna()
    ds = ds.drop("Line")
    ds = ds.drop(("Generator", "Load"), errors="ignore")
    ds = ds / 1e3
    ds.attrs["unit"] = "GW"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.installed_capacity_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.installed_capacity_bar.replace('.pdf', '.html'), colors, title='Installed Capacity')

    fig, ax = plt.subplots()
    ds = n.statistics.optimal_capacity()
    ds = ds.drop("Line")
    ds = ds.drop(("Generator", "Load"), errors="ignore")
    ds = ds / 1e3
    ds.attrs["unit"] = "GW"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.optimal_capacity_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.optimal_capacity_bar.replace('.pdf', '.html'), colors, title='Optimal Capacity')

    fig, ax = plt.subplots()
    ds = n.statistics.capex()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.capital_expenditure_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.capital_expenditure_bar.replace('.pdf', '.html'), colors, title='Capital Expenditure')

    fig, ax = plt.subplots()
    ds = n.statistics.opex()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.operational_expenditure_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.operational_expenditure_bar.replace('.pdf', '.html'), colors, title='Operational Expenditure')

    fig, ax = plt.subplots()
    ds = n.statistics.curtailment()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.curtailment_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.curtailment_bar.replace('.pdf', '.html'), colors, title='Curtailment')

    fig, ax = plt.subplots()
    ds = n.statistics.supply()
    ds = ds.drop("Line")
    ds = ds / 1e6
    ds.attrs["unit"] = "TWh"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.supply_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.supply_bar.replace('.pdf', '.html'), colors, title='Supply')

    fig, ax = plt.subplots()
    ds = n.statistics.withdrawal()
    ds = ds.drop("Line")
    ds = ds / -1e6
    ds.attrs["unit"] = "TWh"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.withdrawal_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.withdrawal_bar.replace('.pdf', '.html'), colors, title='Withdrawal')

    fig, ax = plt.subplots()
    ds = n.statistics.market_value()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.market_value_bar)
    plt.close(fig)
    save_plotly_bar(ds, snakemake.output.market_value_bar.replace('.pdf', '.html'), colors, title='Market Value')

    # touch file
    with open(snakemake.output.barplots_touch, "a"):
        pass
