# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>>
#
# SPDX-License-Identifier: MIT
"""
Creates plots for optimised power network topologies and regional generation,
storage and conversion capacities built for the perfect foresight scenario.
"""

import logging

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from pypsa.plot import add_legend_circles, add_legend_lines
import plotly.graph_objects as go

from scripts._helpers import configure_logging, retry, set_scenario_config
from scripts.make_summary import assign_locations
from scripts.plot_power_network import load_projection, rename_techs_tyndp
from scripts.plot_summary import preferred_order

logger = logging.getLogger(__name__)


@retry
def plot_map_perfect(
    n,
    components=["Link", "Store", "StorageUnit", "Generator"],
    bus_size_factor=2e10,
):
    assign_locations(n)
    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)
    # investment periods
    investments = n.snapshots.levels[0]

    costs = {}
    for comp in components:
        df_c = n.df(comp)
        if df_c.empty:
            continue
        df_c["nice_group"] = df_c.carrier.map(rename_techs_tyndp)

        attr = "e_nom_opt" if comp == "Store" else "p_nom_opt"

        active = pd.concat(
            [n.get_active_assets(comp, inv_p).rename(inv_p) for inv_p in investments],
            axis=1,
        ).astype(int)
        capital_cost = n.df(comp)[attr] * n.df(comp).capital_cost
        capital_cost_t = (
            (active.mul(capital_cost, axis=0))
            .groupby([n.df(comp).location, n.df(comp).nice_group])
            .sum()
        )

        capital_cost_t.drop("load", level=1, inplace=True, errors="ignore")

        costs[comp] = capital_cost_t

    costs = pd.concat(costs).groupby(level=[1, 2]).sum()
    costs.drop(costs[costs.sum(axis=1) == 0].index, inplace=True)

    new_columns = preferred_order.intersection(costs.index.levels[1]).append(
        costs.index.levels[1].difference(preferred_order)
    )
    costs = costs.reindex(new_columns, level=1)

    for item in new_columns:
        if item not in snakemake.config["plotting"]["tech_colors"]:
            print(
                "Warning!",
                item,
                "not in config/plotting/tech_colors, assign random color",
            )
            snakemake.config["plotting"]["tech_colors"] = "pink"

    n.links.drop(
        n.links.index[(n.links.carrier != "DC") & (n.links.carrier != "B2B")],
        inplace=True,
    )

    # drop non-bus
    to_drop = costs.index.levels[0].symmetric_difference(n.buses.index)
    if len(to_drop) != 0:
        print("dropping non-buses", to_drop)
        costs.drop(to_drop, level=0, inplace=True, axis=0, errors="ignore")

    # make sure they are removed from index
    costs.index = pd.MultiIndex.from_tuples(costs.index.values)

    # PDF has minimum width, so set these to zero
    line_lower_threshold = 500.0
    line_upper_threshold = 1e4
    linewidth_factor = 2e3
    ac_color = "gray"
    dc_color = "m"

    line_widths = n.lines.s_nom_opt
    link_widths = n.links.p_nom_opt
    linewidth_factor = 2e3
    line_lower_threshold = 0.0
    title = "Today's transmission"

    line_widths[line_widths < line_lower_threshold] = 0.0
    link_widths[link_widths < line_lower_threshold] = 0.0

    line_widths[line_widths > line_upper_threshold] = line_upper_threshold
    link_widths[link_widths > line_upper_threshold] = line_upper_threshold

    for year in costs.columns:
        fig, ax = plt.subplots(subplot_kw={"projection": proj})
        fig.set_size_inches(7, 6)
        fig.suptitle(year)

        n.plot(
            bus_sizes=costs[year] / bus_size_factor,
            bus_colors=snakemake.config["plotting"]["tech_colors"],
            line_colors=ac_color,
            link_colors=dc_color,
            line_widths=line_widths / linewidth_factor,
            link_widths=link_widths / linewidth_factor,
            ax=ax,
            **map_opts,
        )

        sizes = [20, 10, 5]
        labels = [f"{s} bEUR/a" for s in sizes]
        sizes = [s / bus_size_factor * 1e9 for s in sizes]

        legend_kw = dict(
            loc="upper left",
            bbox_to_anchor=(0.01, 1.06),
            labelspacing=0.8,
            frameon=False,
            handletextpad=0,
            title="system cost",
        )

        add_legend_circles(
            ax,
            sizes,
            labels,
            srid=n.srid,
            patch_kw=dict(facecolor="lightgrey"),
            legend_kw=legend_kw,
        )

        sizes = [10, 5]
        labels = [f"{s} GW" for s in sizes]
        scale = 1e3 / linewidth_factor
        sizes = [s * scale for s in sizes]

        legend_kw = dict(
            loc="upper left",
            bbox_to_anchor=(0.27, 1.06),
            frameon=False,
            labelspacing=0.8,
            handletextpad=1,
            title=title,
        )

        add_legend_lines(
            ax, sizes, labels, patch_kw=dict(color="lightgrey"), legend_kw=legend_kw
        )

        legend_kw = dict(
            bbox_to_anchor=(1.52, 1.04),
            frameon=False,
        )

        fig.savefig(snakemake.output[f"map_{year}"], bbox_inches="tight")
        plt.close(fig)

        html_fn = list(snakemake.output)[0]
        html_fn = html_fn.replace('.pdf', '.html').replace('.svg', '.html')
        save_plotly_map_placeholder(html_fn)


def save_plotly_map_placeholder(filename):
    with open(filename, 'w') as f:
        f.write("<html><body><h2>Interactive map not yet supported in Plotly for this plot. Please use the static PDF/SVG.</h2></body></html>")


def plot_power_network_perfect_plotly(n, html_out, title="Perfect Foresight Power Network"):
    bus_trace = go.Scattergeo(
        lon=n.buses['x'],
        lat=n.buses['y'],
        text=n.buses.index,
        mode='markers',
        marker=dict(size=8, color='blue'),
        name='Buses'
    )
    line_traces = []
    for _, line in n.lines.iterrows():
        bus0 = n.buses.loc[line['bus0']]
        bus1 = n.buses.loc[line['bus1']]
        line_traces.append(go.Scattergeo(
            lon=[bus0['x'], bus1['x']],
            lat=[bus0['y'], bus1['y']],
            mode='lines',
            line=dict(width=2, color='gray'),
            opacity=0.5,
            showlegend=False
        ))
    fig = go.Figure([bus_trace] + line_traces)
    fig.update_layout(
        title=title,
        geo=dict(
            scope='europe',
            projection_type='mercator',
            showland=True,
            landcolor='rgb(243, 243, 243)',
            countrycolor='rgb(204, 204, 204)',
        ),
        margin=dict(l=0, r=0, t=40, b=0)
    )
    fig.write_html(html_out)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_power_network_perfect",
            opts="",
            clusters="37",
            sector_opts="4380H-T-H-B-I-A-dist1",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    regions = gpd.read_file(snakemake.input.regions).set_index("name")

    map_opts = snakemake.params.plotting["map"]

    if map_opts["boundaries"] is None:
        map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

    proj = load_projection(snakemake.params.plotting)

    plot_map_perfect(n)
    plot_power_network_perfect_plotly(
        n,
        snakemake.output[0].replace('.pdf', '.html'),
        title="Perfect Foresight Power Network"
    )
