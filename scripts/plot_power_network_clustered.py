# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Plot clustered electricity transmission network.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from matplotlib.lines import Line2D
from pypsa.plot import add_legend_lines
import plotly.graph_objects as go

from scripts._helpers import set_scenario_config
from scripts.plot_power_network import load_projection

def save_plotly_map_placeholder(filename):
    with open(filename, 'w') as f:
        f.write("<html><body><h2>Interactive map not yet supported in Plotly for this plot. Please use the static PDF/SVG.</h2></body></html>")

def plot_power_network_clustered_plotly(n, html_out, title="Clustered Power Network"):
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
            "plot_power_network_clustered",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )
    set_scenario_config(snakemake)

    lw_factor = 2e3

    n = pypsa.Network(snakemake.input.network)

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    proj = load_projection(snakemake.params.plotting)

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": proj})
    regions.to_crs(proj.proj4_init).plot(
        ax=ax, facecolor="none", edgecolor="lightgray", linewidth=0.75
    )
    n.plot(
        ax=ax,
        margin=0.06,
        line_widths=n.lines.s_nom / lw_factor,
        link_colors=n.links.p_nom.apply(
            lambda x: "darkseagreen" if x > 0 else "skyblue"
        ),
        link_widths=2.0,
    )

    sizes = [10, 20]
    labels = [f"HVAC ({s} GW)" for s in sizes]
    scale = 1e3 / lw_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc=[0.25, 0.9],
        frameon=False,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=13,
    )

    add_legend_lines(
        ax, sizes, labels, patch_kw=dict(color="rosybrown"), legend_kw=legend_kw
    )

    handles = [
        Line2D([0], [0], color="darkseagreen", lw=2),
        Line2D([0], [0], color="skyblue", lw=2),
    ]
    plt.legend(
        handles,
        ["HVDC existing", "HVDC planned"],
        frameon=False,
        loc=[0.0, 0.9],
        fontsize=13,
    )

    plt.savefig(snakemake.output.map, bbox_inches="tight")
    plt.close()

    # After saving the Matplotlib map
    html_fn = list(snakemake.output)[0]
    html_fn = html_fn.replace('.pdf', '.html').replace('.svg', '.html')
    save_plotly_map_placeholder(html_fn)

    plot_power_network_clustered_plotly(
        n,
        snakemake.output[0].replace('.pdf', '.html'),
        title="Clustered Power Network"
    )
