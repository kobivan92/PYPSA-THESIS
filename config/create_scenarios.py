import itertools
import yaml

# Determine output path (works with or without Snakemake)
if "snakemake" in globals():
    output_file = snakemake.output[0]
else:
    output_file = "scenarios.yaml"

# 1. CO2 constraint: off/on
co2_options = [
    (False, "no-co2"),
    (True,  "co2")
]
# 2. Grid expansion: off/on (max_extension for lines and links)
grid_options = [
    (False, "nogrid", 0,     0),
    (True,  "grid",  20000, 30000)
]
# 3. Sector coupling: off/on
sector_options = [
    (False, "nosector", [""]),
    (True,  "sector",   ["TH"])
]
# 4. Storage: none / battery only / H2 only / both (Store list)
storage_options = [
    ([],      "nostorage", []),
    (["battery"], "bat",     ["battery"]),
    (["H2"],      "h2",      ["H2"]),
    (["battery", "H2"], "allstorage", ["battery", "H2"])
]
# 5. Generation mix: solar only, onwind only, or both
gen_options = [
    ("solar",  "solar",  ["solar", "solar-hsat"],            ["solar", "solar-hsat", "OCGT", "CCGT"]),
    ("onwind", "onwind", ["onwind"],                          ["onwind", "OCGT", "CCGT"]),
    ("both",   "both",   ["solar", "solar-hsat", "onwind"], ["solar", "solar-hsat", "onwind", "OCGT", "CCGT"]),
]

scenarios = {}

# Loop over all combinations of the five dimensions
tuple_product = itertools.product(
    co2_options,
    grid_options,
    sector_options,
    storage_options,
    gen_options
)
for co2_flag, grid_flag, sector_flag, storage_flag, gen_flag in tuple_product:
    co2_bool,  co2_name   = co2_flag
    grid_bool, grid_name, line_ext, link_ext = grid_flag
    sect_bool, sect_name, sector_opts = sector_flag
    store_items, stor_name, store_list  = storage_flag
    gen_key,   gen_name, ren_list, gen_list = gen_flag

    # Build a unique scenario name with all parameter parts
    parts = [co2_name, grid_name, sect_name, stor_name, gen_name]
    scenario_name = "-".join(parts)

    # Assemble overrides
    overrides = {
        "run": {
            "name": scenario_name
        },
        "electricity": {
            "co2limit_enable": co2_bool,
            "renewable_carriers": ren_list,
            "extendable_carriers": {
                "Generator": gen_list,
                "StorageUnit": [],
                "Store": store_list,
                "Link": []
            }
        },
        "lines": {
            "max_extension": line_ext
        },
        "links": {
            "max_extension": link_ext
        },
        "scenario": {
            "planning_horizons": [2050],
            "clusters": [25],
            "sector_opts": sector_opts
        },
        "snapshots": {
            "start": "2013-01-01",
            "end":   "2014-01-01",
            "inclusive": "left"
        },
        "foresight": "overnight"
    }

    scenarios[scenario_name] = overrides

# Write scenarios.yaml
with open(output_file, "w") as f:
    yaml.safe_dump(scenarios, f, sort_keys=False, default_flow_style=False)
