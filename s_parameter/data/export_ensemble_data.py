#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np
import pandas as pd

from ..io import read_numpy


SOURCE = "10.5281/zenodo.21297277"
MATCH_KEYS = ["Nf", "Nx", "Ny", "Nz", "Nt", "mass"]


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "base_data",
        help="CSV filename for the data inherited from previous work, into which other data will be appended",
    )
    parser.add_argument("new_data", nargs="+", help="New data to include")
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to place the output",
    )
    return parser.parse_args()


def get_index(data, datum):
    expanded_data = data.copy()
    expanded_data["mass"] = expanded_data["mf"]
    query = " & ".join(f"{key} == {datum[key]}" for key in MATCH_KEYS)
    matches = expanded_data.query(query)
    if len(matches) == 0:
        return None
    elif len(matches) == 1:
        return matches.iloc[0].name
    else:
        raise ValueError(f"Multiple results found matching {query}")


def base_metadata(datum):
    return {key: datum[key] for key in ["Nx", "Ny", "Nz", "Nt", "Nf"]}.update(
        beta=3.8, mf=datum["mass"]
    )


def add_value(data, index, name, values):
    for prefix, value in zip(["value", "error", "syst"], values):
        if not np.isnan(value):
            data.loc[index, f"{prefix}_{name}"] = value


def update_source(data, index):
    source = data.loc[index, "source"]
    if SOURCE not in source:
        data.loc[index, "source"] = f"{source}; {SOURCE}"


def add_S_TM(data, index, datum):
    add_value(data, index, "S_time_moment", datum["S_infinite_t"])
    update_source(data, index)


def add_S_VP(data, index, datum):
    result = datum["pade_fit_result"]["Conserved"]["S"]
    add_value(data, index, "S_vacuum_polarization", result)
    update_source(data, index)


def add_sum_rules(data, index, datum):
    keys = [
        "wsr-i",
        "wsr-ii",
        "wsr-i-normalised",
        "wsr-ii-normalised",
        "ksrf-i",
        "ksrf-ii",
        "dmo",
        "l10-r",
        "dmo-l10-r",
    ]
    for key in keys:
        if key in datum:
            add_value(data, index, key.replace("-", "_"), datum[key])
    update_source(data, index)


def add_Z_A(data, index, datum):
    add_value(data, index, "Z_A", datum["Z_A"])
    update_source(data, index)


def add_meson(data, index, datum):
    channels = {"V": ("rho", "rho"), "A": ("a_1", "a1")}
    something_found = False
    for observable_name, observable_shortname in ("mass", "m"), ("decay_const", "f"):
        channel_name, channel_shortname = channels[datum["channel"]]
        name = f"{channel_name}_{observable_name}"
        shortname = f"{observable_shortname}{channel_shortname}"
        if name in datum["fit_result"]:
            something_found = True
            add_value(data, index, shortname, datum["fit_result"][name])
    if something_found:
        update_source(data, index)


def add_S_TM_infinite_volume(data, datum):
    for systematic in datum["systematics"]:
        index = get_index(data, {**datum, **systematic})
        value = np.concatenate([datum["S_infinite_volume"], [systematic["systematic"]]])
        add_value(data, index, "S_infinite_volume", value)
        update_source(data, index)


def add_data(base_data, new_data):
    data = base_data.copy()
    data["Nt"] = data["T"]
    for key in "Nx", "Ny", "Nz":
        data[key] = data["L"]
    data.drop(columns=["L", "T"], inplace=True)

    for datum in new_data:
        if "S_infinite_volume" in datum:
            add_S_TM_infinite_volume(data, datum)
            continue

        index = get_index(data, datum)
        if index is None:
            data = pd.concat(data, [base_metadata(datum)])
            index = len(data) - 1

        if "S_infinite_t" in datum:
            add_S_TM(data, index, datum)
        elif "pade_fit_result" in datum:
            add_S_VP(data, index, datum)
        elif "ksrf-ii" in datum:
            add_sum_rules(data, index, datum)
        elif "Z_A" in datum:
            add_Z_A(data, index, datum)
        elif "fit_result" in datum:
            add_meson(data, index, datum)
        else:
            breakpoint()

    return data


def sort_key(element):
    prefixes = ["value", "error", "syst_upper", "syst_lower", "syst"]
    for prefix_idx, prefix in enumerate(prefixes):
        if element.startswith(prefix):
            # Trim out the trailing underscore
            observable = element[len(prefix) + 1 :]
            return (observable, prefix_idx)
    else:
        return ("", 0)


def main():
    args = get_args()
    base_data = pd.read_csv(args.base_data, comment="#")
    new_data = [read_numpy(filename) for filename in args.new_data]
    combined_data = add_data(base_data, new_data).sort_index(
        axis="columns",
        key=lambda x: x.map(sort_key),
    )
    combined_data.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
