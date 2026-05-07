#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import pandas as pd

from ..define import define


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_metadata")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def bracket(items):
    return f"\\left\\{{{', '.join(map(str, items))}\\right\\}}"


def get_light_ensembles(metadata):
    subset = metadata.query("plot_light_ensembles").sort_values(by="mf")
    components = [
        f"({datum['Nx']}, {datum['mf']})" for datum in subset.to_dict(orient="records")
    ]
    return define("Light_Ensembles_With_Systematics", bracket(components))


def get_single_multi_volume_ensembles(metadata):
    counts = metadata.value_counts("mf").to_frame().reset_index()
    definitions = []
    for query, label in [("count == 1", "Single"), ("count > 1", "Multi")]:
        masses = sorted(counts.query(query).mf)
        definitions.append(define(f"{label}_Volume_Ensembles", bracket(masses)))
    return "\n".join(definitions)


def get_all_volumes(metadata):
    assert (metadata["Nx"] == metadata["Ny"]).all()
    assert (metadata["Nx"] == metadata["Nz"]).all()
    volumes = set(metadata[["Nt", "Nx"]].itertuples(index=False))
    volume_strings = [
        f"({datum.Nt}, {datum.Nx})" for datum in sorted(volumes, reverse=True)
    ]
    return define("All_Lattice_Volumes", bracket(volume_strings))


def main():
    args = get_args()
    metadata = pd.read_csv(args.input_metadata, comment="#")
    print(get_light_ensembles(metadata), file=args.output_definitions)
    print(get_single_multi_volume_ensembles(metadata), file=args.output_definitions)
    print(get_all_volumes(metadata), file=args.output_definitions)


if __name__ == "__main__":
    main()
