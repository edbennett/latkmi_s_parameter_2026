#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import pandas as pd


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def get_row(metadatum):
    traj_length_slug = {0.5: "*", 1: ""}[metadatum["trajectory_length"]]
    assert metadatum["Nx"] == metadatum["Ny"] == metadatum["Nz"]
    metadatum["configuration_count"] = (
        metadatum["trajectory_count"] // metadatum["configuration_separation"]
    )

    columns = [
        ("Nx", False),
        ("Nt", False),
        ("mf", False),
        ("trajectory_count", True),
        ("configuration_separation", True),
        ("configuration_count", True),
        ("bin_width", True),
        ("source_count", False),
    ]
    elements = [
        "{value}{slug}".format(
            value=metadatum[key], slug=traj_length_slug if use_slug else ""
        )
        for key, use_slug in columns
    ]

    return " & ".join(map(str, elements)) + " \\\\"


def tabulate(metadata):
    header = "\n".join(
        [
            r"\begin{tabular}{llllllll}",
            r"\toprule",
            " & ".join(
                [
                    "$L$",
                    "$T$",
                    "$m_f$",
                    r"$N_{\mathrm{traj}}$",
                    r"$N_{\mathrm{traj}}^{\mathrm{int}}$",
                    r"$N_{\mathrm{conf}}$",
                    r"$N_{\mathrm{traj}}^{\mathrm{bin}}$",
                    r"$N_{\mathrm{src}}$",
                ]
            )
            + r" \\",
            r"\midrule",
        ]
    )
    footer = "\n".join([r"\bottomrule", r"\end{tabular}"])

    content = [
        get_row(metadatum)
        for metadatum in metadata.sort_values(
            by=["Nt", "mf"], ascending=[False, True]
        ).to_dict(orient="records")
    ]
    return "\n".join([header, *content, footer])


def main():
    args = get_args()
    metadata = pd.read_csv(args.input_file, comment="#")
    print(tabulate(metadata), file=args.output_file)


if __name__ == "__main__":
    main()
