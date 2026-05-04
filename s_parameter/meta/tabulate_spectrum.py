#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np
import pandas as pd

from ..tables import formatter


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--input_metadata", required=True)
    parser.add_argument("--input_spectrum", required=True)
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def get_value(spectrum, key, slug):
    value = spectrum[f"value_{key}"]
    if np.isnan(value):
        return "---"
    error = spectrum[f"error_{key}"]
    return (
        formatter((value, error), length_control="smallest", significant_figures=2)
        + slug
    )


def get_row(metadatum, all_spectra):
    matched_spectra = all_spectra.query(
        f"Nf == {metadatum['Nf']} & mf == {metadatum['mf']} & T == {metadatum['Nt']} & L == {metadatum['Nx']}"
    )
    assert len(matched_spectra) == 1
    spectrum = matched_spectra.iloc[0]
    source_slug = {"arXiv:1610.07011": "", "arXiv:2505.08658": r"${}^\dagger$"}[
        spectrum.source
    ]

    row_data = [
        spectrum["L"],
        spectrum["T"],
        spectrum["mf"],
        get_value(spectrum, "fpi", source_slug),
        get_value(spectrum, "mpi", source_slug),
        get_value(spectrum, "mrho", source_slug),
        get_value(spectrum, "ma1", source_slug),
    ]
    return " & ".join(map(str, row_data)) + " \\\\"


def tabulate(metadata, spectrum):
    header = "\n".join(
        [
            r"\begin{tabular}{llllllll}",
            r"\toprule",
            " & ".join(
                [
                    "$L$",
                    "$T$",
                    "$m_f$",
                    r"$aF_\pi$",
                    r"$aM_\pi$",
                    r"$aM_{\rho(\mathrm{PV})}$",
                    r"$aM_{a_1(\mathrm{PV})}$",
                ]
            )
            + r" \\",
            r"\midrule",
        ]
    )
    footer = "\n".join([r"\bottomrule", r"\end{tabular}"])
    content = [
        get_row(metadatum, spectrum)
        for metadatum in metadata.sort_values(
            by=["Nt", "mf"], ascending=[False, True]
        ).to_dict(orient="records")
    ]
    return "\n".join([header, *content, footer])


def main():
    args = get_args()
    metadata = pd.read_csv(args.input_metadata, comment="#")
    spectrum = pd.read_csv(args.input_spectrum, comment="#")
    print(tabulate(metadata, spectrum), file=args.output_file)


if __name__ == "__main__":
    main()
