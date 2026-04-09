"""
Tools for plotting.
"""

from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd

from .io import read_numpy
from .comparison import compute_sum_rules


def save_or_show(fig, plot_target):
    """
    Output the given `fig`,
    and flush it from the buffer.
    `plot_target` may be one of:

        None: Output to the screen.

        "/dev/null": Do not output, only close the fig.

        Any other `str`: Output to the specified filename.

        A Matplotlib backend: Output via the given backend.
    """
    if plot_target:
        if isinstance(plot_target, str):
            if plot_target != "/dev/null":
                fig.savefig(plot_target)
        else:
            plot_target.savefig(fig)

        plt.close(fig)
    else:
        plt.show()


class Props:
    _markers = ["o", "s", "^", "v", "<", "H", ">", "+", "x", "D", "h", "p", "8", "*"]
    _colours = mpl.color_sequences["tab10"] + ["black", "seagreen", "violet"]

    def __init__(self, length_only=False, mass_only=False):
        self._props = {}
        self.length_only = length_only
        self.mass_only = mass_only
        assert not (length_only and mass_only)

    def burn(self):
        return self._colours.pop(0), self._markers.pop(0)

    def get(self, datum):
        lattice_size = datum["Nx"]
        assert lattice_size == datum["Ny"] and lattice_size == datum["Nz"]
        mass = datum["mass"]

        if self.length_only:
            key = lattice_size
            label = f"$L={lattice_size}$"
        elif self.mass_only:
            key = mass
            label = f"$am_f={mass}$"
        else:
            key = mass, lattice_size
            label = f"$am_f={mass}$, $L={lattice_size}$"

        if key not in self._props:
            self._props[key] = (
                self._colours.pop(0),
                self._markers.pop(0),
                label,
            )

        return self._props[key]


def get_args(old_data=False, select_plot=False):
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    if old_data:
        parser.add_argument("--previous_data", required=True)
    if select_plot:
        parser.add_argument(
            "--plot_type",
            required=True,
            choices=["frho-fpi", "frho-fa1", "wsr-i-normalised", "wsr-ii-normalised"],
        )
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def comparison_plot_main(callback, old_data=False, select_plot=False):
    args = get_args(old_data, select_plot)
    plt.style.use(args.plot_styles)

    new_data = [read_numpy(input_file) for input_file in args.input_files]

    params = {}
    if old_data:
        params["old_data"] = pd.read_csv(args.previous_data)
    if select_plot:
        params["plot_type"] = args.plot_type

    fig = callback(new_data, **params)
    save_or_show(fig, args.output_file)


def iterate_attribute(data, key, props=None, sort=True):
    ordering = {"length": reversed, "mass": lambda x: x}
    datum_keys = {"length": "Nx", "mass": "mass"}
    values = ordering[key](sorted(set(datum[datum_keys[key]] for datum in data)))
    if not props:
        props = Props(**{f"{key}_only": True})
    for value in values:
        subset = [datum for datum in data if datum[datum_keys[key]] == value]
        yield subset, props.get(subset[0])


def iterate_lengths(data):
    return iterate_attribute(data, "length")


def plot_new_series(ax, data, key, colour, marker, label):
    masses = [datum["mass"] for datum in data]
    values, errors = zip(*[datum[key] for datum in data])
    ax.errorbar(
        masses,
        values,
        yerr=errors,
        linestyle="none",
        marker=marker,
        color=colour,
        label=label,
    )


def add_qcd_value(ax, numerator, denominator=None):
    # S. Navas et al. (Particle Data Group),
    # Phys. Rev. D 110, 030001 (2024) and 2025 update.
    qcd_values = {
        # https://pdglive.lbl.gov/Particle.action?init=0&node=M009&home=MXXX005
        "rho_mass": 770.26,
        # https://pdglive.lbl.gov/DataBlock.action?node=M010M
        "a_1_mass": 1230,
        # https://arxiv.org/pdf/1507.02541
        "rho_decay_const": 221.1,
        # WHERE CAN THIS NUMBER COME FROM???
        "a_1_decay_const": 300,
        # ??????
        "pi_decay_const": 130.2,
    }

    if numerator in qcd_values:
        value = qcd_values[numerator]
        if denominator is not None:
            value /= qcd_values[denominator]
    elif numerator in compute_sum_rules.rules:
        value = compute_sum_rules.rules[numerator](qcd_values)
    else:
        raise NotImplementedError

    ax.plot(
        [0],
        [value],
        marker="*",
        markersize=10,
        color="black",
        fillstyle="full",
        label="Real-world QCD",
        linestyle="none",
    )
