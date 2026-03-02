"""
Tools for plotting.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt


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

    def __init__(self, length_only=False):
        self._props = {}
        self.length_only = length_only

    def get(self, datum):
        lattice_size = datum["Nx"]
        assert lattice_size == datum["Ny"] and lattice_size == datum["Nz"]
        mass = datum["mass"]

        if self.length_only:
            key = lattice_size
            label = f"$L={lattice_size}$"
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
