"""
Tools for plotting.
"""

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
