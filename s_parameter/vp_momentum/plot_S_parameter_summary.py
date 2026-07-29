#!/usr/bin/env python3

from .plot_S_parameter import get_S, main


if __name__ == "__main__":
    main(get_S, r"$S|_{\textnormal{1-doublet}}$", (0.0, 0.3), flatten=True)
