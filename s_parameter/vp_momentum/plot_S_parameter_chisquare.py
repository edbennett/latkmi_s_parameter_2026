#!/usr/bin/env python3


from .plot_S_parameter import main


def get_chisquare(datum):
    fit_result = datum["pade_fit_result"]["Conserved"]
    return fit_result["chisquare"] / fit_result["dof"]


if __name__ == "__main__":
    main(get_chisquare, r"$\chi^2/\textnormal{dof}$", (1e-4, 5.0), "log")
