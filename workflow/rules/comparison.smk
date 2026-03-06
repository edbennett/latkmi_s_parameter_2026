rule plot_tm_vp:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        secondary_flags=lambda wildcards, input: " ".join(
            f"--secondary_input_file {filename}" for filename in input.vp_data
        ),
    input:
        tm_data=[
            rules.S_parameter_tm_fit.output.data.format(**row)
            for row in metadata.to_dict(orient="records")
        ],
        vp_data=[
            rules.pade_fit_renormalised.output.data.format(**row, upper_bound="final")
            for row in metadata.to_dict(orient="records")
        ],
        plot_styles=config["plot_styles"],
        script="s_parameter/comparison/plot_tm_vp.py",
    output:
        plot="assets/plots/S_tm_vp_comparison.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.tm_data} {params.secondary_flags} "
        "--output_file {output.plot} --plot_styles {input.plot_styles}"


rule tabulate_vp_tm:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        tm_data=[
            rules.S_parameter_tm_fit.output.data.format(**row)
            for row in metadata.to_dict(orient="records")
        ],
        vp_data=[
            rules.pade_fit_renormalised.output.data.format(**row, upper_bound="final")
            for row in metadata.to_dict(orient="records")
        ],
        script="s_parameter/comparison/tabulate_S_parameter.py",
    output:
        plot="assets/tables/S_parameter.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.tm_data} {input.vp_data} "
        "--output_file {output.plot}"


rule plot_rho:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel="V")
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/plot_rho_mass.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/rho_mass_comparison.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.previous_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_a_1_rho_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel=channel)
            for channel in ["V", "A"]
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/plot_a_1_rho_ratio.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/a_1_over_rho_comparison.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.previous_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_rho_a_1_decay_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel=channel)
            for channel in ["V", "A"]
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        script="s_parameter/comparison/plot_rho_a_1_decay_ratio.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/rho_over_a_1_decay_comparison.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data}  "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_rho_pi_decay_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel="V")
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/plot_rho_pi_decay_ratio.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/rho_over_pi_decay_comparison.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.previous_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_ksrf_i_ii:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel="V")
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/plot_ksrf_i_ii.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/ksrf_i_ii.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.previous_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_ksrf_ii_lsd:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel="V")
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        previous_data=config["spectrum_file"],
        lsd_data="external_data/lsd_prd19_spectra_nf08_table_1_3_4.csv",
        script="s_parameter/comparison/plot_ksrf_ii_lsd.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/ksrf_ii_lsd.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.previous_data} "
        "--lsd_data {input.lsd_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"
