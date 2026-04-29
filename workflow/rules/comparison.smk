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


rule compute_sum_rules:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            f"processed_data/{subdir_format}/{channel}_mass_decay.json"
            for channel in ["V", "A"]
        ],
        infinite_volume_data=rules.infinite_volume_extrapolate.output.data,
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/compute_sum_rules.py",
    output:
        data=f"processed_data/{subdir_format}/sum_rules.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} {input.infinite_volume_data} "
        "--previous_data {input.previous_data} "
        "--output_file {output.data}"


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


rule comparison_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule simple_comparison_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.compute_sum_rules.output.data.format(**row)
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        script="s_parameter/comparison/plot_simple_ratio.py",
        plot_styles=config["plot_styles"],
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --plot_type {params.plot_type} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule comparison_plot_with_previous:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.previous_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


use rule comparison_plot_with_previous as plot_a_1_rho_ratio with:
    input:
        data=rules.simple_comparison_plot.input.data,
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/plot_a_1_rho_ratio.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/a_1_over_rho_mass_comparison.pdf",


use rule simple_comparison_plot as plot_rho_a_1_decay_ratio with:
    params:
        plot_type="frho-fa1",
    output:
        plot="assets/plots/rho_over_a_1_decay_comparison.pdf",


use rule simple_comparison_plot as plot_rho_pi_decay_ratio with:
    params:
        plot_type="frho-fpi",
    output:
        plot="assets/plots/rho_over_pi_decay_comparison.pdf",


use rule comparison_plot as plot_ksrf_i_ii with:
    input:
        data=rules.simple_comparison_plot.input.data,
        previous_data=config["spectrum_file"],
        script="s_parameter/comparison/plot_ksrf_i_ii.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/ksrf_i_ii.pdf",


rule plot_ksrf_ii_lsd:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.simple_comparison_plot.input.data,
        lsd_data="external_data/lsd_prd19_spectra_nf08_table_1_3_4.csv",
        script="s_parameter/comparison/plot_ksrf_ii_lsd.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/ksrf_ii_lsd.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --lsd_data {input.lsd_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


use rule simple_comparison_plot as plot_wsr with:
    params:
        plot_type="wsr-{wsr_idx}-normalised",
    output:
        plot="assets/plots/wsr_{wsr_idx}.pdf",


rule tabulate_sum_rules:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.simple_comparison_plot.input.data,
        script="s_parameter/comparison/tabulate_sum_rules.py",
    output:
        table="assets/tables/sum_rules.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule S_parameter_group_by_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        S_data=rules.plot_S_parameter.input.data,
        spectrum_data="previous_data/spectrum.csv",
        script="s_parameter/comparison/plot_S_all_ensembles.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/S_parameter_mpi_L_group_by_mass.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.S_data} "
        "--spectrum_data {input.spectrum_data} "
        "--group_by mass "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule S_parameter_group_by_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        S_data=rules.plot_S_parameter.input.data,
        spectrum_data="previous_data/spectrum.csv",
        script="s_parameter/comparison/plot_S_all_ensembles.py",
        plot_styles=config["plot_styles"],
        lsd_data="external_data/lsd_prd14_spectra_sparameter_table_6.csv",
    output:
        plot="assets/plots/S_parameter_mpi_L_group_by_length.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.S_data} "
        "--spectrum_data {input.spectrum_data} "
        "--lsd_data {input.lsd_data} "
        "--group_by length "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule S_parameter_infinite_volume_comparison:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        time_moment_infinite_volume=[
            rules.infinite_volume_extrapolate.output.data.format(**ensemble)
            for ensemble in metadata.query(
                "plot_large_volume_light_ensembles"
            ).to_dict(orient="records")
        ],
        time_moment_finite_volume=[
            rules.S_parameter_tm_fit.output.data.format(**ensemble)
            for ensemble in metadata.query(
                "plot_large_volume_light_ensembles"
            ).to_dict(orient="records")
        ],
        vp_momentum_finite_volume=[
            rules.pade_fit_renormalised.output.data.format(
                **ensemble, upper_bound="final"
            )
            for ensemble in metadata.query(
                "plot_large_volume_light_ensembles"
            ).to_dict(orient="records")
        ],
        plot_styles=config["plot_styles"],
        script="s_parameter/comparison/plot_infinite_volume_S.py",
    output:
        plot="assets/plots/S_parameter_infinite_volume_comparison.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} "
        "--time_moment_infinite_volume {input.time_moment_infinite_volume} "
        "--time_moment_finite_volume {input.time_moment_finite_volume} "
        "--vp_momentum_finite_volume {input.vp_momentum_finite_volume} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule tabulate_l10_r:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        sum_rule_data=rules.simple_comparison_plot.input.data,
        time_moment_data=rules.S_parameter_infinite_volume_comparison.input.time_moment_infinite_volume,
        script="s_parameter/comparison/tabulate_dmo_l10r.py",
    output:
        table="assets/tables/dmo_l10_r.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.sum_rule_data} {input.time_moment_data} "
        "--output_file {output.table}"


rule plot_l10_r:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.simple_comparison_plot.input.data,
        jlqcd_data="external_data/jlqcd_prl08_l10_r.csv",
        rbc_ukqcd_data="external_data/rbc_ukqcd_prd10_l10_r.csv",
        plot_styles=config["plot_styles"],
        script="s_parameter/comparison/plot_l10_r.py",
    output:
        plot="assets/plots/l10_r.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--jlqcd_data {input.jlqcd_data} --rbc_ukqcd_data {input.rbc_ukqcd_data} "
        "--output_file {output.plot} --plot_styles {input.plot_styles}"
