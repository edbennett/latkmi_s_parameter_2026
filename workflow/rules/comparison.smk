rule get_previous_fit_results:
    output:
        data="previous_data/fit_results_2505.08658.csv",
    conda:
        "../envs/zenodo_get.yml"
    shadow: "minimal"
    shell:
        """
        zenodo_get --doi https://doi.org/10.5281/zenodo.17037868 --glob fit_results.csv
        mv fit_results.csv {output.data}
        """


rule get_previous_spectrum:
    output:
        data="previous_data/spectrum_2505.08658.csv",
    conda:
        "../envs/zenodo_get.yml"
    shadow: "minimal"
    shell:
        """
        zenodo_get --doi https://doi.org/10.5281/zenodo.17037868 --glob spectrum.csv
        mv spectrum.csv {output.data}
        """


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
            rules.pade_fit_systematics.output.data.format(**row, kind="renormalised")
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
            rules.pade_fit_systematics.output.data.format(**row, kind="renormalised")
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
        skip_dmo_flag=lambda wildcards: (
            ""
            if get_metadata("plot_light_ensembles")(wildcards)[0]
            else "--skip_rule dmo"
        ),
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
        "--previous_data {input.previous_data} {params.skip_dmo_flag} "
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
        spectrum_data="previous_data/spectrum_2505.08658.csv",
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
        spectrum_data="previous_data/spectrum_2505.08658.csv",
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
            rules.pade_fit_systematics.output.data.format(
                **ensemble, kind="renormalised"
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
        sum_rule_data=[
            rules.compute_sum_rules.output.data.format(**datum)
            for datum in metadata.to_dict(orient="records")
        ],
        time_moment_data=[
            rules.infinite_volume_extrapolate.output.data.format(**datum)
            for datum in metadata[["Nf", "mf"]].drop_duplicates().to_dict(orient="records")
        ],
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


rule compute_lsd_infinite_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        lsd_s_data="external_data/lsd_prd14_spectra_sparameter_table_6.csv",
        lsd_chiral_data="external_data/lsd_prd14_spectra_chiral_limit_table_3_9.csv",
        fit_result=rules.finite_volume_fit.output.data,
        script="s_parameter/comparison/lsd_infinite_volume.py",
    # This is a slow process due to computing all finite volume factors at once
    priority: 50
    threads: 5
    output:
        data="processed_data/lsd_infinite_volume.csv",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} "
        "--s_parameter_data {input.lsd_s_data} --chiral_data {input.lsd_chiral_data} "
        "--fit_result {input.fit_result} "
        "--output_file {output.data}"


rule plot_lsd_infinite_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        time_moment_data=rules.S_parameter_infinite_volume_comparison.input.time_moment_infinite_volume,
        chiral_data=rules.get_previous_fit_results.output.data,
        spectrum_data=config["spectrum_file"],
        lsd_data=rules.compute_lsd_infinite_volume.output.data,
        plot_styles=config["plot_styles"],
        script="s_parameter/comparison/plot_infinite_volume_S_lsd.py",
    output:
        plot="assets/plots/S_parameter_infinite_volume_comparison_lsd.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} "
        "--time_moment_infinite_volume {input.time_moment_data} "
        "--spectrum_data {input.spectrum_data} --chiral_fit_result {input.chiral_data} "
        "--lsd_data {input.lsd_data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule define_lattice_spacing_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        latkmi_data=rules.get_previous_fit_results.output.data,
        lsd_data="external_data/lsd_prd14_spectra_chiral_limit_table_3_9.csv",
        script="s_parameter/definitions/lattice_spacing_matching.py",
    output:
        definitions="assets/definitions/lattice_spacing_matching.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} --chiral_fit_result_latkmi {input.latkmi_data} "
        "--chiral_spectrum_lsd {input.lsd_data} "
        "--output_definitions {output.definitions}"


rule define_l10_r:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.compute_sum_rules.output.data,
        script="s_parameter/definitions/l10_r.py",
    output:
        definitions=f"processed_data/{subdir_format}/l10_r.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output_definitions {output.definitions}"


rule define_spectrum_observations:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=config["spectrum_file"],
        script="s_parameter/definitions/spectrum_observations.py",
    output:
        definitions=f"assets/definitions/spectrum_observations.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output_definitions {output.definitions}"


rule define_tm_S_infinite_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.infinite_volume_extrapolate.output.data,
        script="s_parameter/definitions/tm_infinite_volume.py",
    output:
        definitions="processed_data/nf{Nf}/mf{mf}/infinite_volume_S_tm.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output_definitions {output.definitions}"
