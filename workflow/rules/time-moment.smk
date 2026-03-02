spectrum = pd.read_csv(config["spectrum_file"])

rule extract_v_a:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        Nt=get_metadata("Nt"),
        Nx=get_metadata("Nx"),
        Ny=get_metadata("Ny"),
        Nz=get_metadata("Nz"),
    input:
        data=f"data/{subdir_format}/correlator.log.zst",
        script="s_parameter/vp_momentum/extract_correlators.py",
    output:
        data=f"processed_data/{subdir_format}/v_a.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output {output.data} "
        + " ".join(
            f"--source OneLinkCurrent{channel}{direction}"
            for channel in ["V", "A"]
            for direction in range(4)
        )
        + " --Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        + " ".join(
            f"--sink {current}Current{channel}{direction}"
            for current in ["OneLink", "Conserved"]
            for channel in ["V", "A"]
            for direction in range(4)
        )


rule renormalise_v_a:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        script="s_parameter/time_moment/renormalise_v_a.py",
        v_a=rules.extract_v_a.output.data,
        Z_A=rules.Z_A.output.data,
    output:
        data=f"processed_data/{subdir_format}/renormalised_v_a.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.v_a} {input.Z_A} --output_file {output.data}"


rule S_parameter_tm_small_t:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.renormalise_v_a.output.data,
        script="s_parameter/time_moment/sparam_tm.py",
    output:
        data=f"processed_data/{subdir_format}/s_parameter_small_t.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data}"


rule meson_v_a:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        plateau_start=lambda wildcards: get_metadata(f"plateau_start_{wildcards.channel}"),
        plateau_end=lambda wildcards: get_metadata(f"plateau_end_{wildcards.channel}"),
    input:
        data=rules.renormalise_v_a.output.data,
        script="s_parameter/time_moment/single_correlator_fit.py",
    output:
        data=f"processed_data/{subdir_format}/{{channel}}_mass_decay.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} "
        "--channel {wildcards.channel} "
        "--min_timeslice {params.plateau_start} --max_timeslice {params.plateau_end}"


rule meson_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel=channel)
            for channel in ["A", "V"]
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        script="s_parameter/time_moment/plot_masses.py",
        plot_styles=config["plot_styles"],
    output:
        plot="processed_data/spectrum.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule meson_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.meson_v_a.output.data.format(**row, channel=channel)
            for channel in ["A", "V"]
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        script="s_parameter/time_moment/tabulate_masses.py",
    output:
        plot="processed_data/spectrum.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.plot}"


rule S_parameter_tm_fit:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        plateau_start=get_metadata("plateau_start_va"),
        plateau_end=get_metadata("plateau_end_va"),
    input:
        data=rules.S_parameter_tm_small_t.output.data,
        m_rho=f"processed_data/{subdir_format}/V_mass_decay.json",
        m_a_1=f"processed_data/{subdir_format}/A_mass_decay.json",
        script="s_parameter/time_moment/sparam_fit.py",
    output:
        data=f"processed_data/{subdir_format}/s_parameter_tm_fit.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--input_m_rho {input.m_rho} --input_m_a_1 {input.m_a_1} "
        "--min_timeslice {params.plateau_start} --max_timeslice {params.plateau_end} "
        "--output_file {output.data}"


rule plot_zero_momentum_correlator:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.renormalise_v_a.output.data,
        plot_styles=config["plot_styles"],
        script="s_parameter/time_moment/plot_correlator.py",
    output:
        plot=f"processed_data/{subdir_format}/zero_momentum_correlator.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_fitted_zero_momentum_correlator:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.renormalise_v_a.output.data,
        fit_result=rules.S_parameter_tm_fit.output.data,
        plot_styles=config["plot_styles"],
        script="s_parameter/time_moment/plot_correlator_with_fit.py",
    output:
        plot=f"processed_data/{subdir_format}/fitted_zero_momentum_correlator.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --fit_result {input.fit_result} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_S_eff_timemoment:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.S_parameter_tm_small_t.output.data,
        fit_result=rules.S_parameter_tm_fit.output.data,
        plot_styles=config["plot_styles"],
        script="s_parameter/time_moment/plot_s_bar.py",
    output:
        plot=f"processed_data/{subdir_format}/S_extrapolation.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --fit_result {input.fit_result} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_S_contributions_timemoment:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.S_parameter_tm_fit.output.data.format(**row)
            for row in metadata.to_dict(orient="records")
            if row["plot_light_ensembles"]
        ],
        plot_styles=config["plot_styles"],
        script="s_parameter/time_moment/plot_s_ir.py",
    output:
        plot="assets/plots/S_timemoment_contributions.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"
