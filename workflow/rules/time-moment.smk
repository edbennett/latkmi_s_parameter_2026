spectrum = pd.read_csv(config["spectrum_file"])

rule extract_v_a:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        Nt=get_metadata("Nt"),
        Nx=get_metadata("Nx"),
        Ny=get_metadata("Ny"),
        Nz=get_metadata("Nz"),
        Nf=get_metadata("Nf"),
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
        + " --Nf {params.Nf} "
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
        plot="assets/tables/spectrum.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.plot}"


rule reconstruct_mass_samples:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        script="s_parameter/time_moment/reconstruct_mass_samples.py",
        datum=lambda wildcards: rules.meson_v_a.output.data.format(channel="V", **wildcards),
        spectrum_data=config["spectrum_file"],
    output:
        data=f"processed_data/{subdir_format}/reconstructed_mass_samples.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.datum} --spectrum_data {input.spectrum_data} "
        "--output_file {output.data}"


def only_if_metadata(column, filename, from_df=metadata):
    def only_if_metadata_inner(wildcards):
        (required,) = lookup(
            within=from_df,
            query=metadata_query.format(**wildcards),
            cols=column,
        )
        if required:
            return filename.format(**wildcards)
        return []

    return only_if_metadata_inner


rule S_parameter_tm_fit:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        plateau_start=get_metadata("plateau_start_va"),
        plateau_end=get_metadata("plateau_end_va"),
        rho_vt_flag=lambda wildcards, input: (
            f"--m_rho_vt_samples {input.rho_vt_samples}"
            if input.rho_vt_samples
            else ""
        ),
        a_1_flag=lambda wildcards, input: (
            f"--m_a_1_samples {input.a_1_samples}"
            if input.a_1_samples
            else ""
        ),
    input:
        data=rules.S_parameter_tm_small_t.output.data,
        rho_vt_samples=only_if_metadata(
            "use_VT_timemoment",
            f"processed_data/{subdir_format}/V_mass_decay.json"
        ),
        a_1_samples=only_if_metadata(
            "use_VT_timemoment",
            f"processed_data/{subdir_format}/A_mass_decay.json"
        ),
        mass_samples=rules.reconstruct_mass_samples.output.data,
        script="s_parameter/time_moment/sparam_fit.py",
    output:
        data=f"processed_data/{subdir_format}/s_parameter_tm_fit.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--input_mass_samples {input.mass_samples} "
        "{params.rho_vt_flag} {params.a_1_flag} "
        "--min_timeslice {params.plateau_start} --max_timeslice {params.plateau_end} "
        "--output_file {output.data}"


rule finite_volume_factor:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.S_parameter_tm_fit.output.data,
        spectrum=config["spectrum_file"],
        script="s_parameter/time_moment/finite_volume.py",
    output:
        data=f"processed_data/{subdir_format}/finite_volume_factor.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --previous_data {input.spectrum} "
        "--output_file {output.data}"


def multi_volume_ensembles():
    metadata_filter = metadata.value_counts("mf") > 1
    masses = list(metadata_filter[metadata_filter].index)
    return metadata.query(f"mf in {masses}")


rule finite_volume_fit:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.finite_volume_factor.output.data.format(**ensemble)
            for ensemble in multi_volume_ensembles().to_dict(orient="records")
        ],
        script="s_parameter/time_moment/finite_volume_fit.py",
    output:
        data="processed_data/finite_volume_fit.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.data}"


rule finite_volume_fit_form:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        script="s_parameter/time_moment/finite_volume_fit_form.py",
    output:
        data="processed_data/finite_volume_fit_form.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} --output_file {output.data}"


rule infinite_volume_extrapolate:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=lambda wildcards: [
            rules.finite_volume_factor.output.data.format(**ensemble)
            for ensemble in metadata.query(
                f"Nf == {wildcards.Nf} & mf == {wildcards.mf}"
            ).to_dict(orient="records")
        ],
        fit=rules.finite_volume_fit.output.data,
        script="s_parameter/time_moment/infinite_volume_single_ensemble.py",
    output:
        data="processed_data/nf{Nf}/mf{mf}/infinite_volume_S.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --fit_result {input.fit} "
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


rule plot_S_infinite_volume_fit:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.finite_volume_fit.input.data,
        fit_result=rules.finite_volume_fit.output.data,
        fit_form=rules.finite_volume_fit_form.output.data,
        plot_styles=config["plot_styles"],
        script="s_parameter/time_moment/plot_finite_volume_fit.py",
    output:
        plot="assets/plots/finite_volume_fit.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--fit_result {input.fit_result} --fit_form {input.fit_form} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"
