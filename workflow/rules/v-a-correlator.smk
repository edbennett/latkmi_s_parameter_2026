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


rule meson_v_a:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        plateau_start=lambda wildcards: get_metadata(f"plateau_start_{wildcards.channel}"),
        plateau_end=lambda wildcards: get_metadata(f"plateau_end_{wildcards.channel}"),
    input:
        data=rules.renormalise_v_a.output.data,
        script="s_parameter/spectrum/single_correlator_fit.py",
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
        script="s_parameter/spectrum/plot_masses.py",
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
        script="s_parameter/spectrum/tabulate_masses.py",
    output:
        plot="assets/tables/spectrum.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.plot}"
