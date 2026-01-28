rule extract_av:
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
        data=f"processed_data/{subdir_format}/av.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output {output.data} --source point "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        + " ".join(
            f"--sink {channel}"
            for channel in ["ConservedA4", "A4_1LINK", "PION_PS"]
        )


rule extract_vv_aa:
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
        data=f"processed_data/{subdir_format}/vv_aa.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        + " ".join(
            f"--source OneLinkCurrent{channel}{index}"
            for channel in ["V", "A"]
            for index in [0, 1, 2]
        ) + " " + " ".join(
            f"--sink {current_type}Current{channel}{index}"
            for channel in ["V", "A"]
            for index in [0, 1, 2]
            for current_type in ["OneLink", "Conserved"]
        )


rule extract_vpf:
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
        data=f"processed_data/{subdir_format}/vpf.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} --source OneLinkCurrent "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        + " ".join(
            f"--source OneLinkCurrent{channel}{index}"
            for channel in ["V", "A"]
            for index in [0, 1, 2]
        ) + " " + " ".join(
            f"--sink {current_type}Current{channel}"
            for channel in ["V", "A"]
            for current_type in ["OneLink", "Conserved"]
        )


rule Z_A:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        plateau_start=get_metadata("plateau_start_Z_A"),
        plateau_end=get_metadata("plateau_end_Z_A"),
        bin_width=get_metadata("bin_width"),
    input:
        data=rules.extract_av.output.data,
        script="s_parameter/vp_momentum/compute_Z_A.py",
    output:
        data=f"processed_data/{subdir_format}/Z_A.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output {output.data} --bin_size {params.bin_width} "
        "--tmin {params.plateau_start} --tmax {params.plateau_end}"


rule Z_A_eff_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.Z_A.output.data.format(**datum)
            for datum in metadata.to_dict(orient="records")
            if datum["plot_Z_A_eff"]
        ],
        script="s_parameter/vp_momentum/plot_Z_A.py",

    output:
        plot="assets/plots/Z_A_eff.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot} "
        "--Z_A_lower_bound 1.025"


rule Z_A_fit:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=[
            rules.Z_A.output.data.format(**datum)
            for datum in metadata.to_dict(orient="records")
            if datum["plot_Z_A_eff"]
        ],
        script="s_parameter/vp_momentum/fit_Z_A.py",
    output:
        data="processed_data/Z_A.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.data}"


rule Z_A_fit_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        extra_args=lambda wildcards, input: (
            " ".join(f"--extra_input_file {filename}" for filename in input.extra_data)
        ),
    input:
        data=[
            rules.Z_A.output.data.format(**datum)
            for datum in metadata.to_dict(orient="records")
            if datum["plot_Z_A_eff"]
        ],
        extra_data=[
            rules.Z_A.output.data.format(**datum)
            for datum in metadata.to_dict(orient="records")
            if datum["plot_Z_A_extra"]
        ],
        fit_result=rules.Z_A_fit.output.data,
        script="s_parameter/vp_momentum/plot_Z_A_fit.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/Z_A_fit.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} {params.extra_args} "
        "--fit_result {input.fit_result} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"
