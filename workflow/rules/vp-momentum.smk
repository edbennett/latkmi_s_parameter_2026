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


rule extract_pbp:
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
        data=f"processed_data/{subdir_format}/pbp.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        "--sink 1LPBP --source OneLinkCurrent"


rule extract_momentum_currents:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        Nt=get_metadata("Nt"),
        Nx=get_metadata("Nx"),
        Ny=get_metadata("Ny"),
        Nz=get_metadata("Nz"),
    input:
        data=f"data/{subdir_format}/correlator.log.zst",
        script="s_parameter/vp_momentum/extract_momentum_correlators.py",
    output:
        data=f"processed_data/{subdir_format}/momentum_currents.json.gz",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz}"


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


rule single_Z_A_eff_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.Z_A.output.data,
        script="s_parameter/vp_momentum/plot_Z_A.py",
        plot_styles=config["plot_styles"],
    output:
        plot=f"processed_data/{subdir_format}/Z_A_eff.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


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
        plot_styles=config["plot_styles"],
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


rule VPF:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        bin_width=get_metadata("bin_width"),
    input:
        data=rules.extract_momentum_currents.output.data,
        script="s_parameter/vp_momentum/compute_vpf.py",
    output:
        data=f"processed_data/{subdir_format}/vpf.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output_file {output.data} --bin_size {params.bin_width}"


rule renormalise_VPF:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        vpf=rules.VPF.output.data,
        Z_A=rules.Z_A.output.data,
        script="s_parameter/vp_momentum/renormalise_vpf.py",
    output:
        data=f"processed_data/{subdir_format}/renormalised_vpf.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.vpf} {input.Z_A} --output_file {output.data}"


rule pade_fit_renormalised:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.renormalise_VPF.output.data,
        script="s_parameter/vp_momentum/fit_pade_vpf.py",
    output:
        data=f"processed_data/{subdir_format}/renormalised_pade_fit.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output_file {output.data} --renormalised"


rule pade_fit_bare:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.VPF.output.data,
        script="s_parameter/vp_momentum/fit_pade_vpf.py",
    output:
        data=f"processed_data/{subdir_format}/bare_pade_fit.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--output_file {output.data}"


rule plot_VPF:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.VPF.output.data,
        script="s_parameter/vp_momentum/plot_vpf.py",
        fit_result=rules.pade_fit_bare.output.data,
        plot_styles=config["plot_styles"],
    output:
        plot=f"processed_data/{subdir_format}/vpf.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --fit_result {input.fit_result} "
        "--plot_styles {input.plot_styles} --output_file {output.plot} "
        "--q_squared_upper_bound 0.3333"


rule plot_VPF_renormalised:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.renormalise_VPF.output.data,
        script="s_parameter/vp_momentum/plot_vpf.py",
        fit_result=rules.pade_fit_renormalised.output.data,
        plot_styles=config["plot_styles"],
    output:
        plot=f"processed_data/{subdir_format}/renormalised_vpf.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --fit_result {input.fit_result} "
        "--plot_styles {input.plot_styles} --output_file {output.plot} "
        "--q_squared_upper_bound 0.3333 --renormalised"


def get_ensemble_data(plot_name, data_filename):
    def get_ensemble_data_inner(wildcards):
        return [
            f"processed_data/{subdir_format}/{data_filename}".format(
                **ensemble_metadata
            )
            for ensemble_metadata in metadata.query(
                    f"plot_{plot_name}"
            ).to_dict(orient="records")
        ]
    return get_ensemble_data_inner


rule plot_VPF_multiple:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        fit_result_flags=lambda wildcards, input: " ".join(
            f"--fit_result {fit_result}" for fit_result in input.fit_results
        ),
        q_squared_upper_bound=lambda wildcards: config["vpf_plot_upper_bounds"][wildcards.mass_range],
    input:
        data=lambda wildcards: get_ensemble_data(f"{wildcards.mass_range}_vpf_data", "renormalised_vpf.json")(wildcards),
        fit_results=lambda wildcards: get_ensemble_data(f"{wildcards.mass_range}_vpf_fit", "renormalised_pade_fit.json")(wildcards),
        script="s_parameter/vp_momentum/plot_vpf_multi.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/vpf_{mass_range}_ensembles.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} {params.fit_result_flags} "
        "--plot_styles {input.plot_styles} --output_file {output.plot} "
        "--q_squared_upper_bound {params.q_squared_upper_bound}"


rule plot_S_parameter:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data={
            f"processed_data/{subdir_format}/renormalised_pade_fit.json".format(
                **metadatum
            )
            for metadatum in metadata.to_dict(orient="records")
        },
        script="s_parameter/vp_momentum/plot_s_parameter.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/S_parameter_mf.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"
