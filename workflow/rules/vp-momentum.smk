rule extract_vv_aa:
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
        data=f"processed_data/{subdir_format}/vv_aa.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        "--Nf {params.Nf} "
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
        Nf=get_metadata("Nf"),
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
        "--Nf {params.Nf} --sink 1LPBP --source OneLinkCurrent"


rule extract_momentum_currents:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        Nt=get_metadata("Nt"),
        Nx=get_metadata("Nx"),
        Ny=get_metadata("Ny"),
        Nz=get_metadata("Nz"),
        Nf=get_metadata("Nf"),
    input:
        data=f"data/{subdir_format}/correlator.log.zst",
        script="s_parameter/vp_momentum/extract_momentum_correlators.py",
    output:
        data=f"processed_data/{subdir_format}/momentum_currents.json.gz",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output {output.data} "
        "--Nt {params.Nt} --Nx {params.Nx} --Ny {params.Ny} --Nz {params.Nz} "
        "--Nf {params.Nf}"


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
        data=f"processed_data/{subdir_format}/renormalised_pade_fit_{{upper_bound}}.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--renormalised --upper_bound {wildcards.upper_bound} "
        "--output_file {output.data}"


rule pade_fit_bare:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.VPF.output.data,
        script="s_parameter/vp_momentum/fit_pade_vpf.py",
    output:
        data=f"processed_data/{subdir_format}/bare_pade_fit_{{upper_bound}}.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--upper_bound {wildcards.upper_bound} "
        "--output_file {output.data}"


rule pade_fit_systematics:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=lambda wildcards: (
            f"processed_data/{subdir_format}/{{kind}}_pade_fit_{fit_range}.json"
            for fit_range in (
                ["max2", "max3"] if int(wildcards.Nx) > 18 else ["max2", "max3", "1"]
            )
        ),
        script="s_parameter/vp_momentum/compute_s_parameter_systematics.py",
    output:
        data=f"processed_data/{subdir_format}/{{kind}}_pade_fit.json",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.data}"


rule plot_VPF:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=rules.VPF.output.data,
        script="s_parameter/vp_momentum/plot_vpf.py",
        fit_result=f"processed_data/{subdir_format}/bare_pade_fit.json",
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
        fit_result=f"processed_data/{subdir_format}/renormalised_pade_fit.json",
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
        data=lambda wildcards: get_ensemble_data(f"{wildcards.mass_range}_vpf_data", "vpf.json")(wildcards),
        fit_results=lambda wildcards: get_ensemble_data(f"{wildcards.mass_range}_vpf_fit", "bare_pade_fit.json")(wildcards),
        script="s_parameter/vp_momentum/plot_vpf_multi.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/vpf_{mass_range}_ensembles.pdf",
        definitions="assets/definitions/vpf_{mass_range}_ensembles.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} {params.fit_result_flags} "
        "--plot_styles {input.plot_styles} --output_file {output.plot} "
        "--q_squared_upper_bound {params.q_squared_upper_bound} "
        "--output_definitions {output.definitions}"


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


rule plot_fit_range_comparison_large_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data={
            f"processed_data/{subdir_format}/renormalised_pade_fit_{upper_bound}.json".format(
                **metadatum
            )
            for metadatum in metadata.to_dict(orient="records")
            for upper_bound in ["max2", "max3"]
            if metadatum["Nx"] >= 24
        },
        script="s_parameter/vp_momentum/plot_s_parameter{target}.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/S_parameter{target}_mf_fit_range_large.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule plot_fit_range_comparison_small_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data={
            f"processed_data/{subdir_format}/renormalised_pade_fit_{upper_bound}.json".format(
                **metadatum
            )
            for metadatum in metadata.to_dict(orient="records")
            for upper_bound in ["1", "max2", "max3"]
            if metadatum["Nx"] == 18
        },
        script="s_parameter/vp_momentum/plot_s_parameter{target}.py",
        plot_styles=config["plot_styles"],
    output:
        plot="assets/plots/S_parameter{target}_mf_fit_range_small.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.data} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"


rule vp_chisquare_definitions:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        fit_results=get_ensemble_data(f"heavy_vpf_fit", "renormalised_pade_fit.json"),
        script="s_parameter/definitions/pade_chisquare.py",
    output:
        definitions="assets/definitions/pade_chisquare.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.fit_results} "
        "--output_definitions {output.definitions}"


rule lightest_S_definition:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        fit_result=f"processed_data/{subdir_format}/renormalised_pade_fit.json",
        script="s_parameter/definitions/pade_single.py",
    output:
        definitions=f"processed_data/{subdir_format}/pade_result.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} {input.fit_result} "
        "--output_definitions {output.definitions}"
