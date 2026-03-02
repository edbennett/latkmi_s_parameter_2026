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
