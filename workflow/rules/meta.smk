rule ensemble_table:
    input:
        metadata=config["metadata_file"],
        script="s_parameter/meta/tabulate_ensembles.py",
    output:
        table="assets/tables/ensembles.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python {input.script} {input.metadata} --output_file {output.table}"


rule pv_spectrum_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        metadata=config["metadata_file"],
        spectrum=config["spectrum_file"],
        script="s_parameter/meta/tabulate_spectrum.py",
    output:
        table="assets/tables/quoted_spectrum.tex",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} "
        "--input_metadata {input.metadata} --input_spectrum {input.spectrum} "
        "--output_file {output.table}"


rule pv_spectrum_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        metadata=config["metadata_file"],
        spectrum=config["spectrum_file"],
        plot_styles=config["plot_styles"],
        script="s_parameter/meta/plot_pv_spectrum.py",
    output:
        plot="assets/plots/pv_spectrum.pdf",
    conda:
        "../envs/python.yml"
    shell:
        "python -m {params.module} "
        "--input_metadata {input.metadata} --input_spectrum {input.spectrum} "
        "--plot_styles {input.plot_styles} --output_file {output.plot}"
