rule dummy_add_text:
    message:
        "Dummy rule combining user inputs and automatic downloads."
    params:
        config_text=config["dummy_text"],
    input:
        user_file="resources/user/user_message.md",
        readme="resources/automatic/dummy_readme.md",
    output:
        combined="results/combined_text.md",
    conda:
        "../envs/shell.yaml"
    script:
        "../scripts/dummy_script.py"
