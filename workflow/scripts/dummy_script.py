"""A simple script to serve as an example.

Should be deleted in real workflows.
"""

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import snakemake

config = snakemake.params.config_text
readme = Path(snakemake.input.readme).read_text()
user = Path(snakemake.input.user_file).read_text()

output_text = "\n\n".join([readme, user, config])

Path(snakemake.output.combined).write_text(output_text)
