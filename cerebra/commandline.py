# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401
from functools import partial

import click

from cerebra.hello import hello
from cerebra.datadump import s3_import
from cerebra.germline_filter_TE import germline_filter_TE
from cerebra.get_mutationcounts_table_TE import get_mutationcounts_table_TE
from cerebra.get_specific_mutations_TE import get_specific_mutations_TE
from cerebra.get_mutationalburden import get_mutationalburden
from cerebra.generate_summary_tables_TE import generate_summary_tables_TE
from cerebra.fusion_search import fusion_search
from cerebra.fusions_x_cell import fusions_x_cell
from cerebra.check_coverage import check_coverage

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=['-h', '--help'])

@click.group(options_metavar='', subcommand_metavar='<command>',
             context_settings=settings)
def cli():
    """
    finds mutants in your scRNA-seq experiment
    """
    pass

cli.add_command(hello, name='hello')
cli.add_command(s3_import, name='s3_import')
cli.add_command(germline_filter_TE, name='germline_filter_TE')
cli.add_command(get_mutationcounts_table_TE, name='get_mutationcounts_table_TE')
cli.add_command(get_specific_mutations_TE, name='get_specific_mutations_TE')
cli.add_command(get_mutationalburden, name='get_mutationalburden')
cli.add_command(generate_summary_tables_TE, name='generate_summary_tables_TE')
cli.add_command(fusion_search, name='fusion_search')
cli.add_command(fusions_x_cell, name='fusions_x_cell')
cli.add_command(check_coverage, name='check_coverage')

if __name__ == "__main__":
    cli()
