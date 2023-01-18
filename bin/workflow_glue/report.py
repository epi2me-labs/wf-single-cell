#!/usr/bin/env python
"""Make report."""
import base64
from collections import defaultdict
from pathlib import Path
import re

from dominate.tags import b, figure, img, li, p, ul
import ezcharts
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.plots import util
from ezcharts.util import get_named_logger
import pandas as pd

from .util import wf_parser  # noqa: ABS101


# Setup simple globals
WORKFLOW_NAME = 'wf-single-cell'
REPORT_TITLE = f'{WORKFLOW_NAME} report'
Colors = util.Colors


def umap_plots(img_dirs):
    """Plot representative UMAp plots."""
    plots = defaultdict(dict)
    annotated_plots = defaultdict(list)

    for imgdir in img_dirs:
        sample_id = imgdir.replace('images_', '')
        sample_dir = Path(imgdir)

        for img_path in sample_dir.iterdir():
            umap_match = re.search(
                r'(genes|gene_annotate|transcript|mito)', str(img_path))
            if umap_match:
                umap_type = umap_match.groups()[0]
                if umap_type == 'gene_annotate':
                    annotated_plots[sample_id].append(img_path)
                else:
                    plots[sample_id][umap_type] = img_path

    tabs = Tabs()
    active = True

    for sample_id in plots.keys():
        with tabs.add_tab(sample_id, active):
            active = False
            with Grid(columns=2):
                for dtype in ['genes', 'transcript', 'mito']:
                    plot_path = plots[sample_id][dtype]
                    with open(plot_path, 'rb') as fh:
                        b64img = base64.b64encode(fh.read()).decode()
                    with figure():
                        img(src=f'data:image/png;base64,{b64img}', width=600)
            if len(annotated_plots[sample_id]) > 0:
                b(
                    """Gene expression umaps overlaid with expression levels of
                    genes of interest""")
                with Grid(columns=2):
                    for plot_path in annotated_plots[sample_id]:
                        with open(plot_path, 'rb') as fh:
                            b64img = base64.b64encode(fh.read()).decode()
                        with figure():
                            img(
                                src=f'data:image/png;base64,{b64img}',
                                width=600)


def diagnostic_plots(img_dirs):
    """Diagnostic plots."""
    tabs = Tabs()
    active = True

    for imgdir in img_dirs:
        sample_id = imgdir.replace('images_', '')
        with tabs.add_tab(sample_id, active):
            active = False
            with Grid(columns=2):
                sample_dir = Path(imgdir)
                knee_plot_path = next(
                    sample_dir.glob('*kneeplot.png'))
                saturation_plot = next(
                    sample_dir.glob('*saturation_curves.png'))
                for img_path, width in \
                        [(knee_plot_path, 300), (saturation_plot, 900)]:
                    with open(img_path, 'rb') as fh:
                        b64img = base64.b64encode(fh.read()).decode()
                    with figure():
                        img(src=f'data:image/png;base64,{b64img}', width=width)


def main(args):
    """wf-single-cell report generation."""
    logger = get_named_logger("Report")

    logger.info('Building plots')

    logger.info('Building report')

    # Create report
    report = LabsReport(
        REPORT_TITLE, WORKFLOW_NAME, args.params, args.versions,
        head_resources=[*LAB_head_resources])

    with report.add_section('Read summaries', 'Read summary'):
        SeqSummary(args.read_stats)

    survival_df = pd.read_csv(args.survival, sep='\t', index_col=0)
    wf_summ_df = pd.read_csv(args.wf_summary, sep='\t', index_col=0)

    with report.add_section('Single cell sample summary', 'sc-seq summary'):
        p(
            """This table summarises the number of input reads,
             and the number of cells, genes and transcripts identified within
             each sample."""
        )
        table = DataTable(
            headers=['sample ID',
                     'reads',
                     'cells',
                     'genes',
                     'transcripts',
                     ])
        for row in wf_summ_df.itertuples():
            table.add_row(
                title=row.sample_id,
                columns=[
                    row.n_reads, row.total_cells, row.total_genes,
                    row.total_transcripts])

    with report.add_section('Read surivival by stage', 'Attrition'):
        p(
            """These plots detail the number of remaining reads at different
            stages of the workflow.""")

        # Stage descriptions
        ul(
            li("""full length: Proportion of reads containing adapters in
            expected configurations."""),
            li("""total tagged: Proportion of reads that have been assigned
            corrected UMIs and barcodes."""),
            li("""gene tagged: Proportion of reads assigned to a gene."""),
            li("""transcripts tagged: Proportion of reads assigned a
            transcript.""")
        )

        x_name = 'Workflow stage'
        y_name = 'Proportion of reads [%]'

        order = [
            'full_length', 'total_tagged', 'gene_tagged',
            'transcript_tagged']
        data = survival_df.rename(
            columns={'class': x_name})
        data = data[data[x_name].isin(order)]

        EZChart(
            ezcharts.barplot(
                data=data, x=x_name, y=y_name, hue='sample', order=order),
            theme='epi2melabs')

    with report.add_section('Primer configuration', 'Primers'):
        p(
            """Full length reads are identified by locating read segments
            flanked by known primers in expected orientations:
            adapter1---full_length_read---adapter2.""")

        p(
            """These full length reads can then be
            oriented in the same way and are used in the next stages of the
            workflow.""")
        p(
            """Every library prep will contain some level of artifact reads
            including mis-primed reads and those without adapters.
            These are identified by non-standard primer
            configurations, and are not used for subsequent stages of the
            workflow. The plots here show the proportions of different primer
            configurations within each sample, which can help diagnosing
            library preparation issues.
            The majority of reads should be full_length.""")

        order = [
            'full_length',
            'double_adapter2',
            'single_adapter2',
            'single_adapter1',
            'double_adapter1',
            'no_adapters',
            'other']
        x_name = 'Primer configuration'
        y_name = 'Proportion of reads [%]'
        data = survival_df.rename(columns={'class': x_name})
        data = data[data[x_name].isin(order)]

        EZChart(
            ezcharts.barplot(
                data=data, x=x_name, y=y_name, hue='sample', order=order),
            theme='epi2melabs')

    with report.add_section('Diagnostic plots', 'Diagnostic plots'):

        b('Knee plot')
        p(
            """The knee plot is a quality control for RNA-seq data and
            illustrates the procedure used to filter invalid cells.
            The X-axis represents cells ranked by number of reads
            and the Y-axis reads per barcode. The
            vertical dashed line shows the cutoff. Cells to the right
            of this are assumed to be invalid cells, including dead cells and
            background from empty droplets."""
        )

        b('Saturation plots')
        p(
            """Sequencing saturation provides a view of the amount of library
            complexity that has been captured in the experiment.
            As read depth increases,
            the number of genes and distinct UMIs identified will increase at a
            rate that is dependent on the complexity of the input library.
            A steep slope indicates that new genes or UMIs could still be
            identified by increasing the read coverage. A slope which
            flattens towards higher read coverage indicates that the full
            library complexity is being well captured."""
        )
        ul(
            li("Gene saturation: Genes per cell as a function of depth."),
            li("UMI saturation: UMIs per cell as a function of read depth."),
            li("""Sequencing saturation:  This metric is a measure of the
            proportion of reads that come from a previously observed UMI,
            and is calculated with the following formula:
             1 - (number of unique UMIs / number of reads).""")
        )

        diagnostic_plots(args.images)

    with report.add_section('UMAP projections', 'UMAP'):
        p(
            """This section presents various UMAP projections
            of the data.
            UMAP is an unsupervised algorithm that projects the
            multidimensional single cell expression data into 2
            dimensions. This can reveal structure in the data representing
            different cell types or cells that share common regulatory
            pathways, for example.

            The UMAP algorithm is stochastic, which means running the same data
            multiple times through UMAP with the same settings can lead to
            different projections.
            In order to to have some confidence in the observed clusters,
            it can be useful to run the projection multiple times.
            A single plot for each resulting expression matrix is shown here,
            but 10 replicated plots can be found in the output folder:
            "{out_dir}/{sample_id}/UMAP""")
        umap_plots(args.images)

    report.write(args.output)
    logger.info('Report writing finished')


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")

    parser.add_argument(
        "--read_stats",
        help="fastcat read stats file, with multiple samples concatenated")
    parser.add_argument(
        "--images", nargs='+',
        help="various images to put in report")
    parser.add_argument(
        "--survival",
        help="Read survival data in TSV format")
    parser.add_argument(
        "--wf_summary",
        help="Workflow summary statistics")
    parser.add_argument(
        "--params", help="Workflow params json file")
    parser.add_argument(
        "--versions", help="Workflow versions file")
    parser.add_argument(
        "--output", help="Output HTML file.")
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
