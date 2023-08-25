#!/usr/bin/env python
"""Make report."""
import base64
from pathlib import Path

from dominate.tags import b, figure, h6, img, li, p, ul
import ezcharts
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.plots import Plot, util
from ezcharts.util import get_named_logger
import pandas as pd

from .util import wf_parser  # noqa: ABS101


# Setup simple globals
WORKFLOW_NAME = 'wf-single-cell'
REPORT_TITLE = f'{WORKFLOW_NAME} report'
Colors = util.Colors


def umap_plots(umaps_dirs, genes_file):
    """Plot UMAPs."""
    # Raad single column gene list file
    genes = pd.read_csv(genes_file, header=None)[0]

    sample_tabs = Tabs()

    def _plot(data, title, hue):
        plt = Plot()
        plt.title = dict(text=title)
        plt.title.textStyle = dict(fontSize=14)
        plt.xAxis = dict(name='D1')
        plt.yAxis = dict(name='D2')
        plt.add_dataset(dict(
            source=data.values))

        plt.visualMap = [
            {
                'show': True,
                'text': [str(round(max(data[hue]), 2)), str(round(min(data[hue]), 2))],
                'left': 'right',
                'type': "continuous",
                'min': min(data[hue]),
                'max': max(data[hue]),
                'inRange': {
                    'color': ['blue', 'green', 'yellow',  'red'],
                    'opacity': 0.7,

                },
                'dimension': 2
            }]
        plt.add_series(dict(
            type='scatter',
            datasetIndex=0,
            symbolSize=2,
            symbol='circle'))
        EZChart(plt, theme='epi2melabs')

    # Get data from each sample folder
    for d in umaps_dirs:

        sample_id = d.replace('_umap', '')
        sample_dir = Path(d)

        with sample_tabs.add_tab(sample_id):

            repl_tabs = Tabs()

            gene_umap_files = sample_dir.glob('*gene_umap*')
            transcript_umap_files = sample_dir.glob('*transcript_umap*')

            for i, (
                gene_umap_file, transcript_umap_file
            ) in enumerate(zip(gene_umap_files, transcript_umap_files)):

                with repl_tabs.add_tab(f'umap #{i}'):

                    gene_umap = pd.read_csv(gene_umap_file, sep='\t', index_col=0)
                    gene_expression_file = sample_dir / 'gene_expression'
                    gene_expression = pd.read_csv(
                        gene_expression_file, sep='\t', index_col=0)
                    gene_mean_expression = pd.DataFrame(
                        gene_expression.mean(axis=0), columns=['gene mean expression'])

                    transcript_umap = pd.read_csv(
                        transcript_umap_file, sep='\t', index_col=0)
                    transcript_expression_file = sample_dir / 'transcript_expression'
                    transcript_expression = pd.read_csv(
                        transcript_expression_file, sep='\t', index_col=0)
                    transcript_mean_expression = pd.DataFrame(
                        transcript_expression.mean(axis=0),
                        columns=['transcript mean expression'])

                    mito_expression_file = sample_dir / "mitochondrial_expression"
                    mito_expression = pd.read_csv(
                        mito_expression_file, sep='\t', index_col=0)

                    # Make the plots
                    with Grid(columns=2):

                        gdata = gene_umap.merge(
                            gene_mean_expression, left_index=True, right_index=True)
                        _plot(
                            gdata, title='Gene UMAP / mean gene expression annotation',
                            hue='gene mean expression')

                        tdata = transcript_umap.merge(
                            transcript_mean_expression, left_index=True,
                            right_index=True)
                        _plot(
                            tdata,
                            title='Transcript UMAP / '
                                  'mean transcript expression annotation',
                            hue='transcript mean expression')

                        mdata = gene_umap.merge(
                            mito_expression, left_index=True, right_index=True)
                        _plot(
                            mdata,
                            title='Gene UMAP / mean mito. expression '
                                  'annotation', hue='mito_pct')

                        # Plot expresion levels from a single gene over gene UMAP.
                        for gene in genes:
                            if gene in gene_expression.index:
                                go_data = gene_umap.merge(
                                    gene_expression.loc[gene],
                                    left_index=True, right_index=True)
                                _plot(
                                    go_data,
                                    title=f'Gene umap / single gene expression '
                                          f'annotation: {gene}',
                                    hue=gene)
                            else:
                                h6(f'{gene} not in dataset / has been filtered out')


def diagnostic_plots(img_dirs):
    """Diagnostic plots."""
    tabs = Tabs()

    for imgdir in img_dirs:
        sample_id = imgdir.replace('images_', '')
        with tabs.add_tab(sample_id):
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

    with report.add_section('Read survival by stage', 'Attrition'):
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

    # Check if umaps were skipped
    if not (Path(args.umap_dirs[0]) / 'OPTIONAL_FILE').is_file():
        with report.add_section('UMAP projections', 'UMAP'):
            p(
                """This section presents various UMAP projections
                of the data.
                UMAP is an unsupervised algorithm that projects the
                multidimensional single cell expression data into 2
                dimensions. This could reveal structure in the data representing
                different cell types or cells that share common regulatory
                pathways, for example.

                The UMAP algorithm is stochastic; analysing the same data
                multiple times with UMAP, using identical parameters, can lead to
                visually different projections.
                In order to have some confidence in the observed results,
                it can be useful to run the projection multiple times and so a series of
                UMAP projections can be viewed below.""")
            umap_plots(args.umap_dirs, args.umap_genes)
    else:
        logger.info('Skipping UMAP plotting')

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
        help="Sample directories containing various images to put in report")
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
        "--output", help="Output HTML file")
    parser.add_argument(
        "--umap_dirs", nargs='+',
        help="Sample directories containing umap and gene expression files")
    parser.add_argument(
        "--gene_expression", help="Paths to UMAP TSVs",
        nargs='+')
    parser.add_argument(
        "--umap_genes", help="File containing list of genes to annnotate UMAPs")
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
