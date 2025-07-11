"""Make report."""
import json
import math
from pathlib import Path

from bokeh.models import ColorBar, Span
from bokeh.models.tickers import BasicTicker
from bokeh.transform import linear_cmap
from dominate.tags import a, b, div, h4, h6, li, p, ul
from dominate.util import raw, text
import ezcharts
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.plots import BokehPlot, Plot, util
from ezcharts.util import get_named_logger
import matplotlib
import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


# Setup simple globals
WORKFLOW_NAME = 'wf-single-cell'
REPORT_TITLE = f'{WORKFLOW_NAME} report'
Colors = util.Colors


def visium_spatial_plots(non_hd_coords, sample_dirs, hd=False):
    """
    Plot visium spatial expression.

    :param coords_file: path to 10x barcode coordinate file
    :param sample_dirs: per-sample directories containing gene expression data.
    """
    raw(
        """These plots display the raw, UMI deduplicated, expression levels of the
        genes listed in the <i>genes_of_interest</i> file.
        Each point represents a spot on the 10x Visium expression slide.""")
    tabs = Tabs()

    # Add visium message

    for d in sample_dirs:
        sample_id = d.replace('_expression', '')
        sample_dir = Path(d)

        goi_file = sample_dir / 'raw_goi_expression.tsv'

        try:
            goi_df = pd.read_csv(
                goi_file, sep='\t',
                index_col=None, usecols=['gene', 'barcode', 'count'])
        except (pd.errors.EmptyDataError, ValueError):
            goi_df = None
        if goi_df is None or goi_df.empty:
            raw("No data available for the selected genes of interest")
            return
        goi_df = goi_df.sort_index()  # Not needed

        if hd:
            # Get a mapping of barcodes to X and Y coordinates from the barcode names
            goi_df[['X', 'Y']] = (
                goi_df['barcode'].str.extract(r'(\d*)_(\d*)$').astype(int)
            )
        else:
            # Merge the pre-defined coordinates with the expression data
            xy_df = pd.read_csv(
                non_hd_coords,
                delimiter="\t",
                names=["barcode", "X", "Y"]
            ).set_index("barcode")

            goi_df = goi_df.merge(
                xy_df, left_on='barcode', right_index=True, how='left')

        no_data = goi_df.loc[~goi_df.any(axis=1)]
        goi_df = goi_df.loc[goi_df.any(axis=1)]

        genes = goi_df['gene'].unique()

        with tabs.add_dropdown_menu(sample_id, change_header=False):
            for gene_pairs in [genes[i:i+2] for i in range(0, len(genes), 2)]:
                genes = "; ".join(gene_pairs)

                with tabs.add_dropdown_tab(genes):
                    with Grid(columns=2):
                        for gene_name in gene_pairs:
                            df = goi_df.query('gene == @gene_name')
                            plt = BokehPlot(toolbar_location='right')

                            mpl_cmap = (
                                matplotlib.colors.LinearSegmentedColormap.from_list(
                                    "visium_rainbow_cmap",
                                    [
                                        '#2d009e',
                                        '#42cef5',
                                        '#42f545',
                                        '#f2f542',
                                        '#f54242',
                                        '#960000'
                                    ])
                            )
                            pallete = [
                                matplotlib.colors.rgb2hex(mpl_cmap(c))
                                for c in range(mpl_cmap.N)]
                            cmap = linear_cmap(
                                field_name='count',
                                palette=pallete,
                                low=0,
                                high=np.nanmax(
                                    df['count']), nan_color=(50, 50, 50, 0.1))
                            plt._fig.scatter(
                                x='X', y='Y', source=df, fill_color=cmap,
                                line_alpha=0.0, radius=1.0)
                            color_bar = ColorBar(
                                color_mapper=cmap['transform'], width=8, title='count',
                                ticker=BasicTicker(min_interval=1))
                            plt._fig.title.text = gene_name
                            plt._fig.add_layout(color_bar, 'right')
                            EZChart(plt, width='500px', height='500px')
        if len(no_data) > 1:
            raw(f"<br>The following genes were not found in the expression matrix: "
                f"{', '.join(no_data.index)}")


def umap_plots(umaps_dirs, genes_file):
    """Plot UMAPs."""
    sample_tabs = Tabs()

    def _plot(data, title, hue):
        plt = Plot()
        plt.title = dict(text=title)
        plt.title.textStyle = dict(fontSize=14)
        plt.xAxis = dict(name='D1')
        plt.yAxis = dict(name='D2')
        plt.add_dataset(dict(
            source=data.values))

        try:
            plt.visualMap = [
                {
                    'show': True,
                    'text': [
                        str(round(max(data[hue]), 2)), str(round(min(data[hue]), 2))],
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
        except Exception:
            raise
        plt.add_series(dict(
            type='scatter',
            datasetIndex=0,
            symbolSize=2,
            symbol='circle'))
        EZChart(plt, theme='epi2melabs')

    with sample_tabs.add_dropdown_menu('sample', change_header=True):
        # Get data from each sample folder
        for d in umaps_dirs:

            sample_id = d.replace('_umap', '')
            sample_dir = Path(d)

            with sample_tabs.add_dropdown_tab(sample_id):

                repl_tabs = Tabs()

                gene_umap_files = sample_dir.glob(
                    '*gene_expression_umap*.tsv')
                transcript_umap_files = sample_dir.glob(
                    '*transcript_expression_umap*.tsv')
                goi_file = sample_dir / 'raw_goi_expression.tsv'
                try:
                    goi_df = pd.read_csv(
                        goi_file, sep='\t', index_col='barcode',
                        usecols=['gene', 'barcode', 'count'])
                except (pd.errors.EmptyDataError, ValueError):
                    goi_df = None
                for i, (
                    gene_umap_file, transcript_umap_file
                ) in enumerate(zip(gene_umap_files, transcript_umap_files)):

                    with repl_tabs.add_tab(f'umap #{i}'):  # Tab for UMAP repeat

                        gene_umap = pd.read_csv(gene_umap_file, sep='\t', index_col=0)

                        gene_mean_expression = pd.read_csv(
                            sample_dir / 'gene_mean_expression.tsv',
                            sep='\t', index_col=0,
                        )

                        transcript_umap = pd.read_csv(
                            transcript_umap_file, sep='\t', index_col=0)

                        transcript_mean_expression = pd.read_csv(
                            sample_dir / 'transcript_mean_expression.tsv',
                            sep='\t', index_col=0,
                        )

                        mito_expression_file = (
                            sample_dir / "mitochondrial_expression.tsv")
                        mito_expression = pd.read_csv(
                            mito_expression_file, sep='\t', index_col=0)

                        # Make the plots mean gene transcript and mito plots
                        # with Grid(columns=2):

                        gdata = gene_umap.merge(
                            gene_mean_expression, left_index=True, right_index=True)
                        _plot(
                            gdata,
                            title='Gene UMAP / mean gene expression annotation',
                            hue='mean_expression')

                        tdata = transcript_umap.merge(
                            transcript_mean_expression, left_index=True,
                            right_index=True)
                        _plot(
                            tdata,
                            title='Transcript UMAP / '
                                  'mean transcript expression annotation',
                            hue='mean_expression')

                        mdata = gene_umap.merge(
                            mito_expression, left_index=True, right_index=True)
                        _plot(
                            mdata,
                            title='Gene UMAP / mean mito. expression '
                                  'annotation', hue='mito_pct')

                        # Start another grid to put in genes of interest and
                        # optional SNV plot
                        with Grid(columns=2):

                            # Annotate with gene of interest expression levels.
                            # Read single column gene list file
                            genes = pd.read_csv(genes_file, header=None)[0]
                            goi_tabs = Tabs()

                            with goi_tabs.add_dropdown_menu(
                                    'Genes of interest', change_header=True):
                                if goi_df is None:
                                    raw("No genes of interest data available")
                                else:
                                    for gene in genes:
                                        with goi_tabs.add_dropdown_tab(gene):
                                            gene_df = goi_df.query('gene == @gene')
                                            if not gene_df.empty:
                                                goi_data = gene_umap.merge(
                                                    gene_df,
                                                    left_index=True, right_index=True)
                                                # raise ValueError(goi_data)
                                                if goi_data.empty:
                                                    continue
                                                _plot(
                                                    goi_data,
                                                    title=f'Gene UMAP / single gene '
                                                    f'expression, annotation: {gene}',
                                                    hue='count')
                                            else:
                                                with div():
                                                    p(f"No data for {gene}")

                            # If there's a dataframe of top SNVs, create a further UMAP
                            # figure in the grid with a drop-down of SNVs for overlay
                            snv_tabs = Tabs()

                            top_snv = sample_dir / 'top_snvs.tsv'
                            if top_snv.is_file():
                                snv_df = pd.read_csv(
                                    top_snv, sep='\t', dtype=str, index_col=0)
                                snv_df = snv_df.replace({
                                    '-1': 'no_data',
                                    '0': 'hom-ref',
                                    '1': 'het',
                                    '2': 'hom-alt'
                                })

                                # SNV annotation
                                with snv_tabs.add_dropdown_menu(
                                        'SNV', change_header=False):
                                    for i, snv in snv_df.iterrows():
                                        with snv_tabs.add_dropdown_tab(snv.name):
                                            gene_snv_data = gene_umap.merge(
                                                snv, left_index=True, right_index=True)
                                            plt = ezcharts.scatterplot(
                                                data=gene_snv_data,
                                                x='D1', y='D2',
                                                hue=str(snv.name), s=0.3,
                                                fill_alpha=.4, marker=True)
                                            plt._fig.title = (
                                                'Gene umap / genotype '
                                                f'annotation: {snv.name}'
                                            )
                                            EZChart(plt, theme='epi2melabs')


def saturation_plots(seq_data, gene_data, type_label):
    """Saturation plots."""
    tabs = Tabs()

    df_seq = pd.read_csv(seq_data, sep='\t', index_col=None)
    df_gene = pd.read_csv(gene_data, sep='\t', index_col=None)
    with tabs.add_dropdown_menu('sample', change_header=True):
        for sample_id, df_seq_sample in df_seq.groupby('sample'):

            df_seq_sample.rename(columns={
                'n_reads': 'Num. reads',
                'saturation': 'Saturation',
            }, inplace=True)
            df_gene_sample = df_gene.query("`sample` == @sample_id")
            df_gene_sample.rename(columns={
                'n_umis': 'Num. UMIs',
                'median_feats': f'Num. genes / {type_label}',
                'median_umis': f'Num. UMIs / {type_label}'
            }, inplace=True)
            with tabs.add_dropdown_tab(sample_id):
                with Grid(columns=3):
                    seq_sat_plt = ezcharts.lineplot(
                        data=df_seq_sample, x='Num. reads',
                        y='Saturation',
                        title='Sequencing saturation', theme='epi2melabs',
                        color='blue', marker=False)
                    seq_sat_plt._fig.y_range.end = 1.0
                    hline = Span(
                        location=df_seq_sample[
                            'Saturation'].max(), dimension='width',
                        line_color='black',
                        line_width=1, line_dash='dashed')
                    seq_sat_plt._fig.renderers.extend([hline])
                    EZChart(seq_sat_plt, height='350px')

                    EZChart(
                        ezcharts.lineplot(
                            data=df_gene_sample, x='Num. UMIs',
                            y=f'Num. genes / {type_label}',
                            title='Gene saturation', theme='epi2melabs',
                            color='orange', marker=False),
                        height='350px')

                    EZChart(
                        ezcharts.lineplot(
                            data=df_gene_sample, x='Num. UMIs',
                            y=f'Num. UMIs / {type_label}',
                            title='UMI saturation', theme='epi2melabs',
                            color='green', marker=False),
                        height='350px')


# once https://nanoporetech.atlassian.net/browse/CW-5808 is released,
# We can add sortable=False and use_header=False to the DataTables in this section
def experiment_summary_section(
        report, wf_df, aln_df, cell_counts_df, visium, type_label):
    """Three columns experiment summary section."""
    div_style = "padding: 20px;background-color: rgb(245, 245, 245)"
    with report.add_section('Summary', 'Summary'):
        tabs = Tabs()
        with tabs.add_dropdown_menu('sample', change_header=True):
            for sample, sample_df in wf_df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    with Grid(columns=3):
                        # Column 1 - Experiment summary
                        with div(style=div_style):
                            h4("Experiment summary")
                            counts_label = f'{type_label.capitalize()}s with data' \
                                if visium else 'Estimated cells'
                            df_col1 = pd.DataFrame.from_dict(
                                {
                                    'Input reads': sample_df.loc['reads', 'count'],
                                    counts_label:
                                        sample_df.loc['cells', 'count'],
                                    f'Reads per {type_label} (mean)':
                                        sample_df.loc['mean_reads_per_cell', 'count'],
                                    f'UMIs per {type_label} (median)':
                                        sample_df.loc['median_umis_per_cell', 'count'],
                                    f'Genes per {type_label} (median)':
                                        sample_df.loc['median_genes_per_cell', 'count']
                                }, orient='index')
                            df_col1 = df_col1.apply(
                                lambda x: x.apply(lambda y: f"{y:,.0f}"))
                            # When CW-5808 is released we can use sortable=False and
                            # use_headers=False, the latter will mean we don't need the
                            # previous two lines.
                            DataTable.from_pandas(
                                df_col1, use_index=True, paging=False,
                                searchable=False, sortable=False, use_headers=False)
                        # column 2 - rank plot
                        with div(style=div_style):  # barcode rankplot
                            h4('Barcode rank plot')
                            df_col2 = cell_counts_df[
                                cell_counts_df['sample'] == sample]
                            n_cells = int(sample_df.at['cells', 'count'])
                            # Barcodes are ordered by ascending read count
                            # Get the rank of each barcode in descending order
                            df_col2[f'{type_label.capitalize()} barcode rank'] \
                                = np.arange(0, len(df_col2))[::-1]
                            df_col2.rename(
                                {'count': 'Read count'}, axis=1, inplace=True)

                            rank_plt = ezcharts.lineplot(
                                data=df_col2,
                                x=f'{type_label.capitalize()} barcode rank',
                                y='Read count',
                                theme='epi2melabs', line_width=0.5, marker=False,
                                line_color='blue',
                                bokeh_kwargs={
                                    'y_axis_type': 'log', 'x_axis_type': 'log'})
                            rank_plt._fig.y_range.start = 1
                            vline = Span(
                                location=n_cells, dimension='height',
                                line_color='black',
                                line_width=1, line_dash='dashed')
                            rank_plt._fig.renderers.extend([vline])
                            EZChart(rank_plt, height='300px')

                        # Column 3 - Alignment summary
                        with div(style=div_style):
                            h4('Alignment / feature summary')
                            aln_sample = aln_df.loc[sample]
                            df_col3 = pd.DataFrame.from_dict(
                                {
                                    "Pass reads": aln_sample['reads_aligned'],
                                    "Mapped": aln_sample['primary'],
                                    "Unmapped": aln_sample['unmapped'],
                                    "Supplementary": aln_sample['supplementary'],
                                    'Unique genes':
                                        sample_df.loc['genes', 'count'],
                                    'Unique isoforms':
                                        sample_df.loc['transcripts', 'count']
                                }, orient='index')
                            df_col3 = df_col3.apply(
                                lambda x: x.apply(lambda y: f"{y:,.0f}"))
                            # CW-5328 sortable=False, use_headers=False
                            DataTable.from_pandas(
                                df_col3, use_index=True, paging=False, searchable=False,
                                sortable=False, use_headers=False)
                    # Descriptiopns for the threee table columns.
                    with Grid(columns=3):
                        # Text for experiment summary
                        with div(style=div_style):
                            with ul():
                                li(
                                    b("Input reads: "),
                                    "The total number of reads in the input data."
                                )
                                li(
                                    (
                                        b(f"{type_label.capitalize()}s with data: "),
                                        f"The number of {type_label}s identified by \
                                            the workflow."
                                    ) if visium else (
                                        b("Estimated cells: "),
                                        "The estimated number of real cells identified \
                                            by the workflow."
                                    )
                                )
                                li(
                                    b(f"Mean reads per {type_label}: "),
                                    f"The average number of reads per \
                                     {type_label if visium else 'real cell'}."
                                )
                                li(
                                    b(f"Median UMI counts per {type_label}: "),
                                    f"""The median number of unique molecular
                                    identifiers (UMIs) per {type_label}."""
                                )
                                li(
                                    b(f"Median genes per {type_label}: "),
                                    f"""The median number of unique genes identified per
                                     {'real cell' if not visium else type_label}."""
                                )
                        with div(style=div_style):
                            # Text for rank plot
                            if visium:
                                raw(f"""{type_label.capitalize()}s are ranked by read
                                    count in descending order on the x-axis, and the
                                    read count for each barcode is displayed on the
                                    y-axis. Unlike with single cell data, no filtering
                                    of the barcodes is applied.
                                    <br><br>""")
                            else:
                                raw("""
                                    Cells are ranked by read count
                                    in descending order on the x-axis, and the read
                                    count for each barcode is displayed on the y-axis.
                                    Only high quality barcodes are used to generate the
                                    rank plot
                                    (min qscore 15 and 100% match to the 10x whitelist)
                                    <br><br>

                                    The dashed line indicates the read count threshold
                                    that was determined by the workflow. Barcodes to
                                    the left of this point are considered "real cells",
                                    and those to the right are considered as
                                    non-cell barcodes and are not included in the
                                    downstream analysis."""
                                    )
                        with div(style=div_style):
                            # Text for alignment summary
                            with ul():
                                li(
                                    b("Pass reads: "),
                                    """The total number of reads that passed the input
                                    filtering stages of the analysis. This number
                                    excludes reads where the expected adapters were not
                                    found."""
                                )
                                li(
                                    b("Mapped: "),
                                    "The number of primary alignments."
                                )
                                li(
                                    b("Unmapped: "),
                                    """The number of reads that were not mapped to the
                                    reference genome."""
                                )
                                li(
                                    b("Supplementary: "),
                                    """The number of supplementary alignments. These can
                                    be indicative of fusion genes or chimeric reads."""
                                )
                                li(
                                    b("Unique genes/isoforms: "),
                                    f"""The total number of features identified across
                                    all {type_label}s."""
                                )


def fusion_section(report, fusion_results_dir, visium=False):
    """Make gene fusion report section."""
    # Fusion detection is currently done on the final tagged BAM. So for visium HD
    # data, barcodes will be the original non-binned spatial barcodes.

    type_label = 'spatial barcode' if visium else 'cell'

    per_fusion_file = fusion_results_dir / 'fusion_summary.tsv'
    summary_file = fusion_results_dir / 'fusion_per_sample_summary.tsv'

    try:
        fusion_summary_df = (
            pd.read_csv(per_fusion_file, sep='\t', index_col=None)
            .rename(
                columns=lambda col: col[0].upper() + col.replace('_', ' ')[1:]
                if col != 'UMIs' else col)
        )
    except pd.errors.EmptyDataError:
        fusion_summary_df = None

    fusion_sample_df = (
        pd.read_csv(summary_file, sep='\t', index_col=None)
        .rename(columns=lambda x:  x[0].upper() + x.replace('_', ' ')[1:])
    )

    digit_group_cols = ['Cells with fusions', 'Reads', 'Unique fusions']

    # Rename columns to match data type
    if visium:
        fusion_sample_df.rename({
            'Cells with fusions': 'Spatial barcodes with fusions',
            'Mean unique fusions per cell': 'Mean unique fusions per spatial barcode',
            'Mean fusion reads per cell': 'Mean fusion reads per spatial barcode'
        }, inplace=True)
        digit_group_cols = [
            'Spatial barcodes with fusions', 'Reads',
            'Mean unique fusions per spatial barcode']

    round_cols = [
        f'Mean unique fusions per {type_label}', f'Mean fusion reads per {type_label}']
    fusion_sample_df[digit_group_cols] = \
        fusion_sample_df[digit_group_cols].applymap(lambda x: f'{int(x):,}')
    fusion_sample_df[round_cols] = \
        fusion_sample_df[round_cols].applymap(lambda x: f'{float(x):.1f}')

    with report.add_section('Fusion transcripts', 'Fusion transcripts'):

        with p():
            text("Fusion transcript detection has been performed using")
            a('ctat-LR-fusion', href='https://github.com/TrinityCTAT/CTAT-LR-fusion')
            text(". Only reads with valid barcodes are included in these results.")

        with ul():
            li(
                b(f"{type_label.capitalize()}s with fusions:"),
                f" Number of valid {type_label}s with at least one assigned fusion")
            li(b("Unique fusions:"), " Total fusions detected")
            li(b("Reads:"), " The total number of fusion-supporting reads")
            li(
                b(f"Mean fusion reads per {type_label}:"),
                f" The mean number of fusion-supporting reads per {type_label}")
            li(
                b(f"Mean unique fusions per {type_label}:"),
                f" The mean number of unique fusion transcripts per {type_label}")

        DataTable.from_pandas(
            fusion_sample_df, use_index=False, paging=False, searchable=False)

        tabs = Tabs()
        if fusion_summary_df is None or fusion_summary_df.empty:
            return
        with tabs.add_dropdown_menu('sample', change_header=True):
            for sample, sample_df in fusion_summary_df.groupby('Sample ID'):

                if sample_df.empty:
                    h6(f"No fusions detected for {sample}")
                    continue
                fusion_summary_df['UMIs'] = \
                    fusion_summary_df['UMIs'].apply(lambda x: f"{float(x):.2f}")

                if fusion_summary_df.empty:
                    p(f"No fusions detected for {sample}")
                    continue

                with tabs.add_dropdown_tab(sample):
                    with div():
                        p(
                            """The following table summarises fusion transcripts
                            identified in each sample. Each fusion is defined by gene
                            pair partners and the breakpoint locations.""")
                        with div():
                            with ul():
                                li(
                                    b("Left breakpoint: "),
                                    "Breakpoint location of the first gene"
                                )
                                li(
                                    b("Right breakpoint: "),
                                    "Breakpoint location of the second gene"
                                )
                                with li():
                                    b("Splice type:")
                                    with ul():
                                        li(
                                            b("ONLY_REF_SPLICE: "),
                                            "Breakpoint occurs at annotated junction"
                                        )
                                        li(
                                            b("INCL_NON_REF_SPLICE: "),
                                            "Breakpoint is not at an annotated junction"
                                        )
                                li(
                                    b(f"{type_label}s: "),
                                    f"Number of {type_label}s containing this fusion"
                                )
                                li(
                                    b("UMIs: "),
                                    """Number of deduplicated reads supporting this
                                    fusion site"""
                                )

                    DataTable.from_pandas(
                        sample_df.drop(columns=['Sample ID']),
                        use_index=False, paging=False, searchable=False)


def main(args):
    """wf-single-cell report generation."""
    logger = get_named_logger("Report")
    logger.info('Building report')

    if args.visium_hd:
        type_label = '8 µm bin'
    elif args.visium_spatial_coords:
        type_label = 'spot'
    else:
        type_label = 'cell'
    visium = args.visium_spatial_coords or args.visium_hd

    report = LabsReport(
        'Workflow Single Cell report', 'wf-single-cell',
        args.params, args.versions, args.wf_version,
        head_resources=[*LAB_head_resources])

    wf_df = pd.read_csv(args.survival, sep='\t').set_index("statistic")
    df_aln = pd.read_csv(
        args.bam_stats, sep='\t', index_col=0,
        dtype={
            "primary": int,
            "secondary": int,
            "supplementary": int,
            "unmapped": int,
            "total_reads": int,
            "sample": str
        })

    df_aln = df_aln.rename(
        columns={'sample': 'sample ID', "reads aligned": "reads aligned"})

    df_counts = pd.read_csv(
        args.knee_plot_counts, sep='\t',
        usecols=['barcode', 'count', 'sample'], index_col='barcode')

    # Subsample barcode counts if there are too many barcodes,
    max_counts = 10000
    if len(df_counts) > max_counts:
        df_counts = df_counts[::len(df_counts) // max_counts]

    experiment_summary_section(report, wf_df, df_aln, df_counts, visium, type_label)

    logger.info('Report writing finished')

    with open(args.metadata) as metadata:
        sample_details = [{
            'sample': d['alias'],
            'type': d['type'],
            'barcode': d['barcode']
        } for d in json.load(metadata)]

    with report.add_section('Read summary', 'Read summary'):
        names = tuple(d['sample'] for d in sample_details)
        stats = tuple(args.stats)
        if len(stats) == 1:
            stats = stats[0]
            names = names[0]
        fastcat.SeqSummary(
            stats, sample_names=names, height='350px', read_length_quantile_xend=0.99)

    # set statistic as index column to allow for easy selection in building table.
    # For barplots we'll pull it back into a column
    survival_df = pd.read_csv(args.survival, sep='\t').set_index("statistic")

    with report.add_section('Read assignment summary', 'Read assignment'):
        tabs = Tabs()
        with tabs.add_dropdown_menu('sample', change_header=True):
            for sample, sample_df in survival_df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    with Grid(columns=2):
                        x_name = 'Workflow stage'

                        names = {
                            'full_length': 'Full length',
                            'valid_barcodes': 'Valid barcode',
                            'gene_tagged': 'Gene assigned',
                            'transcript_tagged': 'Transcript assigned'
                        }

                        sample_df = sample_df.rename(index=names)
                        order = names.values()
                        data = sample_df.reset_index(names=x_name)
                        data = data[data[x_name].isin(order)]
                        data = data.drop(columns=['pct_of_input_reads', 'sample_id'])
                        data = data.set_index('Workflow stage', drop=True)
                        data['count'] = data['count'].apply(lambda x: f"{int(x):,}")
                        data['pct_of_fl_reads'] = data[
                            'pct_of_fl_reads'].apply(lambda x: f"{float(x):.2f}%")
                        data = data.rename(columns={
                            'pct_of_fl_reads': r"% full length reads"})
                        data = data.transpose()
                        data.index.name = " "
                        data = data[order]

                        DataTable.from_pandas(
                            data, use_index=True, paging=False, searchable=False)

                        with div():
                            with ul():
                                li(
                                    b("Full length: "),
                                    """ Proportion of reads containing adapters in the
                                    expected configuration. Full-length reads are
                                    carried forward in the workflow to attempt to
                                    assign barcode/UMI."""
                                )
                                li(
                                    b("Valid barcodes: "),
                                    f""" Proportion of reads that have been assigned
                                    corrected {'spatial' if visium else 'cell' }
                                    barcodes and UMIs. All reads with
                                    valid barcodes are used in the subsequent stages
                                    of the workflow."""
                                )
                                li(
                                    b("Gene assigned: "),
                                    """Proportion of reads unambiguously assigned to a
                                    gene."""
                                )
                                li(
                                    b("Transcript assigned: "),
                                    """Proportion of reads unambiguously assigned a
                                    transcript."""
                                )

    if args.fusion_results_dir:
        fusion_section(report, args.fusion_results_dir)

    with report.add_section('Adapter configuration', 'Adapter configuration'):
        order = [
            'full_length',
            'single_adapter1', 'double_adapter1',
            'single_adapter2', 'double_adapter2',
            'no_adapters', 'other']
        x_name = 'Primer configuration'
        y_name = 'pct_of_input_reads'
        data = survival_df.reset_index(names=x_name)
        data = data[data[x_name].isin(order)]

        with Grid(columns=2):

            plt = ezcharts.barplot(
                data=data, x=x_name, y=y_name, hue='sample_id', order=order)
            plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
            EZChart(plt, theme='epi2melabs')

            with div():
                raw(
                    """<br><br>Full length reads are defined as those flanked
                    by primers/adapters in the expected orientations:
                    <i> adapter1---cDNA---adapter2</i>.""")

                p(
                    """These full length reads can then be
                    oriented in the same way and are used in the next stages of the
                    workflow. If `full_length_only` is set to `false` reads with all
                    primer configurations are analysed. """)
                p(
                    """Every library prep will contain some level of non-standard
                    adapter configuration artifacts. These are not used for subsequent
                    stages of the workflow. These plots show the proportions of
                    different adapter configurations within each sample,
                    which can help in diagnosing library preparation issues.
                    For most applications, the majority of reads should be
                    full_length.""")

                p(
                    """The adapters used to identify read segments vary slightly
                    between the supported kits. They are:
                    """
                )
                p("3prime, multiome and visium kits:")
                with ul():
                    li("Adapter1: Read1")
                    li("Adapter2: TSO")
                p("5prime kit:")
                with ul():
                    li("Adapter1: Read1")
                    li("Adapter2: Non-Poly(dT) RT primer")

    with report.add_section('Saturation', 'Saturation'):

        p("""Sequencing saturation is an indication of how well the diversity of a
        library has been captured in an experiment. As sequencing depth increases,
        the number of detected genes and unique molecular identifiers (UMIs) will also
        increase at a rate that depends on the complexity of the input library.
        The curve gradient indicates the rate at which new genes or UMIs are being
        recovered; as saturation increases the the curve flattens. All metrics are
        calculated through random sampling of the complete dataset.""")
        with ul():
            li(
                b("Sequencing saturation:"),
                """ The total number of unique cDNA molecules observed having sampled
                sequencing reads. Calculated as
                1 - (number of unique UMIs / number of reads)."""
            )
            li(
                b("Gene saturation:"),
                f""" Unique genes observed per {type_label}. Calculated as a median
                across {type_label}s, after sampling the expression matrix.""")
            li(
                b("UMI saturation:"),
                f""" Unique UMIs observed per {type_label}. Calculated as a median
                across {type_label}s, after sampling the expression matrix.""")

        saturation_plots(args.seq_saturation, args.gene_saturation, type_label)

    # A UMAP plot of x barcodes takes up x MB in the report. Skip this for now
    # UMAPs make less sense for spatial data than for single cell data anyway
    if not args.visium_hd:
        with report.add_section('UMAP projections', 'UMAP'):
            p(
                f"""This section presents various UMAP projections
                of the data.
                UMAP is an unsupervised algorithm that projects the
                multidimensional {'spatial' if visium else 'single cell'}
                expression data into two
                dimensions. This could reveal structure in the data representing
                different cell types or cells that share common regulatory
                pathways, for example.

                The UMAP algorithm is stochastic; analysing the same data
                multiple times with UMAP, using identical parameters, can lead to
                visually different projections.
                In order to have some confidence in the observed results,
                it can be useful to run the projection multiple times and so a series of
                UMAP projections can be viewed below.""")
            umap_plots(args.expr_dirs, args.umap_genes)

    # Check if visium data were analysed
    if args.visium_spatial_coords or args.visium_hd:
        with report.add_section('Visium spatial plots', 'Visium'):
            visium_spatial_plots(
                args.visium_spatial_coords, args.expr_dirs, args.visium_hd)

    report.write(args.report)
    logger.info('Report writing finished')


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='+',
        help="Fastcat per-read stats, ordered as per entries in --metadata.")
    parser.add_argument(
        "--survival",
        help="Read survival data in TSV format")
    parser.add_argument(
        "--bam_stats",
        help="Alignment summary statistics in TSV format")
    parser.add_argument(
        "--params", help="Workflow params json file")
    parser.add_argument(
        "--versions", help="Workflow versions file")
    parser.add_argument(
        "--expr_dirs", nargs='+',
        help="Sample directories containing umap and gene expression files")
    parser.add_argument(
        "--umap_genes", help="File containing list of genes to annotate UMAPs")
    parser.add_argument(
        "--metadata", default='metadata.json', required=True,
        help="sample metadata")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    parser.add_argument(
        "--visium_spatial_coords", default=False, type=Path,
        help='Non-HD barcode to coordinate map file. For non-HD Visium data.')
    parser.add_argument(
        "--visium_hd", action='store_true',
        help='Data is visium.')
    parser.add_argument(
        "--q_filtered", action='store_true',
        help="True if the input reads were subject to min read quality filtering.")
    parser.add_argument(
        "--seq_saturation", type=Path)
    parser.add_argument(
        "--gene_saturation", type=Path)
    parser.add_argument(
        "--knee_plot_counts", type=Path)
    parser.add_argument(
        '--fusion_results_dir', type=Path, default=None,
        help="Fusion summary directory")
    return parser
