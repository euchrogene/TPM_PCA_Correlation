# TPM_PCA_Correlation

*EuchroGene TPM_PCA_Correlation v1.0, for EuchroGene members.*

Sample-level ordination and clustering of a gene expression matrix. The pipeline takes the TPM table produced by the EuchroGene STAR-RSEM or Salmon quantification pipelines (or any CSV with genes in rows and samples in columns), filters and log-transforms it, and draws the two sample-level figures a transcriptome manuscript is expected to carry: a principal component analysis, and a sample-to-sample correlation heatmap ordered by hierarchical clustering with its dendrogram. The analyses run on [scikit-learn](https://jmlr.org/papers/v12/pedregosa11a.html) and [SciPy](https://doi.org/10.1038/s41592-019-0686-2), and the figures are rendered with [matplotlib](https://doi.org/10.1109/MCSE.2007.55) at journal resolution in vector formats whose text stays editable in Illustrator.

Replicate groups are found in the sample names rather than in a design file: the names are sorted, consecutive names are compared, and names that differ by a single character are placed in the same group, which then gets its own colour and marker across every panel and its own band in the strip beside the heatmap. The run reports how closely replicates agree, whether each group came out as a single block of the tree, and the cophenetic correlation of the clustering, so the figure works as a quality check as well as a printable panel. A design file can be supplied when the naming scheme is irregular.

---

## Installation

**0. Install EG_tools**

```bash
wget https://github.com/euchrogene/EG_tools/raw/main/EG_tools
chmod 777 EG_tools
sudo mv EG_tools /usr/bin
```

**1. Install**

```bash
sudo EG_tools install \
    -r https://github.com/euchrogene/TPM_PCA_Correlation.git \
    -d TPM_PCA_Correlation \
    -e TPM_PCA_Correlation_v.1.0 \
    -m "PCA ordination and correlation clustering of a gene expression TPM matrix"
```

**2. Display installed software**

```bash
EG_tools
```

**3. Show help contents**

```bash
TPM_PCA_Correlation_v.1.0
```

**4. Uninstall**

```bash
sudo EG_tools uninstall -t TPM_PCA_Correlation_v.1.0 -i managene7/tpm-pca-correlation:v1.0
```

## Quick Start

```bash
TPM_PCA_Correlation_v.1.0 \
    -input gene_TPM.csv \
    -species "Vaccinium virgatum"
```

Results land in `gene_TPM_PCA_Correlation_result` beside the input file, and the ordination points carry no sample names, only the group legend.

Sample names beside every point, which is what you want while checking that the replicates behave:

```bash
TPM_PCA_Correlation_v.1.0 \
    -input gene_TPM.csv \
    -out PCA_correlation_result \
    -label true
```

Spearman correlation, complete linkage, correlation values printed in the cells:

```bash
TPM_PCA_Correlation_v.1.0 \
    -input gene_TPM.csv \
    -out PCA_correlation_result \
    -correlation spearman \
    -linkage complete \
    -annotate true
```

Irregular sample names, replicate assignment supplied by hand:

```bash
TPM_PCA_Correlation_v.1.0 \
    -input gene_TPM.csv \
    -out PCA_correlation_result \
    -group_file sample_groups.tsv \
    -top_genes 0
```

---

## Inputs

| Input | Required | Notes |
|---|---|---|
| `-input` | yes | CSV with the gene ID in the first column and one TPM column per sample. The `Gene_ID` + TPM layout written by `RNA_seq_to_TPM_STAR_v.2.0` and the Salmon pipeline is read directly. Non-numeric annotation columns are dropped with a warning, duplicated gene IDs are collapsed by maximum, missing values become 0. |
| `-out` | no | Results folder, created if it does not exist. Defaults to the input file name with `_PCA_Correlation_result` appended, written next to the input file, so `gene_TPM.csv` gives `gene_TPM_PCA_Correlation_result`. |
| `-group_file` | no | Two-column TSV, `sample<TAB>group`, one line per sample. Overrides name-based grouping. Use it when the sample names do not encode the design, for example `Sample_1` against `Sample_10`, where the number of digits differs. |
| `-species` | no | Organism name. It is written into the Methods paragraph, so set it for anything headed to a manuscript. |

---

## Outputs

```
PCA_correlation_result/
├── report.html                          self-contained report, figures inline, Methods included
├── RUN_SPEC.json                        parameters, tool versions, input MD5, worker SHA-256
├── figures/
│   ├── PCA_PC1_PC2.pdf|png|svg          PCA ordination, one colour per replicate group
│   ├── PCA_scree.pdf|png|svg            explained variance per component
│   ├── Sample_correlation_heatmap.*     correlation heatmap, clustered, with dendrogram
│   ├── Sample_dendrogram.*              the tree alone, leaf labels in the group colours
│   └── PCA_correlation_combined.*       PCA and heatmap side by side, manuscript ready
├── tables/
│   ├── sample_groups.tsv                detected replicate groups with colour and marker
│   ├── pca_coordinates.tsv              sample scores on every principal component
│   ├── pca_explained_variance.tsv       eigenvalues, explained and cumulative variance
│   ├── pca_top_loadings.tsv             genes with the largest loading on the plotted axes
│   ├── sample_correlation_matrix.tsv    the full sample-by-sample correlation matrix
│   ├── sample_cluster_order.tsv         dendrogram leaf order with the group of each sample
│   ├── sample_dendrogram.newick         the tree in Newick form, for iTOL, FigTree or ggtree
│   └── filtered_expression_matrix.tsv   the transformed matrix the ordination ran on
└── .euchrogene/
    ├── run_container.sh                 the generated container command, kept for the record
    ├── analysis_summary.json            what the worker measured, input to the report
    └── tool_versions.json               versions queried inside the container at run time
```

Use `figures/PCA_correlation_combined.pdf` for the manuscript and `tables/sample_groups.tsv` to confirm that replicate detection matched the experimental design.

---

## Key Options

| Option | Default | Purpose |
|---|---|---|
| `-correlation` | `pearson` | Correlation between samples. `spearman` is the rank-based alternative, less sensitive to a handful of very highly expressed genes. |
| `-cor_genes` | `all` | Which genes the correlation matrix runs on. `all` uses every gene that passed the expression filter, which is how published sample correlation heatmaps are computed and what keeps the values comparable with them. `top` uses the same variable subset as the ordination, which lowers the correlations and sharpens the contrast between groups. |
| `-linkage` | `average` | Hierarchical clustering on 1 - r distances. `average` is UPGMA, the usual choice for correlation distances; `complete` gives compact clusters; `single` is prone to chaining; `ward` is popular but defined for Euclidean distances, so treat it as a display choice here. |
| `-annotate` | `auto` | Print the correlation value in every cell. `auto` prints them up to 12 samples and widens the cells so the numbers fit. |
| `-heatmap_palette` | `auto` | `auto` picks `sequential` (cream through red to black) while every correlation is positive, and `diverging` (deep blue through cream at zero to deep red) as soon as the matrix holds a negative. Either can be forced, and `viridis`, `magma` and `blues` remain available. |
| `-cor_min` | `auto` | Lower end of the colour scale. `auto` starts at the lowest correlation observed, which is what makes the block structure visible when every sample correlates above 0.9. Under the diverging ramp the scale stays symmetric about zero, so this sets both ends at once. |
| `-dendrogram` | `true` | Draw the clustering tree above the heatmap. |
| `-group_mode` | `chain` | `chain` groups consecutive sorted names that differ by one character. `position` additionally requires every member of a group to differ at the same character position, which keeps `S1_R2` and `S2_R2` apart when a replicate is missing. `none` treats every sample as its own group. |
| `-allow_indel` | `false` | Also group names of unequal length when a single insertion or deletion explains the difference, for example `Sample_9` against `Sample_10`. |
| `-min_tpm` / `-min_samples` | `1.0` / `1` | A gene is kept when its TPM reaches `-min_tpm` in at least `-min_samples` samples. |
| `-top_genes` | `2000` | Number of most variable genes carried into the ordination. `0` uses every gene that passed the filter. The correlation heatmap is unaffected unless `-cor_genes top`. |
| `-log2` / `-scale` | `true` / `false` | log2(TPM + 1) transformation, then optional gene-wise z-scoring. Scaling gives low-expressed genes the same weight as high-expressed ones. |
| `-label` | `false` | `false` draws points only and leaves the legend to identify groups, `true` writes the sample name next to every point, `group` writes one name per replicate group at the group centre. |
| `-legend` | `true` | Draw the group legend beside the panel. Turn it off when `-label group` already names the groups, and leave it on at the default `-label false`, where it is the only thing identifying a group. |
| `-ellipse` / `-connect` | `false` / `false` | 95 percent confidence ellipse per group, and thin lines joining the replicates of a group. |
| `-pc_x` / `-pc_y` | `1` / `2` | Which components the PCA panel shows. |
| `-formats` / `-dpi` / `-width` | `pdf,png,svg` / `600` / `89` | Figure formats, raster resolution, and panel width in mm (89 is one journal column, 183 is two). |
| `-cores` / `-memory` | `8` / `64g` | Resources handed to the container. The ceilings are 100 cores and 500g. |

---

## Notes

- Replicate detection is a heuristic over sample names, not a design file. Check `tables/sample_groups.tsv` on the first run of a new dataset. The chain rule can merge two conditions when a replicate is missing and the two remaining names happen to differ by one character, which is what `-group_mode position` guards against.
- The two analyses deliberately run on different gene sets. Variance selection sharpens the ordination, but it strips the shared expression baseline that makes a correlation heatmap comparable with the ones printed in other papers, so the correlations use every expressed gene by default. On a typical dataset the same samples come out around r = 0.26 to 0.94 over all genes and around r = -0.04 to 0.96 over the top 2000, which is the same data telling two different stories. `-cor_genes top` restores the old single-gene-set behaviour.
- The heatmap is the quality check most reviewers look for. Replicates should correlate more closely with each other than with anything else, and each group should come out as one block of the tree. When either fails the run says so on the terminal, and the Methods paragraph reports the numbers the run actually produced rather than the ones it expected.
- The cophenetic correlation is reported so the tree can be judged rather than trusted. A value well below about 0.8 means the tree is a poor summary of the underlying distances, whatever it looks like.
- Leaf order comes from optimal leaf ordering, which rotates branches to put similar samples next to each other. It changes the drawing order only, never the topology.
- Sample names are off by default in the ordination, so a run with many samples stays readable. Turn them on with `-label true` while checking replicate behaviour, then turn them off again for the manuscript figure. The heatmap always labels its rows and columns, since an unlabelled correlation matrix is useless.
- Sample labels in the ordination are placed by a deterministic collision solver, so rerunning the same command reproduces the same figure. `-label_repel adjusttext` switches to the adjustText package that ships in the image, `-label_repel none` puts every label exactly on its point.
- The heatmap ramp follows the data rather than a fixed default. All-positive matrices get the sequential ramp, where the low end is a pale cream rather than pure white so a low cell never reads as an empty cell. A matrix holding negative correlations gets the diverging ramp on a scale symmetric about zero, so cream falls exactly on r = 0 and the sign change is visible. Negatives are normal under `-scale true`, where gene-wise centring fixes the mean off-diagonal correlation at -1/(n - 1).
- Numbers printed in the cells switch between near-black and white according to the WCAG relative luminance of the cell behind them, so they stay legible from the cream end to the black end of the ramp.
- PDF and SVG output keeps text as editable TrueType (`fonttype 42`) and draws the heatmap as vector cells rather than an embedded raster, so every panel can be restyled in Illustrator or Inkscape without re-running the pipeline.
- The image ships Liberation Sans, which carries the metrics of Arial. `-font` accepts any face installed in the image, and falls back to the first installed sans-serif when the requested one is absent.

---

## Citation

> Pedregosa F, Varoquaux G, Gramfort A, et al. (2011) Scikit-learn: machine learning in Python. *Journal of Machine Learning Research* 12: 2825-2830.

> Virtanen P, Gommers R, Oliphant TE, et al. (2020) SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17: 261-272.

> Hunter JD (2007) Matplotlib: a 2D graphics environment. *Computing in Science & Engineering* 9: 90-95.

> EuchroGene TPM_PCA_Correlation Pipeline v1.0 (2026). EuchroGene, LLC.

**Support:** bioinformatics@euchrogene.com
