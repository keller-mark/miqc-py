# miqc-py

Python implementation of https://github.com/greenelab/miQC.

## Installation

```sh
pip install miqc_py
```

## Usage

Usage follows the scverse API conventions.
Parameter names follow the R implementation of `miQC`.

```python
import miqc_py

# ...

miqc_py.calculate_miqc(adata)
miqc_py.filter_cells(adata)
```

### Plotting

Optionally, we can plot the results with `altair` (not a dependency of `miqc_py` - may need to install first).

```python
import altair as alt

alt.Chart(adata.obs).mark_circle().encode(
    x="total_genes_by_counts:Q",
    y="pct_counts_mito:Q",
    color="prob_compromised:Q"
)
```


## Development

```sh
conda env create -f environment.yml
conda activate miqc-py
```

```sh
jupyter lab
```