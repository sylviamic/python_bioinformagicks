#!/usr/bin/env python

import pytest
import os

import python_bioinformagicks as bim


@pytest.fixture
def sample_adata():
    """
    
    Loads an ad.AnnData object from disk to use
    in tests. This is data from multiple publications
    that has been integrated and annotated before 
    subsetting to 5% of the cells and the top
    2000 highly variable genes.
    
    """

    import anndata as ad

    return ad.read_h5ad("./tests/data/test_adata.h5ad")


def test_plot_grouped_proportions(
    sample_adata,
    factor_to_plot = "celltype",
    split_by = "age"
):
    """

    Testing bim.pl.plot_grouped_proportions
    
    Will generate two barplots next to one another.

    TODO: implement image-similarity-based testing.
    """

    fig = bim.pl.plot_grouped_proportions(
        sample_adata,
        factor_to_plot = factor_to_plot,
        split_by = split_by
    )

    save_path = "./tests/data/plot_grouped_proportions.png"
    fig.savefig(save_path, dpi=300)

    assert (os.path.getsize(save_path) > 100000)


def test_plot_split_embedding(
    sample_adata,
    color = ["Actb", "celltype"],
    groupby = "paper",
    use_rep = "X_umap"
):
    """

    Testing bim.pl.plot_split_embedding
    
    Will generate 3 sets of 2 plots, because the 
    sample_adata object has three papers and we are providing
    two items to color by.

    TODO: implement image-similarity-based testing.
    """

    fig = bim.pl.plot_split_embedding(
        sample_adata,
        color = color,
        groupby = groupby,
        use_rep = use_rep,
        size = 10,
        frameon = False
    )

    save_path = "./tests/data/plot_split_embedding.png"
    fig.savefig(save_path, dpi=150)

    assert (os.path.getsize(save_path) > 100000)