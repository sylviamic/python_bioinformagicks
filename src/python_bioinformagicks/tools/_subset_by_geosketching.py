import anndata as ad
from geosketch import gs


def subset_by_geosketching(
    adata: ad.AnnData, 
    n_cells_to_keep: int = None, 
    frac_cells_to_keep: float = 0.33, 
    use_rep: str = "X_pca"
):
    """
    Subsamples a single-cell dataset using geometric sketching
    to more fairly represent rare and common cell types
    
    Defaults to keeping a fraction of cells. To use an 
    exact target cell number, specify `n_cells_to_keep`, which
    will take priority over `frac_cells_to_keep`. 

    Parameters
    ----------
    adata: ad.AnnData 
        The original anndata object, with `use_rep` calculated
        and stored in the :code:`adata.obsm[use_rep]` slot.
    
    n_cells_to_keep: int (default: None)
        The number of cells to keep in the subsampled
        adata object. Must be less than :code:`len(adata.obs.index)`.
    
    frac_cells_to_keep: float (0,1)
        Fraction of cells to keep; overridden by 
        `n_cells_to_keep` when not None.
    
    use_rep: str (default: "X_pca")
        The latent space representation to use.
        Must be a key in :code:`adata.obsm`.
    
    Returns
    -------
    mask: Anndata object
        A boolean mask of length len(adata.obs.index) indicating
        which cells to keep after geosketching. 
    
    """
    
    n_cells = 0
    if not (n_cells_to_keep is None):
        if ((n_cells_to_keep < len(adata.obs.index)) and (n_cells_to_keep > 0)):
            n_cells = n_cells_to_keep
        else:
            print("[ERROR] n_cells_to_keep specified, but not in range (0, len(adata.obs.index))")
            print("[WARN] defaulting to 33% of cells")
    if (n_cells == 0):
        if ((frac_cells_to_keep >= 1) or (frac_cells_to_keep <= 0)):
            print("[ERROR] frac_cells_to_keep is not in range (0,1)")
            print("[WARN] defaulting to 33% of cells")
            frac_cells_to_keep = 0.33
        n_cells = int(frac_cells_to_keep * len(adata.obs.index))

    X_orig = adata.obsm[use_rep]
    mask = gs(X_orig, n_cells, replace=False)
    
    return mask