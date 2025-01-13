from gprofiler import GProfiler

def do_gprofiler_analysis(
    genes: list[str], 
    organism: str = "hsapiens",
    max_term_size: int = 10000,
    ordered: bool = False,
    sources: list[str] = ["GO:BP", "GO:MF", "REAC", "KEGG", "TF"]
):
    """
    Given an ordered list of genes,
    performs an over-representation analysis (ORA) 
    using the gProfiler API.
    
    Parameters
    ----------
    
    genes: list of str
        List of genes to query, optionally sorted by user
        based on FDR or another significance metric.
    
    organism: str (default: "hsapiens")
        The organism ID to use for the query.
        See https://biit.cs.ut.ee/gprofiler/page/organism-list
        for a list of organism IDs supported
        by gProfiler. 
        Commonly one of: ["hsapiens", "mmusculus"]
    
    max_term_size: uint (default: 10000)
        Filter results table to only include those
        results with term sizes less than this 
        maximum. Set to None to disable.
    
    ordered: bool (default: False)
        If the list of genes is ordered by
        descending significance, set to True, 
        otherwise set to False. A slightly different 
        ORA is performed on ordered gene lists. 
    
    sources: list of str (default: ["GO\\:BP", "GO\\:MF", "REAC", "KEGG", "TF"])
        The list of source databases to consider for 
        gProfiler ORA. Some may only be available for
        certain organisms; see organism list page in
        the gProfiler documentation for more
        information.
    
    Returns
    -------
    ora_df: pandas.DataFrame
        The resulting gProfiler results dataframe.

    References
    ----------

    * https://biit.cs.ut.ee/gprofiler/gost
    * https://doi.org/10.1093/nar/gkad347

    Usage
    -----

    >>> diff_exp_genes = ["SFTPC", "LAMP3", "NAPSA", "SFTPB", "EPCAM", "COL1A1"]
    >>> ora_df = do_gprofiler_analysis(diff_exp_genes, ordered=False)
    >>> print(ora_df.head(1)[["source", "name", "p_value", "intersections"]])
        source                 name   p_value          intersections
    0   REAC  Surfactant metabolism  0.000033  [SFTPC, NAPSA, SFTPB]

    """
    
    gp = GProfiler(return_dataframe=True)
    
    ora_df = gp.profile(
        organism=organism,
        query=genes,
        no_evidences=False,
        ordered=ordered,
        sources=sources
    )
    
    if (max_term_size) & (max_term_size > 0):
        ora_df = ora_df[ora_df["term_size"] < max_term_size]
    
    return ora_df