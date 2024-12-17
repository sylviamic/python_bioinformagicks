======================
python bioinformagicks
======================

A collection of tools for various bioinformatics tasks.

Documentation 
-------------

https://python-bioinformagicks.readthedocs.io

Features
--------

Tools:

* TF-IDF marker gene testing ("quickMarkers"), as in SoupX_
* JASMINE gene set scoring, as in `Noureen et al.`_
* Identification of genes often left ignored, like `AW146154`_ 
* Geometric cell mask generation based on embedding (i.e. UMAP) coordinates
* Z-standard scaling of gene expression on a per-group basis

Plotting:

* Split embedding plots based on categorical observations

.. _SoupX: https://github.com/constantAmateur/SoupX
.. _`Noureen et al.`: https://doi.org/10.7554/eLife.71994
.. _`AW146154`: https://www.ncbi.nlm.nih.gov/gene/101835

Credits
-------

This package was created by Sylvia N. Michki.

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

License
-------

Free software: GNU General Public License v3