Installation
************

In the Sali Lab
===============

If you are working in the Sali lab, you don't need to build and install
MiST - it is already set up for you as a module. Just run
``module load mist`` to load it.

Dependencies
============

* `Python <https://www.python.org>`_ 2.6 or later (Python 3 should work too).

* The `scikit-learn <https://scikit-learn.org/>`_ Python package. The
  ``sklearn`` package is expected to be found in the standard Python
  search path (e.g. ``PYTHONPATH``).

* `nose <https://nose.readthedocs.io/en/latest/>`_ is also needed to run the
  test suite (recommended but not essential).

In the Sali lab, running 
``module load python/scikit`` will get all of these dependencies.

Building
========

Use ``make test`` to test the library, and ``make install`` to install it.
In most cases you will need to tell ``make`` where to install (if running on
a Linux cluster, MiST will need to be installed on a network-accessible
filesystem) with something like
``make PREFIX=/shared/mist install``.
