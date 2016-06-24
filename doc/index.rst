.. coolvib documentation master file, created by
   sphinx-quickstart on Fri Aug 28 13:36:54 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ./coolvib.png
   :align: center
   :width: 85%

Welcome to coolvib's documentation!
===================================

**coolvib** is a python package that calculates quantities related 
to electronic friction acting on adsorbate atoms and molecules on metal surfaces. 
This friction is induced by non-adiabatic coupling between 
electronic and nuclear motion. Such coupling mediates energy loss of adsorbate vibrations 
and induces vibrational cooling - therefore the name **coolvib**.

Features
========

Current Version
^^^^^^^^^^^^^^^^^^^

    * Individual base routines for calculation of spectral function and friction tensor
    * Uses input from `FHI-Aims <www.fhi-berlin.mpg.de/aims/>`_ and `Siesta <icmab.es/siesta/>`_ codes.
    * Convenient workflow for calculation of cartesian friction tensor
    * Useful analysis and visualization tools

Planned for future versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

    * non-momentum conserving spectral function
    * workflow class for calculation of friction along given mode
    * tool that facilitates finite difference calculations
    * parser and routines for `CASTEP <www.castep.org>`_
    * extension of tutorial
    * test suite for main routines

Contents
=========

.. toctree::
   :maxdepth: 2


    Introduction<intro>
    Theory and Methods<methods>
    Tutorial<tutorial>
    Coolvib: Modules and Classes<coolvib>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

