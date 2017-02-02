.. File Output.rst

.. _output:

Output Catalogues
=================

|SCAMP| is run from the shell with the following syntax:

.. _mergedcat:
Merged Catalogues
-----------------

In addition to astrometric header information, |SCAMP| can write out “merged”
catalogues (one per field group). These catalogues contain the calibrated coordinates
and magnitudes of a union of all detections from input catalogues that passed the |SCAMP|
acceptation criteria (S/N, flags).
Merged coordinates and magnitudes are computed using a weighted average of all overlapping
measurements, and are accompanied by estimates of formal errors and 1-sigma uncertainties
on individual measurements.
The ``MERGEDOUTCAT_TYPE`` configuration parameter sets the format of the merged catalogue;
it is set to ``NONE`` by default, which means that no catalogue is created.
The available formats are ``ASCII`` (pure ASCII table), ``ASCII_HEAD`` (ASCII table with
a small header describing the column content), and ``FITS_LDAC`` (FITS binary table).
``FITS_LDAC`` is the recommended format; ``FITS_LDAC`` files are smaller, carry the data
with full precision, and can be read with popular software such as |fv| and |TOPCAT|.
They can be converted to ASCII format with the ``ldactoasc`` command-line utility,
which is part of the |SExtractor| package.

The columns present in the file are:

:Name: ``BADPIXEL_FILTER``\*
:Content: ``N``
:Unit: *Boolean*

  If true (``Y``), input objects with vignettes containing more than
  ``BADPIXEL_NMAX``  pixels flagged by |SExtractor| as bad or from deblended
  neighbours will be rejected. 

-----------------

.. _fullcat:
Full Catalogues
---------------

Other catalogues which can be produced by |SCAMP| are the “merged” catalogues (one per field
group). These catalogues contain the raw and calibrated coordinates and magnitudes of all
individual detections that passed the |SCAMP| acceptation criteria. Each detection is linked to
a parent source through to the ``SOURCE_NUMBER`` identifier, also present in the merged output
catalogue.
The ``FULLOUTCAT_TYPE`` configuration parameter sets the format of the full catalogue; the choice
of options is the same as for ``MERGEDOUTCAT_TYPE``.

The columns present in the file are:

:Name: ``SOURCE_NUMBER``\*
:Content: Source ID
:Unit: -
-----------------
:Name: ``CATALOG_NUMBER``\*
:Content: Catalogue ID
:Unit: -
-----------------
:Name: ``EXTENSION``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``ASTR_INSTRUM``\*
:Content: Astrometric instrument (context) ID
:Unit: -
-----------------
:Name: ``PHOT_INSTRUM``\*
:Content: Photometric instrument (context) ID
:Unit: -
-----------------
:Name: ``X_IMAGE``\*
:Content: ``x`` pixel coordinate of centroid
:Unit: -
-----------------
:Name: ``Y_IMAGE``\*
:Content: ``y`` pixel coordinate of centroid
:Unit: -
-----------------
:Name: ``ERRA_IMAGE``\*
:Content: Position uncertainty (RMS) along major error ellipse axis
:Unit: -
-----------------
:Name: ``ERRB_IMAGE``\*
:Content: Position uncertainty (RMS) along minor error ellipse axis
:Unit: -
-----------------
:Name: ``ERRTHETA_IMAGE``\*
:Content: Position angle of error ellipse (forced to 0 by current versions of
SCAMP which ``isotropise`` input uncertainties.)
:Unit: -
-----------------
:Name: ``ALPHA_J2000``\*
:Content: Calibrated Right-Ascension of centroid in the ICRS (at epoch of observation)
:Unit: -
-----------------
:Name: ``DELTA_J2000``\*
:Content: Calibrated Declination of centroid in the ICRS (at epoch of observation)
:Unit: -
-----------------
:Name: ``ERRA_WORLD``\*
:Content: Position uncertainty (RMS) along major world axis
 (may include additional uncertainty computed by SCAMP)
:Unit: -
-----------------
:Name: ``ERRB_WORLD``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``ERRTHETA_WORLD``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``EPOCH``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``MAG``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``MAGERR``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``MAG_DISP``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``SPREAD_MODEL``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``SPREADERR_MODEL``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``FLAGS_EXTRACTION``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------
:Name: ``FLAGS_SCAMP``\*
:Content: FITS extension of the parent image
:Unit: -
-----------------

