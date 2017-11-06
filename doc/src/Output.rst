.. File Output.rst

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

:Name: ``SOURCE_NUMBER``
:Content: Source ID, the same as in “full” catalogue
:Unit: 

-----------------
:Name: ``NPOS_TOTAL``
:Content: Total number of overlapping positions
:Unit: 

-----------------
:Name: ``NPOS_OK``
:Content: Number of overlapping positions kept for astrometry
:Unit: 

-----------------
:Name: ``ALPHA_J2000``
:Content: (Weighted-)average Right-Ascension
:Unit: deg

-----------------
:Name: ``DELTA_J2000``
:Content: (Weighted-)average Declination
:Unit: deg

-----------------
:Name: ``ERRA_WORLD``
:Content: Position uncertainty (RMS) along major world axis
 (may include additional uncertainty computed by SCAMP)
:Unit: deg

-----------------
:Name: ``ERRB_WORLD``
:Content: Position uncertainty (RMS) along minor world axis
 (may include additional uncertainty computed by SCAMP)
:Unit: deg

-----------------
:Name: ``ERRTHETA_WORLD``
:Content: Position angle of error ellipse (CCW/world-x)
(The current estimation of error ellipse parameters is still done very crudely)
:Unit: deg

-----------------
:Name: ``DISPALPHA_J2000``
:Content: Measured dispersion (RMS) of position along Right-Ascension
:Unit: deg

-----------------
:Name: ``DISPDELTA_J2000``
:Content: Measured dispersion (RMS) of position along Declination
:Unit: deg

-----------------
:Name: ``PMALPHA_J2000``
:Content: Apparent proper motion along Right-Ascension
:Unit: mas/yr

-----------------
:Name: ``PMDELTA_J2000``
:Content: Apparent proper motion along Declination
:Unit: mas/yr

-----------------
:Name: ``PMALPHAERR_J2000``
:Content: Proper motion uncertainty (RMS) along Right Ascension
:Unit: mas/yr

-----------------
:Name: ``PMDELTAERR_J2000``
:Content: Proper motion uncertainty (RMS) along Declination
:Unit: mas/yr

-----------------
:Name: ``PARALLAX_WORLD``
:Content: Trigonometric parallax
:Unit: mas

-----------------
:Name: ``PARALLAXERR_WORLD``
:Content: Trigonometric parallax uncertainty
:Unit: mas

-----------------
:Name: ``CHI2_ASTROM``
:Content: Reduced chi2 per d.o.f.  of the proper motion/parallax fit
:Unit: -

-----------------
:Name: ``EPOCH``
:Content: (Astrometrically-)weighted-average of observation dates
:Unit: yr

-----------------
:Name: ``EPOCH_MIN``
:Content: Earliest observation date
:Unit: yr

-----------------
:Name: ``EPOCH_MAX``
:Content: Latest observation date
:Unit: yr

-----------------
:Name: ``NMAG``
:Content: Number of magnitude measurements in each \photometric instrument"
:Unit: -

-----------------
:Name: ``MAG``
:Content: Vector of ( ux-weighted-)average magnitudes
:Unit: mag

-----------------
:Name: ``MAGERR``
:Content: Vector of magnitude uncertainties
:Unit: mag

-----------------
:Name: ``MAG_DISP``
:Content: Vector of measured magnitude dispersions (RMS)
:Unit: mag

-----------------
:Name: ``COLOR``
:Content: Composite colour index computed by SCAMP
:Unit: mag

-----------------
:Name: ``SPREAD_MODEL``
:Content: (Weighted-)average of ``SPREAD_MODEL``s
:Unit: -

-----------------
:Name: ``SPREADERR_MODEL``
:Content: ``SPREAD_MODEL`` uncertainty
:Unit: -

-----------------
:Name: ``FLAGS_EXTRACTION``
:Content: Arithmetic OR of SExtractor flags over overlapping detection
:Unit: -

-----------------
:Name: ``FLAGS_SCAMP``
:Content: SCAMP flags for this detection
:Unit: -

-----------------


.. _fullcat:
Full Catalogues
---------------

Other catalogues which can be produced by |SCAMP| are the “full” catalogues (one per field
group). These catalogues contain the raw and calibrated coordinates and magnitudes of all
individual detections that passed the |SCAMP| acceptation criteria. Each detection is linked to
a parent source through to the ``SOURCE_NUMBER`` identifier, also present in the merged output
catalogue. The ``CATALOG_NUMBER`` identifier link the detected source to its catalogue (as |SCAMP|
could receive more than one input catalogue), ``CATALOG_NUMBER`` set to 0 beeing reserved to the
reference catalogue.
The ``FULLOUTCAT_TYPE`` configuration parameter sets the format of the full catalogue; the choice
of options is the same as for ``MERGEDOUTCAT_TYPE``.

The columns present in the file are:

:Name: ``SOURCE_NUMBER``
:Content: Source ID, the same as in “merged” catalogue
:Unit: -
-----------------
:Name: ``CATALOG_NUMBER``
:Content: Catalogue ID, 0 beeing reserved to the reference catalogue
:Unit: -
-----------------
:Name: ``EXTENSION``
:Content: FITS extension of the parent image (always set to 1 for single extension images)
:Unit: -
-----------------
:Name: ``ASTR_INSTRUM``
:Content: Astrometric instrument (context) ID
:Unit: -
-----------------
:Name: ``PHOT_INSTRUM``
:Content: Photometric instrument (context) ID
:Unit: -
-----------------
:Name: ``X_IMAGE``
:Content: ``x`` pixel coordinate of centroid
:Unit: pixel
-----------------
:Name: ``Y_IMAGE``
:Content: ``y`` pixel coordinate of centroid
:Unit: pixel
-----------------
:Name: ``ERRA_IMAGE``
:Content: Position uncertainty (RMS) along major error ellipse axis
:Unit: pixel
-----------------
:Name: ``ERRB_IMAGE``
:Content: Position uncertainty (RMS) along minor error ellipse axis
:Unit: pixel
-----------------
:Name: ``ERRTHETA_IMAGE``
:Content: Position angle of error ellipse (forced to 0 by current versions of
SCAMP which ``isotropise`` input uncertainties.)
:Unit: deg
-----------------
:Name: ``ALPHA_J2000``
:Content: Calibrated Right-Ascension of centroid in the ICRS (at epoch of observation)
:Unit: deg
-----------------
:Name: ``DELTA_J2000``
:Content: Calibrated Declination of centroid in the ICRS (at epoch of observation)
:Unit: deg
-----------------
:Name: ``ERRA_WORLD``
:Content: Position uncertainty (RMS) along major world axis
 (may include additional uncertainty computed by SCAMP)
:Unit: deg
-----------------
:Name: ``ERRB_WORLD``
:Content: Position uncertainty (RMS) along minor world axis
 (may include additional uncertainty computed by SCAMP)
:Unit: deg
-----------------
:Name: ``ERRTHETA_WORLD``
:Content: Position angle of error ellipse (CCW/world-x)
(The current estimation of error ellipse parameters is still done very crudely)
:Unit: deg
-----------------
:Name: ``EPOCH``
:Content: Julian date at start of observation
:Unit: yr
-----------------
:Name: ``MAG``
:Content: Calibrated magnitude
:Unit: mag
-----------------
:Name: ``MAGERR``
:Content: Magnitude uncertainty (may include additional uncertainty computed by
SCAMP.)
:Unit: mag
-----------------
:Name: ``MAG_DISP``
:Content: Vector of measured magnitude dispersions (RMS)
:Unit: mag
-----------------
:Name: ``SPREAD_MODEL``
:Content: ``SPREAD_MODEL`` parameter
:Unit: -
-----------------
:Name: ``SPREADERR_MODEL``
:Content: ``SPREAD_MODEL`` uncertainty
:Unit: -
-----------------
:Name: ``FLAGS_EXTRACTION``
:Content: SExtractor flags
:Unit: -
-----------------
:Name: ``FLAGS_SCAMP``
:Content: SCAMP flags for this detection
:Unit: -
-----------------

.. include:: keys.rst
