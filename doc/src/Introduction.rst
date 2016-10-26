.. File Introduction.rst

Introduction
============

|SCAMP|_ (Software for Calibrating AstroMetry and Photometry) is a computer
program that computes astrometric projection parameters from source catalogues
derived from |FITS|_ images. The computed solution is expressed according
to the |WCS|_ standard. The main features of |SCAMP| are:

* Compatibility with |SExtractor|_ |FITS| or Multi-Extension |FITS| catalogue
  format in input,
* Generation of |WCS|-compliant and |SWarp|_-compatible FITS image headers in
  output,
* Automatic grouping of catalogues on the sky
* Selectable on-line astrometric reference catalogue
* Automatic determination of scale, position angle, flipping and coordinate
  shift using fast pattern-matching
* Various astrometric calibration modes for single detectors and detector
  arrays
* Combined astrometric solutions for multi-channel/instrument surveys
* Highly configurable astrometric distortion polynomials
* Correction for differential chromatic refraction
* Proper motion measurements
* Multi-threaded code that takes advantage of multiple processors.
* |VOTable|_-compliant |XML|_ output of meta-data.
* |XSLT|_ filter sheet provided for convenient access to metadata from a
  regular web browser.

.. include:: keys.rst

