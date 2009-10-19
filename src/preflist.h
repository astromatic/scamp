 /*
 				preflist.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Keywords for the configuration file.
*
*	Last modify:	19/10/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include "key.h"

#ifndef _PREFS_H_
#include "prefs.h"
#endif

#ifndef _CHECK_H_
#include "check.h"
#endif

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _MATCH_H_
#include "match.h"
#endif

#ifndef _XML_H_
#include "xml.h"
#endif

#ifdef	USE_THREADS
#define	THREADS_PREFMAX THREADS_NMAX
#else
#define	THREADS_PREFMAX 65535
#endif

int idummy;

pkeystruct key[] =
 {
  {"AHEADER_GLOBAL", P_STRING, prefs.ahead_global},
  {"AHEADER_SUFFIX", P_STRING, prefs.ahead_suffix},
  {"AIRMASS_KEY", P_STRING, prefs.airmass_key},
  {"ASTRCLIP_NSIGMA", P_FLOAT, &prefs.astrclip_nsig, 0,0, 0.0,1e31},
  {"ASTREF_BAND", P_STRING, prefs.astref_bandname},
  {"ASTREF_CATALOG", P_KEY, &prefs.astrefcat, 0,0, 0.0,0.0,
   {"NONE", "FILE", "USNO-A1", "USNO-A2", "USNO-B1", "GSC-1.3", "GSC-2.2",
    "2MASS", "DENIS-3", "UCAC-1", "UCAC-2", "UCAC-3", "SDSS-R3", "SDSS-R5",
     "SDSS-R6", "SDSS-R7", "NOMAD-1", ""}},
  {"ASTREF_WEIGHT", P_FLOAT, &prefs.astref_weight, 0,0, 1e-6,1e6},
  {"ASTREFCAT_NAME", P_STRINGLIST, prefs.astref_name, 0,0,0.0,0.0,
    {""}, 0, MAXNGROUP, &prefs.nastref_name},
  {"ASTREFCENT_KEYS", P_STRINGLIST, prefs.astrefcent_key, 0,0, 0.0,0.0,
   {""}, 2, NAXIS, &prefs.nastrefcent_key},
  {"ASTREFERR_KEYS", P_STRINGLIST, prefs.astreferr_key, 0,0, 0.0,0.0,
   {""}, 2, (NAXIS*(NAXIS+1))/2, &prefs.nastreferr_key},
  {"ASTREFMAG_KEY", P_STRING, prefs.astrefmag_key},
  {"ASTREFMAG_LIMITS", P_FLOATLIST, prefs.astref_maglim, 0,0, -99.0,99.0,
    {""}, 2, 2, &prefs.nastref_maglim},
  {"ASTRINSTRU_KEY", P_STRINGLIST, prefs.astrinstru_key, 0,0, 0.0,0.0,
   {""}, 0, 35, &prefs.nastrinstru_key},
  {"CDSCLIENT_EXEC", P_STRING, &prefs.cdsclient_path},
  {"CENTROID_KEYS", P_STRINGLIST, prefs.centroid_key, 0,0,0.0,0.0,
    {""}, 2, NAXIS, &prefs.ncentroid_key},
  {"CENTROIDERR_KEYS", P_STRINGLIST, prefs.centroiderr_key, 0,0,0.0,0.0,
    {""}, 2, (NAXIS*(NAXIS+1))/2, &prefs.ncentroiderr_key},
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
    {"NONE", "AS_PAIR", "AS_REFPAIR", "AS_XCORR","LL_PAIR", "LL_REFPAIR",
	"LL_XCORR", ""},
    0, MAXCHECK, &prefs.ncheck_type},
  {"CHECKPLOT_ANTIALIAS", P_BOOL, &prefs.cplot_antialiasflag},
  {"CHECKPLOT_CKEY", P_STRING, prefs.cplot_colourkey},
  {"CHECKPLOT_DEV", P_KEYLIST, prefs.cplot_device, 0,0, 0.0,0.0,
    {"NULL", "XWIN", "TK", "PLMETA", "PS", "PSC", "XFIG", "PNG", "JPEG", "PSTEX", "AQT", "PDF", "SVG",
	""},
    0, MAXCHECK, &prefs.ncplot_device},
  {"CHECKPLOT_NAME", P_STRINGLIST, prefs.cplot_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncplot_name},
  {"CHECKPLOT_RES", P_INTLIST, prefs.cplot_res, 0,16384,0.0,0.0,
    {""}, 1, 2, &prefs.ncplot_res},
   {"CHECKPLOT_TYPE", P_KEYLIST, prefs.cplot_type, 0,0, 0.0,0.0,
    {"NONE", "SKY_ALL", "FGROUPS", "PHOTOM", "ASTR_INTERROR1D",
     "ASTR_REFERROR1D", "ASTR_PIXERROR1D", "ASTR_SUBPIXERROR1D", "DISTORTION",
	"SHEAR_VS_AIRMASS", "PHOT_ERROR", "PHOT_ERRORVSMAG", "PHOT_ZPCORR", "ASTR_CHI2",
	"ASTR_INTERROR2D", "ASTR_REFERROR2D", "PHOT_ZPCORR3D",
	"ASTR_COLSHIFT1D" , "ASTR_REFPROPER", "ASTR_INTSYSMAP",
	"ASTR_REFSYSMAP", ""},
    0, MAXCHECK, &prefs.ncplot_type},
  {"CORRECT_COLOURSHIFTS", P_BOOL, &prefs.colourshift_flag},
  {"CROSSID_RADIUS", P_FLOAT, &prefs.crossid_radius, 0,0, 0.0,1e31},
  {"DISTORT_KEYS", P_STRINGLIST, prefs.context_name, 0,0,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_name},
  {"DISTORT_GROUPS", P_INTLIST, prefs.context_group, 1,MAXCONTEXT,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_group},
  {"DISTORT_DEGREES", P_INTLIST, prefs.group_deg, 1,16,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ngroup_deg},
  {"EXPOTIME_KEY", P_STRING, prefs.expotime_key},
  {"EXTINCT_KEY", P_STRING, prefs.extcoeff_key},
  {"FGROUP_RADIUS", P_FLOAT, &prefs.group_radius, 0,0, 0.0,180.0},
  {"FLAGS_MASK", P_INT, &prefs.flags_mask, 0,0xff, 0.0,0.0},
  {"FIXFOCALPLANE_NMIN", P_INT, &prefs.fixfocalplane_nmin, 0, 65535, 0.0, 0.0},
  {"FWHM_THRESHOLDS", P_FLOATLIST, prefs.fwhm_thresh, 0,0, 0.0,1e30,
    {""}, 2, 2, &prefs.nfwhm_thresh},
  {"HEADER_SUFFIX", P_STRING, prefs.head_suffix},
  {"HEADER_TYPE", P_KEY, &prefs.header_type, 0,0, 0.0,0.0,
   {"NORMAL","FOCAL_PLANE",""}},
  {"IMAFLAGS_MASK", P_INT, &prefs.imaflags_mask, 0,0x7fffffff, 0.0,0.0},
  {"MAGZERO_KEY", P_STRING, prefs.magzero_key},
  {"MAGZERO_OUT", P_FLOATLIST, prefs.magzero_out, 0,0, -1000.0,1000.0,
    {""}, 1, MAXPHOTINSTRU, &prefs.nmagzero_out},
  {"MAGZERO_INTERR", P_FLOATLIST, prefs.magzero_interr, 0,0, 0.0,1000.0,
    {""}, 1, MAXPHOTINSTRU, &prefs.nmagzero_interr},
  {"MAGZERO_REFERR", P_FLOATLIST, prefs.magzero_referr, 0,0, 0.0,1000.0,
    {""}, 1, MAXPHOTINSTRU, &prefs.nmagzero_referr},
  {"MATCH", P_BOOL, &prefs.match_flag},
  {"MATCH_FLIPPED", P_BOOL, &prefs.flip_flag},
  {"MATCH_NMAX", P_INT, &prefs.nmatchmax, 0,1000000, 0.0,0.0},
  {"MATCH_RESOL", P_FLOAT, &prefs.match_resol, 0,0, 0.0,1e31},
  {"MERGEDOUTCAT_NAME", P_STRING, prefs.mergedcat_name},
  {"MERGEDOUTCAT_TYPE", P_KEY, &prefs.mergedcat_type, 0,0, 0.0,0.0,
   {"NONE", "ASCII_HEAD", "ASCII", "ASCII_SKYCAT", "ASCII_VOTABLE",
	"FITS_LDAC", ""}},
  {"MOSAIC_TYPE", P_KEYLIST, prefs.mosaic_type, 0,0, 0.0,0.0,
   {"UNCHANGED", "SAME_CRVAL", "SHARE_PROJAXIS", "FIX_FOCALPLANE","LOOSE", ""},
    0, MAXASTRINSTRU, &prefs.nmosaic_type},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"PHOTCLIP_NSIGMA", P_FLOAT, &prefs.photclip_nsig, 0,0, 0.0,1e31},
  {"PHOTFLUX_KEY", P_STRING, prefs.photflux_key},
  {"PHOTFLUXERR_KEY", P_STRING, prefs.photfluxerr_key},
  {"PHOTINSTRU_KEY", P_STRINGLIST, prefs.photinstru_key, 0,0, 0.0,0.0,
   {""}, 0, 35, &prefs.nphotinstru_key},
  {"PHOTOMFLAG_KEY", P_STRING, prefs.photomflag_key},
  {"PIXSCALE_MAXERR", P_FLOAT, &prefs.pixscale_maxerr, 0,0, 1.0,2.0},
  {"POSANGLE_MAXERR", P_FLOAT, &prefs.posangle_maxerr, 0,0, 0.0,180.0},
  {"POSITION_MAXERR", P_FLOATLIST, prefs.position_maxerr, 0,0, 0.0,10800.0,
    {""}, 0, NAXIS, &prefs.nposition_maxerr},
  {"REFOUT_CATPATH", P_STRING, prefs.outref_path},
  {"REF_SERVER", P_STRINGLIST, prefs.ref_server, 0,0, 0.0,0.0,
   {""}, 0, MAX_SERVER, &prefs.nref_server},
  {"REF_PORT", P_INTLIST, prefs.ref_port, 0,65535, 0.0,0.0,
   {""}, 0, MAX_SERVER, &prefs.nref_port},
  {"SAVE_REFCATALOG", P_BOOL, &prefs.outrefcat_flag},
  {"SN_THRESHOLDS", P_FLOATLIST, prefs.sn_thresh, 0,0, 0.0,1e30,
    {""}, 2, 2, &prefs.nsn_thresh},
  {"SOLVE_ASTROM", P_BOOL, &prefs.solvastrom_flag},
  {"SOLVE_PHOTOM", P_BOOL, &prefs.solvphotom_flag},
  {"STABILITY_TYPE", P_KEYLIST, prefs.stability_type, 0,0, 0.0,0.0,
    {"EXPOSURE", "GROUP", "INSTRUMENT", "ALL",""},
    0, MAXFILE, &prefs.nstability_type},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","LOG", "FULL",""}},
  {"XML_NAME", P_STRING, prefs.xml_name},
  {"XSL_URL", P_STRING, prefs.xsl_name},
  {"WEIGHTFLAGS_MASK", P_INT, &prefs.wflags_mask, 0,0xff, 0.0,0.0},
  {"WRITE_XML", P_BOOL, &prefs.xml_flag},
  {""}
 };

char			keylist[sizeof(key)/sizeof(pkeystruct)][32];
const char		notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
" ",
"#----------------------------- Field grouping ---------------------------------",
" ",
"FGROUP_RADIUS          1.0             # Max dist (deg) between field groups",
" ",
"#---------------------------- Reference catalogs ------------------------------",
" ",
"REF_SERVER         cocat1.u-strasbg.fr # Internet addresses of catalog servers",
"*REF_PORT               1660            # Ports to connect to catalog servers",
"*CDSCLIENT_EXEC         " CDSCLIENT "         # CDSclient executable",
"ASTREF_CATALOG         USNO-B1         # NONE, FILE, USNO-A1,USNO-A2,USNO-B1,",
"                                       # GSC-1.3,GSC-2.2, UCAC-1,UCAC-2,UCAC-3",
"                                       # NOMAD-1, 2MASS, DENIS-3,",
"                                       # SDSS-R3,SDSS-R5,SDSS-R6 or SDSS-R7",
"ASTREF_BAND            DEFAULT         # Photom. band for astr.ref.magnitudes",
"                                       # or DEFAULT, BLUEST, or REDDEST",
"ASTREFCAT_NAME         astrefcat.cat   # Local astrometric reference catalogs",
"ASTREFCENT_KEYS        X_WORLD,Y_WORLD # Local ref.cat.centroid parameters",
"ASTREFERR_KEYS         ERRA_WORLD, ERRB_WORLD, ERRTHETA_WORLD",
"                                       # Local ref.cat.error ellipse parameters",
"ASTREFMAG_KEY          MAG             # Local ref.cat.magnitude parameter",
"ASTREFMAG_LIMITS       -99.0,99.0      # Select magnitude range in ASTREF_BAND",
"SAVE_REFCATALOG        N               # Save ref catalogs in FITS-LDAC format?",
"REFOUT_CATPATH         .               # Save path for reference catalogs",
" ",
"#--------------------------- Merged output catalogs ---------------------------",
" ",
"MERGEDOUTCAT_NAME      scamp.cat       # Merged output catalog filename",
"MERGEDOUTCAT_TYPE      NONE            # NONE, ASCII_HEAD, ASCII, FITS_LDAC",
" ",
"#----------------------------- Pattern matching -------------------------------",
" ",
"MATCH                  Y               # Do pattern-matching (Y/N) ?",
"MATCH_NMAX             0               # Max.number of detections for MATCHing",
"                                       # (0=auto)",
"PIXSCALE_MAXERR        1.2             # Max scale-factor uncertainty",
"POSANGLE_MAXERR        5.0             # Max position-angle uncertainty (deg)",
"POSITION_MAXERR        1.0             # Max positional uncertainty (arcmin)",
"MATCH_RESOL            0               # Matching resolution (arcsec); 0=auto",
"MATCH_FLIPPED          N               # Allow matching with flipped axes?",
"MOSAIC_TYPE            UNCHANGED       # UNCHANGED, SAME_CRVAL, SHARE_PROJAXIS,",
"                                       # FIX_FOCALPLANE or LOOSE",
"*FIXFOCALPLANE_NMIN     1               # Min number of dets for FIX_FOCALPLANE",
" ",
"#---------------------------- Cross-identification ----------------------------",
" ",
"CROSSID_RADIUS         2.0             # Cross-id initial radius (arcsec)",
" ",
"#---------------------------- Astrometric solution ----------------------------",
" ",
"SOLVE_ASTROM           Y               # Compute astrometric solution (Y/N) ?",
"ASTRINSTRU_KEY         FILTER,QRUNID   # FITS keyword(s) defining the astrom",
"STABILITY_TYPE         INSTRUMENT      # EXPOSURE, GROUP, INSTRUMENT or FILE",
"CENTROID_KEYS          XWIN_IMAGE,YWIN_IMAGE # Cat. parameters for centroiding",
"CENTROIDERR_KEYS       ERRAWIN_IMAGE,ERRBWIN_IMAGE,ERRTHETAWIN_IMAGE",
"                                       # Cat. params for centroid err ellipse",
"DISTORT_KEYS           XWIN_IMAGE,YWIN_IMAGE # Cat. parameters or FITS keywords",
"DISTORT_GROUPS         1,1             # Polynom group for each context key",
"DISTORT_DEGREES        3               # Polynom degree for each group",
"*ASTREF_WEIGHT          1.0             # Relative weight of ref.astrom.cat.",
"*ASTRCLIP_NSIGMA        3.0             # Astrom. clipping threshold in sigmas",
"*CORRECT_COLOURSHIFTS   N               # Correct for colour shifts (Y/N)?",
" ",
"#---------------------------- Photometric solution ----------------------------",
" ",
"SOLVE_PHOTOM           Y               # Compute photometric solution (Y/N) ?",
"MAGZERO_OUT            0.0             # Magnitude zero-point(s) in output",
"MAGZERO_INTERR         0.01            # Internal mag.zero-point accuracy",
"MAGZERO_REFERR         0.03            # Photom.field mag.zero-point accuracy",
"PHOTINSTRU_KEY         FILTER          # FITS keyword(s) defining the photom.",
"MAGZERO_KEY            PHOT_C          # FITS keyword for the mag zero-point",
"EXPOTIME_KEY           EXPTIME         # FITS keyword for the exposure time (s)",
"AIRMASS_KEY            AIRMASS         # FITS keyword for the airmass",
"EXTINCT_KEY            PHOT_K          # FITS keyword for the extinction coeff",
"PHOTOMFLAG_KEY         PHOTFLAG        # FITS keyword for the photometry flag",
"PHOTFLUX_KEY           FLUX_AUTO       # Catalog param. for the flux measurement",
"PHOTFLUXERR_KEY        FLUXERR_AUTO    # Catalog parameter for the flux error",
"*PHOTCLIP_NSIGMA        3.0             # Photom.clipping threshold in sigmas",
" ",
"#------------------------------- Check-plots ----------------------------------",
" ",
"*CHECKPLOT_CKEY         SCAMPCOL        # FITS keyword for PLPLOT field colour",
"CHECKPLOT_DEV          PNG             # NULL, XWIN, TK, PS, PSC, XFIG, PNG,",
"                                       # JPEG, AQT, PDF or SVG",
"*CHECKPLOT_RES          0               # Check-plot resolution (0 = default)",
"*CHECKPLOT_ANTIALIAS    Y               # Anti-aliasing using convert (Y/N) ?",
"CHECKPLOT_TYPE         " "FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,"
"ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR",
"CHECKPLOT_NAME         "
"fgroups,distort,astr_interror2d,astr_interror1d,astr_referror2d,"
"astr_referror1d,astr_chi2,psphot_error # Check-plot filename(s)",
" ",
"*#------------------------------- Check-images ---------------------------------",
"* ",
"*CHECKIMAGE_TYPE        NONE            # NONE, AS_PAIR, AS_REFPAIR, or AS_XCORR",
"*CHECKIMAGE_NAME        check.fits      # Check-image filename(s)",
"* ",
"#------------------------------ Miscellaneous ---------------------------------",
" ",
"SN_THRESHOLDS          10.0,100.0      # S/N thresholds (in sigmas) for all and",
"                                       # high-SN sample",
"FWHM_THRESHOLDS        0.0,100.0       # FWHM thresholds (in pixels) for sources",
"*FLAGS_MASK             0x00f0          # Rejection mask on SEx FLAGS",
"*WEIGHTFLAGS_MASK       0x00ff          # Rejection mask on SEx FLAGS_WEIGHT",
"*IMAFLAGS_MASK          0x0             # Rejection mask on SEx IMAFLAGS_ISO",
"*AHEADER_GLOBAL         scamp.ahead     # Filename of the global INPUT header",
"AHEADER_SUFFIX         .ahead          # Filename extension for additional",
"                                       # INPUT headers",
"HEADER_SUFFIX          .head           # Filename extension for OUTPUT headers",
"*HEADER_TYPE            NORMAL          # NORMAL or FOCAL_PLANE",
"VERBOSE_TYPE           NORMAL          # QUIET, NORMAL, LOG or FULL",
"WRITE_XML              Y               # Write XML file (Y/N)?",
"XML_NAME               scamp.xml       # Filename for XML output",
"*XSL_URL                " XSL_URL,
"*                                       # Filename for XSL style-sheet",
#ifdef USE_THREADS
"NTHREADS               0               # Number of simultaneous threads for",
"                                       # the SMP version of " BANNER,
"                                       # 0 = automatic",
#else
"NTHREADS               1               # 1 single thread",
#endif
""
 };

