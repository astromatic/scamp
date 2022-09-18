/*
*				main.c
*
* Command line parsing
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2021 IAP/CNRS/SorbonneU
*
*	License:		GNU General Public License
*
*	SCAMP is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	SCAMP is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		10/02/2021
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_PLPLOT
#include PLPLOT_H
#endif

#include "define.h"
#include "globals.h"
#include "prefs.h"
#include "cplot.h"

#define		SYNTAX \
"scamp catalog1 [catalog2,...][@catalog_list1 [@catalog_list2 ...]]\n" \
"\t\t[-c <config_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file: " BANNER " -d \n" \
"> to dump a default extended configuration file: " BANNER " -dd \n"

extern const char	notokstr[];

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   FILE         *fp;
   double	tdiff, fields, dets;
   char		**argkey, **argval,
                *str,*listbuf;
   int		a, l, narg, nim, ntok, opt,opt2;

#ifdef HAVE_SETLINEBUF
/* flush output buffer at each line */
  setlinebuf(stderr);
#endif

  if (argc<2)
    {
    fprintf(OUTPUT, "\n         %s  version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nWritten by %s\n", AUTHORS);
    fprintf(OUTPUT, "Copyright %s\n", COPYRIGHT);
    fprintf(OUTPUT, "\nvisit %s\n", WEBSITE);
    fprintf(OUTPUT, "\n%s\n", DISCLAIMER);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

#ifdef HAVE_PLPLOT
  if (argc>2)
    plparseopts(&argc, (char **)argv, PL_PARSE_SKIP);
#endif

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/* Default parameters */
  prefs.command_line = argv;
  prefs.ncommand_line = argc;
  prefs.nfile = 1;
  prefs.file_name[0] = "catalog";
  strcpy(prefs.prefs_name, "scamp.conf");
  narg = nim = 0;
  listbuf = (char *)NULL;
for (a=1; a<argc; a++)
  for (a=1; a<argc; a++)
    {
    if (argv[a] && *(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefs.prefs_name, argv[++a]);
            break;
          case 'd':
            dumpprefs(opt2=='d' ? 1 : 0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(EXIT_SUCCESS);
            break;
          case 'h':
            fprintf(OUTPUT, "\nSYNTAX: %s", SYNTAX);
#ifdef HAVE_PLPLOT
            fprintf(OUTPUT, "\nPLPLOT-specific options:\n");
            plparseopts(&argc, (char **)argv, PL_PARSE_SKIP);
#endif
            exit(EXIT_SUCCESS);
            break;
          default:
            error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
          }
        }
      else {
        argkey[narg] = &argv[a][1];
        if (!argv[a+1]	// Manage keywords with missing value
        	|| (*argv[a+1] == '-' && (argv[a+2] && *argv[a+2] != '-'))) {
          argval[narg++] = calloc(1,1);
        } else
          argval[narg++] = argv[++a];
        }       
      }
    else
      {
/*---- The input image filename(s) */
      for(; a<argc; a++)
        {
        if (!argv[a])
          continue;
        if (*argv[a]=='-')
          break;
        str = (*argv[a] == '@'? listbuf=list_to_str(argv[a]+1) : argv[a]);
        for (ntok=0; (str=strtok(ntok?NULL:str, notokstr)); nim++,ntok++)
          if (nim<MAXFILE)
            prefs.file_name[nim] = str;
          else
            error(EXIT_FAILURE, "*Error*: Too many input catalogues: ", str);
        }
      a--;
      }
    }
  prefs.nfile = nim;

  readprefs(prefs.prefs_name, argkey, argval, narg);
  useprefs();

  free(argkey);
  free(argval);

  makeit();

  free(listbuf);

  endprefs();

  NFPRINTF(OUTPUT, "");
  tdiff = prefs.time_diff>0.0? prefs.time_diff : 1.0;
  fields =  (double)prefs.nfile/tdiff;
  dets = (double)prefs.ndets/tdiff;
  NPRINTF(OUTPUT,
	"> All done (in %.0f s: %.*f field%s/s , %.0f detection%s/s)\n",
	prefs.time_diff, fields>0.1? 1: 2, fields, fields>1.0? "s":"",
	dets, dets>1.0? "s":"");

  exit(EXIT_SUCCESS);
  }

