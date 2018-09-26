/*
 *               htmlout.c
 *
 * Generate scamp html output
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *   This file part of:  SCAMP
 *
 *   Copyright:      (C) 2002-2018 Emmanuel Bertin -- IAP/CNRS/UPMC
 *
 *   License:        GNU General Public License
 *
 *   SCAMP is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   SCAMP is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   You should have received a copy of the GNU General Public License
 *   along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 *   Last modified:      13/03/2018
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <stdio.h>
#include "prefs.h"

/* DOC: see documentation in htmlout.h */
void
HtmlOut_write()
{
    char ch;

    FILE *tpl = fopen(prefs.html_tpl, "r");
    if (!tpl) {
        perror(prefs.html_tpl);
        return;
    }

    FILE *json = fopen(prefs.json_name, "r");
    if (!json) {
        perror("scamp.json");
        fclose(tpl);
        return;
    }

    FILE *out = fopen(prefs.html_name, "w");
    if (!out) {
        perror(prefs.html_name);
        fclose(tpl);
        fclose(json);
        return;
    }

    while ((ch = fgetc(tpl)) != EOF)
        fputc(ch, out);
    fclose(tpl);

    fprintf(out, "<script>\n\tvar scamp_data = $.parseJSON('");

    while ((ch = fgetc(json)) != EOF) {
        if (ch == '\'')
            fputc('\\', out);
        fputc(ch, out);
    }
    fclose(json);

    fprintf(out, "');\n</script>");
    fclose(out);

}
