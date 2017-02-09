#ifndef __INPUT_H__
#define __INPUT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mdsys_t.h"

#ifndef BLEN
  #define BLEN 200
#endif

int get_a_line(FILE *fp, char *buf);
int get_mdsys_stdin( mdsys_t *sys,
                     char restfile[BLEN],
                     char trajfile[BLEN],
                     char ergfile[BLEN],
                     int *nprint );

int get_rest_file(const char *restfile, mdsys_t *sys);

#endif

