/***************************************************************************
 *   Copyright (C) 2009 by Nicola Di Mauro                                 *
 *   ndm@di.uniba.it                                                       *
 *                                                                         *
 *   Department of Computer Science                                        *
 *   University of Bari, Italy                                             *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
  int               4 B = 32 b
  unsigned int      4 B = 32 b
  long int          4 B = 32 b
  long long int     8 B = 64 b
  double            8 B = 64 b
  long double       12 B = 96 b

  n agents: array of 2^n-1 double

Agents  Elements        B	        KB	 MB	GB
26	67108864	536870912	524288	 512	0,5
27	134217728	1073741824	1048576	 1024	1
28	268435456	2147483648	2097152	 2048	2
29	536870912	4294967296	4194304	 4096	4
30	1073741824	8589934592	8388608	 8192	8
31	2147483648	17179869184	16777216 16384	16
32	4294967296	34359738368	33554432 32768	32

Quindi è ragionevole rappresentare fino ad un massimo di 27 agenti per 2GB di RAM

*/

#ifndef _CHARACTERISTIC_H
#define _CHARACTERISTIC_H

#include <math.h>
#include <time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <stdio.h>

/* characteristic function */

extern int __verbose;

typedef double *chf;

#define random_number(_distance) ( (int) ((double) (_distance) *rand()/(RAND_MAX+1.0)))

/** Possiamo gestire al più 31 agenti */

/** */
chf chf_new (const int);

/** */
void chf_free (chf);

void chf_init_uniform (const int, chf *);
void chf_init_uniform_scaled (const int, chf *);
void chf_init_normal (const int, chf *);
void chf_init_normal_scaled (const int, chf *);

void chf_init_normal_distributed (const int, chf *);

int chf_length (const int, const int);

void chf_print (const int, const chf *, const char *);
#endif
