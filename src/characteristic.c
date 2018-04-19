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

#include "characteristic.h"

int *coalitionLength;
chf characteristic_function;
unsigned long int *powers;

/** Maskera di bit */
unsigned long int Mask[32] = {
  0x1, 0x2, 0x4, 0x8,
  0x10, 0x20, 0x40, 0x80,
  0x100, 0x200, 0x400, 0x800,
  0x1000, 0x2000, 0x4000, 0x8000,
  0x10000, 0x20000, 0x40000, 0x80000,
  0x100000, 0x200000, 0x400000, 0x800000,
  0x1000000, 0x2000000, 0x4000000, 0x8000000,
  0x10000000, 0x20000000, 0x40000000, 0x80000000
};

/** */
chf
chf_new (const int nAgents)
{
  chf result;
  double i;

  i = pow (2, nAgents) - 1;
  result = (double *) calloc (i, sizeof (double));
  return (result);
}

/** */
void
chf_free (chf CHF)
{
  free (CHF);
}

/**
   1      1000  1
   2      0100  2
   3      0010  4
   4      0001  8
   12     1100  3
   13     1010  5
   14     1001  9
   23     0110  6
   24     0101  10
   34     0011  12
   123    1110  7
   124    1101  11
   134    1011  13
   234    0111  14
   1234   1111  15

   split di 1234 = 15
   1/234 1  14
   2/134 2  13
   3/124 4  11
   4/123 8  7
   12/34 3  12
   13/24 5  10
   14/23 9  6

   split di 123 = 7
   1/23 1  6
   2/13 2  5
   3/12 4  3

   split di 134 = 13
   1/34 1 12
   3/14 4 9
   4/13 8 5
   
   split di 234 = 14
   2/34 2 12
   3/24 4 10
   4/23 8 6

*/

void
chf_init_normal_distributed (const int nAgents, chf * CHF)
{
  double i;
  register int k;

  printf ("\n[Computing Normally Distributed  Characteristic Function]");
  fflush (stdout);
  long long int savedtime = time (NULL);
  srand (savedtime);
  printf (" srand(%lld)", savedtime);
  i = pow (2, nAgents) - 1;

  double fac, r, v1, v2;
  double std, mean;

  coalitionLength = (int *) calloc (i, sizeof (int));
  for (k = 0; k < i; k++)
    coalitionLength[k] = chf_length (nAgents, k + 1);

  for (k = 0; k < i; k++)
    {
  /**
     Marsaglia polar method for N(0,1)
     Choose random points (x, y) in the square −1 < x < 1, −1 < y < 1 until
     s = x^2 + y^2 < 1
     the return
     x * sqrt((-2 * log(s)) / s)

     for N(mu, sigma) = N(0,1) * sqrt(sigma) + mu
   */

      do
	{
	  v1 = 2 * (double) rand () / (RAND_MAX + 1.0) - 1;
	  v2 = 2 * (double) rand () / (RAND_MAX + 1.0) - 1;
	  r = v1 * v1 + v2 * v2;
	}
      while (r > 1.0);
      fac = sqrt ((-2 * log (r)) / r);

      mean = coalitionLength[k];
      std = mean;

      (*CHF)[k] = v1 * fac * sqrt (std) + mean;
    }
}

void
chf_init_normal (const int nAgents, chf * CHF)
{
  double i;
  register int k;

  printf ("\n[Computing Normal(1,0.01) Characteristic Function]");
  fflush (stdout);
  long long int savedtime = time (NULL);
  srand (savedtime);
  printf (" srand(%lld)", savedtime);
  i = pow (2, nAgents) - 1;

  double fac, r, v1, v2;
  double std = 0.1, mean = 1.0;

  for (k = 0; k < i; k++)
    {
  /**
     Marsaglia polar method for N(0,1)
     Choose random points (x, y) in the square −1 < x < 1, −1 < y < 1 until
     s = x^2 + y^2 < 1
     the return
     x * sqrt((-2 * log(s)) / s)

     for N(mu, sigma) = N(0,1) * sqrt(sigma) + mu
   */

      do
	{
	  v1 = 2 * (double) rand () / (RAND_MAX + 1.0) - 1;
	  v2 = 2 * (double) rand () / (RAND_MAX + 1.0) - 1;
	  r = v1 * v1 + v2 * v2;
	}
      while (r > 1.0);
      fac = sqrt ((-2 * log (r)) / r);
      (*CHF)[k] = v1 * fac * std + mean;
    }
}

void
chf_init_normal_scaled (const int nAgents, chf * CHF)
{
  double i;
  register int k;

  printf ("\n[Computing Normal(1,0.01) Characteristic Function]");
  fflush (stdout);
  long long int savedtime = time (NULL);
  srand (savedtime);
  printf (" srand(%lld)", savedtime);
  i = pow (2, nAgents) - 1;

  double fac, r, v1, v2;
  double std = 0.1, mean = 1.0;

  coalitionLength = (int *) calloc (i, sizeof (int));
  for (k = 0; k < i; k++)
    coalitionLength[k] = chf_length (nAgents, k + 1);

  for (k = 0; k < i; k++)
    {
  /**
     Marsaglia polar method for N(0,1)
     Choose random points (x, y) in the square −1 < x < 1, −1 < y < 1 until
     s = x^2 + y^2 < 1
     the return
     x * sqrt((-2 * log(s)) / s)

     for N(mu, sigma) = N(0,1) * sqrt(sigma) + mu
   */

      do
	{
	  v1 = 2 * (double) rand () / (RAND_MAX + 1.0) - 1;
	  v2 = 2 * (double) rand () / (RAND_MAX + 1.0) - 1;
	  r = v1 * v1 + v2 * v2;
	}
      while (r > 1.0);
      fac = sqrt ((-2 * log (r)) / r);
      (*CHF)[k] = (v1 * fac * std + mean) * coalitionLength[k];
    }
}


void
chf_init_uniform (const int nAgents, chf * CHF)
{
  double i;
  register int k;

  printf ("\n[Computing Uniform(0,1) Characteristic Function]");
  fflush (stdout);
  long long int savedtime = time (NULL);
  srand (savedtime);
  printf (" srand(%lld)", savedtime);
  i = pow (2, nAgents) - 1;
  for (k = 0; k < i; k++)
    (*CHF)[k] = (double) rand () / (RAND_MAX + 1.0);
}


void
chf_init_uniform_scaled (const int nAgents, chf * CHF)
{
  double i;
  register int k;

  printf ("\n[Computing |C|*Uniform(0,1) Characteristic Function]");
  long long int savedtime = time (NULL);
  srand (savedtime);
  printf (" srand(%lld)", savedtime);

  i = pow (2, nAgents) - 1;
  coalitionLength = (int *) calloc (i, sizeof (int));
  /*  printf ("\n  [Computing coalitions' length]"); */
  fflush (stdout);
  for (k = 0; k < i; k++)
    coalitionLength[k] = chf_length (nAgents, k + 1);
  for (k = 0; k < i; k++)
    {
      (*CHF)[k] = (double) rand () / (RAND_MAX + 1.0) * coalitionLength[k];
    }
}

void
chf_print (const int nAgents, const chf * CHF, const char * fileName)
{

  unsigned long int i;
  int k;
  long int ncoal = pow (2, nAgents);

  FILE *fp;

  if ((fp = fopen(fileName, "wb")) == NULL) {
    printf("Impossibile aprire il file\n");
    exit(1);
  }

  printf ("\nPrinting characteristic function\n");

  fprintf(fp, "\n");
  fprintf(fp, " %ld\n",ncoal);
  fprintf(fp, " %d\n",nAgents);

  for (i = 0; i < ncoal; i++)
    fprintf(fp, " %.10f", (*CHF)[i]);

  fprintf(fp,"\n");
  /* for (i = 0; i < pow (2, nAgents) - 1; i++)
    {
      printf ("\n[");
      for (k = 0; k < nAgents; k++)
	printf ("%d ", (((i + 1) & Mask[k]) > 0));
      printf ("] %f", (*CHF)[i]);
      }*/
  fclose(fp);
}


/*

*/

int
chf_length (const int nAgents, const int chfi)
{
  register int k;
  int r = 0;
  for (k = 0; k < nAgents; k++)
    if ((chfi & Mask[k]) > 0)
      r++;
  return (r);
}
