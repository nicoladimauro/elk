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

#include "dp.h"

static double dp_computesol (unsigned long int, int *, unsigned long int *, unsigned long int *, double *, int);
static int dp_next_coalition (int *, int);
static int dp_next_split_structure (int *, int);

/* dynamic programming */

/* possiamo rappresentare il tutto in una sola matrice 
   
   [ID][a1,...,an][NC][F11][F12][F2]
   
   ID: codice della coalizione da 1 a 2^n-1
   ai: agente nella coalizione (=0 o =1)
   NC: n.ro colaizioni risultanti

   C    CV  ID  L
   1     1  1   
   2     2  2
   3     4  3
   4     8  4

   12    3  5
   13    5  6
   14    9  7
   23    6  8
   24    10 9
   34    12 10 

   123   7  11
   124   11 12
   134   13 13
   234   14 14

   1234  15 15

   C    CV  ID
   1     1  1 
   12    3  5
   13    5  6
   14    9  7
   123   7  11
   124   11 12
   134   13 13
   1234  15 15
   2     2  2
   23    6  8
   24    10 9
   234   14 14
   3     4  3
   34    12 10 
   4     8  4


   CV è l'indice dove trovo il valore della coalizione nella funzione caratteristica
   ID è l'indice in questo vettore

   [a1,a2,...,an] è la coalizione C
   
   \sum_{i=1 \ldots |C|} 2^(a_i-1)
   
   (|C|-1)*4+

*/


int dp_next_split_structure (int *input, int n) {
  int i = 0;

  if (input[i])  {
		while (input[i] && i < n - 1)	{
			input[i] = 0;
			i++;
		}
		input[i] = 1;
	}
  else    {
		input[i] = 1;
	}
  if (i < n - 1)
    return (1);
  else
    return (0);
}

/*

  11000 ->
   10100
   10010
   10001
   01100
   01010
   01001
   00110
   00101
   00011

  11100 ->

   11010
   11001
   10110
   10101
   10011
   01110
   01101
   00111


  11110000 ->

  11101000
  11100100
  11100010
  11100001
  11011000
  11010100
  11010010
  11010001
  11001100
  11001010
  11001001
  11000110
  11000101
  11000011

*/

int dp_next_coalition (int *input, int n) {
  int i = n - 1, k;

  while (!input[i] && i >= 0)
    i--;

  if (i < n - 1)  {
      /* sposta in avanti l'ultimo 1 */
		input[i] = 0; 
		input[i + 1] =  1;
		return (1);
	}  else     {
		int ones = 1;
		int zeros = 0;
		i--;
		int trovato = 0, pos;
		while (!trovato && i >= 0)	{
			if (!input[i])
				zeros++;
			else	    {
	      if (zeros > 0)		{
					trovato = 1;
					pos = i;
				}
	      ones++;
	    }
			i--;
		}
		if (trovato)	{
			input[pos] = 0;
			for (k = 1; k <= ones; k++)
				input[pos + k] = 1;
			k = pos + ones + 1;
			while (k < n)	    {
	      input[k] = 0;
	      k++;
	    }
			return (1);
		}
	}
  return (0);
}


double
dp_computesol (unsigned long int coal, int *splitted, unsigned long int *f11,
	       unsigned long int *f12, double *f2, int nAgents)
{
  double val = 0.0;
  int k;

  /*printf("\nCoalition index: %d",coal); */
  if (splitted[coal])
    {
      val = dp_computesol (f11[coal], splitted, f11, f12, f2, nAgents) +
	dp_computesol (f12[coal], splitted, f11, f12, f2, nAgents);
    }
  else
    {
      printf ("{");
      for (k = 0; k < nAgents; k++)
	if (((coal + 1) & Mask[k]) > 0)
	  printf ("%d ", k + 1);
      printf ("}");
      return (f2[coal]);
    }
  return (val);

}

double
dp (const int nAgents)
{

  int k, i, cl;
  unsigned long long index;
  unsigned long long nsplits = 0;

  printf ("\n\n== DP: Dynamic Programming ==");
  fflush (stdout);

  int *coal;
  if ((coal = (int *) calloc (nAgents, sizeof (int))) == NULL)
    {
      printf ("\nImpossible to allocate %d sizeof (int)\n", nAgents);
      exit (1);
    }

  int *coalindex;
  if ((coalindex = (int *) calloc (nAgents, sizeof (int))) == NULL)
    {
      printf ("\nImpossible to allocate %d sizeof (int)\n", nAgents);
      exit (1);
    }

  int *coalSplit;
  if ((coalSplit = (int *) calloc (nAgents, sizeof (int))) == NULL)
    {
      printf ("\nImpossible to allocate %d sizeof (int)\n", nAgents);
      exit (1);
    }
  int *splitted;
  if ((splitted =
       (int *) calloc ((pow (2, nAgents) - 1), sizeof (int))) == NULL)
    {
      printf ("\nImpossible to allocate %f sizeof (int)\n",
	      pow (2, nAgents) - 1);
      exit (1);
    }
  unsigned long int *f11;
  if ((f11 =
       (unsigned long int *) calloc ((pow (2, nAgents) - 1),
				     sizeof (unsigned long int))) == NULL)
    {
      printf ("\nImpossible to allocate %f sizeof (unsigned long int)\n",
	      pow (2, nAgents) - 1);
      exit (1);
    }
  unsigned long int *f12;
  if ((f12 =
       (unsigned long int *) calloc ((pow (2, nAgents) - 1),
				     sizeof (unsigned long int))) == NULL)
    {
      printf ("\nImpossible to allocate %f sizeof (unsigned long int)\n",
	      pow (2, nAgents) - 1);
      exit (1);
    }
  double *f2;
  if ((f2 =
       (double *) calloc ((pow (2, nAgents) - 1), sizeof (double))) == NULL)
    {
      printf ("\nImpossible to allocate %f sizeof (double)\n",
	      pow (2, nAgents) - 1);
      exit (1);
    }


  printf ("\n[Computing matrix F1 and F2]");
  fflush (stdout);

  /* inizializzo i valori per le coalizioni di lunghezza 1 */
  for (k = 0; k < nAgents; k++)
    coal[k] = 0;

  for (k = 0; k < nAgents; k++)    {
		coal[k] = 1;
		if (k)
			coal[k - 1] = 0;
		/* compute coalition value */
		index = 0;
		for (i = 0; i < nAgents; i++)
			if (coal[i])
				index += powers[i];
		index--;
		
		f2[index] = characteristic_function[index];
		splitted[index] = 0;
		f11[index] = index;

		if (__verbose)	{
			printf ("\n[");
			for (i = 0; i < nAgents; i++)
				if (coal[i])
					printf (" %d ", i + 1);
			printf ("] %lld %f", index, f2[index]);
		}
	}

  double maxVal;
  int maxSplit;
  unsigned long long int maxf11, maxf12;

  /* calcolo valori per coalizioni di lunghezza cl */
  for (cl = 2; cl <= nAgents; cl++)    {
		//      printf ("\n-> CL %d", cl);fflush(stdout);
		if (__verbose)
			printf ("\n-> CL %d", cl);
		/* imposto la coalizione iniziale */
		for (k = 0; k < cl; k++)
			coal[k] = 1;
		for (k = cl; k < nAgents; k++)
			coal[k] = 0;

     /* computeBestSplitting */

   unsigned long long int index1, index2;

	 do	{
	  /*index = 0;
	     for (i = 0; i < nAgents; i++)
	     if (coal[i])
	     index += powers[i];
	     index--; */
	  k = 0;
	  index = 0;
	  for (i = 0; i < nAgents; i++)	    {
			if (coal[i])		{
				index += powers[i];
				coalindex[k] = i;
				k++;
			}
		}
	  index--;


	  maxVal = characteristic_function[index];
	  maxSplit = 0;
	  maxf11 = index;

	  for (k = 0; k < nAgents; k++)
	    coalSplit[k] = 0;

	  if (__verbose)  {
			printf ("\n[ ");
			for (i = 0; i < nAgents; i++)
				if (coal[i])
					printf ("%d ", i + 1);
			printf ("] %lld %f", index, maxVal);
		}

	  while (dp_next_split_structure (coalSplit, cl))	    {
			/* trovato uno split */
			nsplits++;
			if (__verbose)		{
				k = 0;
				printf ("\n\t[ ");
				for (i = 0; i < nAgents; i++)
					if (coal[i])
						{
							if (coalSplit[k])
								printf ("%d ", i + 1);
							k++;
						}
				printf ("][ ");
				k = 0;
		  for (i = 0; i < nAgents; i++)
		    if (coal[i])
		      {
			if (!coalSplit[k])
			  printf ("%d ", i + 1);
			k++;
		      }
		  printf ("]");
		}

	      index1 = 0;
	      for (i = 0; i < cl; i++)
		index1 += powers[coalindex[i]] * coalSplit[i];


	      /*              index1 = 0;
	         k = 0;
	         for (i = 0; i < nAgents; i++)
	         if (coal[i])
	         {
	         if (coalSplit[k])
	         {
	         index1 += powers[i];
	         }
	         k++;
	         }
	       */
	      index1--;
	      index2 = index - index1 - 1;

	      //printf("\n tt %lld %lld %lld",index1,index2,index);
	      //getchar();

	      if (__verbose)
		printf (" %lld:%f %lld:%f %f", index1, f2[index1], index2,
			f2[index2], f2[index1] + f2[index2]);

	      double val = f2[index1] + f2[index2];
	      if (val > maxVal)
		{
		  maxVal = val;
		  maxSplit = 1;
		  maxf11 = index1;
		  maxf12 = index2;
		}
	    }

	  if (__verbose)
	    {
	      printf ("\n  MAX: %lld", maxf11);
	      if (maxSplit)
		printf (" %lld %f", maxf12, maxVal);
	      else
		printf (" # %f", maxVal);
	    }

	  splitted[index] = maxSplit;
	  f11[index] = maxf11;
	  f12[index] = maxf12;
	  f2[index] = maxVal;

	}
      while (dp_next_coalition (coal, nAgents));

    }

  unsigned long int start = 0;
  for (k = 0; k < nAgents; k++)
    start += powers[k];

  double best;

  printf ("\n>Solution: ");
  best = dp_computesol (start - 1, splitted, f11, f12, f2, nAgents);
  printf (" val: %f (Splits: %lld)", best, nsplits);

  /* trovare la soluzione */

  free (coal);
  free (coalindex);
  free (coalSplit);
  free (splitted);
  free (f11);
  free (f12);
  free (f2);
  return (best);
}
