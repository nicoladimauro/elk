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

#include "sandholm.h"

void
sandholm (int nAgents)
{
  struct rusage earlier, later, interval;


  /* calcolo valore coalizione bottom */

  /** LEVEL 1 */
  cs best_cs = cs_new (nAgents);
  unsigned long long int num_coalitions = pow (2, nAgents) - 2;
  long long int nodes = 0;

  printf ("\n == Sandholm ==");
  getrusage (RUSAGE_SELF, &earlier);
  double best_val = characteristic_function[num_coalitions];
  int i;
  for (i = 0; i < nAgents; i++)
    best_cs[i] = 1;
  best_cs[nAgents] = nAgents;
  best_cs[nAgents * 2] = nAgents;
  best_cs[nAgents * 2 + 1] = 1;
  nodes++;
  getrusage (RUSAGE_SELF, &later);
  printf ("\n #Val: %f, nodes %lld, time %.6f (l: 1)", best_val, nodes,
	  (double) timeval_diff (&interval, &later, &earlier) / 1000000);
  /* calcolo valore coalizioni ottenute splittando la bottom */

  /** LEVEL 2 */
  int *ip = (int *) calloc (nAgents + 1, sizeof (int));
  int *ipm = (int *) calloc (nAgents + 1, sizeof (int));
  int found = 0;

  firstKCoalitionStructure (ip, ipm, nAgents, 2);
  nodes++;
  /* compute the index of the CS */
  long long int index = 0;
  for (i = 0; i < nAgents; i++)
    if (ip[i])
      index += powers[i];
  index--;

  double val = characteristic_function[index] +
    characteristic_function[num_coalitions - index - 1];
  long long int max_index;

  if (val > best_val)
    {
      best_val = val;
      max_index = index;
      found = 1;
      getrusage (RUSAGE_SELF, &later);
      printf ("\n #Val: %f, nodes %lld, time %.6f (l: 2)", best_val, nodes,
	      (double) timeval_diff (&interval, &later, &earlier) / 1000000);
    }

  while (nextKCoalitionStructure (ip, ipm, nAgents, 2))
    {
      nodes++;
      /* compute the index of the next CS */
      index = 0;
      for (i = 0; i < nAgents; i++)
	if (ip[i])
	  index += powers[i];
      index--;

      val =
	characteristic_function[index] +
	characteristic_function[num_coalitions - index - 1];
      if (val > best_val)
	{
	  best_val = val;
	  max_index = index;
	  found = 1;
	  getrusage (RUSAGE_SELF, &later);
	  printf ("\n #Val: %f, nodes %lld, time %.6f (l: 2)", best_val,
		  nodes, (double) timeval_diff (&interval, &later,
						&earlier) / 1000000);
	}
    }

  if (found)
    {
      long long int index_cs = max_index + 1, k = 0;
      while (index_cs > 0)
	{
	  best_cs[k] = index_cs % 2 + 1;
	  best_cs[nAgents + best_cs[k] - 1]++;
	  index_cs = index_cs / 2;
	  k++;
	}
      for (i = k; i < nAgents; i++)
	best_cs[i] = 1;
      for (i = nAgents + 2; i < nAgents * 2; i++)
	best_cs[i] = 0;
      best_cs[nAgents * 2] = nAgents;
      best_cs[nAgents * 2 + 1] = 2;
    }


  /* parto dal top */

  /** LEVEL n */
  val = 0.0;
  for (i = 0; i < nAgents; i++)
    val += characteristic_function[powers[i] - 1];
  if (val > best_val)
    {
      best_val = val;
      for (i = 0; i < nAgents; i++)
	{
	  best_cs[i] = i + 1;
	  best_cs[nAgents + i] = 1;
	}
      best_cs[nAgents * 2] = nAgents;
      best_cs[nAgents * 2 + 1] = nAgents;

      getrusage (RUSAGE_SELF, &later);
      printf ("\n #Val: %f, nodes %lld, time %.6f (l: %d)", best_val, nodes,
	      (double) timeval_diff (&interval, &later, &earlier) / 1000000,
	      nAgents);
    }

  nodes++;
  /* SCENDO: ad ogni livello ne fondo due */
  /* ad ogni livello devo generare tutte le CS tale che |CS|=k */

  int level, k;
  for (level = nAgents - 1; level > 2; level--)
    {
      firstKCoalitionStructure (ip, ipm, nAgents, level);
      nodes++;
      /* compute the index of the CS */
      long long int index;
      double val = 0.0;

      for (k = 0; k < level; k++)
	{
	  index = 0;
	  for (i = 0; i < nAgents; i++)
	    if (ip[i] == k)
	      index += powers[i];
	  index--;

	  val += characteristic_function[index];
	}

      if (val > best_val)
	{
	  best_val = val;
	  for (i = 0; i < nAgents * 2 + 2; i++)
	    best_cs[i] = 0;
	  for (i = 0; i < nAgents; i++)
	    {
	      best_cs[i] = ip[i] + 1;
	      best_cs[nAgents + ip[i]]++;
	    }
	  best_cs[nAgents * 2] = nAgents;
	  best_cs[nAgents * 2 + 1] = level;
	  getrusage (RUSAGE_SELF, &later);
	  printf ("\n #Val: %f, nodes %lld, time %.6f (l: %d)", best_val,
		  nodes, (double) timeval_diff (&interval, &later,
						&earlier) / 1000000, level);
	}
      while (nextKCoalitionStructure (ip, ipm, nAgents, level))
	{
	  nodes++;
	  /* compute the index of the next CS */
	  val = 0.0;
	  for (k = 0; k < level; k++)
	    {
	      index = 0;
	      for (i = 0; i < nAgents; i++)
		if (ip[i] == k)
		  index += powers[i];
	      index--;

	      val += characteristic_function[index];
	    }

	  if (val > best_val)
	    {
	      best_val = val;
	      for (i = 0; i < nAgents * 2 + 2; i++)
		best_cs[i] = 0;
	      for (i = 0; i < nAgents; i++)
		{
		  best_cs[i] = ip[i] + 1;
		  best_cs[nAgents + ip[i]]++;
		}
	      best_cs[nAgents * 2] = nAgents;
	      best_cs[nAgents * 2 + 1] = level;
	      getrusage (RUSAGE_SELF, &later);
	      printf ("\n #Val: %f, nodes %lld, time %.6f (l: %d)", best_val,
		      nodes, (double) timeval_diff (&interval, &later,
						    &earlier) / 1000000,
		      level);
	    }
	}
    }
  getrusage (RUSAGE_SELF, &later);
  printf ("\n #Val: %f, nodes %lld, time %.6f", best_val, nodes,
	  (double) timeval_diff (&interval, &later, &earlier) / 1000000);
  printf ("\n>Solution: ");
  int p;
  for (p = 1; p <= best_cs[nAgents * 2 + 1]; p++)
    {
      printf ("{");
      for (i = 0; i < nAgents; i++)
	if (best_cs[i] == p)
	  printf ("%d ", i + 1);
      printf ("}");
    }

  free (ip);
  free (ipm);
  cs_free (best_cs);
}
