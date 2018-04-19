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

#include "coalition.h"

extern chf characteristic_function;
extern int __verbose;
extern unsigned long int *powers;

/** */
cs
cs_new (const int nAgents)
{
  cs result;
  register int i;
  const int limit = nAgents * 2 + 1;

  result = (unsigned int *) calloc (nAgents * 2 + 2, sizeof (unsigned int));
  for (i = 0; i <= limit; i++)
    result[i] = (unsigned int) 0;
  return (result);
}

/** */
void
cs_free (const cs CS)
{
  free (CS);
}

/**  */
void
cs_print (const unsigned int nAgents, const cs CS)
{
  unsigned int i, k;

  printf ("\n");
  for (i = 1; i <= CS[2 * nAgents + 1]; i++)
    {
      printf ("{");
      for (k = 0; k < nAgents; k++)
	if (CS[k] == i)
	  printf (" %d", k + 1);
      printf ("}");
    }

  printf (": [");
  for (i = 0; i < nAgents; i++)
    printf ("%d ", CS[i]);
  printf ("][");
  for (i = 0; i < nAgents; i++)
    printf ("%d ", CS[nAgents + i]);
  printf ("][%d][%d]", CS[2 * nAgents], CS[2 * nAgents + 1]);
}

/** imposta la coalizione per un agente */
void
cs_add_agent (const int nAgents, const int agent, const int coalition,
	      cs * CS)
{
  (*CS)[agent - 1] = coalition;
  (*CS)[nAgents + coalition - 1]++;
  (*CS)[nAgents * 2]++;
}


/**  */
void
cs_remove_agent (const int nAgents, const int agent, const int coalition,
		 cs * CS)
{
  (*CS)[nAgents + (*CS)[agent - 1] - 1]--;
  if (CS[nAgents + (*CS)[agent - 1] - 1] == 0)
    (*CS)[nAgents * 2 + 1]--;
  (*CS)[nAgents * 2]--;
  (*CS)[agent - 1] = 0;
}



cslistptr
cslist_insert_front (cslistptr List, const cs CS, const int nAgents)
{
  cs CSNew;
  cslistptr Node;

  Node = (cslistptr) malloc (sizeof (struct cslist));
  CSNew = cs_new (nAgents);
  memcpy (CSNew, CS, (nAgents * 2 + 2) * sizeof (int));
  Node->CS = CSNew;
  /*  Node->value = cs_compute_value (nAgents, CS); */
  Node->succ = List;
  return Node;

}

/*
cslistptr
cslist_insert_back (cslistptr List, const cs CS, const int nAgents)
{
  cs CSNew;
  cslistptr Node, iter;

  Node = (cslistptr) malloc (sizeof (struct cslist));
  CSNew = cs_new (nAgents);
  memcpy (CSNew, CS, (nAgents * 2 + 2) * sizeof (int));
  Node->CS = CSNew;
  Node->value = cs_compute_value (nAgents, CS);
  Node->succ = NULL;

  iter = List;
  while (iter->succ != NULL)
    iter = iter->succ;
  iter->succ = Node;

  return List;
}
*/

/*
cslistptr
cslist_insert_back_agent (cslistptr List, const cs CS, const int nAgents,
			  const int agent, const cslistptr last)
{
  cs CSNew;
  cslistptr Node;

  Node = (cslistptr) malloc (sizeof (struct cslist));
  CSNew = cs_new (nAgents);
  memcpy (CSNew, CS, (nAgents * 2 + 2) * sizeof (int));
  CSNew[agent - 1] = CSNew[nAgents * 2 + 1];
  CSNew[nAgents + CSNew[nAgents * 2 + 1] - 1]++;
  CSNew[nAgents * 2]++;


  Node->CS = CSNew;
  Node->value = cs_compute_value (nAgents, CSNew);
  Node->succ = NULL;

  last->succ = Node;

  return List;
}
*/

/*
cslistptr
cslist_remove_front (cslistptr List)
{
  cslistptr iter;
  iter = List->succ;

  cs_free (List->CS);
  free (List);
  return (iter);
}
*/

/*
cslistptr
cslist_remove_back (cslistptr List)
{
  cslistptr iter, preciter;
  iter = List->succ;
  preciter = iter;

  while (iter->succ != NULL)
    {
      preciter = iter;
      iter = iter->succ;
    }

  cs_free (iter->CS);
  free (iter);
  preciter->succ = NULL;
  return (List);
}
*/

 /*
    cslistptr
    cslist_last (const cslistptr List)
    {
    cslistptr iter;

    iter = List;
    while (iter->succ != NULL)
    iter = iter->succ;

    return (iter);
    }
  */

void
cslist_free (cslistptr List)
{
  cslistptr iter, succ;

  iter = List;
  while (iter != NULL)
    {
      succ = iter->succ;
      cs_free (iter->CS);
      free (iter);
      iter = succ;
    }
}


/* DA FARE LO SCAMBIO

   1/23/4/56

   123/4/56
   23/14/56
   23/4/156
   ...
   1/2/3/4/56
   ...
   SCAMBIO
   2/13/4/56

*/

double
cs_compute_value (const unsigned int nAgents, const cs CS)
{
  register unsigned int c;
  register unsigned int k;
  double structureValue = 0.0;
  extern unsigned long int *powers;
  const unsigned int numCoalitions = c_get_num_coalitions (nAgents, CS);

  extern int *indexv;

  for (c = 0; c < numCoalitions; c++)
    indexv[c] = 0;

  for (k = 1; k <= nAgents; k++)
    indexv[c_get_coalition (k, CS) - 1] += powers[k - 1];
  for (c = 0; c < numCoalitions; c++)
    if (indexv[c])
      structureValue += characteristic_function[indexv[c] - 1];

  /*  for (c = 1; c <= numCoalitions; c++)
     {
     coalitionIndex = 0;
     for (k = 1; k <= nAgents; k++)
     if (c_get_coalition (k, CS) == c)
     coalitionIndex += powers[k - 1];
     if (coalitionIndex != 0)
     structureValue += characteristic_function[coalitionIndex - 1];
     }
   */
  return (structureValue);
}
