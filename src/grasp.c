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

#include "grasp.h"

double grasp_nodes_succ;
double grasp_nodes_neigh;
extern chf characteristic_function;
extern unsigned long int *powers;

int *indexv;


static cslistptr grasp_cs_neighbours (const cs, const unsigned int);
static cslistptr grasp_cs_successors (const cs, const double,
				      const unsigned int);
static cslistptr grasp_cs_select_best (const cslistptr, int);
static cslistptr grasp_cs_best_rcl (const cslistptr, const double, int);



static cslistptr
grasp_cs_select_best1 (const cslistptr ListCS, int nAgents, double val)
{
  register cslistptr iter, maxPos;
  double max;

  int select, c = 0;

  iter = ListCS;
  while (iter != NULL)
    {
      if (iter->value > val)
	c++;
      iter = iter->succ;
    }
  
  select = (double) rand () / RAND_MAX  * c + 1;  

  //  printf("\n NN %d/%d",select,c);

  if (c)
    {
      c = 0;
      iter = ListCS;
      while (iter != NULL && c < select)
	{
	  if (iter->value > val)
	    {
	      c++;
	      maxPos = iter;
	    }
	  iter = iter->succ;
	}      
    }
  else
    {
      iter = ListCS;
      maxPos = iter;
      if (iter != NULL)
	{
	  max = iter->value;
	  iter = iter->succ;
	}
      while (iter != NULL)
	{
	  if (iter->value > max)
	    {
	      max = iter->value;
	      maxPos = iter;
	    }
	  iter = iter->succ;
	}
    }

  cs CSNew;
  cslistptr Node;
  Node = (cslistptr) malloc (sizeof (struct cslist));
  CSNew = cs_new (nAgents);
  memcpy (CSNew, maxPos->CS, (nAgents * 2 + 2) * sizeof (int));
  Node->CS = CSNew;
  Node->value = maxPos->value;
  Node->succ = NULL;
  return Node;
}

double
grasp (const int nAgents, const int maxIterations)
{
  cs currentCS;
  int level, iteration;
  cslistptr pbest, best = NULL;
  double newValue, bestValue;
  cslistptr Ne;
  double absolute = 0.0;
  double current_cs_value;
  double alpha;

  indexv = (int *) malloc (nAgents * sizeof (int));

  grasp_nodes_succ = 0.0;
  grasp_nodes_neigh = 0.0;

  for (iteration = 0; iteration < maxIterations; iteration++) 
    {
      if (__verbose)
	printf ("\n== Iteration %d ==", iteration + 1);

      alpha = (double) rand () / (RAND_MAX + 1.0);
      //printf("\n  %f",alpha);

      pbest = NULL;
      currentCS = cs_new (nAgents);
      level = 0;
      int terminated = 0;
      current_cs_value = 0.0;
      while (level < nAgents && !terminated)
	{
	  cslistptr successorsCS =
	    grasp_cs_successors (currentCS, current_cs_value, nAgents);
	  if (level == 0)
	    cs_free (currentCS);
	  if (successorsCS == NULL)
	    terminated = 1;
	  else
	    {
	      /* compute the rcl set and select the best */
	      if (best != pbest)
		cslist_free (pbest);
	      pbest = grasp_cs_best_rcl (successorsCS, alpha, nAgents);
	      cslist_free (successorsCS);
	      /*pbest = grasp_cs_select_best(successorsCS); */
	      currentCS = pbest->CS;
	      current_cs_value = pbest->value;
	      if (__verbose)
		{
		  printf ("\n* Level %d", level);
		  cs_print (nAgents, currentCS);
		  printf (" %f", pbest->value);
		}
	      level++;
	    }
	  /*      printf("\n  >I:%d v:%f n:%lld (s)",
	     iteration+1,
	     pbest->value,
	     (long long int) (grasp_nodes_succ + grasp_nodes_neigh)); */
	}

      best = pbest;
      pbest = NULL;
      bestValue = best->value;
      newValue = 0;


      do
	{
	  if (newValue > bestValue)
	    bestValue = newValue;
	  Ne = grasp_cs_neighbours (currentCS, nAgents);
	  //      cslist_free(pbest);
	  //pbest = grasp_cs_select_best1 (Ne, nAgents, bestValue);
	  pbest = grasp_cs_select_best (Ne, nAgents);
	  cslist_free (Ne);
	  if (__verbose)
	    {
	      printf ("\n-----------------------");
	      cs_print (nAgents, pbest->CS);
	      printf (" %f", pbest->value);
	    }
	  newValue = pbest->value;
	  currentCS = pbest->CS;
	  if (newValue > bestValue)
	    {
	      cslist_free (best);
	      best = pbest;
	    }
	  /*	        printf("\n  >I:%d v:%f n:%lld (n)",
			iteration+1,
			pbest->value,
			(long long int) (grasp_nodes_succ + grasp_nodes_neigh)); 
	  */
	}
      while (newValue > bestValue);
      /*    printf("\nBest solution: "); */
      if (__verbose)
	{
	  cs_print (nAgents, best->CS);
	  printf (" %f", best->value);
	}
      if (best->value > absolute)
	absolute = best->value;

      //      cs_free(currentCS);      
      if (best != pbest)
	{
	  cslist_free (best);
	  cslist_free (pbest);
	}
      else
	cslist_free (best);

    }

  free (indexv);
  return (absolute);
}

/** data una CS costruisce tutti i suoi neighbours */
static cslistptr
grasp_cs_neighbours (const cs CS, const unsigned int nAgents)
{
  register cslistptr listCS;
  register unsigned int c;
  register unsigned int k;
  register unsigned int j;
  register unsigned int ca;
  register unsigned int c1, c2;
  unsigned int tmp;

  listCS = cslist_new ();

  if (__verbose)
    {
      printf ("\n Neighbours of ");
      cs_print (nAgents, CS);
    }
  /* per ogni agente cambio la coalizione di appartenenza mettendolo in 
     una di quelle già esistenti */
  for (k = 1; k <= c_get_num_agents (nAgents, CS); k++)
    {
      for (c = 1; c <= c_get_num_coalitions (nAgents, CS); c++)
	{
	  if ((ca = c_get_coalition (k, CS)) != c)
	    {
	      listCS = cslist_insert_front (listCS, CS, nAgents);
	      grasp_nodes_neigh++;
	      listCS->CS[k - 1] = c;
	      listCS->CS[nAgents + ca - 1]--;
	      listCS->CS[nAgents + c - 1]++;

	      if (CS[nAgents + ca - 1] == 1)
		{
		  /* in questo caso ho eliminato una coalizione dalla struttura e devo controllare
		     se sono necessarie delle rinumerazioni

		     es:  [4 3 2 1 ][1 1 1 1 ][4][4] --> [4 3 4 1 ][1 0 1 2 ][4][3] 
		     non deve essere 
		     [4 3 4 1 ][1 0 1 2 ][4][3]  ma 
		     [3 2 3 1 ][1 1 2 0 ][4][3] 

		   */

		  listCS->CS[nAgents * 2 + 1]--;
		  for (j = 0; j < nAgents; j++)
		    if (listCS->CS[j] > CS[k - 1])
		      listCS->CS[j]--;
		  for (j = CS[k - 1]; j <= nAgents; j++)
		    listCS->CS[nAgents + j - 1] = listCS->CS[nAgents + j];
		  listCS->CS[nAgents * 2 - 1] = 0;

		  /*      cs_print(nAgents,listCS->CS);
		     printf(" %f",listCS->value); */
		}
	      listCS->value = cs_compute_value (nAgents, listCS->CS);
	    }
	}
    }
  /* Scambio due agenti in coalizioni diverse */
  for (k = 1; k <= c_get_num_agents (nAgents, CS); k++)
    {
      for (c = k + 1; c <= c_get_num_agents (nAgents, CS); c++)
	{
	  if (c_get_coalition (k, CS) != c_get_coalition (c, CS))
	    {
	      /* li posso scambiare */
	      listCS = cslist_insert_front (listCS, CS, nAgents);
	      grasp_nodes_neigh++;
	      tmp = listCS->CS[k - 1];
	      listCS->CS[k - 1] = listCS->CS[c - 1];
	      listCS->CS[c - 1] = tmp;
	      listCS->value = cs_compute_value (nAgents, listCS->CS);

	    }
	}
    }
  /* fondo due coalizioni */

  for (c1 = 1; c1 <= c_get_num_coalitions (nAgents, CS); c1++)
    {
      for (c2 = c1 + 1; c2 <= c_get_num_coalitions (nAgents, CS); c2++)
	{
	  listCS = cslist_insert_front (listCS, CS, nAgents);
	  grasp_nodes_neigh++;
	  for (k = 1; k <= c_get_num_agents (nAgents, CS); k++)
	    if (c_get_coalition (k, listCS->CS) == c2)
	      {
		listCS->CS[k - 1] = c1;
		listCS->CS[nAgents + c1 - 1]++;
	      }
	  listCS->CS[nAgents + c2 - 1] = 0;
	  listCS->CS[nAgents * 2 + 1]--;

	  /* ho tolto la coalizione c2 e provvedo a rinumerare */
	  for (j = 0; j < nAgents; j++)
	    if (listCS->CS[j] > c2)
	      listCS->CS[j]--;
	  for (j = c2; j <= nAgents; j++)
	    listCS->CS[nAgents + j - 1] = listCS->CS[nAgents + j];
	  listCS->CS[nAgents * 2 - 1] = 0;
	  listCS->value = cs_compute_value (nAgents, listCS->CS);

	}
    }

  /* fondo tre coalizioni */

  /*
  int c3;

  for (c1 = 1; c1 <= c_get_num_coalitions (nAgents, CS); c1++)
    {
      for (c2 = c1 + 1; c2 <= c_get_num_coalitions (nAgents, CS); c2++)
	{
	  for (c3 = c2 + 1; c3 <= c_get_num_coalitions (nAgents, CS); c3++)
	    {
	      listCS = cslist_insert_front (listCS, CS, nAgents);
	      grasp_nodes_neigh++;
	      for (k = 1; k <= c_get_num_agents (nAgents, CS); k++)
		if (c_get_coalition (k, listCS->CS) == c2 || c_get_coalition (k, listCS->CS) == c2)
		  {
		    listCS->CS[k - 1] = c1;
		    listCS->CS[nAgents + c1 - 1]++;
		  }
	      listCS->CS[nAgents + c2 - 1] = 0;
	      listCS->CS[nAgents * 2 + 1]--;
	      listCS->CS[nAgents + c3 - 1] = 0;
	      listCS->CS[nAgents * 2 + 1]--;
	      
	      // ho tolto la coalizione c2 e c3 e provvedo a rinumerare 
	      for (j = 0; j < nAgents; j++)
		if (listCS->CS[j] > c2)
		  listCS->CS[j]=listCS->CS[j]-1;
	      for (j = c2; j <= nAgents; j++)
		listCS->CS[nAgents + j - 1] = listCS->CS[nAgents + j];
	      for (j = c3; j <= nAgents; j++)
		listCS->CS[nAgents + j - 1] = listCS->CS[nAgents + j];
	      listCS->CS[nAgents * 2 - 1] = 0;
	      listCS->CS[nAgents * 2 - 2] = 0;
	      listCS->value = cs_compute_value (nAgents, listCS->CS);
	    }
	}
    }
  */

  /* per tutte le coalizioni con più di un agente ne tiro a turno fuori uno 
     e lo inserisco in una nuova coalizione solo se non ho raggiunto il numero 
     massimo di coalizioni */
  if (c_get_num_coalitions (nAgents, CS) < nAgents)
    for (k = 1; k <= c_get_num_agents (nAgents, CS); k++)
      {
	//printf("\n tt %d",k);
	if (CS[nAgents + CS[k - 1] - 1] > 1)
	  {
	    //	    printf("a ");
	    int trovato = 0;	    
	    /* per coalizioni di size 2 devo metterne fuori solo 1 */
	    if (CS[nAgents + CS[k - 1] - 1] == 2)
	      {
		int trovato = 0;
		for (j=0;j<k-1 && !trovato;j++)
		  if (CS[j]==CS[k-1]) 
		    trovato = 1;
	      }
	    if (!trovato)
	      {
		//		printf("b ");
		listCS = cslist_insert_front (listCS, CS, nAgents);
		grasp_nodes_neigh++;
		listCS->CS[nAgents * 2 + 1]++;
		listCS->CS[nAgents + CS[k - 1] - 1]--;
		listCS->CS[k - 1] = listCS->CS[nAgents * 2 + 1];
		listCS->CS[nAgents + listCS->CS[k - 1] - 1]++;
		listCS->value = cs_compute_value (nAgents, listCS->CS);
	      }
	    //	    cs_print(nAgents,listCS->CS);
	    //	    printf(" %f",listCS->value); 
	  }
      }

  /* random split di coalizioni con almeno 4 elementi
   Coalizioni con due o tre elementi vengono splitate dallo shift  */
  /*    if (c_get_num_coalitions (nAgents, CS) < nAgents)
    for (k = 1; k <= c_get_num_coalitions (nAgents, CS); k++)
      {
	if (CS[nAgents + k - 1] > 3)
	  {
	    for (j=0;j<5;j++)
	      {
		int out = 0;
		listCS = cslist_insert_front (listCS, CS, nAgents);
		grasp_nodes_neigh++;
		listCS->CS[nAgents * 2 + 1]++;
		for (c=0; c<nAgents; c++)
		  {
		    if (CS[c] == k)
		      {
			double prob = (double) rand () / RAND_MAX;
			if (prob > 0.5 && out < (CS[nAgents + k - 1]-1))
			  {
			    listCS->CS[nAgents + CS[c] - 1]--;		
			    listCS->CS[c] = listCS->CS[nAgents * 2 + 1];
			    listCS->CS[nAgents + listCS->CS[c] - 1]++;			
			    out++;
			  }
		      }
		  }
		listCS->value = cs_compute_value (nAgents, listCS->CS);	    	    
	      }
	  }
      }
  */
  return (listCS);
}

static cslistptr
grasp_cs_best_rcl (const cslistptr ListCS, const double alpha, int nAgents)
{
  register cslistptr iter, maxPos, minPos;
  double max, min, threshold;
  int k, s, t;

  if (__verbose)
    printf ("\n\t-- RCL");

  iter = ListCS;

  maxPos = iter;
  minPos = iter;
  if (iter != NULL)
    {
      max = iter->value;
      min = iter->value;
      iter = iter->succ;
    }
  t = 1;
  while (iter != NULL)
    {
      if (iter->value > max)
	{
	  max = iter->value;
	  maxPos = iter;
	}
      else if (iter->value < min)
	{
	  min = iter->value;
	  minPos = iter;
	}
      iter = iter->succ;
      t++;
    }

  if (alpha >= 1.0)
    threshold = max;
  else
    threshold = min + alpha * (max - min);
  if (__verbose)
    printf (" min %f, max %f, t: %f", min, max, threshold);
  k = 0;
  iter = ListCS;
  while (iter != NULL)
    {
      if (iter->value >= threshold)
	k++;
      iter = iter->succ;
    }
  if (__verbose)
    printf (" k/t: %d/%d", k, t);
  s = ((double) rand () / (double) RAND_MAX) * (double) k + 1;
  if (__verbose)
    printf (" s: %d", s);
  fflush (stdout);
  k = 0;
  iter = ListCS;
  while (iter != NULL && k < s)
    {
      if (iter->value >= threshold)
	k++;
      if (k < s)
	iter = iter->succ;
    }


  cs CSNew;
  cslistptr Node;
  Node = (cslistptr) malloc (sizeof (struct cslist));
  CSNew = cs_new (nAgents);
  memcpy (CSNew, iter->CS, (nAgents * 2 + 2) * sizeof (int));
  Node->CS = CSNew;
  Node->value = iter->value;
  Node->succ = NULL;
  return Node;

}

static cslistptr
grasp_cs_select_best (const cslistptr ListCS, int nAgents)
{
  register cslistptr iter, maxPos;
  double max;

  iter = ListCS;

  maxPos = iter;
  if (iter != NULL)
    {
      //                  printf("\n --> ");cs_print (nAgents, iter->CS);
      //          printf (" %f", iter->value);

      max = iter->value;
      iter = iter->succ;
    }
  while (iter != NULL)
    {
      //          printf("\n --> ");cs_print (nAgents, iter->CS);
      //          printf (" %f", iter->value);

      if (iter->value > max)
	{
	  max = iter->value;
	  maxPos = iter;
	}
      iter = iter->succ;
    }

  cs CSNew;
  cslistptr Node;
  Node = (cslistptr) malloc (sizeof (struct cslist));
  CSNew = cs_new (nAgents);
  memcpy (CSNew, maxPos->CS, (nAgents * 2 + 2) * sizeof (int));
  Node->CS = CSNew;
  Node->value = maxPos->value;
  Node->succ = NULL;
  return Node;
}

/** nella ricerca greedy data una CS costruisce tutti i suoi successori di livello successivo */
static cslistptr
grasp_cs_successors (const cs CS, const double cs_value,
		     const unsigned int nAgents)
{
  register cslistptr listCS;

  listCS = cslist_new ();

  /* per ogni possibile agente che posso aggiungere lo aggiungo o ad una delle
     coalizioni già presenti o ad una nuova */
  unsigned int k, c, i;

  for (k = 1; k <= nAgents; k++)
    {
      if (!c_is_assigned_agent (k, CS))
	{
	  /* aggiungo l'agente ad una coalizione già esistente */
	  for (c = 1; c <= c_get_num_coalitions (nAgents, CS); c++)
	    {
	      listCS = cslist_insert_front (listCS, CS, nAgents);
	      grasp_nodes_succ++;
	      listCS->CS[nAgents + c - 1]++;
	      listCS->CS[k - 1] = c;
	      listCS->CS[nAgents * 2]++;
	      /*
	         listCS->value = cs_compute_value (nAgents, listCS->CS);
	       */
	      /* il nuovo valore corrisponde a sottrarre da cs_value il valore
	         della vecchia coalizione e ad aggiungere il valore di quella nuova 
	       */
	      unsigned long long int coal_index = 0;
	      int count = 0;
	      for (i = 1; i <= nAgents; i++)
		{
		  if (c_get_coalition (i, listCS->CS) == c)
		    {
		      coal_index += powers[i - 1];
		      count++;
		    }
		}
	      if (count == 1)
		listCS->value = characteristic_function[coal_index - 1];
	      else
		listCS->value =
		  cs_value -
		  characteristic_function[coal_index - powers[k - 1] - 1]
		  + characteristic_function[coal_index - 1];

	      /*      cs_print(nAgents,listCS->CS);
	         printf(" %f",listCS->value); */
	    }
	  /* aggiungo l'agente ad una nuova coalizione */
	  listCS = cslist_insert_front (listCS, CS, nAgents);
	  grasp_nodes_succ++;
	  listCS->CS[nAgents * 2 + 1]++;
	  listCS->CS[nAgents * 2]++;
	  listCS->CS[k - 1] = listCS->CS[nAgents * 2 + 1];
	  listCS->CS[nAgents + listCS->CS[k - 1] - 1]++;

	  /*      listCS->value = cs_compute_value (nAgents, listCS->CS); */
	  listCS->value =
	    cs_value + characteristic_function[powers[k - 1] - 1];

	  /*       cs_print(nAgents,listCS->CS);
	     printf(" %f",listCS->value); */

	}
    }

  return (listCS);
}
