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

/*! \file coalition.h
    \brief coalition header
*/

#ifndef _COALITION_H
#define _COALITION_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "characteristic.h"

/** struttura per la memorizzazione di una coalition structure 
    [c1,c2,...,ca][n1,n2,...na][A][C]
    dove 
     - a è il masimo numero di agenti 
     - ci=1..C indica a quale coalizione è assegnato l'agente i 
     - ni=1..a indica quanti agenti sono assegnati alla coalizione i
     - A=1..a indica quanti agenti sono assegnati alla CS
     - C=1..a indica quante coalizioni costituiscono la CS
*/

typedef unsigned int *cs;

#define c_get_coalition(_agent, _CS) (_CS[_agent-1])
#define c_get_num_coalitions(_nAgents, _CS) (_CS[_nAgents*2+1])
#define c_is_assigned_agent(_agent, _CS) (_CS[_agent-1])
#define c_get_num_agents(_nAgents, _CS) (_CS[_nAgents*2])
#define c_incr_num_coalitions(_nAgents, _CS) (_CS[_nAgents*2+1]++)
#define cslist_get_cs(_cslistptr) (_cslistptr->CS)
#define cslist_succ(_cslistptr) (_cslistptr->succ)
#define cslist_new() (NULL)


/** */
cs cs_new (const int);

/** */
void cs_free (const cs);

/**  */
void cs_print (const unsigned int, cs);

/** imposta la coalizione per un agente */
void cs_add_agent (const int, const int, const int, cs *);

/**  */
void cs_remove_agent (const int, const int, const int, cs *);


typedef struct cslist *cslistptr;
struct cslist
{
  cs CS;
  double value;
  cslistptr succ;
};


cslistptr cslist_insert_front (cslistptr, const cs, const int);
//cslistptr cslist_insert_back (cslistptr, const cs, const int);
//cslistptr cslist_insert_back_agent (cslistptr, const cs, const int, const int, const cslistptr);
//cslistptr cslist_insert_back (cslistptr, const cs, const int);
//cslistptr cslist_last (const cslistptr);

//cslistptr cslist_remove_front (cslistptr);
//cslistptr cslist_remove_back (cslistptr);
double cs_compute_value (const unsigned int, const cs);
void cslist_free (cslistptr);

#endif
