#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#define PITZER_EXTERNAL extern
#include "pitzer.h"

static char const svnid[] =
  "$Id: pitzer_structures.c 4 2009-04-21 17:29:29Z delucia $";


static struct pitz_param *pitz_param_alloc (void);
static int pitz_param_init (struct pitz_param *pitz_param_ptr);
static struct pitz_param *pitz_param_duplicate (struct pitz_param *old_ptr);
static int pitz_param_copy (struct pitz_param *old_ptr,
			    struct pitz_param *new_ptr);


/* **********************************************************************
 *
 *   Routines related to structure "pitz_param"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct pitz_param *
pitz_param_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct pitz_param *pitz_param_ptr;
  pitz_param_ptr =
    (struct pitz_param *) PHRQ_malloc (sizeof (struct pitz_param));
  if (pitz_param_ptr == NULL)
    malloc_error ();
  return (pitz_param_ptr);
}

/* ---------------------------------------------------------------------- */
int
pitz_param_init (struct pitz_param *pitz_param_ptr)
/* ---------------------------------------------------------------------- */
{
  int i;
/*
 *   Frees all data associated with pitz_param structure.
 */

  if (pitz_param_ptr == NULL)
    return (ERROR);
  pitz_param_ptr->species[0] = NULL;
  pitz_param_ptr->species[1] = NULL;
  pitz_param_ptr->species[2] = NULL;
  pitz_param_ptr->ispec[0] = -1;
  pitz_param_ptr->ispec[1] = -1;
  pitz_param_ptr->ispec[2] = -1;
  pitz_param_ptr->type = TYPE_Other;
  pitz_param_ptr->p = 0.0;
  pitz_param_ptr->U.b0 = 0.0;
  for (i = 0; i < 5; i++)
  {
    pitz_param_ptr->a[i] = 0.0;
  }
  pitz_param_ptr->alpha = 0.0;
  pitz_param_ptr->thetas = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct pitz_param *
pitz_param_read (char *string, int n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read pitzer parameter info from string
 *   n is number of species (character values)
 *          
 */
  int l, i, j, k;
  char *ptr;
  char token[2 * MAX_LENGTH];
  struct pitz_param pzp, *pzp_ptr;

  if (n != 2 && n != 3)
    return (NULL);
  if (string == NULL)
    return (NULL);

  pitz_param_init (&pzp);
  ptr = string;
  if (copy_token (token, &ptr, &l) == EMPTY)
    return (NULL);
  ptr = string;
  for (i = 0; i < n; i++)
  {
    if (copy_token (token, &ptr, &l) == EMPTY)
      return (NULL);
    pzp.species[i] = string_hsave (token);
  }
  k = 0;
  for (i = 0; i < 5; i++)
  {
    if (copy_token (token, &ptr, &l) == EMPTY)
      break;
    j = sscanf (token, "%le", &pzp.a[i]);
    if (j <= 0)
      break;
    k++;
  }
  if (k <= 0)
    return (NULL);
  pzp_ptr = pitz_param_duplicate (&pzp);
  return (pzp_ptr);
}

/* ---------------------------------------------------------------------- */
struct pitz_param *
pitz_param_duplicate (struct pitz_param *old_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space and makes duplicate copy of pitz_param structure
 */
  struct pitz_param *new_ptr;

  new_ptr = pitz_param_alloc ();
  pitz_param_init (new_ptr);
/*
 *   Copy data
 */
  pitz_param_copy (old_ptr, new_ptr);
  return (new_ptr);
}

/* ---------------------------------------------------------------------- */
int
pitz_param_copy (struct pitz_param *old_ptr, struct pitz_param *new_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies pitz_param data from old_ptr to new location, new_ptr.
 *   Space for the new_ptr structure must already be malloced.
 */
/*
 *   Store data for structure pitz_param
 */
  memcpy (new_ptr, old_ptr, sizeof (struct pitz_param));
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pitz_param_search (struct pitz_param *pzp_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of pitz_params for same type and species
 *  Returns -1 if not found, index number in pitz_params if found
 */
  int i;
  if (pzp_ptr == NULL)
    return -1;
  if (pzp_ptr->type == TYPE_Other)
    return -1;
  for (i = 0; i < count_pitz_param; i++)
  {
    if (pitz_params[i]->type == pzp_ptr->type &&
	pitz_params[i]->species[0] == pzp_ptr->species[0] &&
	pitz_params[i]->species[1] == pzp_ptr->species[1] &&
	pitz_params[i]->species[2] == pzp_ptr->species[2])
    {
      break;
    }
  }
  if (i >= count_pitz_param)
  {
    return -1;
  }
  return i;
}

/* **********************************************************************
 *
 *   Routines related to structure "theta_parm"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct theta_param *
theta_param_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct theta_param *theta_param_ptr;
  theta_param_ptr =
    (struct theta_param *) PHRQ_malloc (sizeof (struct theta_param));
  if (theta_param_ptr == NULL)
    malloc_error ();
  return (theta_param_ptr);
}

/* ---------------------------------------------------------------------- */
int
theta_param_init (struct theta_param *theta_param_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees all data associated with theta_param structure.
 */

  if (theta_param_ptr == NULL)
    return (ERROR);
  theta_param_ptr->zj = 0;
  theta_param_ptr->zk = 0;
  theta_param_ptr->etheta = 0;
  theta_param_ptr->ethetap = 0;
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct theta_param *
theta_param_search (double zj, double zk)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of theta_params for same charge
 *  Returns NULL if not found, index number in theta_params if found
 */
  int i;
  for (i = 0; i < count_theta_param; i++)
  {
    if ((theta_params[i]->zj == zj && theta_params[i]->zk == zk) ||
	(theta_params[i]->zj == zk && theta_params[i]->zk == zj))
    {
      return theta_params[i];
    }
  }
  return NULL;
}
