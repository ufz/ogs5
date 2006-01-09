#define EXTERNAL extern
#include "global.h"
#include "output.h"
/*#include <malloc.h>*/
#include <stdlib.h>
#include <string.h>
#include <assert.h>
static char const svnid[] =
  "$Id: phqalloc.c 4 2009-04-21 17:29:29Z delucia $";

#if defined(PHREEQCI_GUI)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

typedef struct PHRQMemHeader
{
  struct PHRQMemHeader *pNext;	/* memory allocated just after this one */
  struct PHRQMemHeader *pPrev;	/* memory allocated just prior to this one */
  size_t size;			/* memory request + sizeof(PHRQMemHeader) */
#if !defined(NDEBUG)
  char *szFileName;		/* file name */
  int nLine;			/* line number */
  int dummy;			/* alignment */
#endif
} PHRQMemHeader;

/*
 * This pointer points to the last allocated block
 * Note: s_pTail->pNext should always be NULL
 */
static PHRQMemHeader *s_pTail = NULL;


/* ---------------------------------------------------------------------- */
#if !defined(NDEBUG)
void *
PHRQ_malloc (size_t size, const char *szFileName, int nLine)
#else
void *
PHRQ_malloc (size_t size)
#endif
/* ---------------------------------------------------------------------- */
{
  PHRQMemHeader *p;
  if (svnid == NULL)
    fprintf (stderr, " ");

  assert ((s_pTail == NULL) || (s_pTail->pNext == NULL));

  p = (PHRQMemHeader *) malloc (sizeof (PHRQMemHeader) + size);

  if (p == NULL)
    return NULL;

  memset (p, 0, sizeof (PHRQMemHeader) + size);
  p->pNext = NULL;

  if ((p->pPrev = s_pTail) != NULL)
  {
    s_pTail->pNext = p;
  }

  p->size = sizeof (PHRQMemHeader) + size;
#if !defined(NDEBUG)
  p->szFileName = (char *) malloc (strlen (szFileName) + 1);
  if (p->szFileName)
    strcpy (p->szFileName, szFileName);
  p->nLine = nLine;
#endif

  s_pTail = p;
  p++;
  return ((void *) (p));
}

/* ---------------------------------------------------------------------- */
void
PHRQ_free (void *ptr)
/* ---------------------------------------------------------------------- */
{
  PHRQMemHeader *p;

  assert ((s_pTail == NULL) || (s_pTail->pNext == NULL));

  if (ptr == NULL)
    return;

  p = (PHRQMemHeader *) ptr - 1;

  if (p->pNext != NULL)
  {
    p->pNext->pPrev = p->pPrev;
  }
  else
  {
    /* Handle special case when (p == s_pTail) */
    assert (s_pTail != NULL);
    assert (p == s_pTail);
    s_pTail = p->pPrev;
  }

  if (p->pPrev)
  {
    p->pPrev->pNext = p->pNext;
  }

#if !defined(NDEBUG)
  free (p->szFileName);
#endif

  free (p);
}

/* ---------------------------------------------------------------------- */
void
PHRQ_free_all (void)
/* ---------------------------------------------------------------------- */
{
  assert ((s_pTail == NULL) || (s_pTail->pNext == NULL));

  if (s_pTail == NULL)
  {
#if !defined(NDEBUG)
    output_msg (OUTPUT_MESSAGE, "No memory leaks\n");
#endif
    return;
  }
  while (s_pTail->pPrev != NULL)
  {
    s_pTail = s_pTail->pPrev;
#if !defined(NDEBUG)
    output_msg (OUTPUT_MESSAGE, "%s(%d) %p: freed in PHRQ_free_all\n",
		s_pTail->pNext->szFileName, s_pTail->pNext->nLine,
		(void *) (s_pTail->pNext + 1));
    free (s_pTail->pNext->szFileName);
#endif
    free (s_pTail->pNext);
  }

#if !defined(NDEBUG)
  output_msg (OUTPUT_MESSAGE, "%s(%d) %p: freed in PHRQ_free_all\n",
	      s_pTail->szFileName, s_pTail->nLine, (void *) (s_pTail + 1));
  if (phast == TRUE)
  {
    output_msg (OUTPUT_STDERR, "%s(%d) %p: freed in PHRQ_free_all, %d\n",
		s_pTail->szFileName, s_pTail->nLine, (void *) (s_pTail + 1),
		phreeqc_mpi_myself);
  }
  free (s_pTail->szFileName);
#endif
  free (s_pTail);
  s_pTail = NULL;
}

/* ---------------------------------------------------------------------- */
void *
PHRQ_calloc (size_t num, size_t size
#if !defined(NDEBUG)
	     , const char *szFileName, int nLine
#endif
  )
/* ---------------------------------------------------------------------- */
{
  PHRQMemHeader *p;

  assert ((s_pTail == NULL) || (s_pTail->pNext == NULL));

  p = (PHRQMemHeader *) malloc (sizeof (PHRQMemHeader) + size * num);

  if (p == NULL)
    return NULL;

  p->pNext = NULL;

  if ((p->pPrev = s_pTail) != NULL)
  {
    s_pTail->pNext = p;
  }

#if !defined(NDEBUG)
  p->size = sizeof (PHRQMemHeader) + size * num;
  p->szFileName = (char *) malloc (strlen (szFileName) + 1);
  if (p->szFileName)
    strcpy (p->szFileName, szFileName);
  p->nLine = nLine;
#endif

  s_pTail = p;
  p++;
  return memset (p, 0, size * num);
}

/* ---------------------------------------------------------------------- */
void *
PHRQ_realloc (void *ptr, size_t size
#if !defined(NDEBUG)
	      , const char *szFileName, int nLine
#endif
  )
/* ---------------------------------------------------------------------- */
{
  PHRQMemHeader *p;
  size_t new_size;
  size_t old_size;

  if (ptr == NULL)
  {
    return PHRQ_malloc (size
#if !defined(NDEBUG)
			, szFileName, nLine
#endif
      );
  }

  assert ((s_pTail == NULL) || (s_pTail->pNext == NULL));

  p = (PHRQMemHeader *) ptr - 1;

  new_size = sizeof (PHRQMemHeader) + size;

  old_size = p->size;
  p = (PHRQMemHeader *) realloc (p, new_size);
  if (p != NULL)
  {
    p->size = new_size;
    if (new_size > old_size)
    {
      memset ((char *) p + old_size, 0, new_size - old_size);
    }
  }

  if (p == NULL)
    return NULL;

  if (p->pPrev != NULL)
  {
    p->pPrev->pNext = p;
  }

  if (p->pNext != NULL)
  {
    p->pNext->pPrev = p;
  }
  else
  {
    s_pTail = p;
  }

#if !defined(NDEBUG)
  free (p->szFileName);
  p->szFileName = (char *) malloc (strlen (szFileName) + 1);
  if (p->szFileName)
    strcpy (p->szFileName, szFileName);
  p->nLine = nLine;
#endif

  p++;
  return ((void *) (p));
}
