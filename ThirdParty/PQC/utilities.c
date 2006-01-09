#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: utilities.c 4 2009-04-21 17:29:29Z delucia $";

#ifdef PHREEQ98
extern int AutoLoadOutputFile, CreateToC;
extern int ProcessMessages, ShowProgress, ShowProgressWindow, ShowChart;
extern int outputlinenr;
extern int stop_calculations;
void AddToCEntry (char *a, int l, int i);
void ApplicationProcessMessages (void);
/* void check_line_breaks(char *s); */
char err_str98[80];
int copy_title (char *token_ptr, char **ptr, int *length);
extern int clean_up_null (void);
#endif

static int isamong (char c, const char *s_l);

/* ---------------------------------------------------------------------- */
int
add_elt_list (struct elt_list *elt_list_ptr, LDBLE coef)
/* ---------------------------------------------------------------------- */
{
  struct elt_list *elt_list_ptr1;
  if (svnid == NULL)
    fprintf (stderr, " ");

  if (elt_list_ptr == NULL)
    return (OK);

  for (elt_list_ptr1 = elt_list_ptr; elt_list_ptr1->elt != NULL;
       elt_list_ptr1++)
  {
    if (count_elts >= max_elts)
    {
      space ((void **) ((void *) &elt_list), count_elts, &max_elts,
	     sizeof (struct elt_list));
    }
    elt_list[count_elts].elt = elt_list_ptr1->elt;
    elt_list[count_elts].coef = elt_list_ptr1->coef * coef;
    count_elts++;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_alk (struct reaction * rxn_ptr)
/* ---------------------------------------------------------------------- */
{
  int i;
  LDBLE return_value;
  struct master *master_ptr;

  return_value = 0.0;
  for (i = 1; rxn_ptr->token[i].s != NULL; i++)
  {
    master_ptr = rxn_ptr->token[i].s->secondary;
    if (master_ptr == NULL)
    {
      master_ptr = rxn_ptr->token[i].s->primary;
    }
    if (master_ptr == NULL)
    {
      sprintf (error_string, "Non-master species in secondary reaction, %s.",
	       rxn_ptr->token[0].s->name);
      error_msg (error_string, CONTINUE);
      input_error++;
      break;
    }
    return_value += rxn_ptr->token[i].coef * master_ptr->alk;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
compute_gfw (const char *string, LDBLE * gfw)
/* ---------------------------------------------------------------------- */
{
/*
 *    Input:  string contains a chemical formula
 *    Output:  gfw contains the calculated gfw
 */
  int i;
  char token[MAX_LENGTH];
  char *ptr;

  count_elts = 0;
  paren_count = 0;
  strcpy (token, string);
  ptr = token;
  if (get_elts_in_species (&ptr, 1.0) == ERROR)
  {
    return (ERROR);
  }
  *gfw = 0.0;
  for (i = 0; i < count_elts; i++)
  {
    if (elt_list[i].elt->gfw <= 0.0)
    {
      return (ERROR);
    }
    *gfw += elt_list[i].coef * (elt_list[i].elt)->gfw;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
copy_token (char *token_ptr, char **ptr, int *length)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies from **ptr to *token_ptr until first space is encountered.
 *
 *   Arguments:
 *      *token_ptr  output, place to store token
 *
 *     **ptr        input, character string to read token from
 *                  output, next position after token
 *
 *       length     output, length of token
 *
 *   Returns:
 *      UPPER,
 *      LOWER,
 *      DIGIT,
 *      EMPTY,
 *      UNKNOWN.
 */
  int i, return_value;
  char c;

/*
 *   Read to end of whitespace
 */
  while (isspace ((int) (c = **ptr)))
    (*ptr)++;
/*
 *   Check what we have
 */
  if (isupper ((int) c) || c == '[')
  {
    return_value = UPPER;
  }
  else if (islower ((int) c))
  {
    return_value = LOWER;
  }
  else if (isdigit ((int) c) || c == '.' || c == '-')
  {
    return_value = DIGIT;
  }
  else if (c == '\0')
  {
    return_value = EMPTY;
  }
  else
  {
    return_value = UNKNOWN;
  }
/*
 *   Begin copying to token
 */
  i = 0;
  while ((!isspace ((int) (c = **ptr))) &&
	 /*              c != ',' && */
	 c != ';' && c != '\0')
  {
    token_ptr[i] = c;
    (*ptr)++;
    i++;
  }
  token_ptr[i] = '\0';
  *length = i;
#ifdef PHREEQ98
  if ((return_value == DIGIT) && (strstr (token_ptr, ",") != NULL))
  {
    sprintf (error_string, "Commas are not allowed as decimal separator: %s.",
	     token_ptr);
    error_msg (error_string, CONTINUE);
  }
#endif
  return (return_value);
}

#ifdef PHREEQ98
/* ---------------------------------------------------------------------- */
int
copy_title (char *token_ptr, char **ptr, int *length)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies from **ptr to *token_ptr until first space or comma is encountered.
 *
 *   Arguments:
 *      *token_ptr  output, place to store token
 *
 *     **ptr        input, character string to read token from
 *                  output, next position after token
 *
 *       length     output, length of token
 *
 *   Returns:
 *      UPPER,
 *      LOWER,
 *      DIGIT,
 *      EMPTY,
 *      UNKNOWN.
 */
  int i, return_value;
  char c;
  int Quote = FALSE;

/*
 *   Read to end of whitespace
 */
  while (isspace ((int) (c = **ptr)) || (c == ',') || (c == '"'))
  {
    if (c == '"')
      Quote = TRUE;
    (*ptr)++;
  }
/*
 *   Check what we have
 */
  if (isupper ((int) c))
  {
    return_value = UPPER;
  }
  else if (islower ((int) c))
  {
    return_value = LOWER;
  }
  else if (isdigit ((int) c) || c == '.' || c == '-')
  {
    return_value = DIGIT;
  }
  else if (c == '\0')
  {
    return_value = EMPTY;
  }
  else
  {
    return_value = UNKNOWN;
  }
/*
 *   Begin copying to token
 */
  i = 0;
  if (Quote == TRUE)
  {
    while (((int) (c = **ptr) != '"') && c != '\0')
    {
      token_ptr[i] = c;
      (*ptr)++;
      i++;
    }
  }
  else
  {
    while ((!isspace ((int) (c = **ptr))) &&
	   c != ',' && c != ';' && c != '\0')
    {
      token_ptr[i] = c;
      (*ptr)++;
      i++;
    }
  }
  token_ptr[i] = '\0';
  *length = i;
  return (return_value);
}
#endif
/* ---------------------------------------------------------------------- */
int
dup_print (const char *ptr, int emphasis)
/* ---------------------------------------------------------------------- */
{
/*
 *   print character string to output and logfile
 *   if emphasis == TRUE the print is set off by
 *   a row of dashes before and after the character string.
 *
 */
  int l, i;
  char *dash;

  if (pr.headings == FALSE)
    return (OK);
#ifdef PHREEQ98
  if ((CreateToC == TRUE) && (AutoLoadOutputFile == TRUE))
  {
    if (strstr (ptr, "Reading") == ptr)
      AddToCEntry ((char *) ptr, 1, outputlinenr);
    else if (strstr (ptr, "Beginning") == ptr)
      AddToCEntry ((char *) ptr, 2, outputlinenr);
    else if ((strstr (ptr, "TITLE") != ptr) && (strstr (ptr, "End") != ptr))
      AddToCEntry ((char *) ptr, 3, outputlinenr);
  }
#endif
  l = (int) strlen (ptr);
  dash = (char *) PHRQ_malloc ((size_t) (l + 2) * sizeof (char));
  if (dash == NULL)
    malloc_error ();
  if (emphasis == TRUE)
  {
    for (i = 0; i < l; i++)
      dash[i] = '-';
    dash[i] = '\0';
    output_msg (OUTPUT_MESSAGE, "%s\n%s\n%s\n\n", dash, ptr, dash);
    output_msg (OUTPUT_LOG, "%s\n%s\n%s\n\n", dash, ptr, dash);
  }
  else
  {
    output_msg (OUTPUT_MESSAGE, "%s\n\n", ptr);
    output_msg (OUTPUT_LOG, "%s\n\n", ptr);
  }
  dash = (char *) free_check_null (dash);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
equal (LDBLE a, LDBLE b, LDBLE eps)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks equality between two LDBLE precision numbers
 */
  if (fabs (a - b) <= eps)
    return (TRUE);
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
void *
free_check_null (void *ptr)
/* ---------------------------------------------------------------------- */
{
  if (ptr != NULL)
  {
    free (ptr);
  }
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
get_token (char **eqnaddr, char *string, LDBLE * z, int *l)
/* ---------------------------------------------------------------------- */
/*
 *   Function finds next species in equation, coefficient has already
 *   been removed.  Also determines charge on the species.
 *
 *   Arguments:
 *      *eqnaddr   input, pointer to position in eqn to start parsing
 *                 output, pointer to a pointer to next position in eqn to start
 *                         parsing.
 *      *string    input pointer to place to store token
 *      *z         charge on token
 *      *l         length of token
 *
 *   Returns:
 *      ERROR,
 *      OK.
 */
{
  int i, j;
  int ltoken, lcharge;
  char c;
  char *ptr, *ptr1, *rest;
  char charge[MAX_LENGTH];

  rest = *eqnaddr;
  ptr = *eqnaddr;
  i = 0;
/*
 *   Find end of token or begining of charge
 */
  while (((c = *ptr) != '+') && (c != '-') && (c != '=') && (c != '\0'))
  {
    string[i++] = c;
    if (c == '[')
    {
      ptr++;
      while ((c = *ptr) != ']')
      {
	if (c == '\0')
	{
	  sprintf (error_string,
		   "No final bracket \"]\" for element name, %s.", string);
	  error_msg (error_string, CONTINUE);
	  return (ERROR);
	}
	string[i++] = c;
	if (i >= MAX_LENGTH)
	{
	  output_msg (OUTPUT_MESSAGE,
		      "Species name greater than MAX_LENGTH (%d) characters.\n%s\n",
		      MAX_LENGTH, string);
	  return (ERROR);
	}
	ptr++;
      }
      string[i++] = c;
    }

    /* check for overflow of space */
    if (i >= MAX_LENGTH)
    {
      output_msg (OUTPUT_MESSAGE,
		  "Species name greater than MAX_LENGTH (%d) characters.\n%s\n",
		  MAX_LENGTH, string);
      return (ERROR);
    }
    ptr++;
  }
  string[i] = '\0';
  ltoken = i;
/*
 *   Check for an empty string
 */
  if (i == 0)
  {
    sprintf (error_string, "NULL string detected in get_token, %s.", rest);
    error_msg (error_string, CONTINUE);
    return (ERROR);
  }
/*
 *   End of token is = or \0, charge is zero
 */
  if (c == '=' || c == '\0')
  {
    *eqnaddr = ptr;
    lcharge = 0;
    *z = 0.0;
  }
  else
  {
/*
 *   Copy characters into charge until next species or end is detected
 */
    j = 0;
    ptr1 = ptr;
    while ((isalpha ((int) (c = *ptr1)) == FALSE) &&
	   (c != '(') &&
	   (c != ')') &&
	   (c != ']') && (c != '[') && (c != '=') && (c != '\0'))
    {
      charge[j++] = c;
      /* error if no more space */
      if (j >= MAX_LENGTH)
      {
	error_msg
	  ("The charge on a species has exceeded MAX_LENGTH characters.",
	   CONTINUE);
	return (ERROR);
      }
      ptr1++;
    }
/*
 *   Go back to last + or - if not end of side,
 *   everything before the last + or - in charge is part of the charge
 */
    if ((c != '=') && (c != '\0'))
    {
      while (((c = *ptr1) != '+') && (c != '-'))
      {
	j--;
	ptr1--;
      }
    }
    charge[j] = '\0';
    lcharge = j;
    *eqnaddr = ptr1;
/*
 *   Charge has been written, now need to check if charge has legal format
 */
    if (get_charge (charge, z) == OK)
    {
      strcat (string, charge);
    }
    else
    {
      return (ERROR);
    }
  }
  *l = ltoken + lcharge;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
isamong (char c, const char *s_l)
/* ---------------------------------------------------------------------- */
/*
 *   Function checks if c is among the characters in the string s
 *
 *   Arguments:
 *      c     input, character to check
 *     *s     string of characters
 *
 *   Returns:
 *      TRUE  if c is in set,
 *      FALSE if c in not in set.
 */
{
  int i;

  for (i = 0; s_l[i] != '\0'; i++)
  {
    if (c == s_l[i])
    {
      return (TRUE);
    }
  }
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
int
islegit (const char c)
/* ---------------------------------------------------------------------- */
/*
 *   Function checks for legal characters for chemical equations
 *
 *   Argument:
 *      c     input, character to check
 *
 *   Returns:
 *      TRUE  if c is in set,
 *      FALSE if c in not in set.
 */
{
  if (isalpha ((int) c) || isdigit ((int) c) || isamong (c, "+-=().:_[]"))
  {
    return (TRUE);
  }
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
void
malloc_error (void)
/* ---------------------------------------------------------------------- */
{
  error_msg ("NULL pointer returned from malloc or realloc.", CONTINUE);
  error_msg ("Program terminating.", STOP);
  return;
}

/* ---------------------------------------------------------------------- */
int
parse_couple (char *token)
/* ---------------------------------------------------------------------- */
{
/*
 *   Parse couple puts redox couples in standard form
 *   "+" is removed and couples are rewritten in sort
 *    order.
 */
  int e1, e2, p1, p2;
  char *ptr;
  char elt1[MAX_LENGTH], elt2[MAX_LENGTH], paren1[MAX_LENGTH],
    paren2[MAX_LENGTH];

  if (strcmp_nocase_arg1 (token, "pe") == 0)
  {
    str_tolower (token);
    return (OK);
  }
  while (replace ("+", "", token) == TRUE);
  ptr = token;
  get_elt (&ptr, elt1, &e1);
  if (*ptr != '(')
  {
    sprintf (error_string, "Element name must be followed by "
	     "parentheses in redox couple, %s.", token);
    error_msg (error_string, CONTINUE);
    parse_error++;
    return (ERROR);
  }
  paren_count = 1;
  paren1[0] = '(';
  p1 = 1;
  while (ptr != '\0')
  {
    ptr++;
    if (*ptr == '/' || *ptr == '\0')
    {
      sprintf (error_string,
	       "End of line or  " "/"
	       " encountered before end of parentheses, %s.", token);
      error_msg (error_string, CONTINUE);
      return (ERROR);
    }
    paren1[p1++] = *ptr;
    if (*ptr == '(')
      paren_count++;
    if (*ptr == ')')
      paren_count--;
    if (paren_count == 0)
      break;
  }
  paren1[p1] = '\0';
  ptr++;
  if (*ptr != '/')
  {
    sprintf (error_string, " " "/" " must follow parentheses "
	     "ending first half of redox couple, %s.", token);
    error_msg (error_string, CONTINUE);
    parse_error++;
    return (ERROR);
  }
  ptr++;
  get_elt (&ptr, elt2, &e2);
  if (strcmp (elt1, elt2) != 0)
  {
    sprintf (error_string, "Redox couple must be two redox states "
	     "of the same element, %s.", token);
    error_msg (error_string, CONTINUE);
    return (ERROR);
  }
  if (*ptr != '(')
  {
    sprintf (error_string, "Element name must be followed by "
	     "parentheses in redox couple, %s.", token);
    error_msg (error_string, CONTINUE);
    parse_error++;
    return (ERROR);
  }
  paren2[0] = '(';
  paren_count = 1;
  p2 = 1;
  while (*ptr != '\0')
  {
    ptr++;
    if (*ptr == '/' || *ptr == '\0')
    {
      sprintf (error_string, "End of line or " "/" " encountered"
	       " before end of parentheses, %s.", token);
      error_msg (error_string, CONTINUE);
      return (ERROR);
    }

    paren2[p2++] = *ptr;
    if (*ptr == '(')
      paren_count++;
    if (*ptr == ')')
      paren_count--;
    if (paren_count == 0)
      break;
  }
  paren2[p2] = '\0';
  if (strcmp (paren1, paren2) < 0)
  {
    strcpy (token, elt1);
    strcat (token, paren1);
    strcat (token, "/");
    strcat (token, elt2);
    strcat (token, paren2);
  }
  else if (strcmp (paren1, paren2) > 0)
  {
    strcpy (token, elt2);
    strcat (token, paren2);
    strcat (token, "/");
    strcat (token, elt1);
    strcat (token, paren1);
  }
  else
  {
    sprintf (error_string, "Both parts of redox couple are the same, %s.",
	     token);
    error_msg (error_string, CONTINUE);
    return (ERROR);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_centered (const char *string)
/* ---------------------------------------------------------------------- */
{
  int i, l, l1, l2;
  char token[MAX_LENGTH];

#ifdef PHREEQ98
  if ((CreateToC == TRUE) && (AutoLoadOutputFile == TRUE))
    AddToCEntry ((char *) string, 4, outputlinenr);
#endif
  l = (int) strlen (string);
  l1 = (79 - l) / 2;
  l2 = 79 - l - l1;
  for (i = 0; i < l1; i++)
    token[i] = '-';
  token[i] = '\0';
  strcat (token, string);
  for (i = 0; i < l2; i++)
    token[i + l1 + l] = '-';
  token[79] = '\0';
  output_msg (OUTPUT_MESSAGE, "%s\n\n", token);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
replace (const char *str1, const char *str2, char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function replaces str1 with str2 in str
 *
 *   Arguments:
 *      str1     search str for str1
 *      str2     replace str1 if str1 found in str
 *      str      string to be searched
 *
 *   Returns
 *      TRUE     if string was replaced
 *      FALSE    if string was not replaced
 */
  int l, l1, l2;
  char *ptr_start;

  ptr_start = strstr (str, str1);
/*
 *   Str1 not found, return
 */
  if (ptr_start == NULL)
    return (FALSE);
/*
 *   Str1 found, replace Str1 with Str2
 */
  l = (int) strlen (str);
  l1 = (int) strlen (str1);
  l2 = (int) strlen (str2);
/*
 *   Make gap in str long enough for str2
 */
#ifdef SKIP
  if (l2 < l1)
  {
    for (ptr = (ptr_start + l1); ptr < ptr_start + l; ptr++)
    {
      ptr1 = ptr + l2 - l1;
      *ptr1 = *ptr;
      if (*ptr == '\0')
	break;
    }
  }
  else
  {
    for (ptr = (str + l); ptr >= ptr_start + l1; ptr--)
    {
      ptr1 = ptr + l2 - l1;
      *ptr1 = *ptr;
    }
  }
#endif
  /* The plus one includes the terminating NULL */
  memmove (ptr_start + l2, ptr_start + l1, l - (ptr_start - str + l1) + 1);
/*
 *   Copy str2 into str
 */
#ifdef SKIP
  ptr1 = ptr_start;
  for (ptr = (char *) str2; *ptr != '\0'; ptr++)
  {
    *ptr1 = *ptr;
    ptr1++;
  }
#endif
  memcpy (ptr_start, str2, l2);
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
void
space (void **ptr, int i, int *max, int struct_size)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine has 4 functions, allocate space, reallocate space, test to
 *   determine whether space is available, and free space.
 *
 *   Arguments:
 *      ptr          pointer to malloced space
 *      i            value for test
 *         i = INIT,          allocates space
 *         i >= 0 && i < max, space is available, return
 *         i >= max,          reallocate space
 *         i = FREE,          free space.
 *      max          maximum value for i with available space
 *      struct size  size of structure to be allocated
 */
/*
 *   Return if space exists
 */
  if ((i >= 0) && (i + 1 < *max))
  {
    return;
  }
/*
 *   Realloc space
 */
  if (i + 1 >= *max)
  {
    if (*max > 1000)
    {
      *max += 1000;
    }
    else
    {
      *max *= 2;
    }
    if (i + 1 > *max)
      *max = i + 1;
    *ptr = PHRQ_realloc (*ptr, (size_t) (*max) * struct_size);
    if (*ptr == NULL)
      malloc_error ();
    return;
  }
/*
 *   Allocate space
 */
  if (i == INIT)
  {
/*		free(*ptr); */
    *ptr = PHRQ_malloc ((size_t) (*max) * struct_size);
    if (*ptr == NULL)
      malloc_error ();
    return;
  }
/*
 *   Free space
 */
/*
	if ( i == FREE ) {
		free(*ptr);
		return;
	}
 */
/*
 *   Error
 */
  error_msg ("Illegal argument to function space.", CONTINUE);
  error_msg ("Program terminating.", STOP);
  return;
}

/* ---------------------------------------------------------------------- */
void
squeeze_white (char *s_l)
/* ---------------------------------------------------------------------- */
/*
 *   Delete all white space from string s
 *
 *   Argument:
 *      *s_l input, character string, possibly containing white space
 *           output, character string with all white space removed
 *
 *   Return: void
 */
{
  int i, j;

  for (i = j = 0; s_l[i] != '\0'; i++)
  {
    if (!isspace ((int) s_l[i]))
      s_l[j++] = s_l[i];
  }
  s_l[j] = '\0';
}

/* ---------------------------------------------------------------------- */
void
str_tolower (char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Replaces string, str, with same string, lower case
 */
  char *ptr;
  ptr = str;
  while (*ptr != '\0')
  {
    *ptr = (char) tolower (*ptr);
    ptr++;
  }
}

/* ---------------------------------------------------------------------- */
void
str_toupper (char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Replaces string, str, with same string, lower case
 */
  char *ptr;
  ptr = str;
  while (*ptr != '\0')
  {
    *ptr = (char) toupper (*ptr);
    ptr++;
  }
}

/* ---------------------------------------------------------------------- */
int
strcmp_nocase (const char *str1, const char *str2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare two strings disregarding case
 */
  int c1, c2;
  while ((c1 = tolower (*str1++)) == (c2 = tolower (*str2++)))
  {
    if (c1 == '\0')
      return (0);
  }
  if (c1 < c2)
    return (-1);
  return (1);
}

/* ---------------------------------------------------------------------- */
int
strcmp_nocase_arg1 (const char *str1, const char *str2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare two strings disregarding case
 */
  int c1, c2;
  while ((c1 = tolower (*str1++)) == (c2 = *str2++))
  {
    if (c1 == '\0')
      return (0);
  }
  if (c1 < c2)
    return (-1);
  return (1);
}

/* ---------------------------------------------------------------------- */
char *
string_duplicate (const char *token)
/* ---------------------------------------------------------------------- */
{
  int l;
  char *str;

  if (token == NULL)
    return NULL;
  l = (int) strlen (token);
  str = (char *) PHRQ_malloc ((size_t) (l + 1) * sizeof (char));
  if (str == NULL)
    malloc_error ();
  strcpy (str, token);
  return (str);
}

/* ---------------------------------------------------------------------- */
char *
string_hsave (const char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *      Save character string str
 *
 *      Arguments:
 *         str   input string to save.
 *         l     input, length of str
 *
 *      Returns:
 *         starting address of saved string (str)
 */
  char *new_string;
  ENTRY item, *found_item;

  new_string = string_duplicate (str);
  item.key = new_string;
  item.data = new_string;
  found_item = hsearch_multi (strings_hash_table, item, ENTER);
  if (found_item->key == new_string)
  {
    count_strings++;
    return (new_string);
  }
  new_string = (char *) free_check_null (new_string);
  return (found_item->key);
}

/* ---------------------------------------------------------------------- */
LDBLE
under (LDBLE xval)
/* ---------------------------------------------------------------------- */
{
/*
 *   Exponentiate a number, x, but censor large and small numbers
 *   log values less than MIN_LM are set to 0
 *   log values greater than MAX_LM are set to 10**MAX_LM
 */
/*	if (xval < MIN_LM) { */
  if (xval < -40.)
  {
    return (0.0);
  }
/*	if (xval > MAX_LM) { */
  if (xval > 3.)
  {
    return (1.0e3);
/*		return ( pow (10.0, MAX_LM));*/
  }
/*	return (pow (10.0, xval)); */
  return (exp (xval * LOG_10));
}

/* ---------------------------------------------------------------------- */
int
backspace_screen (int spaces)
/* ---------------------------------------------------------------------- */
{
  int i;
  char token[MAX_LENGTH];
  for (i = 0; i < spaces; i++)
  {
    token[i] = '\b';
  }
  token[i] = '\0';
  output_msg (OUTPUT_SCREEN, "%s", token);
  return (OK);
}

#ifndef PHREEQCI_GUI
/* ---------------------------------------------------------------------- */
int
status (int count, const char *str)
/* ---------------------------------------------------------------------- */
{
  static int spinner;
  char sim_str[20];
  char state_str[45];
  char spin_str[2];
/*	char all[MAX_LENGTH]; */
#ifdef PHREEQ98
  if (ProcessMessages)
    ApplicationProcessMessages ();
  if (stop_calculations == TRUE)
    error_msg ("Execution canceled by user.", STOP);
#endif

  if (pr.status == FALSE || phast == TRUE)
    return (OK);
  sprintf (sim_str, "Simulation %d.", simulation);
  sprintf (state_str, " ");
  sprintf (spin_str, " ");

  if (state == INITIALIZE)
  {
    output_msg (OUTPUT_SCREEN, "\n%-80s", "Initializing...");

    status_on = TRUE;
    return (OK);
  }
/*
 *   If str is defined print it
 */
  if (str != NULL)
  {
    if (status_on == TRUE)
    {
      backspace_screen (80);
    }
    else
    {
      status_on = TRUE;
    }
#ifdef DOS
    backspace_screen (80);
    /* if (state == TRANSPORT ) backspace_screen(80); */
    output_msg (OUTPUT_SCREEN, "%-79s", str);
#else
    output_msg (OUTPUT_SCREEN, "%-80s", str);
#endif
  }
  else if (state != TRANSPORT && state != PHAST)
  {
    if (state == INITIAL_SOLUTION)
    {
      sprintf (state_str, "Initial solution %d.", use.solution_ptr->n_user);
    }
    else if (state == INITIAL_EXCHANGE)
    {
      sprintf (state_str, "Initial exchange %d.", use.exchange_ptr->n_user);
    }
    else if (state == INITIAL_SURFACE)
    {
      sprintf (state_str, "Initial surface %d.", use.surface_ptr->n_user);
    }
    else if (state == INVERSE)
    {
      sprintf (state_str, "Inverse %d. Models = %d.",
	       use.inverse_ptr->n_user, count);
    }
    else if (state == REACTION)
    {
      if (use.kinetics_in == TRUE)
      {
	sprintf (state_str, "Kinetic step %d.", reaction_step);
      }
      else
      {
	sprintf (state_str, "Reaction step %d.", reaction_step);
      }
    }
    else if (state == ADVECTION || state == TRANSPORT)
    {
      if (state == ADVECTION)
      {
	sprintf (state_str, "Advection, shift %d.", advection_step);
      }
      else if (state == TRANSPORT)
      {
	sprintf (state_str, "Transport, shift %d.", transport_step);
      }
      spinner++;
      if (spinner == 1)
      {
	spin_str[0] = '/';
      }
      else if (spinner == 2)
      {
	spin_str[0] = '-';
      }
      else
      {
	spin_str[0] = '\\';
	spinner = 0;
      }
    }
    if (status_on == TRUE)
    {
      backspace_screen (80);
    }
    else
    {
      status_on = TRUE;
    }
    if (use.kinetics_in == TRUE)
    {
#ifdef DOS
      backspace_screen (80);
      output_msg (OUTPUT_SCREEN, "%-15s%-27s%37s", sim_str, state_str, " ");
#else
      output_msg (OUTPUT_SCREEN, "%-15s%-27s%38s", sim_str, state_str, " ");
#endif

    }
    else
    {
#ifdef DOS
      backspace_screen (80);
      output_msg (OUTPUT_SCREEN, "%-15s%-27s%1s%36s", sim_str, state_str,
		  spin_str, " ");
#else
      output_msg (OUTPUT_SCREEN, "%-15s%-27s%1s%37s", sim_str, state_str,
		  spin_str, " ");
#endif
    }
  }
  return (OK);
}
#endif /*PHREEQCI_GUI */
/*
** Dynamic hashing, after CACM April 1988 pp 446-457, by Per-Ake Larson.
** Coded into C, with minor code improvements, and with hsearch(3) interface,
** by ejp@ausmelb.oz, Jul 26, 1988: 13:16;
** also, hcreate/hdestroy routines added to simulate hsearch(3).
**
** These routines simulate hsearch(3) and family, with the important
** difference that the hash table is dynamic - can grow indefinitely
** beyond its original size (as supplied to hcreate()).
**
** Performance appears to be comparable to that of hsearch(3).
** The 'source-code' options referred to in hsearch(3)'s 'man' page
** are not implemented; otherwise functionality is identical.
**
** Compilation controls:
** DEBUG controls some informative traces, mainly for debugging.
** HASH_STATISTICS causes HashAccesses and HashCollisions to be maintained;
** when combined with DEBUG, these are displayed by hdestroy().
**
** Problems & fixes to ejp@ausmelb.oz. WARNING: relies on pre-processor
** concatenation property, in probably unnecessary code 'optimisation'.
** Esmond Pitt, Austec (Asia/Pacific) Ltd
** ...!uunet.UU.NET!munnari!ausmelb!ejp,ejp@ausmelb.oz
*/

# include	<assert.h>

/*
** Fast arithmetic, relying on powers of 2,
** and on pre-processor concatenation property
*/

/* rewrote to remove MUL and DIV */
#ifdef SKIP
# define MUL(x,y)		((x) << (y/**/Shift))
# define DIV(x,y)		((x) >> (y/**/Shift))
#endif
# define MOD(x,y)		((x) & ((y)-1))

/*
** local data templates
*/


/*
** external routines
*/

/*
extern char	*calloc();
extern int	free();
 */
/*
** Entry points
*/


/*
** Internal routines
*/

static Address Hash_multi (HashTable * Table, char *Key);
static void ExpandTable_multi (HashTable * Table);

/*
** Local data
*/

#ifdef HASH_STATISTICS
static long HashAccesses, HashCollisions;
#endif

/*
** Code
*/

int
hcreate_multi (unsigned Count, HashTable ** HashTable_ptr)
{
  int i;
  HashTable *Table;
  /*
   ** Adjust Count to be nearest higher power of 2,
   ** minimum SegmentSize, then convert into segments.
   */
  i = SegmentSize;
  while (i < (int) Count)
    i <<= 1;
/*    Count = DIV(i,SegmentSize); */
  Count = ((i) >> (SegmentSizeShift));

  Table = (HashTable *) PHRQ_calloc (sizeof (HashTable), 1);
  *HashTable_ptr = Table;

  if (Table == NULL)
    return (0);
  /*
   ** resets are redundant - done by calloc(3)
   **
   Table->SegmentCount = Table->p = Table->KeyCount = 0;
   */
  /*
   ** Allocate initial 'i' segments of buckets
   */
  for (i = 0; i < (int) Count; i++)
  {
    Table->Directory[i] =
      (Segment *) PHRQ_calloc (sizeof (Segment), SegmentSize);
    if (Table->Directory[i] == NULL)
    {
      hdestroy_multi (Table);
      return (0);
    }
    Table->SegmentCount++;
  }
/*    Table->maxp = MUL(Count,SegmentSize); */
  Table->maxp = (short) ((Count) << (SegmentSizeShift));
  Table->MinLoadFactor = 1;
  Table->MaxLoadFactor = DefaultMaxLoadFactor;
#ifdef DEBUG
  output_msg (OUTPUT_STDERR,
	      "[hcreate] Table %x Count %d maxp %d SegmentCount %d\n",
	      Table, Count, Table->maxp, Table->SegmentCount);
#endif
#ifdef HASH_STATISTICS
  HashAccesses = HashCollisions = 0;
#endif
  return (1);
}

void
hdestroy_multi (HashTable * Table)
{
  int i, j;
  Segment *seg;
  Element *p, *q;

  if (Table != NULL)
  {
    for (i = 0; i < Table->SegmentCount; i++)
    {
      /* test probably unnecessary        */
      if ((seg = Table->Directory[i]) != NULL)
      {
	for (j = 0; j < SegmentSize; j++)
	{
	  p = seg[j];
	  while (p != NULL)
	  {
	    q = p->Next;
	    free ((void *) p);
	    p = q;
	  }
	}
	free (Table->Directory[i]);
      }
    }
    free (Table);
    /*      Table = NULL; */
#if defined(HASH_STATISTICS) && defined(DEBUG)
    output_msg (OUTPUT_STDERR,
		"[hdestroy] Accesses %ld Collisions %ld\n",
		HashAccesses, HashCollisions);
#endif
  }
}

ENTRY *
hsearch_multi (HashTable * Table, ENTRY item, ACTION action)
/* ACTION       FIND/ENTER	*/
{
  Address h;
  Segment *CurrentSegment;
  int SegmentIndex;
  int SegmentDir;
  Segment *p, q;

  assert (Table != NULL);	/* Kinder really than return(NULL);     */
#ifdef HASH_STATISTICS
  HashAccesses++;
#endif
  h = Hash_multi (Table, item.key);
/*    SegmentDir = DIV(h,SegmentSize); */
  SegmentDir = ((h) >> (SegmentSizeShift));
  SegmentIndex = MOD (h, SegmentSize);
  /*
   ** valid segment ensured by Hash()
   */
  CurrentSegment = Table->Directory[SegmentDir];
  assert (CurrentSegment != NULL);	/* bad failure if tripped       */
  p = &CurrentSegment[SegmentIndex];
  q = *p;
  /*
   ** Follow collision chain
   */
  while (q != NULL && strcmp (q->Key, item.key))
  {
    p = &q->Next;
    q = *p;
#ifdef HASH_STATISTICS
    HashCollisions++;
#endif
  }
  if (q != NULL			/* found        */
      || action == FIND		/* not found, search only       */
      || (q = (Element *) PHRQ_calloc (sizeof (Element), 1)) == NULL	/* not found, no room   */
    )
    return ((ENTRY *) q);
  *p = q;			/* link into chain      */
  /*
   ** Initialize new element
   */
  q->Key = item.key;
  q->Data = (char *) item.data;
  /*
   ** cleared by calloc(3)
   q->Next = NULL;
   */
  /*
   ** Table over-full?
   */
/*    if (++Table->KeyCount / MUL(Table->SegmentCount,SegmentSize) > Table->MaxLoadFactor) */
  if (++Table->KeyCount / ((Table->SegmentCount) << (SegmentSizeShift)) >
      Table->MaxLoadFactor)
    ExpandTable_multi (Table);	/* doesn't affect q     */
  return ((ENTRY *) q);
}

/*
** Internal routines
*/

static Address
Hash_multi (HashTable * Table, char *Key)
{
  Address h, address;
  register unsigned char *k = (unsigned char *) Key;

  h = 0;
  /*
   ** Convert string to integer
   */
  while (*k)
    h = h * Prime1 ^ (*k++ - ' ');
  h %= Prime2;
  address = MOD (h, Table->maxp);
  if (address < (unsigned long) Table->p)
    address = MOD (h, (Table->maxp << 1));	/* h % (2*Table->maxp)  */
  return (address);
}

void
ExpandTable_multi (HashTable * Table)
{
  Address NewAddress;
  int OldSegmentIndex, NewSegmentIndex;
  int OldSegmentDir, NewSegmentDir;
  Segment *OldSegment, *NewSegment;
  Element *Current, **Previous, **LastOfNew;

/*    if (Table->maxp + Table->p < MUL(DirectorySize,SegmentSize)) */
  if (Table->maxp + Table->p < ((DirectorySize) << (SegmentSizeShift)))
  {
    /*
     ** Locate the bucket to be split
     */
/*	OldSegmentDir = DIV(Table->p,SegmentSize); */
    OldSegmentDir = ((Table->p) >> (SegmentSizeShift));
    OldSegment = Table->Directory[OldSegmentDir];
    OldSegmentIndex = MOD (Table->p, SegmentSize);
    /*
     ** Expand address space; if necessary create a new segment
     */
    NewAddress = Table->maxp + Table->p;
/*	NewSegmentDir = DIV(NewAddress,SegmentSize); */
    NewSegmentDir = ((NewAddress) >> (SegmentSizeShift));
    NewSegmentIndex = MOD (NewAddress, SegmentSize);
    if (NewSegmentIndex == 0)
      Table->Directory[NewSegmentDir] =
	(Segment *) PHRQ_calloc (sizeof (Segment), SegmentSize);
    NewSegment = Table->Directory[NewSegmentDir];
    /*
     ** Adjust state variables
     */
    Table->p++;
    if (Table->p == Table->maxp)
    {
      Table->maxp <<= 1;	/* Table->maxp *= 2     */
      Table->p = 0;
    }
    Table->SegmentCount++;
    /*
     ** Relocate records to the new bucket
     */
    Previous = &OldSegment[OldSegmentIndex];
    Current = *Previous;
    LastOfNew = &NewSegment[NewSegmentIndex];
    *LastOfNew = NULL;
    while (Current != NULL)
    {
      if (Hash_multi (Table, Current->Key) == NewAddress)
      {
	/*
	 ** Attach it to the end of the new chain
	 */
	*LastOfNew = Current;
	/*
	 ** Remove it from old chain
	 */
	*Previous = Current->Next;
	LastOfNew = &Current->Next;
	Current = Current->Next;
	*LastOfNew = NULL;
      }
      else
      {
	/*
	 ** leave it on the old chain
	 */
	Previous = &Current->Next;
	Current = Current->Next;
      }
    }
  }
}


void
free_hash_strings (HashTable * Table)
{
  int i, j;
  Segment *seg;
  Element *p, *q;

  if (Table != NULL)
  {
    for (i = 0; i < Table->SegmentCount; i++)
    {
      /* test probably unnecessary        */
      if ((seg = Table->Directory[i]) != NULL)
      {
	for (j = 0; j < SegmentSize; j++)
	{
	  p = seg[j];
	  while (p != NULL)
	  {
	    q = p->Next;
	    p->Data = (char *) free_check_null ((void *) p->Data);
	    p = q;
	  }
	}
      }
    }
#if defined(HASH_STATISTICS) && defined(DEBUG)
    output_msg (OUTPUT_STDERR,
		"[hdestroy] Accesses %ld Collisions %ld\n",
		HashAccesses, HashCollisions);
#endif
  }
}

/* ---------------------------------------------------------------------- */
int
string_trim (char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function trims white space from left and right of string
 *
 *   Arguments:
 *      str      string to trime
 *
 *   Returns
 *      TRUE     if string was changed
 *      FALSE    if string was not changed
 *      EMPTY    if string is all whitespace
 */
  int i, l, start, end, length;
  char *ptr_start;

  l = (int) strlen (str);
  /*
   *   leading whitespace
   */
  for (i = 0; i < l; i++)
  {
    if (isspace ((int) str[i]))
      continue;
    break;
  }
  if (i == l)
    return (EMPTY);
  start = i;
  ptr_start = &(str[i]);
  /*
   *   trailing whitespace
   */
  for (i = l - 1; i >= 0; i--)
  {
    if (isspace ((int) str[i]))
      continue;
    break;
  }
  end = i;
  if (start == 0 && end == l)
    return (FALSE);
  length = end - start + 1;
  memmove ((void *) str, (void *) ptr_start, (size_t) length);
  str[length] = '\0';

  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
string_trim_right (char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function trims white space from right of string
 *
 *   Arguments:
 *      str      string to trime
 *
 *   Returns
 *      TRUE     if string was changed
 *      FALSE    if string was not changed
 *      EMPTY    if string is all whitespace
 */
  int i, l, end, length;

  l = (int) strlen (str);
  for (i = l - 1; i >= 0; i--)
  {
    if (isspace ((int) str[i]))
      continue;
    break;
  }
  end = i;
  length = end + 1;
  str[length] = '\0';
  if (end == 0)
    return (EMPTY);
  if (end == l)
    return (FALSE);
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
string_trim_left (char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function trims white space from left of string
 *
 *   Arguments:
 *      str      string to trime
 *
 *   Returns
 *      TRUE     if string was changed
 *      FALSE    if string was not changed
 *      EMPTY    if string is all whitespace
 */
  int i, l, start, end, length;
  char *ptr_start;

  l = (int) strlen (str);
  /*
   *   leading whitespace
   */
  for (i = 0; i < l; i++)
  {
    if (isspace ((int) str[i]))
      continue;
    break;
  }
  if (i == l)
    return (EMPTY);
  start = i;
  ptr_start = &(str[i]);
  end = l;
  if (start == 0 && end == l)
    return (FALSE);
  length = end - start + 1;
  memmove ((void *) str, (void *) ptr_start, (size_t) length);
  str[length] = '\0';

  return (TRUE);
}

/* ---------------------------------------------------------------------- */
char *
string_pad (char *str, int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function returns new string padded to width i
 *       or returns str if longer
 *   Arguments:
 *      str      string to pad
 *
 *   Returns
 *      new string of with i
 */
  int j, l, max;
  char *str_ptr;

  l = (int) strlen (str);
  max = l;
  if (l < i)
    max = i;
  str_ptr = (char *) PHRQ_malloc ((size_t) ((max + 1) * sizeof (char)));
  if (str_ptr == NULL)
    malloc_error ();
  strcpy (str_ptr, str);
  if (i > l)
  {
    for (j = l; j < i; j++)
    {
      str_ptr[j] = ' ';
    }
    str_ptr[i] = '\0';
  }
  return (str_ptr);
}
/* ---------------------------------------------------------------------- */
void 
zero_double (LDBLE *target, int n)
/* ---------------------------------------------------------------------- */
{
  int i;

  if (n > zeros_max)
  {
    zeros = (LDBLE *) PHRQ_realloc(zeros, (size_t) (n * sizeof(LDBLE)));
    if (zeros == NULL) malloc_error();
    for (i = zeros_max; i < n; i++)
    {
      zeros[i] = 0.0;
    }
    zeros_max = n;
  }
  memcpy((void *) target, (void *) zeros, (size_t) (n * sizeof(LDBLE)));
  return;
}
