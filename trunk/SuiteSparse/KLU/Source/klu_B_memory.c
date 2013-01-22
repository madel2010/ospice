/* ========================================================================== */
/* === KLU_memory =========================================================== */
/* ========================================================================== */

/* KLU memory management routines:
 *
 * KLU_malloc                   malloc wrapper
 * KLU_free                     free wrapper
 * KLU_realloc                  realloc wrapper
 */

#include "klu_internal.h"

/* ========================================================================== */
/* === KLU_add_size_t ======================================================= */
/* ========================================================================== */

/* Safely compute a+b, and check for size_t overflow */



/* ========================================================================== */
/* === KLU_malloc =========================================================== */
/* ========================================================================== */

/* Wrapper around malloc routine (mxMalloc for a mexFunction).  Allocates
 * space of size MAX(1,n)*size, where size is normally a sizeof (...).
 *
 * This routine and KLU_realloc do not set Common->status to KLU_OK on success,
 * so that a sequence of KLU_malloc's or KLU_realloc's can be used.  If any of
 * them fails, the Common->status will hold the most recent error status.
 *
 * Usage, for a pointer to Int:
 *
 *      p = KLU_malloc (n, sizeof (Int), Common)
 *
 * Uses a pointer to the malloc routine (or its equivalent) defined in Common.
 */


void *KLU_new        /* returns pointer to the newly malloc'd block */
(
    /* ---- input ---- */
    size_t n,           /* number of items */
    size_t size,        /* size of each item */
    /* --------------- */
    KLU_common *Common
)
{
    void *p ;
    size_t s ;
    Int ok = TRUE ;

    if (Common == NULL)
    {
        p = NULL ;
    }
    else if (size == 0)
    {
        /* size must be > 0 */
        Common->status = KLU_INVALID ;
        p = NULL ;
    }
    else if (n >= Int_MAX)
    {
        /* object is too big to allocate; p[i] where i is an Int will not
         * be enough. */
        Common->status = KLU_TOO_LARGE ;
        p = NULL ;
    }
    else
    {
        /* call malloc, or its equivalent */
        s = KLU_mult_size_t (MAX (1,n), size, &ok) ;
        //p = ok ? ((Common->malloc_memory) (s)) : NULL ;

	p = (void*) new Entry[n];
	
        if (p == NULL)
        {
            /* failure: out of memory */
            Common->status = KLU_OUT_OF_MEMORY ;
        }
        else
        {
            Common->memusage += s ;
            Common->mempeak = MAX (Common->mempeak, Common->memusage) ;
        }
    }
    return (p) ;
}

    

/* ========================================================================== */
/* === KLU_free ============================================================= */
/* ========================================================================== */

/* Wrapper around free routine (mxFree for a mexFunction).  Returns NULL,
 * which can be assigned to the pointer being freed, as in:
 *
 *      p = KLU_free (p, n, sizeof (int), Common) ;
 */


void *KLU_delete          /* always returns NULL */
(
    /* ---- in/out --- */
    void *p,            /* block of memory to free */
    /* ---- input --- */
    size_t n,           /* size of block to free, in # of items */
    size_t size,        /* size of each item */
    /* --------------- */
    KLU_common *Common
)
{
    size_t s ;
    Int ok = TRUE ;
    if (p != NULL && Common != NULL)
    {
        /* only free the object if the pointer is not NULL */
        /* call free, or its equivalent */
        //(Common->free_memory) (p) ;
        delete[] (Entry*) p;
	s = KLU_mult_size_t (MAX (1,n), size, &ok) ;
        Common->memusage -= s ;
    }
    /* return NULL, and the caller should assign this to p.  This avoids
     * freeing the same pointer twice. */
    return (NULL) ;
}

