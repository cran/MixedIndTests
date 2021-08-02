#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
  extern void estdep(void *, void *, void *, void *, void *, void *);
extern void estdep_serial(void *, void *, void *, void *, void *, void *);
extern void Sn_serial0(void *, void *, void *, void *);
extern void stats_nonserial(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void stats_serial(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void stats_serial_bin(void *, void *, void *, void *, void *, void *, void *, void *);
extern void statsim(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"estdep",           (DL_FUNC) &estdep,            6},
  {"estdep_serial",    (DL_FUNC) &estdep_serial,     6},
  {"Sn_serial0",       (DL_FUNC) &Sn_serial0,        4},
  {"stats_nonserial",  (DL_FUNC) &stats_nonserial,  10},
  {"stats_serial",     (DL_FUNC) &stats_serial,     10},
  {"stats_serial_bin", (DL_FUNC) &stats_serial_bin,  8},
  {"statsim",          (DL_FUNC) &statsim,           7},
  {NULL, NULL, 0}
};

void R_init_MixedIndTests(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
