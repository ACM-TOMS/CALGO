static double
minstd_get_double (void *vstate)
{
  return minstd_get (vstate) / 2147483647.0;
}

static void
minstd_set (void *vstate, unsigned long int s)
{
  minstd_state_t *state = (minstd_state_t *) vstate;

  if (s == 0)
    s = 1;  /* default seed is 1 */

  state->x = s;

  return;
}

static const gsl_rng_type minstd_type =
{"minstd",      /* name */
 2147483646,      /* RAND_MAX */
 1,               /* RAND_MIN */
 sizeof (minstd_state_t),
 &minstd_set,
 &minstd_get,
 &minstd_get_double};

const gsl_rng_type *gsl_rng_minstd = &minstd_type;
