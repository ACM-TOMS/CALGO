#include "pampac.h"

void delete_options (options_struct *opts) {

  if (opts->input_filename!=NULL) {
    debug_print (5, opts, __func__, "Erasing opts->input_filename (%s)\n",
                 opts->input_filename);
    free (opts->input_filename);
  }

  if (opts->tree_base_filename!=NULL) {
    debug_print (5, opts, __func__,
                 "Erasing opts->tree_base_filename (%s)\n",
                 opts->tree_base_filename);
    free (opts->tree_base_filename);
  }

  if (opts->scale_factors!=NULL) {
    debug_print (5, opts, __func__,
                 "Erasing opts->scale_factors (%g, %g, ...)\n",
                 opts->scale_factors[0], opts->scale_factors[1]);
    free (opts->scale_factors);
  }
  return;
}
