#include "pampac.h"
#include <string.h>
/**********************************************************************/
/* Uses input text file to determine problem and parameters required  */
/* by the algorithm to control and tune performance.                  */
/**********************************************************************/
bool
parse_options (int argc, char *argv[], options_struct *opts) {
  FILE *param_file;

  debug_print (4, opts, __func__, "Opening parameter file for parsing.\n");
  if (argc<2) {
    debug_print (0, opts, __func__,
                 "Parameter_file required from command-line input.\n");
    return false;
  }
  param_file = fopen (argv[1], "r");
  if (param_file == NULL) {
    debug_print (0, opts, __func__, "Error Opening File.\n");
    return false;
  }

  /* Parse param_file, one line at a time */
  size_t n_bytes = 0;
  char *line = NULL;
  ssize_t bytes_read = getline (&line, &n_bytes, param_file);
  while (bytes_read!=-1) {
    char param[7];
    bool has_parsed;
    char* parameter_name = malloc (bytes_read * sizeof (char));
    char* value = malloc (bytes_read * sizeof (char));
    if ( sscanf (line, "%s %s %s", param, parameter_name, value)==3 &&
         strcmp (param, "@param") == 0 )
      has_parsed = assign_options (parameter_name, value, opts);
    free (parameter_name);
    parameter_name = NULL;
    free (value);
    value = NULL;
    if (!has_parsed) {
      debug_print (0, opts, __func__, "Invalid option parsed.\n");
      fclose (param_file);
      return false;
    }
    bytes_read = getline (&line, &n_bytes, param_file);
  }
  fclose (param_file);

  /* Set opts.lambda_dir using opts.h_init */
  double h = opts->h_init;
  /* yields +1 or -1 (direction of initial step along curve) */
  opts->lambda_dir = (h > 0) ? +1.0 : -1.0;
  opts->h_init = (h > 0) ? h : -h;

  bool has_succeeded = validate_options (opts);
  return has_succeeded;
}
