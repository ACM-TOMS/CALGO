struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

struct gsl_function_struct
{
  double * function(double x, void * params);
  void * params;
};


  double (* function) (double x, void * params);
  void * params;


  double * function(double x, void * params);
  void * params;



int a=5,c[2]={1,2},b=3;
