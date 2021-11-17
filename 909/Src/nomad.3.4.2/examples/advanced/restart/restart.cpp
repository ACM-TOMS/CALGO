/*--------------------------------------------------------------------------*/
/*  example of a program that makes NOMAD restarts after failed iterations  */
/*--------------------------------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

/*----------------------------------------*/
/*               the problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
private:

  double _mesh_update_basis;
  int    _initial_mesh_index;
  int    _mesh_index;
  Point  _initial_mesh_size;

public:
  My_Evaluator  ( const Parameters & p ) :
    Evaluator           ( p                                 ) ,
    _mesh_update_basis  ( p.get_mesh_update_basis().value() ) ,
    _initial_mesh_index ( p.get_initial_mesh_index()        ) ,
    _mesh_index         ( _initial_mesh_index               ) ,
    _initial_mesh_size  ( p.get_initial_mesh_size()         )   {}

  ~My_Evaluator ( void ) {}

  int get_mesh_index ( void ) const { return _mesh_index; }

  void get_mesh_size ( Point & mesh_size ) const
  {
    Mesh::get_delta_m ( mesh_size           ,
			_initial_mesh_size  ,
			_mesh_update_basis  ,
			_initial_mesh_index ,
			_mesh_index           );
  }

  virtual bool eval_x ( Eval_Point   & x          ,
			const Double & h_max      ,
			bool         & count_eval   ) const;

  virtual void update_iteration ( success_type              success      ,
				  const Stats             & stats        ,
				  const Evaluator_Control & ev_control   ,
				  const Barrier           & true_barrier ,
				  const Barrier           & sgte_barrier ,
				  const Pareto_Front      & pareto_front ,
				  bool                    & stop           );
};

/*----------------------------------------*/
/*           user-defined eval_x          */
/*----------------------------------------*/
bool My_Evaluator::eval_x ( Eval_Point   & x          ,
			    const Double & h_max      ,
			    bool         & count_eval   ) const
{
  Double c1 = 0.0 , c2 = 0.0;
  for ( int i = 0 ; i < 5 ; i++ ) {
    c1 += (x[i]-1).pow2();
    c2 += (x[i]+1).pow2();
  }
  x.set_bb_output  ( 0 , x[4]  ); // objective value
  x.set_bb_output  ( 1 , c1-25 ); // constraint 1
  x.set_bb_output  ( 2 , 25-c2 ); // constraint 2
  
  count_eval = true; // count a black-box evaluation
  
  return true;       // the evaluation succeeded
}

/*----------------------------------------*/
/*       updates at each iteration        */
/*----------------------------------------*/
void My_Evaluator::update_iteration ( success_type              success      ,
				      const Stats             & stats        ,
				      const Evaluator_Control & ev_control   ,
				      const Barrier           & true_barrier ,
				      const Barrier           & sgte_barrier ,
				      const Pareto_Front      & pareto_front ,
				      bool                    & stop           )
{
  _mesh_index = Mesh::get_mesh_index();
  if ( success == UNSUCCESSFUL )
    stop = true;
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
  // display:
  Display out ( std::cout );
  out.precision ( DISPLAY_PRECISION_STD );

  try {

    // NOMAD initializations:
    begin ( argc , argv );

    // parameters creation:
    Parameters p ( out );

    p.set_DIMENSION (5);             // number of variables

    vector<bb_output_type> bbot (3); // definition of
    bbot[0] = OBJ;                   // output types
    bbot[1] = EB;
    bbot[2] = EB;
    p.set_BB_OUTPUT_TYPE ( bbot );

    // p.set_DISPLAY_DEGREE ( FULL_DISPLAY );

    p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );

    p.set_X0 ( Point ( 5 , 0.0 ) );  // starting point

    p.set_LOWER_BOUND ( Point ( 5 , -6.0 ) ); // all var. >= -6
    Point ub ( 5 );                           // x_4 and x_5 have no bounds
    ub[0] = 5.0;                              // x_1 <= 5
    ub[1] = 6.0;                              // x_2 <= 6
    ub[2] = 7.0;                              // x_3 <= 7
    p.set_UPPER_BOUND ( ub );

    p.set_MAX_BB_EVAL (100); // the algorithm terminates after
                             // 100 black-box evaluations

    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev ( p );

    // best solutions:
    const Point * bf = NULL , * bi = NULL;

    // algorithm creation:
    Mads mads ( p , &ev );

    // successive runs:
    for ( int i = 0 ; i < 5 ; ++i ) {

      out << endl << open_block ( "MADS run #" + NOMAD::itos(i) );

      // not for the first run:
      if ( i > 0 ) {
	
	// new starting points:
	p.reset_X0();
	if ( bf )
	  p.set_X0 ( *bf );
	if ( bi )
	  p.set_X0 ( *bi );

	// initial mesh:
	p.set_INITIAL_MESH_INDEX ( ev.get_mesh_index() );
	Point initial_mesh_size;
	ev.get_mesh_size ( initial_mesh_size );
	p.set_INITIAL_MESH_SIZE ( initial_mesh_size );

	// parameters validation:
	p.check();

	// reset the Mads object:
	mads.reset ( true , true );
      }

      // the run:
      mads.run();
      bf = mads.get_best_feasible();
      bi = mads.get_best_infeasible();

      out.close_block();
    }
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  Slave::stop_slaves ( out );
  end();

  return EXIT_SUCCESS;
}
