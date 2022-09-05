#include "FLAME.h"
#include "rank_sv.h"

// ============================================================================
int rank_sv( FLA_Obj v, double rcond ) {
  FLA_Obj vT,              v0,
          vB,              nu1,
                           v2;
  int rank;

  FLA_Part_2x1( v,    &vT,
                      &vB,            0, FLA_TOP );

  //// printf( " rcond: %le \n", rcond );
  //// FLA_Obj_show( " sv = [ ", v, "%le", " ] ; " );

  rank = 0;
  while ( FLA_Obj_length( vT ) != FLA_Obj_length( v ) ){
    FLA_Repart_2x1_to_3x1( vT,                &v0,
                        /* ** */            /* *** */
                                              &nu1,
                           vB,                &v2,        1, FLA_BOTTOM );
    /*-----------------------------------------------------------------------*/

    //// printf( " sv(1): %le   ", *((double *) FLA_Obj_buffer( v ) ) );
    //// printf( " sv(i): %le ", *((double *) FLA_Obj_buffer( nu1 ) ) );
    //// printf( "\n" );

    if( ( rcond * *((double *)FLA_Obj_buffer_at_view( v ) ) ) <
        *((double *) FLA_Obj_buffer_at_view( nu1 ) ) ) {
      rank++;
    } else {
      break;
    }

    /*-----------------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &vT,                v0,
                                                  nu1,
                            /* ** */           /* *** */
                              &vB,                v2,     FLA_TOP );
  }
  //// printf( " Rank de sv: %d \n", rank );
  return rank;
}

