// ============================================================================
// MyFLA_Utils:
//   Version:  0.32
//   Date:     2018-04-20
// ============================================================================
#include <math.h>
#include "FLAME.h"
#include "MyFLA_Utils.h"


// ============================================================================
// Declaration of local macros and prototypes.

#define my_abs( a )  ( (a) > 0.0 ? (a) : -(a) )


// ============================================================================
// Set contents of object A to zero.
FLA_Error MyFLA_Obj_set_to_zero( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Obj_set_to_zero:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set contents of object A to one.
FLA_Error MyFLA_Obj_set_to_one( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 1.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 1.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_set_to_one:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Extract the upper triangular part.
FLA_Error MyFLA_Triu( FLA_Obj A, int num_diags ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = max( 0, j - num_diags + 1 ); i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = max( 0, j - num_diags + 1 ); i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Triu:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Extract the lower triangular part.
FLA_Error MyFLA_Tril( FLA_Obj A, int num_diags ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < min( m_A, j - num_diags ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < min( m_A, j - num_diags ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Tril:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set contents of main diagonal of object A to one.
FLA_Error MyFLA_Set_main_diagonal_to_one( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the main diagonal.
      for ( j = 0; j < min( m_A, n_A ); j++ ) {
        buff_A[ j * rs_A + j * cs_A ] = 1.0f;
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the main diagonal.
      for ( j = 0; j < min( m_A, n_A ); j++ ) {
        buff_A[ j * rs_A + j * cs_A ] = 1.0;
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_set_to_one:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set contents of object A to the identity matrix.
FLA_Error MyFLA_Set_to_identity( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0f;
        }
      }
      // Set the main diagonal.
      for ( j = 0; j < min( m_A, n_A ); j++ ) {
        buff_A[ j * rs_A + j * cs_A ] = 1.0f;
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
      // Set the main diagonal.
      for ( j = 0; j < min( m_A, n_A ); j++ ) {
        buff_A[ j * rs_A + j * cs_A ] = 1.0;
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Set_to_identity:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Zero the strictly lower triangular part of a matrix, excluding the
// diagonal.
FLA_Error MyFLA_Zero_strict_lower_triangular( FLA_Obj A ) {

  MyFLA_Triu( A, 0 );
  return FLA_SUCCESS;
}


// ============================================================================
// Zero the strictly upper triangular part of a matrix, excluding the
// diagonal.
FLA_Error MyFLA_Zero_strict_upper_triangular( FLA_Obj A ) {

  MyFLA_Tril( A, 0 );
  return FLA_SUCCESS;
}


// ============================================================================
// Set the strictly lower triangular part of a matrix excluding the diagonal
// to a given float value.
FLA_Error MyFLA_Set_strict_lower_triangular_to_float( FLA_Obj A,
    float value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = j+1; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr,
             "+++ ERROR in MyFLA_Set_strict_lower_triangular_to_float:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set the strictly lower triangular part of a matrix excluding the diagonal
// to a given double value.
FLA_Error MyFLA_Set_strict_lower_triangular_to_double( FLA_Obj A,
    double value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = j+1; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr,
             "+++ ERROR in MyFLA_Set_strict_lower_triangular_to_double:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set the strictly upper triangular part of a matrix excluding the diagonal
// to a given float value.
FLA_Error MyFLA_Set_strict_upper_triangular_to_float( FLA_Obj A,
    float value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < min( j, m_A ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr,
             "+++ ERROR in MyFLA_Set_strict_upper_triangular_to_float:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set the strictly upper triangular part of a matrix excluding the diagonal
// to a given double value.
FLA_Error MyFLA_Set_strict_upper_triangular_to_double( FLA_Obj A,
    double value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < min( j, m_A ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr,
             "+++ ERROR in MyFLA_Set_strict_upper_triangular_to_double:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Sets the full contents of the object to a int value.
FLA_Error MyFLA_Obj_set_to_int( FLA_Obj A, int value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_INT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      int     * buff_A;

      buff_A  = ( int * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Obj_set_to_int:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Sets the full contents of the object to a float value.
FLA_Error MyFLA_Obj_set_to_float( FLA_Obj A, float value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Obj_set_to_float:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Sets the full contents of the object to a float value.
FLA_Error MyFLA_Obj_set_to_double( FLA_Obj A, double value ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = value;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Obj_set_to_double:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error MyFLA_Nrm1( FLA_Obj A, FLA_Obj nrm ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A, max, sum_col;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      max = 0.0f;
      for ( j = 0; j < n_A; j++ ) {
        sum_col = 0.0f;
        for ( i = 0; i < m_A; i++ ) {
          sum_col += my_abs( buff_A[ i * rs_A + j * cs_A ] );
        }
        if( sum_col > max )
          max = sum_col;
        //// printf( " i,j: %d %d  sum_c: %g  max: %g \n", i, j, sum_col, max );
      }
      *( ( float * ) FLA_Obj_buffer_at_view( nrm ) ) = max;
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A, max, sum_col;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      max = 0.0;
      for ( j = 0; j < n_A; j++ ) {
        sum_col = 0.0;
        for ( i = 0; i < m_A; i++ ) {
          sum_col += my_abs( buff_A[ i * rs_A + j * cs_A ] );
        }
        if( sum_col > max )
          max = sum_col;
        //// printf( " i,j: %d %d  sum_c: %g  max: %g \n", i, j, sum_col, max );
      }
      *( ( double * ) FLA_Obj_buffer_at_view( nrm ) ) = max;
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Nrm1:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Returns the frobenius norm of A into nrm.
FLA_Error MyFLA_Frob_norm( FLA_Obj A, FLA_Obj nrm ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A, sum, val;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      sum = 0.0f;
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          val = buff_A[ i * rs_A + j * cs_A ];
          sum += ( val * val );
        }
      }
      *( (float *) FLA_Obj_buffer_at_view( nrm ) ) = sqrtf( sum );
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A, sum, val;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      sum = 0.0;
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          val = buff_A[ i * rs_A + j * cs_A ];
          sum += ( val * val );
        }
      }
      *( (double *) FLA_Obj_buffer_at_view( nrm ) ) = sqrt( sum );
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Frob_norm:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Copy triu( A ) into triu( B ).
FLA_Error MyFLA_Copy_triu( FLA_Obj A, FLA_Obj B ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, rs_B, cs_B, i, j;
      float   * buff_A, * buff_B;

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );
      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      buff_B  = ( float * ) FLA_Obj_buffer_at_view( B );
      rs_B    = FLA_Obj_row_stride( B );
      cs_B    = FLA_Obj_col_stride( B );

      // Copy upper triangular part.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i <= min( j, m_A - 1 ); i++ ) {
          buff_B[ i * rs_B + j * cs_B ] = buff_A[ i * rs_A + j * cs_A ];
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, rs_B, cs_B, i, j;
      double  * buff_A, * buff_B;

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      buff_B  = ( double * ) FLA_Obj_buffer_at_view( B );
      rs_B    = FLA_Obj_row_stride( B );
      cs_B    = FLA_Obj_col_stride( B );

      // Copy upper triangular part.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i <= min( j, m_A - 1 ); i++ ) {
          buff_B[ i * rs_B + j * cs_B ] = buff_A[ i * rs_A + j * cs_A ];
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Copy_triu:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Copy A into B.
FLA_Error MyFLA_Copy( FLA_Obj A, FLA_Obj B ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, rs_B, cs_B, i, j;
      float   * buff_A, * buff_B;

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );
      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      buff_B  = ( float * ) FLA_Obj_buffer_at_view( B );
      rs_B    = FLA_Obj_row_stride( B );
      cs_B    = FLA_Obj_col_stride( B );

      // Copy from A to B.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_B[ i * rs_B + j * cs_B ] = buff_A[ i * rs_A + j * cs_A ];
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, rs_B, cs_B, i, j;
      double  * buff_A, * buff_B;

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      buff_B  = ( double * ) FLA_Obj_buffer_at_view( B );
      rs_B    = FLA_Obj_row_stride( B );
      cs_B    = FLA_Obj_col_stride( B );

      // Copy from A to B.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_B[ i * rs_B + j * cs_B ] = buff_A[ i * rs_A + j * cs_A ];
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Copy:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Set: A := abs( A ).
FLA_Error MyFLA_Abs( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] =
              my_abs( buff_A[ i * rs_A + j * cs_A ]);
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] =
              my_abs( buff_A[ i * rs_A + j * cs_A ] );
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Abs:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Symmetrize from lower matrix: Copy lower part to upper part.
FLA_Error MyFLA_Symmetrize_from_lower_matrix( FLA_Obj A ) {
  int  m_A, n_A;

  // Check matrix is square.
  m_A  = FLA_Obj_length( A );
  n_A  = FLA_Obj_width( A );
  if( m_A != n_A ) {
    fprintf( stderr, "+++ ERROR in MyFLA_Symmetrize_from_lower_matrix:  " );
    fprintf( stderr, "Non-square matrix. \n" );
    return -1;
  }

  // Symmetrize matrix: triu( A, 1 ) := tril( A, -1 ).
  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = j+1; i < min( m_A, n_A ); i++ ) {
          buff_A[ j * rs_A + i * cs_A ] = buff_A[ i * rs_A + j * cs_A ];
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      int     rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = j+1; i < min( m_A, n_A ); i++ ) {
          buff_A[ j * rs_A + i * cs_A ] = buff_A[ i * rs_A + j * cs_A ];
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Symmetrize_from_lower_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Create a matrix of FLA_Obj objects.
FLA_Error MyFLASH_Obj_create( FLA_Datatype datatype, int m, int n,
    FLA_Obj * A ) {
  FLA_Obj_create_ext( datatype, FLA_MATRIX, m, n, m, n, 0, 0, A );

  return FLA_SUCCESS;
}


// ============================================================================
// Copy A into B.
FLA_Error NoFLA_Copy_matrix_d( int m, int n, double * buff_A, int ldim_A,
    double * buff_B, int ldim_B ) {
  int  i, j;

  // Copy from A to B.
  for ( j = 0; j < n; j++ ) {
    for ( i = 0; i < m; i++ ) {
      buff_B[ i + j * ldim_B ] = buff_A[ i + j * ldim_A ];
    }
  }
  return FLA_SUCCESS;
}

// ============================================================================
// Generate random matrix, using rand.
FLA_Error MyFLA_Generate_random_matrix( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      float   * buff_A;

      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = ( float ) rand() / RAND_MAX;
        }
      }
    }
    break;


  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = ( double ) rand() / RAND_MAX;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Generate_random_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
// Generates SPD A.
FLA_Error MyFLA_Generate_spd_matrix( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      FLA_Obj  norm1, A_aux;

      FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & norm1 );
      FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & A_aux );

      FLA_Random_matrix( A_aux );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
                FLA_ONE, A_aux, A_aux, FLA_ZERO, A );
      MyFLA_Nrm1( A, norm1 );
      FLA_Shift_diag( FLA_NO_CONJUGATE, norm1, A );

      FLA_Obj_free( & norm1 );
      FLA_Obj_free( & A_aux );
    }
    break;

  case FLA_DOUBLE:
    {
      FLA_Obj  norm1, A_aux;

      FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & norm1 );
      FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & A_aux );

      FLA_Random_matrix( A_aux );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
                FLA_ONE, A_aux, A_aux, FLA_ZERO, A );
      MyFLA_Nrm1( A, norm1 );
      FLA_Shift_diag( FLA_NO_CONJUGATE, norm1, A );

      FLA_Obj_free( & norm1 );
      FLA_Obj_free( & A_aux );
    }
    break;


  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Generate_int_lower_triangular:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
// Generates lower triangular matrix A with ints starting at "num".
FLA_Error MyFLA_Generate_int_lower_triangular( int num, FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      float  * buff_A, fnum;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Initialize matrix A with integer values starting from "num".
      fnum = ( float ) num;
      for( j = 0; j < n_A; j++ ) {
        for( i = j; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      double * buff_A, fnum;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Initialize matrix A with integer values starting from "num".
      fnum = ( double ) num;
      for( j = 0; j < n_A; j++ ) {
        for( i = j; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0;
        }
      }
    }
    break;


  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Generate_int_lower_triangular:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
// Generates upper triangular matrix A with ints starting at "num".
FLA_Error MyFLA_Generate_int_upper_triangular( int num, FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      float  * buff_A, fnum;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Initialize matrix A with integer values starting from "num".
      fnum = ( float ) num;
      for( j = 0; j < n_A; j++ ) {
        for( i = 0; i < min( j, m_A ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      double * buff_A, fnum;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Initialize matrix A with integer values starting from "num".
      fnum = ( double ) num;
      for( j = 0; j < n_A; j++ ) {
        for( i = 0; i < min( j, m_A ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Generate_int_upper_triangular:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
// Generates dense matrix A with ints starting at "num".
FLA_Error MyFLA_Generate_int_matrix( int num, FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      float  * buff_A, fnum;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Initialize matrix A with integer values starting from "num".
      fnum = ( float ) num;
      for( j = 0; j < n_A; j++ ) {
        for( i = j; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0f;
        }
        for( i = 0; i < min( j, m_A ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0f;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      double * buff_A, fnum;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Initialize matrix A with integer values starting from "num".
      fnum = ( double ) num;
      for( j = 0; j < n_A; j++ ) {
        for( i = j; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0;
        }
        for( i = 0; i < min( j, m_A ); i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = fnum;
          fnum += 1.0;
        }
      }
    }
    break;


  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Generate_int_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
// Scales matrix down.
FLA_Error MyFLA_Scale_matrix_down( FLA_Obj A ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_FLOAT:
    {
      float  * buff_A, maxi, inv_maxi;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( float * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Compute maximum.
      maxi = buff_A[ 0 * rs_A + + 0 * cs_A ];
      for( j = 0; j < n_A; j++ ) {
        for( i = 0; i < m_A; i++ ) {
          if( buff_A[ i * rs_A + j * cs_A ] ) {
            maxi = buff_A[ i * rs_A + j * cs_A ];
          }
        }
      }
      //// printf( "Maxi: %f\n", maxi );

      // Scale matrix down.
      if( maxi != 0.0f ) {
        inv_maxi = 1.0f / maxi;
      } else {
        inv_maxi = 1.0f;
      }
      for( j = 0; j < n_A; j++ ) {
        for( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] *= inv_maxi;
        }
      }
    }
    break;

  case FLA_DOUBLE:
    {
      double * buff_A, maxi, inv_maxi;
      int    m_A, n_A, rs_A, cs_A, i, j;

      // Some initializations.
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Compute maximum.
      maxi = buff_A[ 0 * rs_A + 0 * cs_A ];
      for( j = 0; j < n_A; j++ ) {
        for( i = 0; i < m_A; i++ ) {
          if( buff_A[ i * rs_A + j * cs_A ] ) {
            maxi = buff_A[ i * rs_A + j * cs_A ];
          }
        }
      }
      //// printf( "Maxi: %f\n", maxi );

      // Scale matrix down.
      if( maxi != 0.0 ) {
        inv_maxi = 1.0 / maxi;
      } else {
        inv_maxi = 1.0;
      }
      for( j = 0; j < n_A; j++ ) {
        for( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] *= inv_maxi;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Scale_matrix_down:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
FLA_Error MyFLA_Set_matrix_main_diag_from_vector( FLA_Obj v, FLA_Obj A ) {
// Set main diagonal of matrix A from vector v.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      double  * buff_A, * buff_v;
      int     mn_A, ldim_A, j;

      // Some initializations.
      mn_A    = FLA_Obj_min_dim( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      buff_v  = ( double * ) FLA_Obj_buffer_at_view( v );
      ldim_A  = FLA_Obj_col_stride( A );

      // Main loop.
      for ( j = 0; j < mn_A; j++ ) {
        buff_A[ j + j * ldim_A ] = buff_v[ j ];
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Set_matrix_main_diag_from_vector:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
FLA_Error MyFLA_Set_vector_from_matrix_main_diag( FLA_Obj A, FLA_Obj v ) {
// Set vector v from main diagonal of matrix A.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      double  * buff_A, * buff_v;
      int     mn_A, ldim_A, j;

      // Some initializations.
      mn_A    = FLA_Obj_min_dim( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      buff_v  = ( double * ) FLA_Obj_buffer_at_view( v );
      ldim_A  = FLA_Obj_col_stride( A );

      // Main loop.
      for ( j = 0; j < mn_A; j++ ) {
        buff_v[ j ] = buff_A[ j + j * ldim_A ];
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Set_vector_from_matrix_main_diag:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


