module qmessycall_m
#define CMESSY_TY_ qmessy_ty
#define CALLMESSY_ callqmessy
#define GET_CMESSY_DEFAULTS_ get_qmessy_defaults
#define OPEN_CMESSY_FILES_ open_qmessy_files
#define CLOSE_CMESSY_FILES_ close_qmessy_files
#define ALLOCATE_CMESSY_INTERFACE_ allocate_qmessy_interface
#define DEALLOCATE_CMESSY_INTERFACE_ deallocate_qmessy_interface
#undef pdeft_
#define pdeft_ character, parameter :: T_="q"  
#define QCALLMESSY_ "callqmessy"
#define QGET_CMESSY_DEFAULTS_ "get_qmessy_defaults"
#define QOPEN_CMESSY_FILES_ "open_qmessy_files"
#define QLOSE_CMESSY_FILES_ "close_qmessy_files"
#define QALLOCATE_CMESSY_INTERFACE_ "allocate_qmessy_interface"
#define QDEALLOCATE_CMESSY_INTERFACE_ "deallocate_qmessy_interface"
#define Qmessy_ "qmessy"
#define CTYP_ C_LONG_DOUBLE
  use  messy_q_m, only: messy, messy_ty, rk, OUTPUT_UNIT, ERROR_UNIT
#include "cmessycall_bod.F90"
end module qmessycall_m
