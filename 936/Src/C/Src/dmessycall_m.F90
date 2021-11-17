module dmessycall_m
#define CMESSY_TY_ dmessy_ty
#define CALLMESSY_ calldmessy
#define GET_CMESSY_DEFAULTS_ get_dmessy_defaults
#define OPEN_CMESSY_FILES_ open_dmessy_files
#define CLOSE_CMESSY_FILES_ close_dmessy_files
#define ALLOCATE_CMESSY_INTERFACE_ allocate_dmessy_interface
#define DEALLOCATE_CMESSY_INTERFACE_ deallocate_dmessy_interface
#undef pdeft_
#define pdeft_ character, parameter :: T_="d"  
#define QCALLMESSY_ "calldmessy"
#define QGET_CMESSY_DEFAULTS_ "get_dmessy_defaults"
#define QOPEN_CMESSY_FILES_ "open_dmessy_files"
#define QLOSE_CMESSY_FILES_ "close_dmessy_files"
#define QALLOCATE_CMESSY_INTERFACE_ "allocate_dmessy_interface"
#define QDEALLOCATE_CMESSY_INTERFACE_ "deallocate_dmessy_interface"
#define Qmessy_ "dmessy"
#define CTYP_ C_DOUBLE
  use  messy_d_m, only: messy, messy_ty, rk, OUTPUT_UNIT, ERROR_UNIT
#include "cmessycall_bod.F90"
end module dmessycall_m
