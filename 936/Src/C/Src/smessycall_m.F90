module smessycall_m
#define CMESSY_TY_ smessy_ty
#define CALLMESSY_ callsmessy
#define GET_CMESSY_DEFAULTS_ get_smessy_defaults
#define OPEN_CMESSY_FILES_ open_smessy_files
#define CLOSE_CMESSY_FILES_ close_smessy_files
#define ALLOCATE_CMESSY_INTERFACE_ allocate_smessy_interface
#define DEALLOCATE_CMESSY_INTERFACE_ deallocate_smessy_interface
#undef pdeft_
#define pdeft_ character, parameter :: T_="s"  
#define QCALLMESSY_ "callsmessy"
#define QGET_CMESSY_DEFAULTS_ "get_smessy_defaults"
#define QOPEN_CMESSY_FILES_ "open_smessy_files"
#define QLOSE_CMESSY_FILES_ "close_smessy_files"
#define QALLOCATE_CMESSY_INTERFACE_ "allocate_smessy_interface"
#define QDEALLOCATE_CMESSY_INTERFACE_ "deallocate_smessy_interface"
#define Qmessy_ "qmessy"
#define CTYP_ C_FLOAT
  use  messy_s_m, only: messy, messy_ty, rk, OUTPUT_UNIT, ERROR_UNIT
#include "cmessycall_bod.F90"
end module smessycall_m
