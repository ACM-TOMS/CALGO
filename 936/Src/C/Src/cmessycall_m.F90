#define CMESSYCALL_ plet_ ## messycall_m
#define MESSYUSE_ messy_ ## plet_ ## _m

module CMESSYCALL_
  use  MESSYUSE_, only: messy, messy_ty, rk, OUTPUT_UNIT, ERROR_UNIT
#include "cmessycall_bod.F90"
end module CMESSYCALL_
