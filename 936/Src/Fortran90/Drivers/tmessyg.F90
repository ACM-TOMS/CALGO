program tmessyg ! A test driver illustrating us of the generic capability
  use  messy_m, only: messy, messy_ty_s, messy_ty_d, messy_ty_q, &
                              rk_s, rk_d, rk_q
  type(messy_ty_s) :: e_s
  type(messy_ty_d) :: e_d
  type(messy_ty_q) :: e_q
  call messy(e_s, "$NSingle:$Npi=$R$N e=$R",&
    & rdat=[4.0_rk_s * atan(1.0_rk_s), exp(1.0_rk_s)])
  call messy(e_d, "$NDouble:$Npi=$R$N e=$R",&
    & rdat=[4.0_rk_d * atan(1.0_rk_d), exp(1.0_rk_d)])
  call messy(e_q, "$NQuad:$Npi=$R$N e=$R",&
    & rdat=[4.0_rk_q * atan(1.0_rk_q), exp(1.0_rk_q)])
end program tmessyg
