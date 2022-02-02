module Dose_Factors_m

  ! Dose factors from ICRP Publication 60, Table F_1.

  type :: Dose_Factor_t
    character(7) :: Sym
    real :: HalfLife
    character :: HalfLifeUnit
    real :: F_1_Infant
    real :: Dose_Infant     ! Sv/Bq
    real :: F_1             ! Ages 1-15
    real :: Doses(4)        ! Sv/Bq ages 1, 5, 10, 15
    real :: F1_Adult
    real :: Dose_Adult      ! Sv/Bq Adult
  end type Dose_Factor_t

  ! Dose factors for ingestion
  type(dose_factor_t), parameter :: H_3     = dose_factor_t("H-3", 12.35, "y", 1.0, 1.2E-10, 1.0, &
    & [ 1.2E-10, 7.3E-11, 5.7E-11, 4.2E-11], 1.0, 4.2E-11 )
  type(dose_factor_t), parameter :: Be_7    = dose_factor_t("Be-7", 53.3, "d", 0.02, 1.8E-10, 0.005, &
    & [ 1.3E-10, 7.7E-11, 5.3E-11, 3.5E-11], 0.005, 2.8E-11 )
  type(dose_factor_t), parameter :: Be_10   = dose_factor_t("Be-10", 1.6E6, "y", 0.02, 1.4E-08, 0.005, &
    & [ 8.0E-09, 4.1E-09, 2.4E-09, 1.4E-09], 0.005, 1.1E-09 )
  type(dose_factor_t), parameter :: C_11    = dose_factor_t("C-11", 20.38, "m", 1.0, 2.6E-10, 1.0, &
    & [ 1.5E-10, 7.3E-11, 4.3E-11, 3.0E-11], 1.0, 2.4E-11 )
  type(dose_factor_t), parameter :: C_14    = dose_factor_t("C-14", 5730, "y", 1.0, 1.4E-09, 1.0, &
    & [ 1.6E-09, 9.9E-10, 8.0E-10, 5.7E-10], 1.0, 5.8E-10 )
  type(dose_factor_t), parameter :: F_18    = dose_factor_t("F-18", 109.77, "m", 1.0, 5.2E-10, 1.0, &
    & [ 3.0E-10, 1.5E-10, 9.1E-11, 6.2E-11], 1.0, 4.9E-11 )
  type(dose_factor_t), parameter :: Na_22   = dose_factor_t("Na-22", 2.602, "y", 1.0, 2.1E-08, 1.0, &
    & [ 1.5E-08, 8.4E-09, 5.5E-09, 3.7E-09], 1.0, 3.2E-09 )
  type(dose_factor_t), parameter :: Na_24   = dose_factor_t("Na-24", 15.00, "h", 1.0, 3.5E-09, 1.0, &
    & [ 2.3E-09, 1.2E-09, 7.7E-10, 5.2E-10], 1.0, 4.3E-10 )
  type(dose_factor_t), parameter :: Mg_28   = dose_factor_t("Mg-28", 20.91, "h", 1.0, 1.2E-08, 0.5, &
    & [ 1.4E-08, 7.4E-09, 4.5E-09, 2.7E-09], 0.5, 2.2E-09 )
  type(dose_factor_t), parameter :: Al_26   = dose_factor_t("Al-26", 7.16E5, "y", 0.02, 3.4E-08, 0.01, &
    & [ 2.1E-08, 1.1E-08, 7.1E-09, 4.3E-09], 0.01, 3.5E-09 )
  type(dose_factor_t), parameter :: Si_31   = dose_factor_t("Si-31", 157.3, "m", 0.02, 1.9E-09, 0.01, &
    & [ 1.0E-09, 5.1E-10, 3.0E-10, 1.8E-10], 0.01, 1.6E-10 )
  type(dose_factor_t), parameter :: Si_32   = dose_factor_t("Si-32", 450, "y", 0.02, 7.3E-09, 0.01, &
    & [ 4.1E-09, 2.0E-09, 1.2E-09, 7.0E-10], 0.01, 5.6E-10 )
  type(dose_factor_t), parameter :: P_32    = dose_factor_t("P-32", 14.29, "d", 1.0, 3.1E-08, 0.8, &
    & [ 1.9E-08, 9.4E-09, 5.3E-09, 3.1E-09], 0.8, 2.4E-09 )
  type(dose_factor_t), parameter :: P_33    = dose_factor_t("P-33", 25.4, "d", 1.0, 2.7E-09, 0.8, &
    & [ 1.8E-09, 9.1E-10, 5.3E-10, 3.1E-10], 0.8, 2.4E-10 )
  type(dose_factor_t), parameter :: S_35    = dose_factor_t("S-35", 87.44, "d", 1.0, 7.7E-09, 1.0, &
    & [ 5.4E-09, 2.7E-09, 1.6E-09, 9.5E-10], 1.0, 7.7E-10 )
  type(dose_factor_t), parameter :: Cl_36   = dose_factor_t("Cl-36", 3.01E5, "y", 1.0, 9.8E-09, 1.0, &
    & [ 6.3E-09, 3.2E-09, 1.9E-09, 1.2E-09], 1.0, 9.3E-10 )
  type(dose_factor_t), parameter :: Cl_38   = dose_factor_t("Cl-38", 37.21, "m", 1.0, 1.4E-09, 1.0, &
    & [ 7.7E-10, 3.8E-10, 2.2E-10, 1.5E-10], 1.0, 1.2E-10 )
  type(dose_factor_t), parameter :: Cl_39   = dose_factor_t("Cl-39", 55.6, "m", 1.0, 9.7E-10, 1.0, &
    & [ 5.5E-10, 2.7E-10, 1.6E-10, 1.1E-10], 1.0, 8.5E-11 )
  type(dose_factor_t), parameter :: K_40    = dose_factor_t("K-40", 1.28E9, "y", 1.0, 6.2E-08, 1.0, &
    & [ 4.2E-08, 2.1E-08, 1.3E-08, 7.6E-09], 1.0, 6.2E-09 )
  type(dose_factor_t), parameter :: K_42    = dose_factor_t("K-42", 12.36, "h", 1.0, 5.1E-09, 1.0, &
    & [ 3.0E-09, 1.5E-09, 8.6E-10, 5.4E-10], 1.0, 4.3E-10 )
  type(dose_factor_t), parameter :: K_43    = dose_factor_t("K-43", 22.6, "h", 1.0, 2.3E-09, 1.0, &
    & [ 1.4E-09, 7.6E-10, 4.7E-10, 3.0E-10], 1.0, 2.5E-10 )
  type(dose_factor_t), parameter :: K_44    = dose_factor_t("K-44", 22.13, "m", 1.0, 1.0E-09, 1.0, &
    & [ 5.5E-10, 2.7E-10, 1.6E-10, 1.1E-10], 1.0, 8.4E-11 )
  type(dose_factor_t), parameter :: K_45    = dose_factor_t("K-45", 20, "m", 1.0, 6.2E-10, 1.0, &
    & [ 3.5E-10, 1.7E-10, 9.9E-11, 6.8E-11], 1.0, 5.4E-11 )
  type(dose_factor_t), parameter :: Ca_41   = dose_factor_t("Ca-41", 1.4E5, "y", 0.6, 1.2E-09, 0.4, &
    & [ 5.2E-10, 3.9E-10, 4.8E-10, 5.0E-10], 0.3, 1.9E-10 )
  type(dose_factor_t), parameter :: Ca_45   = dose_factor_t("Ca-45", 163, "d", 0.6, 1.1E-08, 0.4, &
    & [ 4.9E-09, 2.6E-09, 1.8E-09, 1.3E-09], 0.3, 7.1E-10 )
  type(dose_factor_t), parameter :: Ca_47   = dose_factor_t("Ca-47", 4.53, "d", 0.6, 1.3E-08, 0.4, &
    & [ 9.3E-09, 4.9E-09, 3.0E-09, 1.8E-09], 0.3, 1.6E-09 )
  type(dose_factor_t), parameter :: Sc_43   = dose_factor_t("Sc-43", 3.891, "h", 0.001, 1.8E-09, 0.0001, &
    & [ 1.2E-09, 6.1E-10, 3.7E-10, 2.3E-10], 0.0001, 1.9E-10 )
  type(dose_factor_t), parameter :: Sc_44   = dose_factor_t("Sc-44", 3.927, "h", 0.001, 3.5E-09, 0.0001, &
    & [ 2.2E-09, 1.2E-09, 7.1E-10, 4.4E-10], 0.0001, 3.5E-10 )
  type(dose_factor_t), parameter :: Sc_44m  = dose_factor_t("Sc-44m", 58.6, "h", 0.001, 2.4E-08, 0.0001, &
    & [ 1.6E-08, 8.3E-09, 5.1E-09, 3.1E-09], 0.0001, 2.4E-09 )
  type(dose_factor_t), parameter :: Sc_46   = dose_factor_t("Sc-46", 83.83, "d", 0.001, 1.1E-08, 0.0001, &
    & [ 7.9E-09, 4.4E-09, 2.9E-09, 1.8E-09], 0.0001, 1.5E-09 )
  type(dose_factor_t), parameter :: Sc_47   = dose_factor_t("Sc-47", 3.351, "d", 0.001, 6.1E-09, 0.0001, &
    & [ 3.9E-09, 2.0E-09, 1.2E-09, 6.8E-10], 0.0001, 5.4E-10 )
  type(dose_factor_t), parameter :: Sc_48   = dose_factor_t("Sc-48", 43.7, "h", 0.001, 1.3E-08, 0.0001, &
    & [ 9.3E-09, 5.1E-09, 3.3E-09, 2.1E-09], 0.0001, 1.7E-09 )
  type(dose_factor_t), parameter :: Sc_49   = dose_factor_t("Sc-49", 57.4, "m", 0.001, 1.0E-09, 0.0001, &
    & [ 5.7E-10, 2.8E-10, 1.6E-10, 1.0E-10], 0.0001, 8.2E-11 )
  type(dose_factor_t), parameter :: Ti_44   = dose_factor_t("Ti-44", 47.3, "y", 0.02, 5.5E-08, 0.01, &
    & [ 3.1E-08, 1.7E-08, 1.1E-08, 6.9E-09], 0.01, 5.8E-09 )
  type(dose_factor_t), parameter :: Ti_45   = dose_factor_t("Ti-45", 3.08, "h", 0.02, 1.6E-09, 0.01, &
    & [ 9.8E-10, 5.0E-10, 3.1E-10, 1.9E-10], 0.01, 1.5E-10 )
  type(dose_factor_t), parameter :: V_47    = dose_factor_t("V-47", 32.6, "m", 0.02, 7.3E-10, 0.01, &
    & [ 4.1E-10, 2.0E-10, 1.2E-10, 8.0E-11], 0.01, 6.3E-11 )
  type(dose_factor_t), parameter :: V_48    = dose_factor_t("V-48", 16.238, "d", 0.02, 1.5E-08, 0.01, &
    & [ 1.1E-08, 5.9E-09, 3.9E-09, 2.5E-09], 0.01, 2.0E-09 )
  type(dose_factor_t), parameter :: V_49    = dose_factor_t("V-49", 330, "d", 0.02, 2.2E-10, 0.01, &
    & [ 1.4E-10, 6.9E-11, 4.0E-11, 2.3E-11], 0.01, 1.8E-11 )
  type(dose_factor_t), parameter :: Cr_48   = dose_factor_t("Cr-48", 22.96, "h", 0.2, 1.4E-09, 0.1, &
    & [ 9.9E-10, 5.7E-10, 3.8E-10, 2.5E-10], 0.1, 2.0E-10 )
  type(dose_factor_t), parameter :: Cr_49   = dose_factor_t("Cr-49", 42.09, "m", 0.2, 6.8E-10, 0.1, &
    & [ 3.9E-10, 2.0E-10, 1.1E-10, 7.7E-11], 0.1, 6.1E-11 )
  type(dose_factor_t), parameter :: Cr_51   = dose_factor_t("Cr-51", 27.704, "d", 0.2, 3.5E-10, 0.1, &
    & [ 2.3E-10, 1.2E-10, 7.8E-11, 4.8E-11], 0.1, 3.8E-11 )
  type(dose_factor_t), parameter :: Mn_51   = dose_factor_t("Mn-51", 46.2, "m", 0.2, 1.1E-09, 0.1, &
    & [ 6.1E-10, 3.0E-10, 1.8E-10, 1.2E-10], 0.1, 9.3E-11 )
  type(dose_factor_t), parameter :: Mn_52   = dose_factor_t("Mn-52", 5.591, "d", 0.2, 1.2E-08, 0.1, &
    & [ 8.8E-09, 5.1E-09, 3.4E-09, 2.2E-09], 0.1, 1.8E-09 )
  type(dose_factor_t), parameter :: Mn_52m  = dose_factor_t("Mn-52m", 21.1, "m", 0.2, 7.8E-10, 0.1, &
    & [ 4.4E-10, 2.2E-10, 1.3E-10, 8.8E-11], 0.1, 6.9E-11 )
  type(dose_factor_t), parameter :: Mn_53   = dose_factor_t("Mn-53", 3.7E6, "y", 0.2, 4.1E-10, 0.1, &
    & [ 2.2E-10, 1.1E-10, 6.5E-11, 3.7E-11], 0.1, 3.0E-11 )
  type(dose_factor_t), parameter :: Mn_54   = dose_factor_t("Mn-54", 312.5, "d", 0.2, 5.4E-09, 0.1, &
    & [ 3.1E-09, 1.9E-09, 1.3E-09, 8.7E-10], 0.1, 7.1E-10 )
  type(dose_factor_t), parameter :: Mn_56   = dose_factor_t("Mn-56", 2.5785, "h", 0.2, 2.7E-09, 0.1, &
    & [ 1.7E-09, 8.5E-10, 5.1E-10, 3.2E-10], 0.1, 2.5E-10 )
  type(dose_factor_t), parameter :: Fe_52   = dose_factor_t("Fe-52", 8.275, "h", 0.6, 1.3E-08, 0.2, &
    & [ 9.1E-09, 4.6E-09, 2.8E-09, 1.7E-09], 0.3, 1.4E-09 )
  type(dose_factor_t), parameter :: Fe_55   = dose_factor_t("Fe-55", 2.7, "y", 0.6, 7.6E-09, 0.2, &
    & [ 2.4E-09, 1.7E-09, 1.1E-09, 7.7E-10], 0.3, 3.3E-10 )
  type(dose_factor_t), parameter :: Fe_59   = dose_factor_t("Fe-59", 44.529, "d", 0.6, 3.9E-08, 0.2, &
    & [ 1.3E-08, 7.5E-09, 4.7E-09, 3.1E-09], 0.3, 1.8E-09 )
  type(dose_factor_t), parameter :: Fe_60   = dose_factor_t("Fe-60", 1E5, "y", 0.6, 7.9E-07, 0.2, &
    & [ 2.7E-07, 2.7E-07, 2.5E-07, 2.3E-07], 0.3, 1.1E-07 )
  type(dose_factor_t), parameter :: Co_55   = dose_factor_t("Co-55", 17.54, "h", 0.6, 6.0E-09, 0.3, &
    & [ 5.5E-09, 2.9E-09, 1.8E-09, 1.1E-09], 0.3, 1.0E-09 )
  type(dose_factor_t), parameter :: Co_56   = dose_factor_t("Co-56", 78.76, "d", 0.6, 2.5E-08, 0.3, &
    & [ 1.5E-08, 8.8E-09, 5.8E-09, 3.8E-09], 0.3, 2.5E-09 )
  type(dose_factor_t), parameter :: Co_57   = dose_factor_t("Co-57", 270.9, "d", 0.6, 2.9E-09, 0.3, &
    & [ 1.6E-09, 8.9E-10, 5.8E-10, 3.7E-10], 0.3, 2.1E-10 )
  type(dose_factor_t), parameter :: Co_58   = dose_factor_t("Co-58", 70.80, "d", 0.6, 7.3E-09, 0.3, &
    & [ 4.4E-09, 2.6E-09, 1.7E-09, 1.1E-09], 0.3, 7.4E-10 )
  type(dose_factor_t), parameter :: Co_58m  = dose_factor_t("Co-58m", 9.15, "h", 0.6, 2.0E-10, 0.3, &
    & [ 1.5E-10, 7.8E-11, 4.7E-11, 2.8E-11], 0.3, 2.4E-11 )
  type(dose_factor_t), parameter :: Co_60   = dose_factor_t("Co-60", 5.271, "y", 0.6, 5.4E-08, 0.3, &
    & [ 2.7E-08, 1.7E-08, 1.1E-08, 7.9E-09], 0.3, 3.4E-09 )
  type(dose_factor_t), parameter :: Co_60m  = dose_factor_t("Co-60m", 10.47, "m", 0.6, 2.2E-11, 0.3, &
    & [ 1.2E-11, 5.7E-12, 3.2E-12, 2.2E-12], 0.3, 1.7E-12 )
  type(dose_factor_t), parameter :: Co_61   = dose_factor_t("Co-61", 1.65, "h", 0.6, 8.2E-10, 0.3, &
    & [ 5.1E-10, 2.5E-10, 1.4E-10, 9.2E-11], 0.3, 7.4E-11 )
  type(dose_factor_t), parameter :: Co_62m  = dose_factor_t("Co-62m", 13.91, "m", 0.6, 5.3E-10, 0.3, &
    & [ 3.0E-10, 1.5E-10, 8.7E-11, 6.0E-11], 0.3, 4.7E-11 )
  type(dose_factor_t), parameter :: Ni_56   = dose_factor_t("Ni-56", 6.10, "d", 0.1, 5.3E-09, 0.05, &
    & [ 4.0E-09, 2.3E-09, 1.6E-09, 1.1E-09], 0.05, 8.6E-10 )
  type(dose_factor_t), parameter :: Ni_57   = dose_factor_t("Ni-57", 36.08, "h", 0.1, 6.8E-09, 0.05, &
    & [ 4.9E-09, 2.7E-09, 1.7E-09, 1.1E-09], 0.05, 8.7E-10 )
  type(dose_factor_t), parameter :: Ni_59   = dose_factor_t("Ni-59", 7.5E4, "y", 0.1, 6.4E-10, 0.05, &
    & [ 3.4E-10, 1.9E-10, 1.1E-10, 7.3E-11], 0.05, 6.3E-11 )
  type(dose_factor_t), parameter :: Ni_63   = dose_factor_t("Ni-63", 96, "y", 0.1, 1.6E-09, 0.05, &
    & [ 8.4E-10, 4.6E-10, 2.8E-10, 1.8E-10], 0.05, 1.5E-10 )
  type(dose_factor_t), parameter :: Ni_65   = dose_factor_t("Ni-65", 2.520, "h", 0.1, 2.1E-09, 0.05, &
    & [ 1.3E-09, 6.3E-10, 3.8E-10, 2.3E-10], 0.05, 1.8E-10 )
  type(dose_factor_t), parameter :: Ni_66   = dose_factor_t("Ni-66", 54.6, "h", 0.1, 3.3E-08, 0.05, &
    & [ 2.2E-08, 1.1E-08, 6.6E-09, 3.7E-09], 0.05, 3.0E-09 )
  type(dose_factor_t), parameter :: Cu_60   = dose_factor_t("Cu-60", 23.2, "m", 1.0, 7.0E-10, 0.5, &
    & [ 4.2E-10, 2.2E-10, 1.3E-10, 8.9E-11], 0.5, 7.0E-11 )
  type(dose_factor_t), parameter :: Cu_61   = dose_factor_t("Cu-61", 3.408, "h", 1.0, 7.1E-10, 0.5, &
    & [ 7.5E-10, 3.9E-10, 2.3E-10, 1.5E-10], 0.5, 1.2E-10 )
  type(dose_factor_t), parameter :: Cu_64   = dose_factor_t("Cu-64", 12.701, "h", 1.0, 5.2E-10, 0.5, &
    & [ 8.3E-10, 4.2E-10, 2.5E-10, 1.5E-10], 0.5, 1.2E-10 )
  type(dose_factor_t), parameter :: Cu_67   = dose_factor_t("Cu-67", 61.86, "h", 1.0, 2.1E-09, 0.5, &
    & [ 2.4E-09, 1.2E-09, 7.2E-10, 4.2E-10], 0.5, 3.4E-10 )
  type(dose_factor_t), parameter :: Zn_62   = dose_factor_t("Zn-62", 9.26, "h", 1.0, 4.2E-09, 0.5, &
    & [ 6.5E-09, 3.3E-09, 2.0E-09, 1.2E-09], 0.5, 9.4E-10 )
  type(dose_factor_t), parameter :: Zn_63   = dose_factor_t("Zn-63", 38.1, "m", 1.0, 8.7E-10, 0.5, &
    & [ 5.2E-10, 2.6E-10, 1.5E-10, 1.0E-10], 0.5, 7.9E-11 )
  type(dose_factor_t), parameter :: Zn_65   = dose_factor_t("Zn-65", 243.9, "d", 1.0, 3.6E-08, 0.5, &
    & [ 1.6E-08, 9.7E-09, 6.4E-09, 4.5E-09], 0.5, 3.9E-09 )
  type(dose_factor_t), parameter :: Zn_69   = dose_factor_t("Zn-69", 57, "m", 1.0, 3.5E-10, 0.5, &
    & [ 2.2E-10, 1.1E-10, 6.0E-11, 3.9E-11], 0.5, 3.1E-11 )
  type(dose_factor_t), parameter :: Zn_69m  = dose_factor_t("Zn-69m", 13.76, "h", 1.0, 1.3E-09, 0.5, &
    & [ 2.3E-09, 1.2E-09, 7.0E-10, 4.1E-10], 0.5, 3.3E-10 )
  type(dose_factor_t), parameter :: Zn_71m  = dose_factor_t("Zn-71m", 3.92, "h", 1.0, 1.4E-09, 0.5, &
    & [ 1.5E-09, 7.8E-10, 4.8E-10, 3.0E-10], 0.5, 2.4E-10 )
  type(dose_factor_t), parameter :: Zn_72   = dose_factor_t("Zn-72", 46.5, "h", 1.0, 8.7E-09, 0.5, &
    & [ 8.6E-09, 4.5E-09, 2.8E-09, 1.7E-09], 0.5, 1.4E-09 )
  type(dose_factor_t), parameter :: Ga_65   = dose_factor_t("Ga-65", 15.2, "m", 0.01, 4.3E-10, 0.001, &
    & [ 2.4E-10, 1.2E-10, 6.9E-11, 4.7E-11], 0.001, 3.7E-11 )
  type(dose_factor_t), parameter :: Ga_66   = dose_factor_t("Ga-66", 9.40, "h", 0.01, 1.2E-08, 0.001, &
    & [ 7.9E-09, 4.0E-09, 2.5E-09, 1.5E-09], 0.001, 1.2E-09 )
  type(dose_factor_t), parameter :: Ga_67   = dose_factor_t("Ga-67", 78.26, "h", 0.01, 1.8E-09, 0.001, &
    & [ 1.2E-09, 6.4E-10, 4.0E-10, 2.4E-10], 0.001, 1.9E-10 )
  type(dose_factor_t), parameter :: Ga_68   = dose_factor_t("Ga-68", 68.0, "m", 0.01, 1.2E-09, 0.001, &
    & [ 6.7E-10, 3.4E-10, 2.0E-10, 1.3E-10], 0.001, 1.0E-10 )
  type(dose_factor_t), parameter :: Ga_70   = dose_factor_t("Ga-70", 21.15, "m", 0.01, 3.9E-10, 0.001, &
    & [ 2.2E-10, 1.0E-10, 5.9E-11, 4.0E-11], 0.001, 3.1E-11 )
  type(dose_factor_t), parameter :: Ga_72   = dose_factor_t("Ga-72", 14.1, "h", 0.01, 1.0E-08, 0.001, &
    & [ 6.8E-09, 3.6E-09, 2.2E-09, 1.4E-09], 0.001, 1.1E-09 )
  type(dose_factor_t), parameter :: Ga_73   = dose_factor_t("Ga-73", 4.91, "h", 0.01, 3.0E-09, 0.001, &
    & [ 1.9E-09, 9.3E-10, 5.5E-10, 3.3E-10], 0.001, 2.6E-10 )
  type(dose_factor_t), parameter :: Ge_66   = dose_factor_t("Ge-66", 2.27, "h", 1.0, 8.3E-10, 1.0, &
    & [ 5.3E-10, 2.9E-10, 1.9E-10, 1.3E-10], 1.0, 1.0E-10 )
  type(dose_factor_t), parameter :: Ge_67   = dose_factor_t("Ge-67", 18.7, "m", 1.0, 7.7E-10, 1.0, &
    & [ 4.2E-10, 2.1E-10, 1.2E-10, 8.2E-11], 1.0, 6.5E-11 )
  type(dose_factor_t), parameter :: Ge_68   = dose_factor_t("Ge-68", 288, "d", 1.0, 1.2E-08, 1.0, &
    & [ 8.0E-09, 4.2E-09, 2.6E-09, 1.6E-09], 1.0, 1.3E-09 )
  type(dose_factor_t), parameter :: Ge_69   = dose_factor_t("Ge-69", 39.05, "h", 1.0, 2.0E-09, 1.0, &
    & [ 1.3E-09, 7.1E-10, 4.6E-10, 3.0E-10], 1.0, 2.4E-10 )
  type(dose_factor_t), parameter :: Ge_71   = dose_factor_t("Ge-71", 11.8, "d", 1.0, 1.2E-10, 1.0, &
    & [ 7.8E-11, 4.0E-11, 2.4E-11, 1.5E-11], 1.0, 1.2E-11 )
  type(dose_factor_t), parameter :: Ge_75   = dose_factor_t("Ge-75", 82.78, "m", 1.0, 5.5E-10, 1.0, &
    & [ 3.1E-10, 1.5E-10, 8.7E-11, 5.9E-11], 1.0, 4.6E-11 )
  type(dose_factor_t), parameter :: Ge_77   = dose_factor_t("Ge-77", 11.30, "h", 1.0, 3.0E-09, 1.0, &
    & [ 1.8E-09, 9.9E-10, 6.2E-10, 4.1E-10], 1.0, 3.3E-10 )
  type(dose_factor_t), parameter :: Ge_78   = dose_factor_t("Ge-78", 87, "m", 1.0, 1.2E-09, 1.0, &
    & [ 7.0E-10, 3.6E-10, 2.2E-10, 1.5E-10], 1.0, 1.2E-10 )
  type(dose_factor_t), parameter :: As_69   = dose_factor_t("As-69", 15.2, "m", 1.0, 6.6E-10, 0.5, &
    & [ 3.7E-10, 1.8E-10, 1.1E-10, 7.2E-11], 0.5, 5.7E-11 )
  type(dose_factor_t), parameter :: As_70   = dose_factor_t("As-70", 52.6, "m", 1.0, 1.2E-09, 0.5, &
    & [ 7.8E-10, 4.1E-10, 2.5E-10, 1.7E-10], 0.5, 1.3E-10 )
  type(dose_factor_t), parameter :: As_71   = dose_factor_t("As-71", 64.8, "h", 1.0, 2.8E-09, 0.5, &
    & [ 2.8E-09, 1.5E-09, 9.3E-10, 5.7E-10], 0.5, 4.6E-10 )
  type(dose_factor_t), parameter :: As_72   = dose_factor_t("As-72", 26.0, "h", 1.0, 1.1E-08, 0.5, &
    & [ 1.2E-08, 6.3E-09, 3.8E-09, 2.3E-09], 0.5, 1.8E-09 )
  type(dose_factor_t), parameter :: As_73   = dose_factor_t("As-73", 80.30, "d", 1.0, 2.6E-09, 0.5, &
    & [ 1.9E-09, 9.3E-10, 5.6E-10, 3.2E-10], 0.5, 2.6E-10 )
  type(dose_factor_t), parameter :: As_74   = dose_factor_t("As-74", 17.76, "d", 1.0, 1.0E-08, 0.5, &
    & [ 8.2E-09, 4.3E-09, 2.6E-09, 1.6E-09], 0.5, 1.3E-09 )
  type(dose_factor_t), parameter :: As_76   = dose_factor_t("As-76", 26.32, "h", 1.0, 1.0E-08, 0.5, &
    & [ 1.1E-08, 5.8E-09, 3.4E-09, 2.0E-09], 0.5, 1.6E-09 )
  type(dose_factor_t), parameter :: As_77   = dose_factor_t("As-77", 38.8, "h", 1.0, 2.7E-09, 0.5, &
    & [ 2.9E-09, 1.5E-09, 8.7E-10, 5.0E-10], 0.5, 4.0E-10 )
  type(dose_factor_t), parameter :: As_78   = dose_factor_t("As-78", 90.7, "m", 1.0, 2.0E-09, 0.5, &
    & [ 1.4E-09, 7.0E-10, 4.1E-10, 2.7E-10], 0.5, 2.1E-10 )
  type(dose_factor_t), parameter :: Se_70   = dose_factor_t("Se-70", 41.0, "m", 1.0, 1.0E-09, 0.8, &
    & [ 7.1E-10, 3.6E-10, 2.2E-10, 1.5E-10], 0.8, 1.2E-10 )
  type(dose_factor_t), parameter :: Se_73   = dose_factor_t("Se-73", 7.15, "h", 1.0, 1.6E-09, 0.8, &
    & [ 1.4E-09, 7.4E-10, 4.8E-10, 2.5E-10], 0.8, 2.1E-10 )
  type(dose_factor_t), parameter :: Se_73m  = dose_factor_t("Se-73m", 39, "m", 1.0, 2.6E-10, 0.8, &
    & [ 1.8E-10, 9.5E-11, 5.9E-11, 3.5E-11], 0.8, 2.8E-11 )
  type(dose_factor_t), parameter :: Se_75   = dose_factor_t("Se-75", 119.8, "d", 1.0, 2.0E-08, 0.8, &
    & [ 1.3E-08, 8.3E-09, 6.0E-09, 3.1E-09], 0.8, 2.6E-09 )
  type(dose_factor_t), parameter :: Se_79   = dose_factor_t("Se-79", 65000, "y", 1.0, 4.1E-08, 0.8, &
    & [ 2.8E-08, 1.9E-08, 1.4E-08, 4.1E-09], 0.8, 2.9E-09 )
  type(dose_factor_t), parameter :: Se_81   = dose_factor_t("Se-81", 18.5, "m", 1.0, 3.4E-10, 0.8, &
    & [ 1.9E-10, 9.0E-11, 5.1E-11, 3.4E-11], 0.8, 2.7E-11 )
  type(dose_factor_t), parameter :: Se_81m  = dose_factor_t("Se-81m", 57.25, "m", 1.0, 6.0E-10, 0.8, &
    & [ 3.7E-10, 1.8E-10, 1.1E-10, 6.7E-11], 0.8, 5.3E-11 )
  type(dose_factor_t), parameter :: Se_83   = dose_factor_t("Se-83", 22.5, "m", 1.0, 4.6E-10, 0.8, &
    & [ 2.9E-10, 1.5E-10, 8.7E-11, 5.9E-11], 0.8, 4.7E-11 )
  type(dose_factor_t), parameter :: Br_74   = dose_factor_t("Br-74", 25.3, "m", 1.0, 9.0E-10, 1.0, &
    & [ 5.2E-10, 2.6E-10, 1.5E-10, 1.1E-10], 1.0, 8.4E-11 )
  type(dose_factor_t), parameter :: Br_74m  = dose_factor_t("Br-74m", 41.5, "m", 1.0, 1.5E-09, 1.0, &
    & [ 8.5E-10, 4.3E-10, 2.5E-10, 1.7E-10], 1.0, 1.4E-10 )
  type(dose_factor_t), parameter :: Br_75   = dose_factor_t("Br-75", 98, "m", 1.0, 8.5E-10, 1.0, &
    & [ 4.9E-10, 2.5E-10, 1.5E-10, 9.9E-11], 1.0, 7.9E-11 )
  type(dose_factor_t), parameter :: Br_76   = dose_factor_t("Br-76", 16.2, "h", 1.0, 4.2E-09, 1.0, &
    & [ 2.7E-09, 1.4E-09, 8.7E-10, 5.6E-10], 1.0, 4.6E-10 )
  type(dose_factor_t), parameter :: Br_77   = dose_factor_t("Br-77", 56, "h", 1.0, 6.3E-10, 1.0, &
    & [ 4.4E-10, 2.5E-10, 1.7E-10, 1.1E-10], 1.0, 9.6E-11 )
  type(dose_factor_t), parameter :: Br_80   = dose_factor_t("Br-80", 17.4, "m", 1.0, 3.9E-10, 1.0, &
    & [ 2.1E-10, 1.0E-10, 5.8E-11, 3.9E-11], 1.0, 3.1E-11 )
  type(dose_factor_t), parameter :: Br_80m  = dose_factor_t("Br-80m", 4.42, "h", 1.0, 1.4E-09, 1.0, &
    & [ 8.0E-10, 3.9E-10, 2.3E-10, 1.4E-10], 1.0, 1.1E-10 )
  type(dose_factor_t), parameter :: Br_82   = dose_factor_t("Br-82", 35.30, "h", 1.0, 3.7E-09, 1.0, &
    & [ 2.6E-09, 1.5E-09, 9.5E-10, 6.4E-10], 1.0, 5.4E-10 )
  type(dose_factor_t), parameter :: Br_83   = dose_factor_t("Br-83", 2.39, "h", 1.0, 5.3E-10, 1.0, &
    & [ 3.0E-10, 1.4E-10, 8.3E-11, 5.5E-11], 1.0, 4.3E-11 )
  type(dose_factor_t), parameter :: Br_84   = dose_factor_t("Br-84", 31.80, "m", 1.0, 1.0E-09, 1.0, &
    & [ 5.8E-10, 2.8E-10, 1.6E-10, 1.1E-10], 1.0, 8.8E-11 )
  type(dose_factor_t), parameter :: Rb_79   = dose_factor_t("Rb-79", 22.9, "m", 1.0, 5.7E-10, 1.0, &
    & [ 3.2E-10, 1.6E-10, 9.2E-11, 6.3E-11], 1.0, 5.0E-11 )
  type(dose_factor_t), parameter :: Rb_81   = dose_factor_t("Rb-81", 4.58, "h", 1.0, 5.4E-10, 1.0, &
    & [ 3.2E-10, 1.6E-10, 1.0E-10, 6.7E-11], 1.0, 5.4E-11 )
  type(dose_factor_t), parameter :: Rb_81m  = dose_factor_t("Rb-81m", 32, "m", 1.0, 1.1E-10, 1.0, &
    & [ 6.2E-11, 3.1E-11, 1.8E-11, 1.2E-11], 1.0, 9.7E-12 )
  type(dose_factor_t), parameter :: Rb_82m  = dose_factor_t("Rb-82m", 6.2, "h", 1.0, 8.7E-10, 1.0, &
    & [ 5.9E-10, 3.4E-10, 2.2E-10, 1.5E-10], 1.0, 1.3E-10 )
  type(dose_factor_t), parameter :: Rb_83   = dose_factor_t("Rb-83", 86.2, "d", 1.0, 1.1E-08, 1.0, &
    & [ 8.4E-09, 4.9E-09, 3.2E-09, 2.2E-09], 1.0, 1.9E-09 )
  type(dose_factor_t), parameter :: Rb_84   = dose_factor_t("Rb-84", 32.77, "d", 1.0, 2.0E-08, 1.0, &
    & [ 1.4E-08, 7.9E-09, 5.0E-09, 3.3E-09], 1.0, 2.8E-09 )
  type(dose_factor_t), parameter :: Rb_86   = dose_factor_t("Rb-86", 18.66, "d", 1.0, 3.1E-08, 1.0, &
    & [ 2.0E-08, 9.9E-09, 5.9E-09, 3.5E-09], 1.0, 2.8E-09 )
  type(dose_factor_t), parameter :: Rb_87   = dose_factor_t("Rb-87", 4.7E10, "y", 1.0, 1.5E-08, 1.0, &
    & [ 1.0E-08, 5.2E-09, 3.1E-09, 1.8E-09], 1.0, 1.5E-09 )
  type(dose_factor_t), parameter :: Rb_88   = dose_factor_t("Rb-88", 17.8, "m", 1.0, 1.1E-09, 1.0, &
    & [ 6.2E-10, 3.0E-10, 1.7E-10, 1.2E-10], 1.0, 9.0E-11 )
  type(dose_factor_t), parameter :: Rb_89   = dose_factor_t("Rb-89", 15.2, "m", 1.0, 5.4E-10, 1.0, &
    & [ 3.0E-10, 1.5E-10, 8.6E-11, 5.9E-11], 1.0, 4.7E-11 )
  type(dose_factor_t), parameter :: Sr_80   = dose_factor_t("Sr-80", 100, "m", 0.6, 3.7E-09, 0.4, &
    & [ 2.3E-09, 1.1E-09, 6.5E-10, 4.2E-10], 0.3, 3.4E-10 )
  type(dose_factor_t), parameter :: Sr_81   = dose_factor_t("Sr-81", 25.5, "m", 0.6, 8.4E-10, 0.4, &
    & [ 4.9E-10, 2.4E-10, 1.4E-10, 9.6E-11], 0.3, 7.7E-11 )
  type(dose_factor_t), parameter :: Sr_82   = dose_factor_t("Sr-82", 25.0, "d", 0.6, 7.2E-08, 0.4, &
    & [ 4.1E-08, 2.1E-08, 1.3E-08, 8.7E-09], 0.3, 6.1E-09 )
  type(dose_factor_t), parameter :: Sr_83   = dose_factor_t("Sr-83", 32.4, "h", 0.6, 3.4E-09, 0.4, &
    & [ 2.7E-09, 1.4E-09, 9.1E-10, 5.7E-10], 0.3, 4.9E-10 )
  type(dose_factor_t), parameter :: Sr_85   = dose_factor_t("Sr-85", 64.84, "d", 0.6, 7.7E-09, 0.4, &
    & [ 3.1E-09, 1.7E-09, 1.5E-09, 1.3E-09], 0.3, 5.6E-10 )
  type(dose_factor_t), parameter :: Sr_85m  = dose_factor_t("Sr-85m", 69.5, "m", 0.6, 4.5E-11, 0.4, &
    & [ 3.0E-11, 1.7E-11, 1.1E-11, 7.8E-12], 0.3, 6.1E-12 )
  type(dose_factor_t), parameter :: Sr_87m  = dose_factor_t("Sr-87m", 2.805, "h", 0.6, 2.4E-10, 0.4, &
    & [ 1.7E-10, 9.0E-11, 5.6E-11, 3.6E-11], 0.3, 3.0E-11 )
  type(dose_factor_t), parameter :: Sr_89   = dose_factor_t("Sr-89", 50.5, "d", 0.6, 3.6E-08, 0.4, &
    & [ 1.8E-08, 8.9E-09, 5.8E-09, 4.0E-09], 0.3, 2.6E-09 )
  type(dose_factor_t), parameter :: Sr_90   = dose_factor_t("Sr-90", 29.12, "y", 0.6, 2.3E-07, 0.4, &
    & [ 7.3E-08, 4.7E-08, 6.0E-08, 8.0E-08], 0.3, 2.8E-08 )
  type(dose_factor_t), parameter :: Sr_91   = dose_factor_t("Sr-91", 9.5, "h", 0.6, 5.2E-09, 0.4, &
    & [ 4.0E-09, 2.1E-09, 1.2E-09, 7.4E-10], 0.3, 6.5E-10 )
  type(dose_factor_t), parameter :: Sr_92   = dose_factor_t("Sr-92", 2.71, "h", 0.6, 3.4E-09, 0.4, &
    & [ 2.7E-09, 1.4E-09, 8.2E-10, 4.8E-10], 0.3, 4.3E-10 )
  type(dose_factor_t), parameter :: Y_86    = dose_factor_t("Y-86", 14.74, "h", 0.001, 7.6E-09, 0.0001, &
    & [ 5.2E-09, 2.9E-09, 1.9E-09, 1.2E-09], 0.0001, 9.6E-10 )
  type(dose_factor_t), parameter :: Y_86m   = dose_factor_t("Y-86m", 48, "m", 0.001, 4.5E-10, 0.0001, &
    & [ 3.1E-10, 1.7E-10, 1.1E-10, 7.1E-11], 0.0001, 5.6E-11 )
  type(dose_factor_t), parameter :: Y_87    = dose_factor_t("Y-87", 80.3, "h", 0.001, 4.6E-09, 0.0001, &
    & [ 3.2E-09, 1.8E-09, 1.1E-09, 7.0E-10], 0.0001, 5.5E-10 )
  type(dose_factor_t), parameter :: Y_88    = dose_factor_t("Y-88", 106.64, "d", 0.001, 8.1E-09, 0.0001, &
    & [ 6.0E-09, 3.5E-09, 2.4E-09, 1.6E-09], 0.0001, 1.3E-09 )
  type(dose_factor_t), parameter :: Y_90    = dose_factor_t("Y-90", 64.0, "h", 0.001, 3.1E-08, 0.0001, &
    & [ 2.0E-08, 1.0E-08, 5.9E-09, 3.3E-09], 0.0001, 2.7E-09 )
  type(dose_factor_t), parameter :: Y_90m   = dose_factor_t("Y-90m", 3.19, "h", 0.001, 1.8E-09, 0.0001, &
    & [ 1.2E-09, 6.1E-10, 3.7E-10, 2.2E-10], 0.0001, 1.7E-10 )
  type(dose_factor_t), parameter :: Y_91    = dose_factor_t("Y-91", 58.51, "d", 0.001, 2.8E-08, 0.0001, &
    & [ 1.8E-08, 8.8E-09, 5.2E-09, 2.9E-09], 0.0001, 2.4E-09 )
  type(dose_factor_t), parameter :: Y_91m   = dose_factor_t("Y-91m", 49.71, "m", 0.001, 9.2E-11, 0.0001, &
    & [ 6.0E-11, 3.3E-11, 2.1E-11, 1.4E-11], 0.0001, 1.1E-11 )
  type(dose_factor_t), parameter :: Y_92    = dose_factor_t("Y-92", 3.54, "h", 0.001, 5.9E-09, 0.0001, &
    & [ 3.6E-09, 1.8E-09, 1.0E-09, 6.2E-10], 0.0001, 4.9E-10 )
  type(dose_factor_t), parameter :: Y_93    = dose_factor_t("Y-93", 10.1, "h", 0.001, 1.4E-08, 0.0001, &
    & [ 8.5E-09, 4.3E-09, 2.5E-09, 1.4E-09], 0.0001, 1.2E-09 )
  type(dose_factor_t), parameter :: Y_94    = dose_factor_t("Y-94", 19.1, "m", 0.001, 9.9E-10, 0.0001, &
    & [ 5.5E-10, 2.7E-10, 1.5E-10, 1.0E-10], 0.0001, 8.1E-11 )
  type(dose_factor_t), parameter :: Y_95    = dose_factor_t("Y-95", 10.7, "m", 0.001, 5.7E-10, 0.0001, &
    & [ 3.1E-10, 1.5E-10, 8.7E-11, 5.9E-11], 0.0001, 4.6E-11 )
  type(dose_factor_t), parameter :: Zr_86   = dose_factor_t("Zr-86", 16.5, "h", 0.02, 6.9E-09, 0.01, &
    & [ 4.8E-09, 2.7E-09, 1.7E-09, 1.1E-09], 0.01, 8.6E-10 )
  type(dose_factor_t), parameter :: Zr_88   = dose_factor_t("Zr-88", 83.4, "d", 0.02, 2.8E-09, 0.01, &
    & [ 2.0E-09, 1.2E-09, 8.0E-10, 5.4E-10], 0.01, 4.5E-10 )
  type(dose_factor_t), parameter :: Zr_89   = dose_factor_t("Zr-89", 78.43, "h", 0.02, 6.5E-09, 0.01, &
    & [ 4.5E-09, 2.5E-09, 1.6E-09, 9.9E-10], 0.01, 7.9E-10 )
  type(dose_factor_t), parameter :: Zr_93   = dose_factor_t("Zr-93", 1.53E6, "y", 0.02, 1.2E-09, 0.01, &
    & [ 7.6E-10, 5.1E-10, 5.8E-10, 8.6E-10], 0.01, 1.1E-09 )
  type(dose_factor_t), parameter :: Zr_95   = dose_factor_t("Zr-95", 63.98, "d", 0.02, 8.5E-09, 0.01, &
    & [ 5.6E-09, 3.0E-09, 1.9E-09, 1.2E-09], 0.01, 9.5E-10 )
  type(dose_factor_t), parameter :: Zr_97   = dose_factor_t("Zr-97", 16.90, "h", 0.02, 2.2E-08, 0.01, &
    & [ 1.4E-08, 7.3E-09, 4.4E-09, 2.6E-09], 0.01, 2.1E-09 )
  type(dose_factor_t), parameter :: Nb_88   = dose_factor_t("Nb-88", 14.3, "m", 0.02, 6.7E-10, 0.01, &
    & [ 3.8E-10, 1.9E-10, 1.1E-10, 7.9E-11], 0.01, 6.3E-11 )
  type(dose_factor_t), parameter :: Nb_89   = dose_factor_t("Nb-89", 122, "m", 0.02, 3.0E-09, 0.01, &
    & [ 2.0E-09, 1.0E-09, 6.0E-10, 3.4E-10], 0.01, 2.7E-10 )
  type(dose_factor_t), parameter :: Nb_89m  = dose_factor_t("Nb-89m", 66, "m", 0.02, 1.5E-09, 0.01, &
    & [ 8.7E-10, 4.4E-10, 2.7E-10, 1.8E-10], 0.01, 1.4E-10 )
  type(dose_factor_t), parameter :: Nb_90   = dose_factor_t("Nb-90", 14.60, "h", 0.02, 1.1E-08, 0.01, &
    & [ 7.2E-09, 3.9E-09, 2.5E-09, 1.6E-09], 0.01, 1.2E-09 )
  type(dose_factor_t), parameter :: Nb_93m  = dose_factor_t("Nb-93m", 13.6, "y", 0.02, 1.5E-09, 0.01, &
    & [ 9.1E-10, 4.6E-10, 2.7E-10, 1.5E-10], 0.01, 1.2E-10 )
  type(dose_factor_t), parameter :: Nb_94   = dose_factor_t("Nb-94", 2.03E4, "y", 0.02, 1.5E-08, 0.01, &
    & [ 9.7E-09, 5.3E-09, 3.4E-09, 2.1E-09], 0.01, 1.7E-09 )
  type(dose_factor_t), parameter :: Nb_95   = dose_factor_t("Nb-95", 35.15, "d", 0.02, 4.6E-09, 0.01, &
    & [ 3.2E-09, 1.8E-09, 1.1E-09, 7.4E-10], 0.01, 5.8E-10 )
  type(dose_factor_t), parameter :: Nb_95m  = dose_factor_t("Nb-95m", 86.6, "h", 0.02, 6.4E-09, 0.01, &
    & [ 4.1E-09, 2.1E-09, 1.2E-09, 7.1E-10], 0.01, 5.6E-10 )
  type(dose_factor_t), parameter :: Nb_96   = dose_factor_t("Nb-96", 23.35, "h", 0.02, 9.2E-09, 0.01, &
    & [ 6.3E-09, 3.4E-09, 2.2E-09, 1.4E-09], 0.01, 1.1E-09 )
  type(dose_factor_t), parameter :: Nb_97   = dose_factor_t("Nb-97", 72.1, "m", 0.02, 7.7E-10, 0.01, &
    & [ 4.5E-10, 2.3E-10, 1.3E-10, 8.7E-11], 0.01, 6.8E-11 )
  type(dose_factor_t), parameter :: Nb_98m  = dose_factor_t("Nb-98m", 51.5, "m", 0.02, 1.2E-09, 0.01, &
    & [ 7.1E-10, 3.6E-10, 2.2E-10, 1.4E-10], 0.01, 1.1E-10 )
  type(dose_factor_t), parameter :: Mo_90   = dose_factor_t("Mo-90", 5.67, "h", 1.0, 1.7E-09, 1.0, &
    & [ 1.2E-09, 6.3E-10, 4.0E-10, 2.7E-10], 1.0, 2.2E-10 )
  type(dose_factor_t), parameter :: Mo_93   = dose_factor_t("Mo-93", 3.5E3, "y", 1.0, 7.9E-09, 1.0, &
    & [ 6.9E-09, 5.0E-09, 4.0E-09, 3.4E-09], 1.0, 3.1E-09 )
  type(dose_factor_t), parameter :: Mo_93m  = dose_factor_t("Mo-93m", 6.85, "h", 1.0, 8.0E-10, 1.0, &
    & [ 5.4E-10, 3.1E-10, 2.0E-10, 1.4E-10], 1.0, 1.1E-10 )
  type(dose_factor_t), parameter :: Mo_99   = dose_factor_t("Mo-99", 66.0, "h", 1.0, 5.5E-09, 1.0, &
    & [ 3.5E-09, 1.8E-09, 1.1E-09, 7.6E-10], 1.0, 6.0E-10 )
  type(dose_factor_t), parameter :: Mo_101  = dose_factor_t("Mo-101", 14.62, "m", 1.0, 4.8E-10, 1.0, &
    & [ 2.7E-10, 1.3E-10, 7.6E-11, 5.2E-11], 1.0, 4.1E-11 )
  type(dose_factor_t), parameter :: Tc_93   = dose_factor_t("Tc-93", 2.75, "h", 1.0, 2.7E-10, 0.5, &
    & [ 2.5E-10, 1.5E-10, 9.8E-11, 6.8E-11], 0.5, 5.5E-11 )
  type(dose_factor_t), parameter :: Tc_93m  = dose_factor_t("Tc-93m", 43.5, "m", 1.0, 2.0E-10, 0.5, &
    & [ 1.3E-10, 7.3E-11, 4.6E-11, 3.2E-11], 0.5, 2.5E-11 )
  type(dose_factor_t), parameter :: Tc_94   = dose_factor_t("Tc-94", 293, "m", 1.0, 1.2E-09, 0.5, &
    & [ 1.0E-09, 5.8E-10, 3.7E-10, 2.5E-10], 0.5, 2.0E-10 )
  type(dose_factor_t), parameter :: Tc_94m  = dose_factor_t("Tc-94m", 52, "m", 1.0, 1.3E-09, 0.5, &
    & [ 6.5E-10, 3.3E-10, 1.9E-10, 1.3E-10], 0.5, 1.0E-10 )
  type(dose_factor_t), parameter :: Tc_95   = dose_factor_t("Tc-95", 20.0, "h", 1.0, 9.9E-10, 0.5, &
    & [ 8.7E-10, 5.0E-10, 3.3E-10, 2.3E-10], 0.5, 1.8E-10 )
  type(dose_factor_t), parameter :: Tc_95m  = dose_factor_t("Tc-95m", 61, "d", 1.0, 4.7E-09, 0.5, &
    & [ 2.8E-09, 1.6E-09, 1.0E-09, 7.0E-10], 0.5, 5.6E-10 )
  type(dose_factor_t), parameter :: Tc_96   = dose_factor_t("Tc-96", 4.28, "d", 1.0, 6.7E-09, 0.5, &
    & [ 5.1E-09, 3.0E-09, 2.0E-09, 1.4E-09], 0.5, 1.1E-09 )
  type(dose_factor_t), parameter :: Tc_96m  = dose_factor_t("Tc-96m", 51.5, "m", 1.0, 1.0E-10, 0.5, &
    & [ 6.5E-11, 3.6E-11, 2.3E-11, 1.6E-11], 0.5, 1.2E-11 )
  type(dose_factor_t), parameter :: Tc_97   = dose_factor_t("Tc-97", 2.6E6, "y", 1.0, 9.9E-10, 0.5, &
    & [ 4.9E-10, 2.4E-10, 1.4E-10, 8.8E-11], 0.5, 6.8E-11 )
  type(dose_factor_t), parameter :: Tc_97m  = dose_factor_t("Tc-97m", 87, "d", 1.0, 8.7E-09, 0.5, &
    & [ 4.1E-09, 2.0E-09, 1.1E-09, 7.0E-10], 0.5, 5.5E-10 )
  type(dose_factor_t), parameter :: Tc_98   = dose_factor_t("Tc-98", 4.2E6, "y", 1.0, 2.3E-08, 0.5, &
    & [ 1.2E-08, 6.1E-09, 3.7E-09, 2.5E-09], 0.5, 2.0E-09 )
  type(dose_factor_t), parameter :: Tc_99   = dose_factor_t("Tc-99", 2.13E5, "y", 1.0, 1.0E-08, 0.5, &
    & [ 4.8E-09, 2.3E-09, 1.3E-09, 8.2E-10], 0.5, 6.4E-10 )
  type(dose_factor_t), parameter :: Tc_99m  = dose_factor_t("Tc-99m", 6.02, "h", 1.0, 2.0E-10, 0.5, &
    & [ 1.3E-10, 7.2E-11, 4.3E-11, 2.8E-11], 0.5, 2.2E-11 )
  type(dose_factor_t), parameter :: Tc_101  = dose_factor_t("Tc-101", 14.2, "m", 1.0, 2.4E-10, 0.5, &
    & [ 1.3E-10, 6.1E-11, 3.5E-11, 2.4E-11], 0.5, 1.9E-11 )
  type(dose_factor_t), parameter :: Tc_104  = dose_factor_t("Tc-104", 18.2, "m", 1.0, 1.0E-09, 0.5, &
    & [ 5.3E-10, 2.6E-10, 1.5E-10, 1.0E-10], 0.5, 8.0E-11 )
  type(dose_factor_t), parameter :: Ru_94   = dose_factor_t("Ru-94", 51.8, "m", 0.1, 9.3E-10, 0.05, &
    & [ 5.9E-10, 3.1E-10, 1.9E-10, 1.2E-10], 0.05, 9.4E-11 )
  type(dose_factor_t), parameter :: Ru_97   = dose_factor_t("Ru-97", 2.9, "d", 0.1, 1.2E-09, 0.05, &
    & [ 8.5E-10, 4.7E-10, 3.0E-10, 1.9E-10], 0.05, 1.5E-10 )
  type(dose_factor_t), parameter :: Ru_103  = dose_factor_t("Ru-103", 39.28, "d", 0.1, 7.1E-09, 0.05, &
    & [ 4.6E-09, 2.4E-09, 1.5E-09, 9.2E-10], 0.05, 7.3E-10 )
  type(dose_factor_t), parameter :: Ru_105  = dose_factor_t("Ru-105", 4.44, "h", 0.1, 2.7E-09, 0.05, &
    & [ 1.8E-09, 9.1E-10, 5.5E-10, 3.3E-10], 0.05, 2.6E-10 )
  type(dose_factor_t), parameter :: Ru_106  = dose_factor_t("Ru-106", 368.2, "d", 0.1, 8.4E-08, 0.05, &
    & [ 4.9E-08, 2.5E-08, 1.5E-08, 8.6E-09], 0.05, 7.0E-09 )
  type(dose_factor_t), parameter :: Rh_99   = dose_factor_t("Rh-99", 16, "d", 0.1, 4.2E-09, 0.05, &
    & [ 2.9E-09, 1.6E-09, 1.0E-09, 6.5E-10], 0.05, 5.1E-10 )
  type(dose_factor_t), parameter :: Rh_99m  = dose_factor_t("Rh-99m", 4.7, "h", 0.1, 4.9E-10, 0.05, &
    & [ 3.5E-10, 2.0E-10, 1.3E-10, 8.3E-11], 0.05, 6.6E-11 )
  type(dose_factor_t), parameter :: Rh_100  = dose_factor_t("Rh-100", 20.8, "h", 0.1, 4.9E-09, 0.05, &
    & [ 3.6E-09, 2.0E-09, 1.4E-09, 8.8E-10], 0.05, 7.1E-10 )
  type(dose_factor_t), parameter :: Rh_101  = dose_factor_t("Rh-101", 3.2, "y", 0.1, 4.9E-09, 0.05, &
    & [ 2.8E-09, 1.6E-09, 1.0E-09, 6.7E-10], 0.05, 5.5E-10 )
  type(dose_factor_t), parameter :: Rh_101m = dose_factor_t("Rh-101m", 4.34, "d", 0.1, 1.7E-09, 0.05, &
    & [ 1.2E-09, 6.8E-10, 4.4E-10, 2.8E-10], 0.05, 2.2E-10 )
  type(dose_factor_t), parameter :: Rh_102m = dose_factor_t("Rh-102m", 2.9, "y", 0.1, 1.9E-08, 0.05, &
    & [ 1.0E-08, 6.4E-09, 4.3E-09, 3.0E-09], 0.05, 2.6E-09 )
  type(dose_factor_t), parameter :: Rh_102  = dose_factor_t("Rh-102", 207, "d", 0.1, 1.2E-08, 0.05, &
    & [ 7.4E-09, 3.9E-09, 2.4E-09, 1.4E-09], 0.05, 1.2E-09 )
  type(dose_factor_t), parameter :: Rh_103m = dose_factor_t("Rh-103m", 56.12, "m", 0.1, 4.7E-11, 0.05, &
    & [ 2.7E-11, 1.3E-11, 7.4E-12, 4.8E-12], 0.05, 3.8E-12 )
  type(dose_factor_t), parameter :: Rh_105  = dose_factor_t("Rh-105", 35.36, "h", 0.1, 4.0E-09, 0.05, &
    & [ 2.7E-09, 1.3E-09, 8.0E-10, 4.6E-10], 0.05, 3.7E-10 )
  type(dose_factor_t), parameter :: Rh_106m = dose_factor_t("Rh-106m", 132, "m", 0.1, 1.4E-09, 0.05, &
    & [ 9.7E-10, 5.3E-10, 3.3E-10, 2.0E-10], 0.05, 1.6E-10 )
  type(dose_factor_t), parameter :: Rh_107  = dose_factor_t("Rh-107", 21.7, "m", 0.1, 2.9E-10, 0.05, &
    & [ 1.6E-10, 7.9E-11, 4.5E-11, 3.1E-11], 0.05, 2.4E-11 )
  type(dose_factor_t), parameter :: Pd_100  = dose_factor_t("Pd-100", 3.63, "d", 0.05, 7.4E-09, 0.005, &
    & [ 5.2E-09, 2.9E-09, 1.9E-09, 1.2E-09], 0.005, 9.4E-10 )
  type(dose_factor_t), parameter :: Pd_101  = dose_factor_t("Pd-101", 8.27, "h", 0.05, 8.2E-10, 0.005, &
    & [ 5.7E-10, 3.1E-10, 1.9E-10, 1.2E-10], 0.005, 9.4E-11 )
  type(dose_factor_t), parameter :: Pd_103  = dose_factor_t("Pd-103", 16.96, "d", 0.05, 2.2E-09, 0.005, &
    & [ 1.4E-09, 7.2E-10, 4.3E-10, 2.4E-10], 0.005, 1.9E-10 )
  type(dose_factor_t), parameter :: Pd_107  = dose_factor_t("Pd-107", 6.5E6, "y", 0.05, 4.4E-10, 0.005, &
    & [ 2.8E-10, 1.4E-10, 8.1E-11, 4.6E-11], 0.005, 3.7E-11 )
  type(dose_factor_t), parameter :: Pd_109  = dose_factor_t("Pd-109", 13.427, "h", 0.05, 6.3E-09, 0.005, &
    & [ 4.1E-09, 2.0E-09, 1.2E-09, 6.8E-10], 0.005, 5.5E-10 )
  type(dose_factor_t), parameter :: Ag_102  = dose_factor_t("Ag-102", 12.9, "m", 0.1, 4.2E-10, 0.05, &
    & [ 2.4E-10, 1.2E-10, 7.3E-11, 5.0E-11], 0.05, 4.0E-11 )
  type(dose_factor_t), parameter :: Ag_103  = dose_factor_t("Ag-103", 65.7, "m", 0.1, 4.5E-10, 0.05, &
    & [ 2.7E-10, 1.4E-10, 8.3E-11, 5.5E-11], 0.05, 4.3E-11 )
  type(dose_factor_t), parameter :: Ag_104  = dose_factor_t("Ag-104", 69.2, "m", 0.1, 4.3E-10, 0.05, &
    & [ 2.9E-10, 1.7E-10, 1.1E-10, 7.5E-11], 0.05, 6.0E-11 )
  type(dose_factor_t), parameter :: Ag_104m = dose_factor_t("Ag-104m", 33.5, "m", 0.1, 5.6E-10, 0.05, &
    & [ 3.3E-10, 1.7E-10, 1.0E-10, 6.8E-11], 0.05, 5.4E-11 )
  type(dose_factor_t), parameter :: Ag_105  = dose_factor_t("Ag-105", 41.0, "d", 0.1, 3.9E-09, 0.05, &
    & [ 2.5E-09, 1.4E-09, 9.1E-10, 5.9E-10], 0.05, 4.7E-10 )
  type(dose_factor_t), parameter :: Ag_106  = dose_factor_t("Ag-106", 23.96, "m", 0.1, 3.7E-10, 0.05, &
    & [ 2.1E-10, 1.0E-10, 6.0E-11, 4.1E-11], 0.05, 3.2E-11 )
  type(dose_factor_t), parameter :: Ag_106m = dose_factor_t("Ag-106m", 8.41, "d", 0.1, 9.7E-09, 0.05, &
    & [ 6.9E-09, 4.1E-09, 2.8E-09, 1.8E-09], 0.05, 1.5E-09 )
  type(dose_factor_t), parameter :: Ag_108m = dose_factor_t("Ag-108m", 127, "y", 0.1, 2.1E-08, 0.05, &
    & [ 1.1E-08, 6.5E-09, 4.3E-09, 2.8E-09], 0.05, 2.3E-09 )
  type(dose_factor_t), parameter :: Ag_110m = dose_factor_t("Ag-110m", 249.9, "d", 0.1, 2.4E-08, 0.05, &
    & [ 1.4E-08, 7.8E-09, 5.2E-09, 3.4E-09], 0.05, 2.8E-09 )
  type(dose_factor_t), parameter :: Ag_111  = dose_factor_t("Ag-111", 7.45, "d", 0.1, 1.4E-08, 0.05, &
    & [ 9.3E-09, 4.6E-09, 2.7E-09, 1.6E-09], 0.05, 1.3E-09 )
  type(dose_factor_t), parameter :: Ag_112  = dose_factor_t("Ag-112", 3.12, "h", 0.1, 4.9E-09, 0.05, &
    & [ 3.0E-09, 1.5E-09, 8.9E-10, 5.4E-10], 0.05, 4.3E-10 )
  type(dose_factor_t), parameter :: Ag_115  = dose_factor_t("Ag-115", 20.0, "m", 0.1, 7.2E-10, 0.05, &
    & [ 4.1E-10, 2.0E-10, 1.2E-10, 7.7E-11], 0.05, 6.0E-11 )
  type(dose_factor_t), parameter :: Cd_104  = dose_factor_t("Cd-104", 57.7, "m", 0.1, 4.2E-10, 0.05, &
    & [ 2.9E-10, 1.7E-10, 1.1E-10, 7.2E-11], 0.05, 5.4E-11 )
  type(dose_factor_t), parameter :: Cd_107  = dose_factor_t("Cd-107", 6.49, "h", 0.1, 7.1E-10, 0.05, &
    & [ 4.6E-10, 2.3E-10, 1.3E-10, 7.8E-11], 0.05, 6.2E-11 )
  type(dose_factor_t), parameter :: Cd_109  = dose_factor_t("Cd-109", 464, "d", 0.1, 2.1E-08, 0.05, &
    & [ 9.5E-09, 5.5E-09, 3.5E-09, 2.4E-09], 0.05, 2.0E-09 )
  type(dose_factor_t), parameter :: Cd_113  = dose_factor_t("Cd-113", 9.3E15, "y", 0.1, 1.0E-07, 0.05, &
    & [ 4.8E-08, 3.7E-08, 3.0E-08, 2.6E-08], 0.05, 2.5E-08 )
  type(dose_factor_t), parameter :: Cd_113m = dose_factor_t("Cd-113m", 13.6, "y", 0.1, 1.2E-07, 0.05, &
    & [ 5.6E-08, 3.9E-08, 2.9E-08, 2.4E-08], 0.05, 2.3E-08 )
  type(dose_factor_t), parameter :: Cd_115  = dose_factor_t("Cd-115", 53.46, "h", 0.1, 1.4E-08, 0.05, &
    & [ 9.7E-09, 4.9E-09, 2.9E-09, 1.7E-09], 0.05, 1.4E-09 )
  type(dose_factor_t), parameter :: Cd_115m = dose_factor_t("Cd-115m", 44.6, "d", 0.1, 4.1E-08, 0.05, &
    & [ 1.9E-08, 9.7E-09, 6.9E-09, 4.1E-09], 0.05, 3.3E-09 )
  type(dose_factor_t), parameter :: Cd_117  = dose_factor_t("Cd-117", 2.49, "h", 0.1, 2.9E-09, 0.05, &
    & [ 1.9E-09, 9.5E-10, 5.7E-10, 3.5E-10], 0.05, 2.8E-10 )
  type(dose_factor_t), parameter :: Cd_117m = dose_factor_t("Cd-117m", 3.36, "h", 0.1, 2.6E-09, 0.05, &
    & [ 1.7E-09, 9.0E-10, 5.6E-10, 3.5E-10], 0.05, 2.8E-10 )
  type(dose_factor_t), parameter :: In_109  = dose_factor_t("In-109", 4.2, "h", 0.04, 5.2E-10, 0.02, &
    & [ 3.6E-10, 2.0E-10, 1.3E-10, 8.2E-11], 0.02, 6.6E-11 )
  type(dose_factor_t), parameter :: In_110  = dose_factor_t("In-110", 4.9, "h", 0.04, 1.5E-09, 0.02, &
    & [ 1.1E-09, 6.5E-10, 4.4E-10, 3.0E-10], 0.02, 2.4E-10 )
  type(dose_factor_t), parameter :: In_110m = dose_factor_t("In-110m", 69.1, "m", 0.04, 1.1E-09, 0.02, &
    & [ 6.4E-10, 3.2E-10, 1.9E-10, 1.3E-10], 0.02, 1.0E-10 )
  type(dose_factor_t), parameter :: In_111  = dose_factor_t("In-111", 2.83, "d", 0.04, 2.4E-09, 0.02, &
    & [ 1.7E-09, 9.1E-10, 5.9E-10, 3.7E-10], 0.02, 2.9E-10 )
  type(dose_factor_t), parameter :: In_112  = dose_factor_t("In-112", 14.4, "m", 0.04, 1.2E-10, 0.02, &
    & [ 6.7E-11, 3.3E-11, 1.9E-11, 1.3E-11], 0.02, 1.0E-11 )
  type(dose_factor_t), parameter :: In_113m = dose_factor_t("In-113m", 1.658, "h", 0.04, 3.0E-10, 0.02, &
    & [ 1.8E-10, 9.3E-11, 6.2E-11, 3.6E-11], 0.02, 2.8E-11 )
  type(dose_factor_t), parameter :: In_114m = dose_factor_t("In-114m", 49.51, "d", 0.04, 5.6E-08, 0.02, &
    & [ 3.1E-08, 1.5E-08, 9.0E-09, 5.2E-09], 0.02, 4.1E-09 )
  type(dose_factor_t), parameter :: In_115  = dose_factor_t("In-115", 5.1E15, "y", 0.04, 1.3E-07, 0.02, &
    & [ 6.4E-08, 4.8E-08, 4.3E-08, 3.6E-08], 0.02, 3.2E-08 )
  type(dose_factor_t), parameter :: In_115m = dose_factor_t("In-115m", 4.486, "h", 0.04, 9.6E-10, 0.02, &
    & [ 6.0E-10, 3.0E-10, 1.8E-10, 1.1E-10], 0.02, 8.6E-11 )
  type(dose_factor_t), parameter :: In_116m = dose_factor_t("In-116m", 54.15, "m", 0.04, 5.8E-10, 0.02, &
    & [ 3.6E-10, 1.9E-10, 1.2E-10, 8.0E-11], 0.02, 6.4E-11 )
  type(dose_factor_t), parameter :: In_117  = dose_factor_t("In-117", 43.8, "m", 0.04, 3.3E-10, 0.02, &
    & [ 1.9E-10, 9.7E-11, 5.8E-11, 3.9E-11], 0.02, 3.1E-11 )
  type(dose_factor_t), parameter :: In_117m = dose_factor_t("In-117m", 116.5, "m", 0.04, 1.4E-09, 0.02, &
    & [ 8.6E-10, 4.3E-10, 2.5E-10, 1.6E-10], 0.02, 1.2E-10 )
  type(dose_factor_t), parameter :: In_119m = dose_factor_t("In-119m", 18.0, "m", 0.04, 5.9E-10, 0.02, &
    & [ 3.2E-10, 1.6E-10, 8.8E-11, 6.0E-11], 0.02, 4.7E-11 )
  type(dose_factor_t), parameter :: Sn_110  = dose_factor_t("Sn-110", 4.0, "h", 0.04, 3.5E-09, 0.02, &
    & [ 2.3E-09, 1.2E-09, 7.4E-10, 4.4E-10], 0.02, 3.5E-10 )
  type(dose_factor_t), parameter :: Sn_111  = dose_factor_t("Sn-111", 35.3, "m", 0.04, 2.5E-10, 0.02, &
    & [ 1.5E-10, 7.4E-11, 4.4E-11, 3.0E-11], 0.02, 2.3E-11 )
  type(dose_factor_t), parameter :: Sn_113  = dose_factor_t("Sn-113", 115.1, "d", 0.04, 7.8E-09, 0.02, &
    & [ 5.0E-09, 2.6E-09, 1.6E-09, 9.2E-10], 0.02, 7.3E-10 )
  type(dose_factor_t), parameter :: Sn_117m = dose_factor_t("Sn-117m", 13.61, "d", 0.04, 7.7E-09, 0.02, &
    & [ 5.0E-09, 2.5E-09, 1.5E-09, 8.8E-10], 0.02, 7.1E-10 )
  type(dose_factor_t), parameter :: Sn_119m = dose_factor_t("Sn-119m", 293.0, "d", 0.04, 4.1E-09, 0.02, &
    & [ 2.5E-09, 1.3E-09, 7.5E-10, 4.3E-10], 0.02, 3.4E-10 )
  type(dose_factor_t), parameter :: Sn_121  = dose_factor_t("Sn-121", 27.06, "h", 0.04, 2.6E-09, 0.02, &
    & [ 1.7E-09, 8.4E-10, 5.0E-10, 2.8E-10], 0.02, 2.3E-10 )
  type(dose_factor_t), parameter :: Sn_121m = dose_factor_t("Sn-121m", 55, "y", 0.04, 4.6E-09, 0.02, &
    & [ 2.7E-09, 1.4E-09, 8.2E-10, 4.7E-10], 0.02, 3.8E-10 )
  type(dose_factor_t), parameter :: Sn_123  = dose_factor_t("Sn-123", 129.2, "d", 0.04, 2.5E-08, 0.02, &
    & [ 1.6E-08, 7.8E-09, 4.6E-09, 2.6E-09], 0.02, 2.1E-09 )
  type(dose_factor_t), parameter :: Sn_123m = dose_factor_t("Sn-123m", 40.08, "m", 0.04, 4.7E-10, 0.02, &
    & [ 2.6E-10, 1.3E-10, 7.3E-11, 4.9E-11], 0.02, 3.8E-11 )
  type(dose_factor_t), parameter :: Sn_125  = dose_factor_t("Sn-125", 9.64, "d", 0.04, 3.5E-08, 0.02, &
    & [ 2.2E-08, 1.1E-08, 6.7E-09, 3.8E-09], 0.02, 3.1E-09 )
  type(dose_factor_t), parameter :: Sn_126  = dose_factor_t("Sn-126", 1.0E5, "y", 0.04, 5.0E-08, 0.02, &
    & [ 3.0E-08, 1.6E-08, 9.8E-09, 5.9E-09], 0.02, 4.7E-09 )
  type(dose_factor_t), parameter :: Sn_127  = dose_factor_t("Sn-127", 2.10, "h", 0.04, 2.0E-09, 0.02, &
    & [ 1.3E-09, 6.6E-10, 4.0E-10, 2.5E-10], 0.02, 2.0E-10 )
  type(dose_factor_t), parameter :: Sn_128  = dose_factor_t("Sn-128", 59.1, "m", 0.04, 1.6E-09, 0.02, &
    & [ 9.7E-10, 4.9E-10, 3.0E-10, 1.9E-10], 0.02, 1.5E-10 )
  type(dose_factor_t), parameter :: Sb_115  = dose_factor_t("Sb-115", 31.8, "m", 0.2, 2.5E-10, 0.1, &
    & [ 1.5E-10, 7.5E-11, 4.5E-11, 3.1E-11], 0.1, 2.4E-11 )
  type(dose_factor_t), parameter :: Sb_116  = dose_factor_t("Sb-116", 15.8, "m", 0.2, 2.7E-10, 0.1, &
    & [ 1.6E-10, 8.0E-11, 4.8E-11, 3.3E-11], 0.1, 2.6E-11 )
  type(dose_factor_t), parameter :: Sb_116m = dose_factor_t("Sb-116m", 60.3, "m", 0.2, 5.0E-10, 0.1, &
    & [ 3.3E-10, 1.9E-10, 1.2E-10, 8.3E-11], 0.1, 6.7E-11 )
  type(dose_factor_t), parameter :: Sb_117  = dose_factor_t("Sb-117", 2.80, "h", 0.2, 1.6E-10, 0.1, &
    & [ 1.0E-10, 5.6E-11, 3.5E-11, 2.2E-11], 0.1, 1.8E-11 )
  type(dose_factor_t), parameter :: Sb_118m = dose_factor_t("Sb-118m", 5.00, "h", 0.2, 1.3E-09, 0.1, &
    & [ 1.0E-09, 5.8E-10, 3.9E-10, 2.6E-10], 0.1, 2.1E-10 )
  type(dose_factor_t), parameter :: Sb_119  = dose_factor_t("Sb-119", 38.1, "h", 0.2, 8.4E-10, 0.1, &
    & [ 5.8E-10, 3.0E-10, 1.8E-10, 1.0E-10], 0.1, 8.0E-11 )
  type(dose_factor_t), parameter :: Sb_120m = dose_factor_t("Sb-120m", 5.76, "d", 0.2, 8.1E-09, 0.1, &
    & [ 6.0E-09, 3.5E-09, 2.3E-09, 1.6E-09], 0.1, 1.2E-09 )
  type(dose_factor_t), parameter :: Sb_120  = dose_factor_t("Sb-120", 15.89, "m", 0.2, 1.7E-10, 0.1, &
    & [ 9.4E-11, 4.6E-11, 2.7E-11, 1.8E-11], 0.1, 1.4E-11 )
  type(dose_factor_t), parameter :: Sb_122  = dose_factor_t("Sb-122", 2.70, "d", 0.2, 1.8E-08, 0.1, &
    & [ 1.2E-08, 6.1E-09, 3.7E-09, 2.1E-09], 0.1, 1.7E-09 )
  type(dose_factor_t), parameter :: Sb_124  = dose_factor_t("Sb-124", 60.20, "d", 0.2, 2.5E-08, 0.1, &
    & [ 1.6E-08, 8.4E-09, 5.2E-09, 3.2E-09], 0.1, 2.5E-09 )
  type(dose_factor_t), parameter :: Sb_124n = dose_factor_t("Sb-124n", 20.2, "m", 0.2, 8.5E-11, 0.1, &
    & [ 4.9E-11, 2.5E-11, 1.5E-11, 1.0E-11], 0.1, 8.0E-12 )
  type(dose_factor_t), parameter :: Sb_125  = dose_factor_t("Sb-125", 2.77, "y", 0.2, 1.1E-08, 0.1, &
    & [ 6.1E-09, 3.4E-09, 2.1E-09, 1.4E-09], 0.1, 1.1E-09 )
  type(dose_factor_t), parameter :: Sb_126  = dose_factor_t("Sb-126", 12.4, "d", 0.2, 2.0E-08, 0.1, &
    & [ 1.4E-08, 7.6E-09, 4.9E-09, 3.1E-09], 0.1, 2.4E-09 )
  type(dose_factor_t), parameter :: Sb_126m = dose_factor_t("Sb-126m", 19.0, "m", 0.2, 3.9E-10, 0.1, &
    & [ 2.2E-10, 1.1E-10, 6.6E-11, 4.5E-11], 0.1, 3.6E-11 )
  type(dose_factor_t), parameter :: Sb_127  = dose_factor_t("Sb-127", 3.85, "d", 0.2, 1.7E-08, 0.1, &
    & [ 1.2E-08, 5.9E-09, 3.6E-09, 2.1E-09], 0.1, 1.7E-09 )
  type(dose_factor_t), parameter :: Sb_128  = dose_factor_t("Sb-128", 9.01, "h", 0.2, 6.3E-09, 0.1, &
    & [ 4.5E-09, 2.4E-09, 1.5E-09, 9.5E-10], 0.1, 7.6E-10 )
  type(dose_factor_t), parameter :: Sb_128m = dose_factor_t("Sb-128m", 10.4, "m", 0.2, 3.7E-10, 0.1, &
    & [ 2.1E-10, 1.0E-10, 6.0E-11, 4.1E-11], 0.1, 3.3E-11 )
  type(dose_factor_t), parameter :: Sb_129  = dose_factor_t("Sb-129", 4.32, "h", 0.2, 4.3E-09, 0.1, &
    & [ 2.8E-09, 1.5E-09, 8.8E-10, 5.3E-10], 0.1, 4.2E-10 )
  type(dose_factor_t), parameter :: Sb_130  = dose_factor_t("Sb-130", 40, "m", 0.2, 9.1E-10, 0.1, &
    & [ 5.4E-10, 2.8E-10, 1.7E-10, 1.2E-10], 0.1, 9.1E-11 )
  type(dose_factor_t), parameter :: Sb_131  = dose_factor_t("Sb-131", 23, "m", 0.2, 1.1E-09, 0.1, &
    & [ 7.3E-10, 3.9E-10, 2.1E-10, 1.4E-10], 0.1, 1.0E-10 )
  type(dose_factor_t), parameter :: Te_116  = dose_factor_t("Te-116", 2.49, "h", 0.6, 1.4E-09, 0.3, &
    & [ 1.0E-09, 5.5E-10, 3.4E-10, 2.1E-10], 0.3, 1.7E-10 )
  type(dose_factor_t), parameter :: Te_121  = dose_factor_t("Te-121", 17, "d", 0.6, 3.1E-09, 0.3, &
    & [ 2.0E-09, 1.2E-09, 8.0E-10, 5.4E-10], 0.3, 4.3E-10 )
  type(dose_factor_t), parameter :: Te_121m = dose_factor_t("Te-121m", 154, "d", 0.6, 2.7E-08, 0.3, &
    & [ 1.2E-08, 6.9E-09, 4.2E-09, 2.8E-09], 0.3, 2.3E-09 )
  type(dose_factor_t), parameter :: Te_123  = dose_factor_t("Te-123", 1E13, "y", 0.6, 2.0E-08, 0.3, &
    & [ 9.3E-09, 6.9E-09, 5.4E-09, 4.7E-09], 0.3, 4.4E-09 )
  type(dose_factor_t), parameter :: Te_123m = dose_factor_t("Te-123m", 119.7, "d", 0.6, 1.9E-08, 0.3, &
    & [ 8.8E-09, 4.9E-09, 2.8E-09, 1.7E-09], 0.3, 1.4E-09 )
  type(dose_factor_t), parameter :: Te_125m = dose_factor_t("Te-125m", 58, "d", 0.6, 1.3E-08, 0.3, &
    & [ 6.3E-09, 3.3E-09, 1.9E-09, 1.1E-09], 0.3, 8.7E-10 )
  type(dose_factor_t), parameter :: Te_127  = dose_factor_t("Te-127", 9.35, "h", 0.6, 1.5E-09, 0.3, &
    & [ 1.2E-09, 6.2E-10, 3.6E-10, 2.1E-10], 0.3, 1.7E-10 )
  type(dose_factor_t), parameter :: Te_127m = dose_factor_t("Te-127m", 109, "d", 0.6, 4.1E-08, 0.3, &
    & [ 1.8E-08, 9.5E-09, 5.2E-09, 3.0E-09], 0.3, 2.3E-09 )
  type(dose_factor_t), parameter :: Te_129  = dose_factor_t("Te-129", 69.6, "m", 0.6, 7.5E-10, 0.3, &
    & [ 4.4E-10, 2.1E-10, 1.2E-10, 8.0E-11], 0.3, 6.3E-11 )
  type(dose_factor_t), parameter :: Te_129m = dose_factor_t("Te-129m", 33.6, "d", 0.6, 4.4E-08, 0.3, &
    & [ 2.4E-08, 1.2E-08, 6.6E-09, 3.9E-09], 0.3, 3.0E-09 )
  type(dose_factor_t), parameter :: Te_131  = dose_factor_t("Te-131", 25.0, "m", 0.6, 9.0E-10, 0.3, &
    & [ 6.6E-10, 3.5E-10, 1.9E-10, 1.2E-10], 0.3, 8.7E-11 )
  type(dose_factor_t), parameter :: Te_131m = dose_factor_t("Te-131m", 30, "h", 0.6, 2.0E-08, 0.3, &
    & [ 1.4E-08, 7.8E-09, 4.3E-09, 2.7E-09], 0.3, 1.9E-09 )
  type(dose_factor_t), parameter :: Te_132  = dose_factor_t("Te-132", 78.2, "h", 0.6, 4.8E-08, 0.3, &
    & [ 3.0E-08, 1.6E-08, 8.3E-09, 5.3E-09], 0.3, 3.8E-09 )
  type(dose_factor_t), parameter :: Te_133  = dose_factor_t("Te-133", 12.45, "m", 0.6, 8.4E-10, 0.3, &
    & [ 6.3E-10, 3.3E-10, 1.6E-10, 1.1E-10], 0.3, 7.2E-11 )
  type(dose_factor_t), parameter :: Te_133m = dose_factor_t("Te-133m", 55.4, "m", 0.6, 3.1E-09, 0.3, &
    & [ 2.4E-09, 1.3E-09, 6.3E-10, 4.1E-10], 0.3, 2.8E-10 )
  type(dose_factor_t), parameter :: Te_134  = dose_factor_t("Te-134", 41.8, "m", 0.6, 1.1E-09, 0.3, &
    & [ 7.5E-10, 3.9E-10, 2.2E-10, 1.4E-10], 0.3, 1.1E-10 )
  type(dose_factor_t), parameter :: I_120   = dose_factor_t("I-120", 81.0, "m", 1.0, 3.9E-09, 1.0, &
    & [ 2.8E-09, 1.4E-09, 7.2E-10, 4.8E-10], 1.0, 3.4E-10 )
  type(dose_factor_t), parameter :: I_120m  = dose_factor_t("I-120m", 53, "m", 1.0, 2.3E-09, 1.0, &
    & [ 1.5E-09, 7.8E-10, 4.2E-10, 2.9E-10], 1.0, 2.1E-10 )
  type(dose_factor_t), parameter :: I_121   = dose_factor_t("I-121", 2.12, "h", 1.0, 6.2E-10, 1.0, &
    & [ 5.3E-10, 3.1E-10, 1.7E-10, 1.2E-10], 1.0, 8.2E-11 )
  type(dose_factor_t), parameter :: I_123   = dose_factor_t("I-123", 13.2, "h", 1.0, 2.2E-09, 1.0, &
    & [ 1.9E-09, 1.1E-09, 4.9E-10, 3.3E-10], 1.0, 2.1E-10 )
  type(dose_factor_t), parameter :: I_124   = dose_factor_t("I-124", 4.18, "d", 1.0, 1.2E-07, 1.0, &
    & [ 1.1E-07, 6.3E-08, 3.1E-08, 2.0E-08], 1.0, 1.3E-08 )
  type(dose_factor_t), parameter :: I_125   = dose_factor_t("I-125", 60.14, "d", 1.0, 5.2E-08, 1.0, &
    & [ 5.7E-08, 4.1E-08, 3.1E-08, 2.2E-08], 1.0, 1.5E-08 )
  type(dose_factor_t), parameter :: I_126   = dose_factor_t("I-126", 13.02, "d", 1.0, 2.1E-07, 1.0, &
    & [ 2.1E-07, 1.3E-07, 6.8E-08, 4.5E-08], 1.0, 2.9E-08 )
  type(dose_factor_t), parameter :: I_128   = dose_factor_t("I-128", 24.99, "m", 1.0, 5.7E-10, 1.0, &
    & [ 3.3E-10, 1.6E-10, 8.9E-11, 6.0E-11], 1.0, 4.6E-11 )
  type(dose_factor_t), parameter :: I_129   = dose_factor_t("I-129", 1.57E7, "y", 1.0, 1.8E-07, 1.0, &
    & [ 2.2E-07, 1.7E-07, 1.9E-07, 1.4E-07], 1.0, 1.1E-07 )
  type(dose_factor_t), parameter :: I_130   = dose_factor_t("I-130", 12.36, "h", 1.0, 2.1E-08, 1.0, &
    & [ 1.8E-08, 9.8E-09, 4.6E-09, 3.0E-09], 1.0, 2.0E-09 )
  type(dose_factor_t), parameter :: I_131   = dose_factor_t("I-131", 8.04, "d", 1.0, 1.8E-07, 1.0, &
    & [ 1.8E-07, 1.0E-07, 5.2E-08, 3.4E-08], 1.0, 2.2E-08 )
  type(dose_factor_t), parameter :: I_132   = dose_factor_t("I-132", 2.30, "h", 1.0, 3.0E-09, 1.0, &
    & [ 2.4E-09, 1.3E-09, 6.2E-10, 4.1E-10], 1.0, 2.9E-10 )
  type(dose_factor_t), parameter :: I_132m  = dose_factor_t("I-132m", 83.6, "m", 1.0, 2.4E-09, 1.0, &
    & [ 2.0E-09, 1.1E-09, 5.0E-10, 3.3E-10], 1.0, 2.2E-10 )
  type(dose_factor_t), parameter :: I_133   = dose_factor_t("I-133", 20.8, "h", 1.0, 4.9E-08, 1.0, &
    & [ 4.4E-08, 2.3E-08, 1.0E-08, 6.8E-09], 1.0, 4.3E-09 )
  type(dose_factor_t), parameter :: I_134   = dose_factor_t("I-134", 52.6, "m", 1.0, 1.1E-09, 1.0, &
    & [ 7.5E-10, 3.9E-10, 2.1E-10, 1.4E-10], 1.0, 1.1E-10 )
  type(dose_factor_t), parameter :: I_135   = dose_factor_t("I-135", 6.61, "h", 1.0, 1.0E-08, 1.0, &
    & [ 8.9E-09, 4.7E-09, 2.2E-09, 1.4E-09], 1.0, 9.3E-10 )
  type(dose_factor_t), parameter :: Cs_125  = dose_factor_t("Cs-125", 45, "m", 1.0, 3.9E-10, 1.0, &
    & [ 2.2E-10, 1.1E-10, 6.5E-11, 4.4E-11], 1.0, 3.5E-11 )
  type(dose_factor_t), parameter :: Cs_127  = dose_factor_t("Cs-127", 6.25, "h", 1.0, 1.8E-10, 1.0, &
    & [ 1.2E-10, 6.6E-11, 4.2E-11, 2.9E-11], 1.0, 2.4E-11 )
  type(dose_factor_t), parameter :: Cs_129  = dose_factor_t("Cs-129", 32.06, "h", 1.0, 4.4E-10, 1.0, &
    & [ 3.0E-10, 1.7E-10, 1.1E-10, 7.2E-11], 1.0, 6.0E-11 )
  type(dose_factor_t), parameter :: Cs_130  = dose_factor_t("Cs-130", 29.9, "m", 1.0, 3.3E-10, 1.0, &
    & [ 1.8E-10, 9.0E-11, 5.2E-11, 3.6E-11], 1.0, 2.8E-11 )
  type(dose_factor_t), parameter :: Cs_131  = dose_factor_t("Cs-131", 9.69, "d", 1.0, 4.6E-10, 1.0, &
    & [ 2.9E-10, 1.6E-10, 1.0E-10, 6.9E-11], 1.0, 5.8E-11 )
  type(dose_factor_t), parameter :: Cs_132  = dose_factor_t("Cs-132", 6.475, "d", 1.0, 2.7E-09, 1.0, &
    & [ 1.8E-09, 1.1E-09, 7.7E-10, 5.7E-10], 1.0, 5.0E-10 )
  type(dose_factor_t), parameter :: Cs_134  = dose_factor_t("Cs-134", 2.062, "y", 1.0, 2.6E-08, 1.0, &
    & [ 1.6E-08, 1.3E-08, 1.4E-08, 1.9E-08], 1.0, 1.9E-08 )
  type(dose_factor_t), parameter :: Cs_134m = dose_factor_t("Cs-134m", 2.90, "h", 1.0, 2.1E-10, 1.0, &
    & [ 1.2E-10, 5.9E-11, 3.5E-11, 2.5E-11], 1.0, 2.0E-11 )
  type(dose_factor_t), parameter :: Cs_135  = dose_factor_t("Cs-135", 2.3E6, "y", 1.0, 4.1E-09, 1.0, &
    & [ 2.3E-09, 1.7E-09, 1.7E-09, 2.0E-09], 1.0, 2.0E-09 )
  type(dose_factor_t), parameter :: Cs_135m = dose_factor_t("Cs-135m", 53, "m", 1.0, 1.3E-10, 1.0, &
    & [ 8.6E-11, 4.9E-11, 3.2E-11, 2.3E-11], 1.0, 1.9E-11 )
  type(dose_factor_t), parameter :: Cs_136  = dose_factor_t("Cs-136", 13.1, "d", 1.0, 1.5E-08, 1.0, &
    & [ 9.5E-09, 6.1E-09, 4.4E-09, 3.4E-09], 1.0, 3.0E-09 )
  type(dose_factor_t), parameter :: Cs_137  = dose_factor_t("Cs-137", 30.0, "y", 1.0, 2.1E-08, 1.0, &
    & [ 1.2E-08, 9.6E-09, 1.0E-08, 1.3E-08], 1.0, 1.3E-08 )
  type(dose_factor_t), parameter :: Cs_138  = dose_factor_t("Cs-138", 32.2, "m", 1.0, 1.1E-09, 1.0, &
    & [ 5.9E-10, 2.9E-10, 1.7E-10, 1.2E-10], 1.0, 9.2E-11 )
  type(dose_factor_t), parameter :: Ba_126  = dose_factor_t("Ba-126", 96.5, "m", 0.6, 2.7E-09, 0.3, &
    & [ 1.7E-09, 8.5E-10, 5.0E-10, 3.1E-10], 0.3, 2.6E-10 )
  type(dose_factor_t), parameter :: Ba_128  = dose_factor_t("Ba-128", 2.43, "d", 0.6, 2.0E-08, 0.3, &
    & [ 1.7E-08, 9.0E-09, 5.2E-09, 3.0E-09], 0.3, 2.7E-09 )
  type(dose_factor_t), parameter :: Ba_131  = dose_factor_t("Ba-131", 11.8, "d", 0.6, 4.2E-09, 0.3, &
    & [ 2.6E-09, 1.4E-09, 9.4E-10, 6.2E-10], 0.3, 4.5E-10 )
  type(dose_factor_t), parameter :: Ba_131m = dose_factor_t("Ba-131m", 14.6, "m", 0.6, 5.8E-11, 0.3, &
    & [ 3.2E-11, 1.6E-11, 9.3E-12, 6.3E-12], 0.3, 4.9E-12 )
  type(dose_factor_t), parameter :: Ba_133  = dose_factor_t("Ba-133", 10.74, "y", 0.6, 2.2E-08, 0.3, &
    & [ 6.2E-09, 3.9E-09, 4.6E-09, 7.3E-09], 0.3, 1.5E-09 )
  type(dose_factor_t), parameter :: Ba_133m = dose_factor_t("Ba-133m", 38.9, "h", 0.6, 4.2E-09, 0.3, &
    & [ 3.6E-09, 1.8E-09, 1.1E-09, 5.9E-10], 0.3, 5.4E-10 )
  type(dose_factor_t), parameter :: Ba_135m = dose_factor_t("Ba-135m", 28.7, "h", 0.6, 3.3E-09, 0.3, &
    & [ 2.9E-09, 1.5E-09, 8.5E-10, 4.7E-10], 0.3, 4.3E-10 )
  type(dose_factor_t), parameter :: Ba_139  = dose_factor_t("Ba-139", 82.7, "m", 0.6, 1.4E-09, 0.3, &
    & [ 8.4E-10, 4.1E-10, 2.4E-10, 1.5E-10], 0.3, 1.2E-10 )
  type(dose_factor_t), parameter :: Ba_140  = dose_factor_t("Ba-140", 12.74, "d", 0.6, 3.2E-08, 0.3, &
    & [ 1.8E-08, 9.2E-09, 5.8E-09, 3.7E-09], 0.3, 2.6E-09 )
  type(dose_factor_t), parameter :: Ba_141  = dose_factor_t("Ba-141", 18.27, "m", 0.6, 7.6E-10, 0.3, &
    & [ 4.7E-10, 2.3E-10, 1.3E-10, 8.6E-11], 0.3, 7.0E-11 )
  type(dose_factor_t), parameter :: Ba_142  = dose_factor_t("Ba-142", 10.6, "m", 0.6, 3.6E-10, 0.3, &
    & [ 2.2E-10, 1.1E-10, 6.6E-11, 4.3E-11], 0.3, 3.5E-11 )
  type(dose_factor_t), parameter :: La_131  = dose_factor_t("La-131", 59, "m", 0.005, 3.5E-10, 0.0005, &
    & [ 2.1E-10, 1.1E-10, 6.6E-11, 4.4E-11], 0.0005, 3.5E-11 )
  type(dose_factor_t), parameter :: La_132  = dose_factor_t("La-132", 4.8, "h", 0.005, 3.8E-09, 0.0005, &
    & [ 2.4E-09, 1.3E-09, 7.8E-10, 4.8E-10], 0.0005, 3.9E-10 )
  type(dose_factor_t), parameter :: La_135  = dose_factor_t("La-135", 19.5, "h", 0.005, 2.8E-10, 0.0005, &
    & [ 1.9E-10, 1.0E-10, 6.4E-11, 3.9E-11], 0.0005, 3.0E-11 )
  type(dose_factor_t), parameter :: La_137  = dose_factor_t("La-137", 6E4, "y", 0.005, 1.1E-09, 0.0005, &
    & [ 4.5E-10, 2.5E-10, 1.6E-10, 1.0E-10], 0.0005, 8.1E-11 )
  type(dose_factor_t), parameter :: La_138  = dose_factor_t("La-138", 1.35E11, "y", 0.005, 1.3E-08, 0.0005, &
    & [ 4.6E-09, 2.7E-09, 1.9E-09, 1.3E-09], 0.0005, 1.1E-09 )
  type(dose_factor_t), parameter :: La_140  = dose_factor_t("La-140", 40.272, "h", 0.005, 2.0E-08, 0.0005, &
    & [ 1.3E-08, 6.8E-09, 4.2E-09, 2.5E-09], 0.0005, 2.0E-09 )
  type(dose_factor_t), parameter :: La_141  = dose_factor_t("La-141", 3.93, "h", 0.005, 4.3E-09, 0.0005, &
    & [ 2.6E-09, 1.3E-09, 7.6E-10, 4.5E-10], 0.0005, 3.6E-10 )
  type(dose_factor_t), parameter :: La_142  = dose_factor_t("La-142", 92.5, "m", 0.005, 1.9E-09, 0.0005, &
    & [ 1.1E-09, 5.8E-10, 3.5E-10, 2.3E-10], 0.0005, 1.8E-10 )
  type(dose_factor_t), parameter :: La_143  = dose_factor_t("La-143", 14.23, "m", 0.005, 6.9E-10, 0.0005, &
    & [ 3.9E-10, 1.9E-10, 1.1E-10, 7.1E-11], 0.0005, 5.6E-11 )
  type(dose_factor_t), parameter :: Ce_134  = dose_factor_t("Ce-134", 72.0, "h", 0.005, 2.8E-08, 0.0005, &
    & [ 1.8E-08, 9.1E-09, 5.5E-09, 3.2E-09], 0.0005, 2.5E-09 )
  type(dose_factor_t), parameter :: Ce_135  = dose_factor_t("Ce-135", 17.6, "h", 0.005, 7.0E-09, 0.0005, &
    & [ 4.7E-09, 2.6E-09, 1.6E-09, 1.0E-09], 0.0005, 7.9E-10 )
  type(dose_factor_t), parameter :: Ce_137  = dose_factor_t("Ce-137", 9.0, "h", 0.005, 2.6E-10, 0.0005, &
    & [ 1.7E-10, 8.8E-11, 5.4E-11, 3.2E-11], 0.0005, 2.5E-11 )
  type(dose_factor_t), parameter :: Ce_137m = dose_factor_t("Ce-137m", 34.4, "h", 0.005, 6.1E-09, 0.0005, &
    & [ 3.9E-09, 2.0E-09, 1.2E-09, 6.8E-10], 0.0005, 5.4E-10 )
  type(dose_factor_t), parameter :: Ce_139  = dose_factor_t("Ce-139", 137.66, "d", 0.005, 2.6E-09, 0.0005, &
    & [ 1.6E-09, 8.6E-10, 5.4E-10, 3.3E-10], 0.0005, 2.6E-10 )
  type(dose_factor_t), parameter :: Ce_141  = dose_factor_t("Ce-141", 32.501, "d", 0.005, 8.1E-09, 0.0005, &
    & [ 5.1E-09, 2.6E-09, 1.5E-09, 8.8E-10], 0.0005, 7.1E-10 )
  type(dose_factor_t), parameter :: Ce_143  = dose_factor_t("Ce-143", 33.0, "h", 0.005, 1.2E-08, 0.0005, &
    & [ 8.0E-09, 4.1E-09, 2.4E-09, 1.4E-09], 0.0005, 1.1E-09 )
  type(dose_factor_t), parameter :: Ce_144  = dose_factor_t("Ce-144", 284.3, "d", 0.005, 6.6E-08, 0.0005, &
    & [ 3.9E-08, 1.9E-08, 1.1E-08, 6.5E-09], 0.0005, 5.2E-09 )
  type(dose_factor_t), parameter :: Pr_136  = dose_factor_t("Pr-136", 13.1, "m", 0.005, 3.7E-10, 0.0005, &
    & [ 2.1E-10, 1.0E-10, 6.1E-11, 4.2E-11], 0.0005, 3.3E-11 )
  type(dose_factor_t), parameter :: Pr_137  = dose_factor_t("Pr-137", 76.6, "m", 0.005, 4.1E-10, 0.0005, &
    & [ 2.5E-10, 1.3E-10, 7.7E-11, 5.0E-11], 0.0005, 4.0E-11 )
  type(dose_factor_t), parameter :: Pr_138m = dose_factor_t("Pr-138m", 2.1, "h", 0.005, 1.0E-09, 0.0005, &
    & [ 7.4E-10, 4.1E-10, 2.6E-10, 1.6E-10], 0.0005, 1.3E-10 )
  type(dose_factor_t), parameter :: Pr_139  = dose_factor_t("Pr-139", 4.51, "h", 0.005, 3.2E-10, 0.0005, &
    & [ 2.0E-10, 1.1E-10, 6.5E-11, 4.0E-11], 0.0005, 3.1E-11 )
  type(dose_factor_t), parameter :: Pr_142  = dose_factor_t("Pr-142", 19.13, "h", 0.005, 1.5E-08, 0.0005, &
    & [ 9.8E-09, 4.9E-09, 2.9E-09, 1.6E-09], 0.0005, 1.3E-09 )
  type(dose_factor_t), parameter :: Pr_142m = dose_factor_t("Pr-142m", 14.6, "m", 0.005, 2.0E-10, 0.0005, &
    & [ 1.2E-10, 6.2E-11, 3.7E-11, 2.1E-11], 0.0005, 1.7E-11 )
  type(dose_factor_t), parameter :: Pr_143  = dose_factor_t("Pr-143", 13.56, "d", 0.005, 1.4E-08, 0.0005, &
    & [ 8.7E-09, 4.3E-09, 2.6E-09, 1.5E-09], 0.0005, 1.2E-09 )
  type(dose_factor_t), parameter :: Pr_144  = dose_factor_t("Pr-144", 17.28, "m", 0.005, 6.4E-10, 0.0005, &
    & [ 3.5E-10, 1.7E-10, 9.5E-11, 6.5E-11], 0.0005, 5.0E-11 )
  type(dose_factor_t), parameter :: Pr_145  = dose_factor_t("Pr-145", 5.98, "h", 0.005, 4.7E-09, 0.0005, &
    & [ 2.9E-09, 1.4E-09, 8.5E-10, 4.9E-10], 0.0005, 3.9E-10 )
  type(dose_factor_t), parameter :: Pr_147  = dose_factor_t("Pr-147", 13.6, "m", 0.005, 3.9E-10, 0.0005, &
    & [ 2.2E-10, 1.1E-10, 6.1E-11, 4.2E-11], 0.0005, 3.3E-11 )
  type(dose_factor_t), parameter :: Nd_136  = dose_factor_t("Nd-136", 50.65, "m", 0.005, 1.0E-09, 0.0005, &
    & [ 6.1E-10, 3.1E-10, 1.9E-10, 1.2E-10], 0.0005, 9.9E-11 )
  type(dose_factor_t), parameter :: Nd_138  = dose_factor_t("Nd-138", 5.04, "h", 0.005, 7.2E-09, 0.0005, &
    & [ 4.5E-09, 2.3E-09, 1.3E-09, 8.0E-10], 0.0005, 6.4E-10 )
  type(dose_factor_t), parameter :: Nd_139  = dose_factor_t("Nd-139", 29.7, "m", 0.005, 2.1E-10, 0.0005, &
    & [ 1.2E-10, 6.3E-11, 3.7E-11, 2.5E-11], 0.0005, 2.0E-11 )
  type(dose_factor_t), parameter :: Nd_139m = dose_factor_t("Nd-139m", 5.5, "h", 0.005, 2.1E-09, 0.0005, &
    & [ 1.4E-09, 7.8E-10, 5.0E-10, 3.1E-10], 0.0005, 2.5E-10 )
  type(dose_factor_t), parameter :: Nd_141  = dose_factor_t("Nd-141", 2.49, "h", 0.005, 7.8E-11, 0.0005, &
    & [ 5.0E-11, 2.7E-11, 1.6E-11, 1.0E-11], 0.0005, 8.3E-12 )
  type(dose_factor_t), parameter :: Nd_147  = dose_factor_t("Nd-147", 10.98, "d", 0.005, 1.2E-08, 0.0005, &
    & [ 7.8E-09, 3.9E-09, 2.3E-09, 1.3E-09], 0.0005, 1.1E-09 )
  type(dose_factor_t), parameter :: Nd_149  = dose_factor_t("Nd-149", 1.73, "h", 0.005, 1.4E-09, 0.0005, &
    & [ 8.7E-10, 4.3E-10, 2.6E-10, 1.6E-10], 0.0005, 1.2E-10 )
  type(dose_factor_t), parameter :: Nd_151  = dose_factor_t("Nd-151", 12.44, "m", 0.005, 3.4E-10, 0.0005, &
    & [ 2.0E-10, 9.7E-11, 5.7E-11, 3.8E-11], 0.0005, 3.0E-11 )
  type(dose_factor_t), parameter :: Pm_141  = dose_factor_t("Pm-141", 20.90, "m", 0.005, 4.2E-10, 0.0005, &
    & [ 2.4E-10, 1.2E-10, 6.8E-11, 4.6E-11], 0.0005, 3.6E-11 )
  type(dose_factor_t), parameter :: Pm_143  = dose_factor_t("Pm-143", 265, "d", 0.005, 1.9E-09, 0.0005, &
    & [ 1.2E-09, 6.7E-10, 4.4E-10, 2.9E-10], 0.0005, 2.3E-10 )
  type(dose_factor_t), parameter :: Pm_144  = dose_factor_t("Pm-144", 363, "d", 0.005, 7.6E-09, 0.0005, &
    & [ 4.7E-09, 2.7E-09, 1.8E-09, 1.2E-09], 0.0005, 9.7E-10 )
  type(dose_factor_t), parameter :: Pm_145  = dose_factor_t("Pm-145", 17.7, "y", 0.005, 1.5E-09, 0.0005, &
    & [ 6.8E-10, 3.7E-10, 2.3E-10, 1.4E-10], 0.0005, 1.1E-10 )
  type(dose_factor_t), parameter :: Pm_146  = dose_factor_t("Pm-146", 2020, "d", 0.005, 1.0E-08, 0.0005, &
    & [ 5.1E-09, 2.8E-09, 1.8E-09, 1.1E-09], 0.0005, 9.0E-10 )
  type(dose_factor_t), parameter :: Pm_147  = dose_factor_t("Pm-147", 2.6234, "y", 0.005, 3.6E-09, 0.0005, &
    & [ 1.9E-09, 9.6E-10, 5.7E-10, 3.2E-10], 0.0005, 2.6E-10 )
  type(dose_factor_t), parameter :: Pm_148  = dose_factor_t("Pm-148", 5.37, "d", 0.005, 3.0E-08, 0.0005, &
    & [ 1.9E-08, 9.7E-09, 5.8E-09, 3.3E-09], 0.0005, 2.7E-09 )
  type(dose_factor_t), parameter :: Pm_148m = dose_factor_t("Pm-148m", 41.3, "d", 0.005, 1.5E-08, 0.0005, &
    & [ 1.0E-08, 5.5E-09, 3.5E-09, 2.2E-09], 0.0005, 1.7E-09 )
  type(dose_factor_t), parameter :: Pm_149  = dose_factor_t("Pm-149", 53.08, "h", 0.005, 1.2E-08, 0.0005, &
    & [ 7.4E-09, 3.7E-09, 2.2E-09, 1.2E-09], 0.0005, 9.9E-10 )
  type(dose_factor_t), parameter :: Pm_150  = dose_factor_t("Pm-150", 2.68, "h", 0.005, 2.8E-09, 0.0005, &
    & [ 1.7E-09, 8.7E-10, 5.2E-10, 3.2E-10], 0.0005, 2.6E-10 )
  type(dose_factor_t), parameter :: Pm_151  = dose_factor_t("Pm-151", 28.40, "h", 0.005, 8.0E-09, 0.0005, &
    & [ 5.1E-09, 2.6E-09, 1.6E-09, 9.1E-10], 0.0005, 7.3E-10 )
  type(dose_factor_t), parameter :: Sm_141  = dose_factor_t("Sm-141", 10.2, "m", 0.005, 4.5E-10, 0.0005, &
    & [ 2.5E-10, 1.3E-10, 7.3E-11, 5.0E-11], 0.0005, 3.9E-11 )
  type(dose_factor_t), parameter :: Sm_141m = dose_factor_t("Sm-141m", 22.6, "m", 0.005, 7.0E-10, 0.0005, &
    & [ 4.0E-10, 2.0E-10, 1.2E-10, 8.2E-11], 0.0005, 6.5E-11 )
  type(dose_factor_t), parameter :: Sm_142  = dose_factor_t("Sm-142", 72.49, "m", 0.005, 2.2E-09, 0.0005, &
    & [ 1.3E-09, 6.2E-10, 3.6E-10, 2.4E-10], 0.0005, 1.9E-10 )
  type(dose_factor_t), parameter :: Sm_145  = dose_factor_t("Sm-145", 340, "d", 0.005, 2.4E-09, 0.0005, &
    & [ 1.4E-09, 7.3E-10, 4.5E-10, 2.7E-10], 0.0005, 2.1E-10 )
  type(dose_factor_t), parameter :: Sm_146  = dose_factor_t("Sm-146", 1.03E8, "y", 0.005, 1.5E-06, 0.0005, &
    & [ 1.5E-07, 1.0E-07, 7.0E-08, 5.8E-08], 0.0005, 5.4E-08 )
  type(dose_factor_t), parameter :: Sm_147  = dose_factor_t("Sm-147", 1.06E11, "y", 0.005, 1.4E-06, 0.0005, &
    & [ 1.4E-07, 9.2E-08, 6.4E-08, 5.2E-08], 0.0005, 4.9E-08 )
  type(dose_factor_t), parameter :: Sm_151  = dose_factor_t("Sm-151", 90, "y", 0.005, 1.5E-09, 0.0005, &
    & [ 6.4E-10, 3.3E-10, 2.0E-10, 1.2E-10], 0.0005, 9.8E-11 )
  type(dose_factor_t), parameter :: Sm_153  = dose_factor_t("Sm-153", 46.7, "h", 0.005, 8.4E-09, 0.0005, &
    & [ 5.4E-09, 2.7E-09, 1.6E-09, 9.2E-10], 0.0005, 7.4E-10 )
  type(dose_factor_t), parameter :: Sm_155  = dose_factor_t("Sm-155", 22.1, "m", 0.005, 3.6E-10, 0.0005, &
    & [ 2.0E-10, 9.7E-11, 5.5E-11, 3.7E-11], 0.0005, 2.9E-11 )
  type(dose_factor_t), parameter :: Sm_156  = dose_factor_t("Sm-156", 9.4, "h", 0.005, 2.8E-09, 0.0005, &
    & [ 1.8E-09, 9.0E-10, 5.4E-10, 3.1E-10], 0.0005, 2.5E-10 )
  type(dose_factor_t), parameter :: Eu_145  = dose_factor_t("Eu-145", 5.94, "d", 0.005, 5.1E-09, 0.0005, &
    & [ 3.7E-09, 2.1E-09, 1.4E-09, 9.4E-10], 0.0005, 7.5E-10 )
  type(dose_factor_t), parameter :: Eu_146  = dose_factor_t("Eu-146", 4.61, "d", 0.005, 8.5E-09, 0.0005, &
    & [ 6.2E-09, 3.6E-09, 2.4E-09, 1.6E-09], 0.0005, 1.3E-09 )
  type(dose_factor_t), parameter :: Eu_147  = dose_factor_t("Eu-147", 24, "d", 0.005, 3.7E-09, 0.0005, &
    & [ 2.5E-09, 1.4E-09, 8.9E-10, 5.6E-10], 0.0005, 4.4E-10 )
  type(dose_factor_t), parameter :: Eu_148  = dose_factor_t("Eu-148", 54.5, "d", 0.005, 8.5E-09, 0.0005, &
    & [ 6.0E-09, 3.5E-09, 2.4E-09, 1.6E-09], 0.0005, 1.3E-09 )
  type(dose_factor_t), parameter :: Eu_149  = dose_factor_t("Eu-149", 93.1, "d", 0.005, 9.7E-10, 0.0005, &
    & [ 6.3E-10, 3.4E-10, 2.1E-10, 1.3E-10], 0.0005, 1.0E-10 )
  type(dose_factor_t), parameter :: Eu_150  = dose_factor_t("Eu-150", 34.2, "y", 0.005, 1.3E-08, 0.0005, &
    & [ 5.7E-09, 3.4E-09, 2.3E-09, 1.5E-09], 0.0005, 1.3E-09 )
  type(dose_factor_t), parameter :: Eu_150m = dose_factor_t("Eu-150m", 12.62, "h", 0.005, 4.4E-09, 0.0005, &
    & [ 2.8E-09, 1.4E-09, 8.2E-10, 4.7E-10], 0.0005, 3.8E-10 )
  type(dose_factor_t), parameter :: Eu_152  = dose_factor_t("Eu-152", 13.33, "y", 0.005, 1.6E-08, 0.0005, &
    & [ 7.4E-09, 4.1E-09, 2.6E-09, 1.7E-09], 0.0005, 1.4E-09 )
  type(dose_factor_t), parameter :: Eu_152m = dose_factor_t("Eu-152m", 9.32, "h", 0.005, 5.7E-09, 0.0005, &
    & [ 3.6E-09, 1.8E-09, 1.1E-09, 6.2E-10], 0.0005, 5.0E-10 )
  type(dose_factor_t), parameter :: Eu_154  = dose_factor_t("Eu-154", 8.8, "y", 0.005, 2.5E-08, 0.0005, &
    & [ 1.2E-08, 6.5E-09, 4.1E-09, 2.5E-09], 0.0005, 2.0E-09 )
  type(dose_factor_t), parameter :: Eu_155  = dose_factor_t("Eu-155", 4.96, "y", 0.005, 4.3E-09, 0.0005, &
    & [ 2.2E-09, 1.1E-09, 6.8E-10, 4.0E-10], 0.0005, 3.2E-10 )
  type(dose_factor_t), parameter :: Eu_156  = dose_factor_t("Eu-156", 15.19, "d", 0.005, 2.2E-08, 0.0005, &
    & [ 1.5E-08, 7.5E-09, 4.6E-09, 2.7E-09], 0.0005, 2.2E-09 )
  type(dose_factor_t), parameter :: Eu_157  = dose_factor_t("Eu-157", 15.15, "h", 0.005, 6.7E-09, 0.0005, &
    & [ 4.3E-09, 2.2E-09, 1.3E-09, 7.5E-10], 0.0005, 6.0E-10 )
  type(dose_factor_t), parameter :: Eu_158  = dose_factor_t("Eu-158", 45.9, "m", 0.005, 1.1E-09, 0.0005, &
    & [ 6.2E-10, 3.1E-10, 1.8E-10, 1.2E-10], 0.0005, 9.4E-11 )
  type(dose_factor_t), parameter :: Gd_145  = dose_factor_t("Gd-145", 22.9, "m", 0.005, 4.5E-10, 0.0005, &
    & [ 2.6E-10, 1.3E-10, 8.1E-11, 5.6E-11], 0.0005, 4.4E-11 )
  type(dose_factor_t), parameter :: Gd_146  = dose_factor_t("Gd-146", 48.3, "d", 0.005, 9.4E-09, 0.0005, &
    & [ 6.0E-09, 3.2E-09, 2.0E-09, 1.2E-09], 0.0005, 9.6E-10 )
  type(dose_factor_t), parameter :: Gd_147  = dose_factor_t("Gd-147", 38.1, "h", 0.005, 4.5E-09, 0.0005, &
    & [ 3.2E-09, 1.8E-09, 1.2E-09, 7.7E-10], 0.0005, 6.1E-10 )
  type(dose_factor_t), parameter :: Gd_148  = dose_factor_t("Gd-148", 93, "y", 0.005, 1.7E-06, 0.0005, &
    & [ 1.6E-07, 1.1E-07, 7.3E-08, 5.9E-08], 0.0005, 5.6E-08 )
  type(dose_factor_t), parameter :: Gd_149  = dose_factor_t("Gd-149", 9.4, "d", 0.005, 4.0E-09, 0.0005, &
    & [ 2.7E-09, 1.5E-09, 9.3E-10, 5.7E-10], 0.0005, 4.5E-10 )
  type(dose_factor_t), parameter :: Gd_151  = dose_factor_t("Gd-151", 120, "d", 0.005, 2.1E-09, 0.0005, &
    & [ 1.3E-09, 6.8E-10, 4.2E-10, 2.4E-10], 0.0005, 2.0E-10 )
  type(dose_factor_t), parameter :: Gd_152  = dose_factor_t("Gd-152", 1.08E14, "y", 0.005, 1.2E-06, 0.0005, &
    & [ 1.2E-07, 7.7E-08, 5.3E-08, 4.3E-08], 0.0005, 4.1E-08 )
  type(dose_factor_t), parameter :: Gd_153  = dose_factor_t("Gd-153", 242, "d", 0.005, 2.9E-09, 0.0005, &
    & [ 1.8E-09, 9.4E-10, 5.8E-10, 3.4E-10], 0.0005, 2.7E-10 )
  type(dose_factor_t), parameter :: Gd_159  = dose_factor_t("Gd-159", 18.56, "h", 0.005, 5.7E-09, 0.0005, &
    & [ 3.6E-09, 1.8E-09, 1.1E-09, 6.2E-10], 0.0005, 4.9E-10 )
  type(dose_factor_t), parameter :: Tb_147  = dose_factor_t("Tb-147", 1.65, "h", 0.005, 1.5E-09, 0.0005, &
    & [ 1.0E-09, 5.4E-10, 3.3E-10, 2.0E-10], 0.0005, 1.6E-10 )
  type(dose_factor_t), parameter :: Tb_149  = dose_factor_t("Tb-149", 4.15, "h", 0.005, 2.4E-09, 0.0005, &
    & [ 1.5E-09, 8.0E-10, 5.0E-10, 3.1E-10], 0.0005, 2.5E-10 )
  type(dose_factor_t), parameter :: Tb_150  = dose_factor_t("Tb-150", 3.27, "h", 0.005, 2.5E-09, 0.0005, &
    & [ 1.6E-09, 8.3E-10, 5.1E-10, 3.2E-10], 0.0005, 2.5E-10 )
  type(dose_factor_t), parameter :: Tb_151  = dose_factor_t("Tb-151", 17.6, "h", 0.005, 2.7E-09, 0.0005, &
    & [ 1.9E-09, 1.0E-09, 6.7E-10, 4.2E-10], 0.0005, 3.4E-10 )
  type(dose_factor_t), parameter :: Tb_153  = dose_factor_t("Tb-153", 2.34, "d", 0.005, 2.3E-09, 0.0005, &
    & [ 1.5E-09, 8.2E-10, 5.1E-10, 3.1E-10], 0.0005, 2.5E-10 )
  type(dose_factor_t), parameter :: Tb_154  = dose_factor_t("Tb-154", 21.4, "h", 0.005, 4.7E-09, 0.0005, &
    & [ 3.4E-09, 1.9E-09, 1.3E-09, 8.1E-10], 0.0005, 6.5E-10 )
  type(dose_factor_t), parameter :: Tb_155  = dose_factor_t("Tb-155", 5.32, "d", 0.005, 1.9E-09, 0.0005, &
    & [ 1.3E-09, 6.8E-10, 4.3E-10, 2.6E-10], 0.0005, 2.1E-10 )
  type(dose_factor_t), parameter :: Tb_156  = dose_factor_t("Tb-156", 5.34, "d", 0.005, 9.0E-09, 0.0005, &
    & [ 6.3E-09, 3.5E-09, 2.3E-09, 1.5E-09], 0.0005, 1.2E-09 )
  type(dose_factor_t), parameter :: Tb_156n = dose_factor_t("Tb-156n", 5.0, "h", 0.005, 8.0E-10, 0.0005, &
    & [ 5.2E-10, 2.7E-10, 1.7E-10, 1.0E-10], 0.0005, 8.1E-11 )
  type(dose_factor_t), parameter :: Tb_156m = dose_factor_t("Tb-156m", 24.4, "h", 0.005, 1.5E-09, 0.0005, &
    & [ 1.0E-09, 5.6E-10, 3.5E-10, 2.2E-10], 0.0005, 1.7E-10 )
  type(dose_factor_t), parameter :: Tb_157  = dose_factor_t("Tb-157", 150, "y", 0.005, 4.9E-10, 0.0005, &
    & [ 2.2E-10, 1.1E-10, 6.8E-11, 4.1E-11], 0.0005, 3.4E-11 )
  type(dose_factor_t), parameter :: Tb_158  = dose_factor_t("Tb-158", 150, "y", 0.005, 1.3E-08, 0.0005, &
    & [ 5.9E-09, 3.3E-09, 2.1E-09, 1.4E-09], 0.0005, 1.1E-09 )
  type(dose_factor_t), parameter :: Tb_160  = dose_factor_t("Tb-160", 72.3, "d", 0.005, 1.6E-08, 0.0005, &
    & [ 1.0E-08, 5.4E-09, 3.3E-09, 2.0E-09], 0.0005, 1.6E-09 )
  type(dose_factor_t), parameter :: Tb_161  = dose_factor_t("Tb-161", 6.91, "d", 0.005, 8.3E-09, 0.0005, &
    & [ 5.3E-09, 2.7E-09, 1.6E-09, 9.0E-10], 0.0005, 7.2E-10 )
  type(dose_factor_t), parameter :: Dy_155  = dose_factor_t("Dy-155", 10.0, "h", 0.005, 9.7E-10, 0.0005, &
    & [ 6.8E-10, 3.8E-10, 2.5E-10, 1.6E-10], 0.0005, 1.3E-10 )
  type(dose_factor_t), parameter :: Dy_157  = dose_factor_t("Dy-157", 8.1, "h", 0.005, 4.4E-10, 0.0005, &
    & [ 3.1E-10, 1.8E-10, 1.2E-10, 7.7E-11], 0.0005, 6.1E-11 )
  type(dose_factor_t), parameter :: Dy_159  = dose_factor_t("Dy-159", 144.4, "d", 0.005, 1.0E-09, 0.0005, &
    & [ 6.4E-10, 3.4E-10, 2.1E-10, 1.3E-10], 0.0005, 1.0E-10 )
  type(dose_factor_t), parameter :: Dy_165  = dose_factor_t("Dy-165", 2.334, "h", 0.005, 1.3E-09, 0.0005, &
    & [ 7.9E-10, 3.9E-10, 2.3E-10, 1.4E-10], 0.0005, 1.1E-10 )
  type(dose_factor_t), parameter :: Dy_166  = dose_factor_t("Dy-166", 81.6, "h", 0.005, 1.9E-08, 0.0005, &
    & [ 1.2E-08, 6.0E-09, 3.6E-09, 2.0E-09], 0.0005, 1.6E-09 )
  type(dose_factor_t), parameter :: Ho_155  = dose_factor_t("Ho-155", 48, "m", 0.005, 3.8E-10, 0.0005, &
    & [ 2.3E-10, 1.2E-10, 7.1E-11, 4.7E-11], 0.0005, 3.7E-11 )
  type(dose_factor_t), parameter :: Ho_157  = dose_factor_t("Ho-157", 12.6, "m", 0.005, 5.8E-11, 0.0005, &
    & [ 3.6E-11, 1.9E-11, 1.2E-11, 8.1E-12], 0.0005, 6.5E-12 )
  type(dose_factor_t), parameter :: Ho_159  = dose_factor_t("Ho-159", 33, "m", 0.005, 7.1E-11, 0.0005, &
    & [ 4.3E-11, 2.3E-11, 1.4E-11, 9.9E-12], 0.0005, 7.9E-12 )
  type(dose_factor_t), parameter :: Ho_161  = dose_factor_t("Ho-161", 2.5, "h", 0.005, 1.4E-10, 0.0005, &
    & [ 8.1E-11, 4.2E-11, 2.5E-11, 1.6E-11], 0.0005, 1.3E-11 )
  type(dose_factor_t), parameter :: Ho_162  = dose_factor_t("Ho-162", 15, "m", 0.005, 3.5E-11, 0.0005, &
    & [ 2.0E-11, 1.0E-11, 6.0E-12, 4.2E-12], 0.0005, 3.3E-12 )
  type(dose_factor_t), parameter :: Ho_162m = dose_factor_t("Ho-162m", 68, "m", 0.005, 2.4E-10, 0.0005, &
    & [ 1.5E-10, 7.9E-11, 4.9E-11, 3.3E-11], 0.0005, 2.6E-11 )
  type(dose_factor_t), parameter :: Ho_164  = dose_factor_t("Ho-164", 29, "m", 0.005, 1.2E-10, 0.0005, &
    & [ 6.5E-11, 3.2E-11, 1.8E-11, 1.2E-11], 0.0005, 9.5E-12 )
  type(dose_factor_t), parameter :: Ho_164m = dose_factor_t("Ho-164m", 37.5, "m", 0.005, 2.0E-10, 0.0005, &
    & [ 1.1E-10, 5.5E-11, 3.2E-11, 2.1E-11], 0.0005, 1.6E-11 )
  type(dose_factor_t), parameter :: Ho_166  = dose_factor_t("Ho-166", 26.80, "h", 0.005, 1.6E-08, 0.0005, &
    & [ 1.0E-08, 5.2E-09, 3.1E-09, 1.7E-09], 0.0005, 1.4E-09 )
  type(dose_factor_t), parameter :: Ho_166m = dose_factor_t("Ho-166m", 1.20E3, "y", 0.005, 2.6E-08, 0.0005, &
    & [ 9.3E-09, 5.3E-09, 3.5E-09, 2.4E-09], 0.0005, 2.0E-09 )
  type(dose_factor_t), parameter :: Ho_167  = dose_factor_t("Ho-167", 3.1, "h", 0.005, 8.8E-10, 0.0005, &
    & [ 5.5E-10, 2.8E-10, 1.7E-10, 1.0E-10], 0.0005, 8.3E-11 )
  type(dose_factor_t), parameter :: Er_161  = dose_factor_t("Er-161", 3.24, "h", 0.005, 6.5E-10, 0.0005, &
    & [ 4.4E-10, 2.4E-10, 1.6E-10, 1.0E-10], 0.0005, 8.0E-11 )
  type(dose_factor_t), parameter :: Er_165  = dose_factor_t("Er-165", 10.36, "h", 0.005, 1.7E-10, 0.0005, &
    & [ 1.1E-10, 6.2E-11, 3.9E-11, 2.4E-11], 0.0005, 1.9E-11 )
  type(dose_factor_t), parameter :: Er_169  = dose_factor_t("Er-169", 9.3, "d", 0.005, 4.4E-09, 0.0005, &
    & [ 2.8E-09, 1.4E-09, 8.2E-10, 4.7E-10], 0.0005, 3.7E-10 )
  type(dose_factor_t), parameter :: Er_171  = dose_factor_t("Er-171", 7.52, "h", 0.005, 4.0E-09, 0.0005, &
    & [ 2.5E-09, 1.3E-09, 7.6E-10, 4.5E-10], 0.0005, 3.6E-10 )
  type(dose_factor_t), parameter :: Er_172  = dose_factor_t("Er-172", 49.3, "h", 0.005, 1.0E-08, 0.0005, &
    & [ 6.8E-09, 3.5E-09, 2.1E-09, 1.3E-09], 0.0005, 1.0E-09 )
  type(dose_factor_t), parameter :: Tm_162  = dose_factor_t("Tm-162", 21.7, "m", 0.005, 2.9E-10, 0.0005, &
    & [ 1.7E-10, 8.7E-11, 5.2E-11, 3.6E-11], 0.0005, 2.9E-11 )
  type(dose_factor_t), parameter :: Tm_166  = dose_factor_t("Tm-166", 7.70, "h", 0.005, 2.1E-09, 0.0005, &
    & [ 1.5E-09, 8.3E-10, 5.5E-10, 3.5E-10], 0.0005, 2.8E-10 )
  type(dose_factor_t), parameter :: Tm_167  = dose_factor_t("Tm-167", 9.24, "d", 0.005, 6.0E-09, 0.0005, &
    & [ 3.9E-09, 2.0E-09, 1.2E-09, 7.0E-10], 0.0005, 5.6E-10 )
  type(dose_factor_t), parameter :: Tm_170  = dose_factor_t("Tm-170", 128.6, "d", 0.005, 1.6E-08, 0.0005, &
    & [ 9.8E-09, 4.9E-09, 2.9E-09, 1.6E-09], 0.0005, 1.3E-09 )
  type(dose_factor_t), parameter :: Tm_171  = dose_factor_t("Tm-171", 1.92, "y", 0.005, 1.5E-09, 0.0005, &
    & [ 7.8E-10, 3.9E-10, 2.3E-10, 1.3E-10], 0.0005, 1.1E-10 )
  type(dose_factor_t), parameter :: Tm_172  = dose_factor_t("Tm-172", 63.6, "h", 0.005, 1.9E-08, 0.0005, &
    & [ 1.2E-08, 6.1E-09, 3.7E-09, 2.1E-09], 0.0005, 1.7E-09 )
  type(dose_factor_t), parameter :: Tm_173  = dose_factor_t("Tm-173", 8.24, "h", 0.005, 3.3E-09, 0.0005, &
    & [ 2.1E-09, 1.1E-09, 6.5E-10, 3.8E-10], 0.0005, 3.1E-10 )
  type(dose_factor_t), parameter :: Tm_175  = dose_factor_t("Tm-175", 15.2, "m", 0.005, 3.1E-10, 0.0005, &
    & [ 1.7E-10, 8.6E-11, 5.0E-11, 3.4E-11], 0.0005, 2.7E-11 )
  type(dose_factor_t), parameter :: Yb_162  = dose_factor_t("Yb-162", 18.9, "m", 0.005, 2.2E-10, 0.0005, &
    & [ 1.3E-10, 6.9E-11, 4.2E-11, 2.9E-11], 0.0005, 2.3E-11 )
  type(dose_factor_t), parameter :: Yb_166  = dose_factor_t("Yb-166", 56.7, "h", 0.005, 7.7E-09, 0.0005, &
    & [ 5.4E-09, 2.9E-09, 1.9E-09, 1.2E-09], 0.0005, 9.5E-10 )
  type(dose_factor_t), parameter :: Yb_167  = dose_factor_t("Yb-167", 17.5, "m", 0.005, 7.0E-11, 0.0005, &
    & [ 4.1E-11, 2.1E-11, 1.2E-11, 8.4E-12], 0.0005, 6.7E-12 )
  type(dose_factor_t), parameter :: Yb_169  = dose_factor_t("Yb-169", 32.01, "d", 0.005, 7.1E-09, 0.0005, &
    & [ 4.6E-09, 2.4E-09, 1.5E-09, 8.8E-10], 0.0005, 7.1E-10 )
  type(dose_factor_t), parameter :: Yb_175  = dose_factor_t("Yb-175", 4.19, "d", 0.005, 5.0E-09, 0.0005, &
    & [ 3.2E-09, 1.6E-09, 9.5E-10, 5.4E-10], 0.0005, 4.4E-10 )
  type(dose_factor_t), parameter :: Yb_177  = dose_factor_t("Yb-177", 1.9, "h", 0.005, 1.0E-09, 0.0005, &
    & [ 6.8E-10, 3.4E-10, 2.0E-10, 1.1E-10], 0.0005, 8.8E-11 )
  type(dose_factor_t), parameter :: Yb_178  = dose_factor_t("Yb-178", 74, "m", 0.005, 1.4E-09, 0.0005, &
    & [ 8.4E-10, 4.2E-10, 2.4E-10, 1.5E-10], 0.0005, 1.2E-10 )
  type(dose_factor_t), parameter :: Lu_169  = dose_factor_t("Lu-169", 34.06, "h", 0.005, 3.5E-09, 0.0005, &
    & [ 2.4E-09, 1.4E-09, 8.9E-10, 5.7E-10], 0.0005, 4.6E-10 )
  type(dose_factor_t), parameter :: Lu_170  = dose_factor_t("Lu-170", 2.00, "d", 0.005, 7.4E-09, 0.0005, &
    & [ 5.2E-09, 2.9E-09, 1.9E-09, 1.2E-09], 0.0005, 9.9E-10 )
  type(dose_factor_t), parameter :: Lu_171  = dose_factor_t("Lu-171", 8.22, "d", 0.005, 5.9E-09, 0.0005, &
    & [ 4.0E-09, 2.2E-09, 1.4E-09, 8.5E-10], 0.0005, 6.7E-10 )
  type(dose_factor_t), parameter :: Lu_172  = dose_factor_t("Lu-172", 6.70, "d", 0.005, 1.0E-08, 0.0005, &
    & [ 7.0E-09, 3.9E-09, 2.5E-09, 1.6E-09], 0.0005, 1.3E-09 )
  type(dose_factor_t), parameter :: Lu_173  = dose_factor_t("Lu-173", 1.37, "y", 0.005, 2.7E-09, 0.0005, &
    & [ 1.6E-09, 8.6E-10, 5.3E-10, 3.2E-10], 0.0005, 2.6E-10 )
  type(dose_factor_t), parameter :: Lu_174  = dose_factor_t("Lu-174", 3.31, "y", 0.005, 3.2E-09, 0.0005, &
    & [ 1.7E-09, 9.1E-10, 5.6E-10, 3.3E-10], 0.0005, 2.7E-10 )
  type(dose_factor_t), parameter :: Lu_174m = dose_factor_t("Lu-174m", 142, "d", 0.005, 6.2E-09, 0.0005, &
    & [ 3.8E-09, 1.9E-09, 1.1E-09, 6.6E-10], 0.0005, 5.3E-10 )
  type(dose_factor_t), parameter :: Lu_176  = dose_factor_t("Lu-176", 3.60E10, "y", 0.005, 2.4E-08, 0.0005, &
    & [ 1.1E-08, 5.7E-09, 3.5E-09, 2.2E-09], 0.0005, 1.8E-09 )
  type(dose_factor_t), parameter :: Lu_176m = dose_factor_t("Lu-176m", 3.68, "h", 0.005, 2.0E-09, 0.0005, &
    & [ 1.2E-09, 6.0E-10, 3.5E-10, 2.1E-10], 0.0005, 1.7E-10 )
  type(dose_factor_t), parameter :: Lu_177  = dose_factor_t("Lu-177", 6.71, "d", 0.005, 6.1E-09, 0.0005, &
    & [ 3.9E-09, 2.0E-09, 1.2E-09, 6.6E-10], 0.0005, 5.3E-10 )
  type(dose_factor_t), parameter :: Lu_177m = dose_factor_t("Lu-177m", 160.9, "d", 0.005, 1.7E-08, 0.0005, &
    & [ 1.1E-08, 5.8E-09, 3.6E-09, 2.1E-09], 0.0005, 1.7E-09 )
  type(dose_factor_t), parameter :: Lu_178  = dose_factor_t("Lu-178", 28.4, "m", 0.005, 5.9E-10, 0.0005, &
    & [ 3.3E-10, 1.6E-10, 9.0E-11, 6.1E-11], 0.0005, 4.7E-11 )
  type(dose_factor_t), parameter :: Lu_178m = dose_factor_t("Lu-178m", 22.7, "m", 0.005, 4.3E-10, 0.0005, &
    & [ 2.4E-10, 1.2E-10, 7.1E-11, 4.9E-11], 0.0005, 3.8E-11 )
  type(dose_factor_t), parameter :: Lu_179  = dose_factor_t("Lu-179", 4.59, "h", 0.005, 2.4E-09, 0.0005, &
    & [ 1.5E-09, 7.5E-10, 4.4E-10, 2.6E-10], 0.0005, 2.1E-10 )
  type(dose_factor_t), parameter :: Hf_170  = dose_factor_t("Hf-170", 16.01, "h", 0.02, 3.9E-09, 0.002, &
    & [ 2.7E-09, 1.5E-09, 9.5E-10, 6.0E-10], 0.002, 4.8E-10 )
  type(dose_factor_t), parameter :: Hf_172  = dose_factor_t("Hf-172", 1.87, "y", 0.02, 1.9E-08, 0.002, &
    & [ 6.1E-09, 3.3E-09, 2.0E-09, 1.3E-09], 0.002, 1.0E-09 )
  type(dose_factor_t), parameter :: Hf_173  = dose_factor_t("Hf-173", 24.0, "h", 0.02, 1.9E-09, 0.002, &
    & [ 1.3E-09, 7.2E-10, 4.6E-10, 2.8E-10], 0.002, 2.3E-10 )
  type(dose_factor_t), parameter :: Hf_175  = dose_factor_t("Hf-175", 70, "d", 0.02, 3.8E-09, 0.002, &
    & [ 2.4E-09, 1.3E-09, 8.4E-10, 5.2E-10], 0.002, 4.1E-10 )
  type(dose_factor_t), parameter :: Hf_177m = dose_factor_t("Hf-177m", 51.4, "m", 0.02, 7.8E-10, 0.002, &
    & [ 4.7E-10, 2.5E-10, 1.5E-10, 1.0E-10], 0.002, 8.1E-11 )
  type(dose_factor_t), parameter :: Hf_178m = dose_factor_t("Hf-178m", 31, "y", 0.02, 7.0E-08, 0.002, &
    & [ 1.9E-08, 1.1E-08, 7.8E-09, 5.5E-09], 0.002, 4.7E-09 )
  type(dose_factor_t), parameter :: Hf_179m = dose_factor_t("Hf-179m", 25.1, "d", 0.02, 1.2E-08, 0.002, &
    & [ 7.8E-09, 4.1E-09, 2.6E-09, 1.6E-09], 0.002, 1.2E-09 )
  type(dose_factor_t), parameter :: Hf_180m = dose_factor_t("Hf-180m", 5.5, "h", 0.02, 1.4E-09, 0.002, &
    & [ 9.7E-10, 5.3E-10, 3.3E-10, 2.1E-10], 0.002, 1.7E-10 )
  type(dose_factor_t), parameter :: Hf_181  = dose_factor_t("Hf-181", 42.4, "d", 0.02, 1.2E-08, 0.002, &
    & [ 7.4E-09, 3.8E-09, 2.3E-09, 1.4E-09], 0.002, 1.1E-09 )
  type(dose_factor_t), parameter :: Hf_182  = dose_factor_t("Hf-182", 9E6, "y", 0.02, 5.6E-08, 0.002, &
    & [ 7.9E-09, 5.4E-09, 4.0E-09, 3.3E-09], 0.002, 3.0E-09 )
  type(dose_factor_t), parameter :: Hf_182m = dose_factor_t("Hf-182m", 61.5, "m", 0.02, 4.1E-10, 0.002, &
    & [ 2.5E-10, 1.3E-10, 7.8E-11, 5.2E-11], 0.002, 4.2E-11 )
  type(dose_factor_t), parameter :: Hf_183  = dose_factor_t("Hf-183", 64, "m", 0.02, 8.1E-10, 0.002, &
    & [ 4.8E-10, 2.4E-10, 1.4E-10, 9.3E-11], 0.002, 7.3E-11 )
  type(dose_factor_t), parameter :: Hf_184  = dose_factor_t("Hf-184", 4.12, "h", 0.02, 5.5E-09, 0.002, &
    & [ 3.6E-09, 1.8E-09, 1.1E-09, 6.6E-10], 0.002, 5.2E-10 )
  type(dose_factor_t), parameter :: Ta_172  = dose_factor_t("Ta-172", 36.8, "m", 0.01, 5.5E-10, 0.001, &
    & [ 3.2E-10, 1.6E-10, 9.8E-11, 6.6E-11], 0.001, 5.3E-11 )
  type(dose_factor_t), parameter :: Ta_173  = dose_factor_t("Ta-173", 3.65, "h", 0.01, 2.0E-09, 0.001, &
    & [ 1.3E-09, 6.5E-10, 3.9E-10, 2.4E-10], 0.001, 1.9E-10 )
  type(dose_factor_t), parameter :: Ta_174  = dose_factor_t("Ta-174", 1.2, "h", 0.01, 6.2E-10, 0.001, &
    & [ 3.7E-10, 1.9E-10, 1.1E-10, 7.2E-11], 0.001, 5.7E-11 )
  type(dose_factor_t), parameter :: Ta_175  = dose_factor_t("Ta-175", 10.5, "h", 0.01, 1.6E-09, 0.001, &
    & [ 1.1E-09, 6.2E-10, 4.0E-10, 2.6E-10], 0.001, 2.1E-10 )
  type(dose_factor_t), parameter :: Ta_176  = dose_factor_t("Ta-176", 8.08, "h", 0.01, 2.4E-09, 0.001, &
    & [ 1.7E-09, 9.2E-10, 6.1E-10, 3.9E-10], 0.001, 3.1E-10 )
  type(dose_factor_t), parameter :: Ta_177  = dose_factor_t("Ta-177", 56.6, "h", 0.01, 1.0E-09, 0.001, &
    & [ 6.9E-10, 3.6E-10, 2.2E-10, 1.3E-10], 0.001, 1.1E-10 )
  type(dose_factor_t), parameter :: Ta_178m = dose_factor_t("Ta-178m", 2.2, "h", 0.01, 6.3E-10, 0.001, &
    & [ 4.5E-10, 2.4E-10, 1.5E-10, 9.1E-11], 0.001, 7.2E-11 )
  type(dose_factor_t), parameter :: Ta_179  = dose_factor_t("Ta-179", 664.9, "d", 0.01, 6.2E-10, 0.001, &
    & [ 4.1E-10, 2.2E-10, 1.3E-10, 8.1E-11], 0.001, 6.5E-11 )
  type(dose_factor_t), parameter :: Ta_180  = dose_factor_t("Ta-180", 8.1, "h", 0.01, 5.8E-10, 0.001, &
    & [ 3.7E-10, 1.9E-10, 1.1E-10, 6.7E-11], 0.001, 5.4E-11 )
  type(dose_factor_t), parameter :: Ta_182  = dose_factor_t("Ta-182", 115.0, "d", 0.01, 1.4E-08, 0.001, &
    & [ 9.4E-09, 5.0E-09, 3.1E-09, 1.9E-09], 0.001, 1.5E-09 )
  type(dose_factor_t), parameter :: Ta_182m = dose_factor_t("Ta-182m", 15.84, "m", 0.01, 1.4E-10, 0.001, &
    & [ 7.5E-11, 3.7E-11, 2.1E-11, 1.5E-11], 0.001, 1.2E-11 )
  type(dose_factor_t), parameter :: Ta_183  = dose_factor_t("Ta-183", 5.1, "d", 0.01, 1.4E-08, 0.001, &
    & [ 9.3E-09, 4.7E-09, 2.8E-09, 1.6E-09], 0.001, 1.3E-09 )
  type(dose_factor_t), parameter :: Ta_184  = dose_factor_t("Ta-184", 8.7, "h", 0.01, 6.7E-09, 0.001, &
    & [ 4.4E-09, 2.3E-09, 1.4E-09, 8.5E-10], 0.001, 6.8E-10 )
  type(dose_factor_t), parameter :: Ta_185  = dose_factor_t("Ta-185", 49, "m", 0.01, 8.3E-10, 0.001, &
    & [ 4.6E-10, 2.3E-10, 1.3E-10, 8.6E-11], 0.001, 6.8E-11 )
  type(dose_factor_t), parameter :: Ta_186  = dose_factor_t("Ta-186", 10.5, "m", 0.01, 3.8E-10, 0.001, &
    & [ 2.1E-10, 1.1E-10, 6.1E-11, 4.2E-11], 0.001, 3.3E-11 )
  type(dose_factor_t), parameter :: W_176   = dose_factor_t("W-176", 2.3, "h", 0.6, 6.8E-10, 0.3, &
    & [ 5.5E-10, 3.0E-10, 2.0E-10, 1.3E-10], 0.3, 1.0E-10 )
  type(dose_factor_t), parameter :: W_177   = dose_factor_t("W-177", 135, "m", 0.6, 4.4E-10, 0.3, &
    & [ 3.2E-10, 1.7E-10, 1.1E-10, 7.2E-11], 0.3, 5.8E-11 )
  type(dose_factor_t), parameter :: W_178   = dose_factor_t("W-178", 21.7, "d", 0.6, 1.8E-09, 0.3, &
    & [ 1.4E-09, 7.3E-10, 4.5E-10, 2.7E-10], 0.3, 2.2E-10 )
  type(dose_factor_t), parameter :: W_179   = dose_factor_t("W-179", 37.5, "m", 0.6, 3.4E-11, 0.3, &
    & [ 2.0E-11, 1.0E-11, 6.2E-12, 4.2E-12], 0.3, 3.3E-12 )
  type(dose_factor_t), parameter :: W_181   = dose_factor_t("W-181", 121.2, "d", 0.6, 6.3E-10, 0.3, &
    & [ 4.7E-10, 2.5E-10, 1.6E-10, 9.5E-11], 0.3, 7.6E-11 )
  type(dose_factor_t), parameter :: W_185   = dose_factor_t("W-185", 75.1, "d", 0.6, 4.4E-09, 0.3, &
    & [ 3.3E-09, 1.6E-09, 9.7E-10, 5.5E-10], 0.3, 4.4E-10 )
  type(dose_factor_t), parameter :: W_187   = dose_factor_t("W-187", 23.9, "h", 0.6, 5.5E-09, 0.3, &
    & [ 4.3E-09, 2.2E-09, 1.3E-09, 7.8E-10], 0.3, 6.3E-10 )
  type(dose_factor_t), parameter :: W_188   = dose_factor_t("W-188", 69.4, "d", 0.6, 2.1E-08, 0.3, &
    & [ 1.5E-08, 7.7E-09, 4.6E-09, 2.6E-09], 0.3, 2.1E-09 )
  type(dose_factor_t), parameter :: Re_177  = dose_factor_t("Re-177", 14.0, "m", 1.0, 2.5E-10, 0.8, &
    & [ 1.4E-10, 7.2E-11, 4.1E-11, 2.8E-11], 0.8, 2.2E-11 )
  type(dose_factor_t), parameter :: Re_178  = dose_factor_t("Re-178", 13.2, "m", 1.0, 2.9E-10, 0.8, &
    & [ 1.6E-10, 7.9E-11, 4.6E-11, 3.1E-11], 0.8, 2.5E-11 )
  type(dose_factor_t), parameter :: Re_181  = dose_factor_t("Re-181", 20, "h", 1.0, 4.2E-09, 0.8, &
    & [ 2.8E-09, 1.4E-09, 8.2E-10, 5.4E-10], 0.8, 4.2E-10 )
  type(dose_factor_t), parameter :: Re_182  = dose_factor_t("Re-182", 64.0, "h", 1.0, 1.4E-08, 0.8, &
    & [ 8.9E-09, 4.7E-09, 2.8E-09, 1.8E-09], 0.8, 1.4E-09 )
  type(dose_factor_t), parameter :: Re_182m = dose_factor_t("Re-182m", 12.7, "h", 1.0, 2.4E-09, 0.8, &
    & [ 1.7E-09, 8.9E-10, 5.2E-10, 3.5E-10], 0.8, 2.7E-10 )
  type(dose_factor_t), parameter :: Re_184  = dose_factor_t("Re-184", 38.0, "d", 1.0, 8.9E-09, 0.8, &
    & [ 5.6E-09, 3.0E-09, 1.8E-09, 1.3E-09], 0.8, 1.0E-09 )
  type(dose_factor_t), parameter :: Re_184m = dose_factor_t("Re-184m", 165, "d", 1.0, 1.7E-08, 0.8, &
    & [ 9.8E-09, 4.9E-09, 2.8E-09, 1.9E-09], 0.8, 1.5E-09 )
  type(dose_factor_t), parameter :: Re_186  = dose_factor_t("Re-186", 90.64, "h", 1.0, 1.9E-08, 0.8, &
    & [ 1.1E-08, 5.5E-09, 3.0E-09, 1.9E-09], 0.8, 1.5E-09 )
  type(dose_factor_t), parameter :: Re_186m = dose_factor_t("Re-186m", 2.0E5, "y", 1.0, 3.0E-08, 0.8, &
    & [ 1.6E-08, 7.6E-09, 4.4E-09, 2.8E-09], 0.8, 2.2E-09 )
  type(dose_factor_t), parameter :: Re_187  = dose_factor_t("Re-187", 5E10, "y", 1.0, 6.8E-11, 0.8, &
    & [ 3.8E-11, 1.8E-11, 1.0E-11, 6.6E-12], 0.8, 5.1E-12 )
  type(dose_factor_t), parameter :: Re_188  = dose_factor_t("Re-188", 16.98, "h", 1.0, 1.7E-08, 0.8, &
    & [ 1.1E-08, 5.4E-09, 2.9E-09, 1.8E-09], 0.8, 1.4E-09 )
  type(dose_factor_t), parameter :: Re_188m = dose_factor_t("Re-188m", 18.6, "m", 1.0, 3.8E-10, 0.8, &
    & [ 2.3E-10, 1.1E-10, 6.1E-11, 4.0E-11], 0.8, 3.0E-11 )
  type(dose_factor_t), parameter :: Re_189  = dose_factor_t("Re-189", 24.3, "h", 1.0, 9.8E-09, 0.8, &
    & [ 6.2E-09, 3.0E-09, 1.6E-09, 1.0E-09], 0.8, 7.8E-10 )
  type(dose_factor_t), parameter :: Os_180  = dose_factor_t("Os-180", 22, "m", 0.02, 1.6E-10, 0.01, &
    & [ 9.8E-11, 5.1E-11, 3.2E-11, 2.2E-11], 0.01, 1.7E-11 )
  type(dose_factor_t), parameter :: Os_181  = dose_factor_t("Os-181", 105, "m", 0.02, 7.6E-10, 0.01, &
    & [ 5.0E-10, 2.7E-10, 1.7E-10, 1.1E-10], 0.01, 8.9E-11 )
  type(dose_factor_t), parameter :: Os_182  = dose_factor_t("Os-182", 22, "h", 0.02, 4.6E-09, 0.01, &
    & [ 3.2E-09, 1.7E-09, 1.1E-09, 7.0E-10], 0.01, 5.6E-10 )
  type(dose_factor_t), parameter :: Os_185  = dose_factor_t("Os-185", 94, "d", 0.02, 3.8E-09, 0.01, &
    & [ 2.6E-09, 1.5E-09, 9.8E-10, 6.5E-10], 0.01, 5.1E-10 )
  type(dose_factor_t), parameter :: Os_189m = dose_factor_t("Os-189m", 6.0, "h", 0.02, 2.1E-10, 0.01, &
    & [ 1.3E-10, 6.5E-11, 3.8E-11, 2.2E-11], 0.01, 1.8E-11 )
  type(dose_factor_t), parameter :: Os_191  = dose_factor_t("Os-191", 15.4, "d", 0.02, 6.3E-09, 0.01, &
    & [ 4.1E-09, 2.1E-09, 1.2E-09, 7.0E-10], 0.01, 5.7E-10 )
  type(dose_factor_t), parameter :: Os_191m = dose_factor_t("Os-191m", 13.03, "h", 0.02, 1.1E-09, 0.01, &
    & [ 7.1E-10, 3.5E-10, 2.1E-10, 1.2E-10], 0.01, 9.6E-11 )
  type(dose_factor_t), parameter :: Os_193  = dose_factor_t("Os-193", 30.0, "h", 0.02, 9.3E-09, 0.01, &
    & [ 6.0E-09, 3.0E-09, 1.8E-09, 1.0E-09], 0.01, 8.1E-10 )
  type(dose_factor_t), parameter :: Os_194  = dose_factor_t("Os-194", 6.0, "y", 0.02, 2.9E-08, 0.01, &
    & [ 1.7E-08, 8.8E-09, 5.2E-09, 3.0E-09], 0.01, 2.4E-09 )
  type(dose_factor_t), parameter :: Ir_182  = dose_factor_t("Ir-182", 15, "m", 0.02, 5.3E-10, 0.01, &
    & [ 3.0E-10, 1.5E-10, 8.9E-11, 6.0E-11], 0.01, 4.8E-11 )
  type(dose_factor_t), parameter :: Ir_184  = dose_factor_t("Ir-184", 3.02, "h", 0.02, 1.5E-09, 0.01, &
    & [ 9.7E-10, 5.2E-10, 3.3E-10, 2.1E-10], 0.01, 1.7E-10 )
  type(dose_factor_t), parameter :: Ir_185  = dose_factor_t("Ir-185", 14.0, "h", 0.02, 2.4E-09, 0.01, &
    & [ 1.6E-09, 8.6E-10, 5.3E-10, 3.3E-10], 0.01, 2.6E-10 )
  type(dose_factor_t), parameter :: Ir_186  = dose_factor_t("Ir-186", 15.8, "h", 0.02, 3.8E-09, 0.01, &
    & [ 2.7E-09, 1.5E-09, 9.6E-10, 6.1E-10], 0.01, 4.9E-10 )
  type(dose_factor_t), parameter :: Ir_186m = dose_factor_t("Ir-186m", 1.75, "h", 0.02, 5.8E-10, 0.01, &
    & [ 3.6E-10, 2.1E-10, 1.3E-10, 7.7E-11], 0.01, 6.1E-11 )
  type(dose_factor_t), parameter :: Ir_187  = dose_factor_t("Ir-187", 10.5, "h", 0.02, 1.1E-09, 0.01, &
    & [ 7.3E-10, 3.9E-10, 2.5E-10, 1.5E-10], 0.01, 1.2E-10 )
  type(dose_factor_t), parameter :: Ir_188  = dose_factor_t("Ir-188", 41.5, "h", 0.02, 4.6E-09, 0.01, &
    & [ 3.3E-09, 1.8E-09, 1.2E-09, 7.9E-10], 0.01, 6.3E-10 )
  type(dose_factor_t), parameter :: Ir_189  = dose_factor_t("Ir-189", 13.3, "d", 0.02, 2.5E-09, 0.01, &
    & [ 1.7E-09, 8.6E-10, 5.2E-10, 3.0E-10], 0.01, 2.4E-10 )
  type(dose_factor_t), parameter :: Ir_190  = dose_factor_t("Ir-190", 12.1, "d", 0.02, 1.0E-08, 0.01, &
    & [ 7.1E-09, 3.9E-09, 2.5E-09, 1.6E-09], 0.01, 1.2E-09 )
  type(dose_factor_t), parameter :: Ir_190n = dose_factor_t("Ir-190n", 3.1, "h", 0.02, 9.4E-10, 0.01, &
    & [ 6.4E-10, 3.5E-10, 2.3E-10, 1.5E-10], 0.01, 1.2E-10 )
  type(dose_factor_t), parameter :: Ir_190m = dose_factor_t("Ir-190m", 1.2, "h", 0.02, 7.9E-11, 0.01, &
    & [ 5.0E-11, 2.6E-11, 1.6E-11, 1.0E-11], 0.01, 8.0E-12 )
  type(dose_factor_t), parameter :: Ir_192  = dose_factor_t("Ir-192", 74.02, "d", 0.02, 1.3E-08, 0.01, &
    & [ 8.7E-09, 4.6E-09, 2.8E-09, 1.7E-09], 0.01, 1.4E-09 )
  type(dose_factor_t), parameter :: Ir_192n = dose_factor_t("Ir-192n", 241, "y", 0.02, 2.8E-09, 0.01, &
    & [ 1.4E-09, 8.3E-10, 5.5E-10, 3.7E-10], 0.01, 3.1E-10 )
  type(dose_factor_t), parameter :: Ir_193m = dose_factor_t("Ir-193m", 11.9, "d", 0.02, 3.2E-09, 0.01, &
    & [ 2.0E-09, 1.0E-09, 6.0E-10, 3.4E-10], 0.01, 2.7E-10 )
  type(dose_factor_t), parameter :: Ir_194  = dose_factor_t("Ir-194", 19.15, "h", 0.02, 1.5E-08, 0.01, &
    & [ 9.8E-09, 4.9E-09, 2.9E-09, 1.7E-09], 0.01, 1.3E-09 )
  type(dose_factor_t), parameter :: Ir_194m = dose_factor_t("Ir-194m", 171, "d", 0.02, 1.7E-08, 0.01, &
    & [ 1.1E-08, 6.4E-09, 4.1E-09, 2.6E-09], 0.01, 2.1E-09 )
  type(dose_factor_t), parameter :: Ir_195  = dose_factor_t("Ir-195", 2.5, "h", 0.02, 1.2E-09, 0.01, &
    & [ 7.3E-10, 3.6E-10, 2.1E-10, 1.3E-10], 0.01, 1.0E-10 )
  type(dose_factor_t), parameter :: Ir_195m = dose_factor_t("Ir-195m", 3.8, "h", 0.02, 2.3E-09, 0.01, &
    & [ 1.5E-09, 7.3E-10, 4.3E-10, 2.6E-10], 0.01, 2.1E-10 )
  type(dose_factor_t), parameter :: Pt_186  = dose_factor_t("Pt-186", 2.0, "h", 0.02, 7.8E-10, 0.01, &
    & [ 5.3E-10, 2.9E-10, 1.8E-10, 1.2E-10], 0.01, 9.3E-11 )
  type(dose_factor_t), parameter :: Pt_188  = dose_factor_t("Pt-188", 10.2, "d", 0.02, 6.7E-09, 0.01, &
    & [ 4.5E-09, 2.4E-09, 1.5E-09, 9.5E-10], 0.01, 7.6E-10 )
  type(dose_factor_t), parameter :: Pt_189  = dose_factor_t("Pt-189", 10.87, "h", 0.02, 1.1E-09, 0.01, &
    & [ 7.4E-10, 3.9E-10, 2.5E-10, 1.5E-10], 0.01, 1.2E-10 )
  type(dose_factor_t), parameter :: Pt_191  = dose_factor_t("Pt-191", 2.8, "d", 0.02, 3.1E-09, 0.01, &
    & [ 2.1E-09, 1.1E-09, 6.9E-10, 4.2E-10], 0.01, 3.4E-10 )
  type(dose_factor_t), parameter :: Pt_193  = dose_factor_t("Pt-193", 50, "y", 0.02, 3.7E-10, 0.01, &
    & [ 2.4E-10, 1.2E-10, 6.9E-11, 3.9E-11], 0.01, 3.1E-11 )
  type(dose_factor_t), parameter :: Pt_193m = dose_factor_t("Pt-193m", 4.33, "d", 0.02, 5.2E-09, 0.01, &
    & [ 3.4E-09, 1.7E-09, 9.9E-10, 5.6E-10], 0.01, 4.5E-10 )
  type(dose_factor_t), parameter :: Pt_195m = dose_factor_t("Pt-195m", 4.02, "d", 0.02, 7.1E-09, 0.01, &
    & [ 4.6E-09, 2.3E-09, 1.4E-09, 7.9E-10], 0.01, 6.3E-10 )
  type(dose_factor_t), parameter :: Pt_197  = dose_factor_t("Pt-197", 18.3, "h", 0.02, 4.7E-09, 0.01, &
    & [ 3.0E-09, 1.5E-09, 8.8E-10, 5.1E-10], 0.01, 4.0E-10 )
  type(dose_factor_t), parameter :: Pt_197m = dose_factor_t("Pt-197m", 94.4, "m", 0.02, 1.0E-09, 0.01, &
    & [ 6.1E-10, 3.0E-10, 1.8E-10, 1.1E-10], 0.01, 8.4E-11 )
  type(dose_factor_t), parameter :: Pt_199  = dose_factor_t("Pt-199", 30.8, "m", 0.02, 4.7E-10, 0.01, &
    & [ 2.7E-10, 1.3E-10, 7.5E-11, 5.0E-11], 0.01, 3.9E-11 )
  type(dose_factor_t), parameter :: Pt_200  = dose_factor_t("Pt-200", 12.5, "h", 0.02, 1.4E-08, 0.01, &
    & [ 8.8E-09, 4.4E-09, 2.6E-09, 1.5E-09], 0.01, 1.2E-09 )
  type(dose_factor_t), parameter :: Au_193  = dose_factor_t("Au-193", 17.65, "h", 0.2, 1.2E-09, 0.1, &
    & [ 8.8E-10, 4.6E-10, 2.8E-10, 1.7E-10], 0.1, 1.3E-10 )
  type(dose_factor_t), parameter :: Au_194  = dose_factor_t("Au-194", 39.5, "h", 0.2, 2.9E-09, 0.1, &
    & [ 2.2E-09, 1.2E-09, 8.1E-10, 5.3E-10], 0.1, 4.2E-10 )
  type(dose_factor_t), parameter :: Au_195  = dose_factor_t("Au-195", 183, "d", 0.2, 2.4E-09, 0.1, &
    & [ 1.7E-09, 8.9E-10, 5.4E-10, 3.2E-10], 0.1, 2.5E-10 )
  type(dose_factor_t), parameter :: Au_198  = dose_factor_t("Au-198", 2.696, "d", 0.2, 1.0E-08, 0.1, &
    & [ 7.2E-09, 3.7E-09, 2.2E-09, 1.3E-09], 0.1, 1.0E-09 )
  type(dose_factor_t), parameter :: Au_198m = dose_factor_t("Au-198m", 2.30, "d", 0.2, 1.2E-08, 0.1, &
    & [ 8.5E-09, 4.4E-09, 2.7E-09, 1.6E-09], 0.1, 1.3E-09 )
  type(dose_factor_t), parameter :: Au_199  = dose_factor_t("Au-199", 3.139, "d", 0.2, 4.5E-09, 0.1, &
    & [ 3.1E-09, 1.6E-09, 9.5E-10, 5.5E-10], 0.1, 4.4E-10 )
  type(dose_factor_t), parameter :: Au_200  = dose_factor_t("Au-200", 48.4, "m", 0.2, 8.3E-10, 0.1, &
    & [ 4.7E-10, 2.3E-10, 1.3E-10, 8.7E-11], 0.1, 6.8E-11 )
  type(dose_factor_t), parameter :: Au_200m = dose_factor_t("Au-200m", 18.7, "h", 0.2, 9.2E-09, 0.1, &
    & [ 6.6E-09, 3.5E-09, 2.2E-09, 1.3E-09], 0.1, 1.1E-09 )
  type(dose_factor_t), parameter :: Au_201  = dose_factor_t("Au-201", 26.4, "m", 0.2, 3.1E-10, 0.1, &
    & [ 1.7E-10, 8.2E-11, 4.6E-11, 3.1E-11], 0.1, 2.4E-11 )
  type(dose_factor_t), parameter :: Hg_193  = dose_factor_t("Hg-193", 3.5, "h", 0.04, 8.5E-10, 0.02, &
    & [ 5.5E-10, 2.8E-10, 1.7E-10, 1.0E-10], 0.02, 8.2E-11 )
  type(dose_factor_t), parameter :: Hg_193m = dose_factor_t("Hg-193m", 11.1, "h", 0.04, 3.6E-09, 0.02, &
    & [ 2.4E-09, 1.3E-09, 8.1E-10, 5.0E-10], 0.02, 4.0E-10 )
  type(dose_factor_t), parameter :: Hg_194  = dose_factor_t("Hg-194", 260, "y", 0.04, 7.2E-09, 0.02, &
    & [ 3.6E-09, 2.6E-09, 1.9E-09, 1.5E-09], 0.02, 1.4E-09 )
  type(dose_factor_t), parameter :: Hg_195  = dose_factor_t("Hg-195", 9.9, "h", 0.04, 9.5E-10, 0.02, &
    & [ 6.3E-10, 3.3E-10, 2.0E-10, 1.2E-10], 0.02, 9.7E-11 )
  type(dose_factor_t), parameter :: Hg_195m = dose_factor_t("Hg-195m", 41.6, "h", 0.04, 5.8E-09, 0.02, &
    & [ 3.8E-09, 2.0E-09, 1.2E-09, 7.0E-10], 0.02, 5.6E-10 )
  type(dose_factor_t), parameter :: Hg_197  = dose_factor_t("Hg-197", 64.1, "h", 0.04, 2.5E-09, 0.02, &
    & [ 1.6E-09, 8.3E-10, 5.0E-10, 2.9E-10], 0.02, 2.3E-10 )
  type(dose_factor_t), parameter :: Hg_197m = dose_factor_t("Hg-197m", 23.8, "h", 0.04, 5.2E-09, 0.02, &
    & [ 3.4E-09, 1.7E-09, 1.0E-09, 5.9E-10], 0.02, 4.7E-10 )
  type(dose_factor_t), parameter :: Hg_199m = dose_factor_t("Hg-199m", 42.6, "m", 0.04, 3.7E-10, 0.02, &
    & [ 2.1E-10, 1.0E-10, 5.9E-11, 3.9E-11], 0.02, 3.1E-11 )
  type(dose_factor_t), parameter :: Hg_203  = dose_factor_t("Hg-203", 46.60, "d", 0.04, 5.5E-09, 0.02, &
    & [ 3.6E-09, 1.8E-09, 1.1E-09, 6.7E-10], 0.02, 5.4E-10 )
  type(dose_factor_t), parameter :: Tl_194  = dose_factor_t("Tl-194", 33, "m", 1.0, 6.1E-11, 1.0, &
    & [ 3.9E-11, 2.2E-11, 1.4E-11, 1.0E-11], 1.0, 8.1E-12 )
  type(dose_factor_t), parameter :: Tl_194m = dose_factor_t("Tl-194m", 32.8, "m", 1.0, 3.8E-10, 1.0, &
    & [ 2.2E-10, 1.2E-10, 7.0E-11, 4.9E-11], 1.0, 4.0E-11 )
  type(dose_factor_t), parameter :: Tl_195  = dose_factor_t("Tl-195", 1.16, "h", 1.0, 2.3E-10, 1.0, &
    & [ 1.4E-10, 7.5E-11, 4.7E-11, 3.3E-11], 1.0, 2.7E-11 )
  type(dose_factor_t), parameter :: Tl_197  = dose_factor_t("Tl-197", 2.84, "h", 1.0, 2.1E-10, 1.0, &
    & [ 1.3E-10, 6.7E-11, 4.2E-11, 2.8E-11], 1.0, 2.3E-11 )
  type(dose_factor_t), parameter :: Tl_198  = dose_factor_t("Tl-198", 5.3, "h", 1.0, 4.7E-10, 1.0, &
    & [ 3.3E-10, 1.9E-10, 1.2E-10, 8.7E-11], 1.0, 7.3E-11 )
  type(dose_factor_t), parameter :: Tl_198m = dose_factor_t("Tl-198m", 1.87, "h", 1.0, 4.8E-10, 1.0, &
    & [ 3.0E-10, 1.6E-10, 9.7E-11, 6.7E-11], 1.0, 5.4E-11 )
  type(dose_factor_t), parameter :: Tl_199  = dose_factor_t("Tl-199", 7.42, "h", 1.0, 2.3E-10, 1.0, &
    & [ 1.5E-10, 7.7E-11, 4.8E-11, 3.2E-11], 1.0, 2.6E-11 )
  type(dose_factor_t), parameter :: Tl_200  = dose_factor_t("Tl-200", 26.1, "h", 1.0, 1.3E-09, 1.0, &
    & [ 9.1E-10, 5.3E-10, 3.5E-10, 2.4E-10], 1.0, 2.0E-10 )
  type(dose_factor_t), parameter :: Tl_201  = dose_factor_t("Tl-201", 3.044, "d", 1.0, 8.4E-10, 1.0, &
    & [ 5.5E-10, 2.9E-10, 1.8E-10, 1.2E-10], 1.0, 9.5E-11 )
  type(dose_factor_t), parameter :: Tl_202  = dose_factor_t("Tl-202", 12.23, "d", 1.0, 2.9E-09, 1.0, &
    & [ 2.1E-09, 1.2E-09, 7.9E-10, 5.4E-10], 1.0, 4.5E-10 )
  type(dose_factor_t), parameter :: Tl_204  = dose_factor_t("Tl-204", 3.779, "y", 1.0, 1.3E-08, 1.0, &
    & [ 8.5E-09, 4.2E-09, 2.5E-09, 1.5E-09], 1.0, 1.2E-09 )
  type(dose_factor_t), parameter :: Pb_195m = dose_factor_t("Pb-195m", 15.8, "m", 0.6, 2.6E-10, 0.4, &
    & [ 1.6E-10, 8.4E-11, 5.2E-11, 3.5E-11], 0.3, 2.9E-11 )
  type(dose_factor_t), parameter :: Pb_198  = dose_factor_t("Pb-198", 2.4, "h", 0.6, 5.9E-10, 0.4, &
    & [ 4.8E-10, 2.7E-10, 1.7E-10, 1.1E-10], 0.3, 1.0E-10 )
  type(dose_factor_t), parameter :: Pb_199  = dose_factor_t("Pb-199", 90, "m", 0.6, 3.5E-10, 0.4, &
    & [ 2.6E-10, 1.5E-10, 9.4E-11, 6.3E-11], 0.3, 5.4E-11 )
  type(dose_factor_t), parameter :: Pb_200  = dose_factor_t("Pb-200", 21.5, "h", 0.6, 2.5E-09, 0.4, &
    & [ 2.0E-09, 1.1E-09, 7.0E-10, 4.4E-10], 0.3, 4.0E-10 )
  type(dose_factor_t), parameter :: Pb_201  = dose_factor_t("Pb-201", 9.4, "h", 0.6, 9.4E-10, 0.4, &
    & [ 7.8E-10, 4.3E-10, 2.7E-10, 1.8E-10], 0.3, 1.6E-10 )
  type(dose_factor_t), parameter :: Pb_202  = dose_factor_t("Pb-202", 3E5, "y", 0.6, 3.4E-08, 0.4, &
    & [ 1.6E-08, 1.3E-08, 1.9E-08, 2.7E-08], 0.3, 8.8E-09 )
  type(dose_factor_t), parameter :: Pb_202m = dose_factor_t("Pb-202m", 3.62, "h", 0.6, 7.6E-10, 0.4, &
    & [ 6.1E-10, 3.5E-10, 2.3E-10, 1.5E-10], 0.3, 1.3E-10 )
  type(dose_factor_t), parameter :: Pb_203  = dose_factor_t("Pb-203", 52.05, "h", 0.6, 1.6E-09, 0.4, &
    & [ 1.3E-09, 6.8E-10, 4.3E-10, 2.7E-10], 0.3, 2.4E-10 )
  type(dose_factor_t), parameter :: Pb_205  = dose_factor_t("Pb-205", 1.43E7, "y", 0.6, 2.1E-09, 0.4, &
    & [ 9.9E-10, 6.2E-10, 6.1E-10, 6.5E-10], 0.3, 2.8E-10 )
  type(dose_factor_t), parameter :: Pb_209  = dose_factor_t("Pb-209", 3.253, "h", 0.6, 5.7E-10, 0.4, &
    & [ 3.8E-10, 1.9E-10, 1.1E-10, 6.6E-11], 0.3, 5.7E-11 )
  type(dose_factor_t), parameter :: Pb_210  = dose_factor_t("Pb-210", 22.3, "y", 0.6, 8.4E-06, 0.4, &
    & [ 3.6E-06, 2.2E-06, 1.9E-06, 1.9E-06], 0.3, 6.9E-07 )
  type(dose_factor_t), parameter :: Pb_211  = dose_factor_t("Pb-211", 36.1, "m", 0.6, 3.1E-09, 0.4, &
    & [ 1.4E-09, 7.1E-10, 4.1E-10, 2.7E-10], 0.3, 1.8E-10 )
  type(dose_factor_t), parameter :: Pb_212  = dose_factor_t("Pb-212", 10.64, "h", 0.6, 1.5E-07, 0.4, &
    & [ 6.3E-08, 3.3E-08, 2.0E-08, 1.3E-08], 0.3, 6.0E-09 )
  type(dose_factor_t), parameter :: Pb_214  = dose_factor_t("Pb-214", 26.8, "m", 0.6, 2.7E-09, 0.4, &
    & [ 1.0E-09, 5.2E-10, 3.1E-10, 2.0E-10], 0.3, 1.4E-10 )
  type(dose_factor_t), parameter :: Bi_200  = dose_factor_t("Bi-200", 36.4, "m", 0.1, 4.2E-10, 0.05, &
    & [ 2.7E-10, 1.5E-10, 9.5E-11, 6.4E-11], 0.05, 5.1E-11 )
  type(dose_factor_t), parameter :: Bi_201  = dose_factor_t("Bi-201", 108, "m", 0.1, 1.0E-09, 0.05, &
    & [ 6.7E-10, 3.6E-10, 2.2E-10, 1.4E-10], 0.05, 1.2E-10 )
  type(dose_factor_t), parameter :: Bi_202  = dose_factor_t("Bi-202", 1.67, "h", 0.1, 6.4E-10, 0.05, &
    & [ 4.4E-10, 2.5E-10, 1.6E-10, 1.1E-10], 0.05, 8.9E-11 )
  type(dose_factor_t), parameter :: Bi_203  = dose_factor_t("Bi-203", 11.76, "h", 0.1, 3.5E-09, 0.05, &
    & [ 2.5E-09, 1.4E-09, 9.3E-10, 6.0E-10], 0.05, 4.8E-10 )
  type(dose_factor_t), parameter :: Bi_205  = dose_factor_t("Bi-205", 15.31, "d", 0.1, 6.1E-09, 0.05, &
    & [ 4.5E-09, 2.6E-09, 1.7E-09, 1.1E-09], 0.05, 9.0E-10 )
  type(dose_factor_t), parameter :: Bi_206  = dose_factor_t("Bi-206", 6.243, "d", 0.1, 1.4E-08, 0.05, &
    & [ 1.0E-08, 5.7E-09, 3.7E-09, 2.4E-09], 0.05, 1.9E-09 )
  type(dose_factor_t), parameter :: Bi_207  = dose_factor_t("Bi-207", 38, "y", 0.1, 1.0E-08, 0.05, &
    & [ 7.1E-09, 3.9E-09, 2.5E-09, 1.6E-09], 0.05, 1.3E-09 )
  type(dose_factor_t), parameter :: Bi_210  = dose_factor_t("Bi-210", 5.012, "d", 0.1, 1.5E-08, 0.05, &
    & [ 9.7E-09, 4.8E-09, 2.9E-09, 1.6E-09], 0.05, 1.3E-09 )
  type(dose_factor_t), parameter :: Bi_210m = dose_factor_t("Bi-210m", 3.0E6, "y", 0.1, 2.1E-07, 0.05, &
    & [ 9.1E-08, 4.7E-08, 3.0E-08, 1.9E-08], 0.05, 1.5E-08 )
  type(dose_factor_t), parameter :: Bi_212  = dose_factor_t("Bi-212", 60.55, "m", 0.1, 3.2E-09, 0.05, &
    & [ 1.8E-09, 8.7E-10, 5.0E-10, 3.3E-10], 0.05, 2.6E-10 )
  type(dose_factor_t), parameter :: Bi_213  = dose_factor_t("Bi-213", 45.65, "m", 0.1, 2.5E-09, 0.05, &
    & [ 1.4E-09, 6.7E-10, 3.9E-10, 2.5E-10], 0.05, 2.0E-10 )
  type(dose_factor_t), parameter :: Bi_214  = dose_factor_t("Bi-214", 19.9, "m", 0.1, 1.4E-09, 0.05, &
    & [ 7.4E-10, 3.6E-10, 2.1E-10, 1.4E-10], 0.05, 1.1E-10 )
  type(dose_factor_t), parameter :: Po_203  = dose_factor_t("Po-203", 36.7, "m", 1.0, 2.9E-10, 0.5, &
    & [ 2.4E-10, 1.3E-10, 8.5E-11, 5.8E-11], 0.5, 4.6E-11 )
  type(dose_factor_t), parameter :: Po_205  = dose_factor_t("Po-205", 1.80, "h", 1.0, 3.5E-10, 0.5, &
    & [ 2.8E-10, 1.6E-10, 1.1E-10, 7.2E-11], 0.5, 5.8E-11 )
  type(dose_factor_t), parameter :: Po_207  = dose_factor_t("Po-207", 350, "m", 1.0, 4.4E-10, 0.5, &
    & [ 5.7E-10, 3.2E-10, 2.1E-10, 1.4E-10], 0.5, 1.1E-10 )
  type(dose_factor_t), parameter :: Po_210  = dose_factor_t("Po-210", 138.38, "d", 1.0, 2.6E-05, 0.5, &
    & [ 8.8E-06, 4.4E-06, 2.6E-06, 1.6E-06], 0.5, 1.2E-06 )
  type(dose_factor_t), parameter :: At_207  = dose_factor_t("At-207", 1.80, "h", 1.0, 2.5E-09, 1.0, &
    & [ 1.6E-09, 8.0E-10, 4.8E-10, 2.9E-10], 1.0, 2.4E-10 )
  type(dose_factor_t), parameter :: At_211  = dose_factor_t("At-211", 7.214, "h", 1.0, 1.2E-07, 1.0, &
    & [ 7.8E-08, 3.8E-08, 2.3E-08, 1.3E-08], 1.0, 1.1E-08 )
  type(dose_factor_t), parameter :: Fr_222  = dose_factor_t("Fr-222", 14.4, "m", 1.0, 6.2E-09, 1.0, &
    & [ 3.9E-09, 2.0E-09, 1.3E-09, 8.5E-10], 1.0, 7.2E-10 )
  type(dose_factor_t), parameter :: Fr_223  = dose_factor_t("Fr-223", 21.8, "m", 1.0, 2.6E-08, 1.0, &
    & [ 1.7E-08, 8.3E-09, 5.0E-09, 2.9E-09], 1.0, 2.4E-09 )
  type(dose_factor_t), parameter :: Ra_223  = dose_factor_t("Ra-223", 11.434, "d", 0.6, 5.3E-06, 0.3, &
    & [ 1.1E-06, 5.7E-07, 4.5E-07, 3.7E-07], 0.3, 1.0E-07 )
  type(dose_factor_t), parameter :: Ra_224  = dose_factor_t("Ra-224", 3.66, "d", 0.6, 2.7E-06, 0.3, &
    & [ 6.6E-07, 3.5E-07, 2.6E-07, 2.0E-07], 0.3, 6.5E-08 )
  type(dose_factor_t), parameter :: Ra_225  = dose_factor_t("Ra-225", 14.8, "d", 0.6, 7.1E-06, 0.3, &
    & [ 1.2E-06, 6.1E-07, 5.0E-07, 4.4E-07], 0.3, 9.9E-08 )
  type(dose_factor_t), parameter :: Ra_226  = dose_factor_t("Ra-226", 1600, "y", 0.6, 4.7E-06, 0.3, &
    & [ 9.6E-07, 6.2E-07, 8.0E-07, 1.5E-06], 0.3, 2.8E-07 )
  type(dose_factor_t), parameter :: Ra_227  = dose_factor_t("Ra-227", 42.2, "m", 0.6, 1.1E-09, 0.3, &
    & [ 4.3E-10, 2.5E-10, 1.7E-10, 1.3E-10], 0.3, 8.1E-11 )
  type(dose_factor_t), parameter :: Ra_228  = dose_factor_t("Ra-228", 5.75, "y", 0.6, 3.0E-05, 0.3, &
    & [ 5.7E-06, 3.4E-06, 3.9E-06, 5.3E-06], 0.3, 6.9E-07 )
  type(dose_factor_t), parameter :: Ac_224  = dose_factor_t("Ac-224", 2.9, "h", 0.005, 1.0E-08, 0.0005, &
    & [ 5.2E-09, 2.6E-09, 1.5E-09, 8.8E-10], 0.0005, 7.0E-10 )
  type(dose_factor_t), parameter :: Ac_225  = dose_factor_t("Ac-225", 10.0, "d", 0.005, 4.6E-07, 0.0005, &
    & [ 1.8E-07, 9.1E-08, 5.4E-08, 3.0E-08], 0.0005, 2.4E-08 )
  type(dose_factor_t), parameter :: Ac_226  = dose_factor_t("Ac-226", 29, "h", 0.005, 1.4E-07, 0.0005, &
    & [ 7.6E-08, 3.8E-08, 2.3E-08, 1.3E-08], 0.0005, 1.0E-08 )
  type(dose_factor_t), parameter :: Ac_227  = dose_factor_t("Ac-227", 21.773, "y", 0.005, 3.3E-05, 0.0005, &
    & [ 3.1E-06, 2.2E-06, 1.5E-06, 1.2E-06], 0.0005, 1.1E-06 )
  type(dose_factor_t), parameter :: Ac_228  = dose_factor_t("Ac-228", 6.13, "h", 0.005, 7.4E-09, 0.0005, &
    & [ 2.8E-09, 1.4E-09, 8.7E-10, 5.3E-10], 0.0005, 4.3E-10 )
  type(dose_factor_t), parameter :: Th_226  = dose_factor_t("Th-226", 30.9, "m", 0.005, 4.4E-09, 0.0005, &
    & [ 2.4E-09, 1.2E-09, 6.7E-10, 4.5E-10], 0.0005, 3.5E-10 )
  type(dose_factor_t), parameter :: Th_227  = dose_factor_t("Th-227", 18.718, "d", 0.005, 3.0E-07, 0.0005, &
    & [ 7.0E-08, 3.6E-08, 2.3E-08, 1.5E-08], 0.0005, 8.8E-09 )
  type(dose_factor_t), parameter :: Th_228  = dose_factor_t("Th-228", 1.9131, "y", 0.005, 3.7E-06, 0.0005, &
    & [ 3.7E-07, 2.2E-07, 1.4E-07, 9.4E-08], 0.0005, 7.2E-08 )
  type(dose_factor_t), parameter :: Th_229  = dose_factor_t("Th-229", 7340, "y", 0.005, 1.1E-05, 0.0005, &
    & [ 1.0E-06, 7.8E-07, 6.2E-07, 5.3E-07], 0.0005, 4.9E-07 )
  type(dose_factor_t), parameter :: Th_230  = dose_factor_t("Th-230", 7.7E4, "y", 0.005, 4.1E-06, 0.0005, &
    & [ 4.1E-07, 3.1E-07, 2.4E-07, 2.2E-07], 0.0005, 2.1E-07 )
  type(dose_factor_t), parameter :: Th_231  = dose_factor_t("Th-231", 25.52, "h", 0.005, 3.9E-09, 0.0005, &
    & [ 2.5E-09, 1.2E-09, 7.4E-10, 4.2E-10], 0.0005, 3.4E-10 )
  type(dose_factor_t), parameter :: Th_232  = dose_factor_t("Th-232", 1.405E10, "y", 0.005, 4.6E-06, 0.0005, &
    & [ 4.5E-07, 3.5E-07, 2.9E-07, 2.5E-07], 0.0005, 2.3E-07 )
  type(dose_factor_t), parameter :: Th_234  = dose_factor_t("Th-234", 24.10, "d", 0.005, 4.0E-08, 0.0005, &
    & [ 2.5E-08, 1.3E-08, 7.4E-09, 4.2E-09], 0.0005, 3.4E-09 )
  type(dose_factor_t), parameter :: Pa_227  = dose_factor_t("Pa-227", 38.3, "m", 0.005, 5.8E-09, 0.0005, &
    & [ 3.2E-09, 1.5E-09, 8.7E-10, 5.8E-10], 0.0005, 4.5E-10 )
  type(dose_factor_t), parameter :: Pa_228  = dose_factor_t("Pa-228", 22, "h", 0.005, 1.2E-08, 0.0005, &
    & [ 4.8E-09, 2.6E-09, 1.6E-09, 9.7E-10], 0.0005, 7.8E-10 )
  type(dose_factor_t), parameter :: Pa_230  = dose_factor_t("Pa-230", 17.4, "d", 0.005, 2.6E-08, 0.0005, &
    & [ 5.7E-09, 3.1E-09, 1.9E-09, 1.1E-09], 0.0005, 9.2E-10 )
  type(dose_factor_t), parameter :: Pa_231  = dose_factor_t("Pa-231", 3.276E4, "y", 0.005, 1.3E-05, 0.0005, &
    & [ 1.3E-06, 1.1E-06, 9.2E-07, 8.0E-07], 0.0005, 7.1E-07 )
  type(dose_factor_t), parameter :: Pa_232  = dose_factor_t("Pa-232", 1.31, "d", 0.005, 7.2E-09, 0.0005, &
    & [ 4.3E-09, 2.3E-09, 1.4E-09, 8.9E-10], 0.0005, 7.2E-10 )
  type(dose_factor_t), parameter :: Pa_233  = dose_factor_t("Pa-233", 27.0, "d", 0.005, 9.7E-09, 0.0005, &
    & [ 6.2E-09, 3.2E-09, 1.9E-09, 1.1E-09], 0.0005, 8.7E-10 )
  type(dose_factor_t), parameter :: Pa_234  = dose_factor_t("Pa-234", 6.70, "h", 0.005, 5.0E-09, 0.0005, &
    & [ 3.2E-09, 1.7E-09, 1.0E-09, 6.4E-10], 0.0005, 5.1E-10 )
  type(dose_factor_t), parameter :: U_230   = dose_factor_t("U-230", 20.8, "d", 0.04, 7.9E-07, 0.02, &
    & [ 3.0E-07, 1.5E-07, 1.0E-07, 6.6E-08], 0.02, 5.6E-08 )
  type(dose_factor_t), parameter :: U_231   = dose_factor_t("U-231", 4.2, "d", 0.04, 3.1E-09, 0.02, &
    & [ 2.0E-09, 1.0E-09, 6.1E-10, 3.6E-10], 0.02, 2.8E-10 )
  type(dose_factor_t), parameter :: U_232   = dose_factor_t("U-232", 72, "y", 0.04, 2.5E-06, 0.02, &
    & [ 8.2E-07, 5.8E-07, 5.7E-07, 6.4E-07], 0.02, 3.3E-07 )
  type(dose_factor_t), parameter :: U_233   = dose_factor_t("U-233", 1.585E5, "y", 0.04, 3.8E-07, 0.02, &
    & [ 1.4E-07, 9.2E-08, 7.8E-08, 7.8E-08], 0.02, 5.1E-08 )
  type(dose_factor_t), parameter :: U_234   = dose_factor_t("U-234", 2.445E5, "y", 0.04, 3.7E-07, 0.02, &
    & [ 1.3E-07, 8.8E-08, 7.4E-08, 7.4E-08], 0.02, 4.9E-08 )
  type(dose_factor_t), parameter :: U_235   = dose_factor_t("U-235", 703.8E6, "y", 0.04, 3.5E-07, 0.02, &
    & [ 1.3E-07, 8.5E-08, 7.1E-08, 7.0E-08], 0.02, 4.7E-08 )
  type(dose_factor_t), parameter :: U_236   = dose_factor_t("U-236", 2.3415E7, "y", 0.04, 3.5E-07, 0.02, &
    & [ 1.3E-07, 8.4E-08, 7.0E-08, 7.0E-08], 0.02, 4.7E-08 )
  type(dose_factor_t), parameter :: U_237   = dose_factor_t("U-237", 6.75, "d", 0.04, 8.3E-09, 0.02, &
    & [ 5.4E-09, 2.8E-09, 1.6E-09, 9.5E-10], 0.02, 7.6E-10 )
  type(dose_factor_t), parameter :: U_238   = dose_factor_t("U-238", 4.468E9, "y", 0.04, 3.4E-07, 0.02, &
    & [ 1.2E-07, 8.0E-08, 6.8E-08, 6.7E-08], 0.02, 4.5E-08 )
  type(dose_factor_t), parameter :: U_239   = dose_factor_t("U-239", 23.54, "m", 0.04, 3.4E-10, 0.02, &
    & [ 1.9E-10, 9.3E-11, 5.4E-11, 3.5E-11], 0.02, 2.7E-11 )
  type(dose_factor_t), parameter :: U_240   = dose_factor_t("U-240", 14.1, "h", 0.04, 1.3E-08, 0.02, &
    & [ 8.1E-09, 4.1E-09, 2.4E-09, 1.4E-09], 0.02, 1.1E-09 )
  type(dose_factor_t), parameter :: Np_232  = dose_factor_t("Np-232", 14.7, "m", 0.005, 8.7E-11, 0.0005, &
    & [ 5.1E-11, 2.7E-11, 1.7E-11, 1.2E-11], 0.0005, 9.7E-12 )
  type(dose_factor_t), parameter :: Np_233  = dose_factor_t("Np-233", 36.2, "m", 0.005, 2.1E-11, 0.0005, &
    & [ 1.3E-11, 6.6E-12, 4.0E-12, 2.8E-12], 0.0005, 2.2E-12 )
  type(dose_factor_t), parameter :: Np_234  = dose_factor_t("Np-234", 4.4, "d", 0.005, 6.2E-09, 0.0005, &
    & [ 4.4E-09, 2.4E-09, 1.6E-09, 1.0E-09], 0.0005, 8.1E-10 )
  type(dose_factor_t), parameter :: Np_235  = dose_factor_t("Np-235", 396.1, "d", 0.005, 7.1E-10, 0.0005, &
    & [ 4.1E-10, 2.0E-10, 1.2E-10, 6.8E-11], 0.0005, 5.3E-11 )
  type(dose_factor_t), parameter :: Np_236  = dose_factor_t("Np-236", 115E3, "y", 0.005, 1.9E-07, 0.0005, &
    & [ 2.4E-08, 1.8E-08, 1.8E-08, 1.8E-08], 0.0005, 1.7E-08 )
  type(dose_factor_t), parameter :: Np_236m = dose_factor_t("Np-236m", 22.5, "h", 0.005, 2.5E-09, 0.0005, &
    & [ 1.3E-09, 6.6E-10, 4.0E-10, 2.4E-10], 0.0005, 1.9E-10 )
  type(dose_factor_t), parameter :: Np_237  = dose_factor_t("Np-237", 2.14E6, "y", 0.005, 2.0E-06, 0.0005, &
    & [ 2.1E-07, 1.4E-07, 1.1E-07, 1.1E-07], 0.0005, 1.1E-07 )
  type(dose_factor_t), parameter :: Np_238  = dose_factor_t("Np-238", 2.117, "d", 0.005, 9.5E-09, 0.0005, &
    & [ 6.2E-09, 3.2E-09, 1.9E-09, 1.1E-09], 0.0005, 9.1E-10 )
  type(dose_factor_t), parameter :: Np_239  = dose_factor_t("Np-239", 2.355, "d", 0.005, 8.9E-09, 0.0005, &
    & [ 5.7E-09, 2.9E-09, 1.7E-09, 1.0E-09], 0.0005, 8.0E-10 )
  type(dose_factor_t), parameter :: Np_240  = dose_factor_t("Np-240", 65, "m", 0.005, 8.7E-10, 0.0005, &
    & [ 5.2E-10, 2.6E-10, 1.6E-10, 1.0E-10], 0.0005, 8.2E-11 )
  type(dose_factor_t), parameter :: Pu_234  = dose_factor_t("Pu-234", 8.8, "h", 0.005, 2.1E-09, 0.0005, &
    & [ 1.1E-09, 5.5E-10, 3.3E-10, 2.0E-10], 0.0005, 1.6E-10 )
  type(dose_factor_t), parameter :: Pu_235  = dose_factor_t("Pu-235", 25.3, "m", 0.005, 2.2E-11, 0.0005, &
    & [ 1.3E-11, 6.5E-12, 3.9E-12, 2.7E-12], 0.0005, 2.1E-12 )
  type(dose_factor_t), parameter :: Pu_236  = dose_factor_t("Pu-236", 2.851, "y", 0.005, 2.1E-06, 0.0005, &
    & [ 2.2E-07, 1.4E-07, 1.0E-07, 8.5E-08], 0.0005, 8.7E-08 )
  type(dose_factor_t), parameter :: Pu_237  = dose_factor_t("Pu-237", 45.3, "d", 0.005, 1.1E-09, 0.0005, &
    & [ 6.9E-10, 3.6E-10, 2.2E-10, 1.3E-10], 0.0005, 1.0E-10 )
  type(dose_factor_t), parameter :: Pu_238  = dose_factor_t("Pu-238", 87.74, "y", 0.005, 4.0E-06, 0.0005, &
    & [ 4.0E-07, 3.1E-07, 2.4E-07, 2.2E-07], 0.0005, 2.3E-07 )
  type(dose_factor_t), parameter :: Pu_239  = dose_factor_t("Pu-239", 24065, "y", 0.005, 4.2E-06, 0.0005, &
    & [ 4.2E-07, 3.3E-07, 2.7E-07, 2.4E-07], 0.0005, 2.5E-07 )
  type(dose_factor_t), parameter :: Pu_240  = dose_factor_t("Pu-240", 6537, "y", 0.005, 4.2E-06, 0.0005, &
    & [ 4.2E-07, 3.3E-07, 2.7E-07, 2.4E-07], 0.0005, 2.5E-07 )
  type(dose_factor_t), parameter :: Pu_241  = dose_factor_t("Pu-241", 14.4, "y", 0.005, 5.6E-08, 0.0005, &
    & [ 5.7E-09, 5.5E-09, 5.1E-09, 4.8E-09], 0.0005, 4.8E-09 )
  type(dose_factor_t), parameter :: Pu_242  = dose_factor_t("Pu-242", 3.763E5, "y", 0.005, 4.0E-06, 0.0005, &
    & [ 4.0E-07, 3.2E-07, 2.6E-07, 2.3E-07], 0.0005, 2.4E-07 )
  type(dose_factor_t), parameter :: Pu_243  = dose_factor_t("Pu-243", 4.956, "h", 0.005, 1.0E-09, 0.0005, &
    & [ 6.2E-10, 3.1E-10, 1.8E-10, 1.1E-10], 0.0005, 8.5E-11 )
  type(dose_factor_t), parameter :: Pu_244  = dose_factor_t("Pu-244", 8.26E7, "y", 0.005, 4.0E-06, 0.0005, &
    & [ 4.1E-07, 3.2E-07, 2.6E-07, 2.3E-07], 0.0005, 2.4E-07 )
  type(dose_factor_t), parameter :: Pu_245  = dose_factor_t("Pu-245", 10.5, "h", 0.005, 8.0E-09, 0.0005, &
    & [ 5.1E-09, 2.6E-09, 1.5E-09, 8.9E-10], 0.0005, 7.2E-10 )
  type(dose_factor_t), parameter :: Pu_246  = dose_factor_t("Pu-246", 10.85, "d", 0.005, 3.6E-08, 0.0005, &
    & [ 2.3E-08, 1.2E-08, 7.1E-09, 4.1E-09], 0.0005, 3.3E-09 )
  type(dose_factor_t), parameter :: Am_237  = dose_factor_t("Am-237", 73.0, "m", 0.005, 1.7E-10, 0.0005, &
    & [ 1.0E-10, 5.5E-11, 3.3E-11, 2.2E-11], 0.0005, 1.8E-11 )
  type(dose_factor_t), parameter :: Am_238  = dose_factor_t("Am-238", 98, "m", 0.005, 2.5E-10, 0.0005, &
    & [ 1.6E-10, 9.1E-11, 5.9E-11, 4.0E-11], 0.0005, 3.2E-11 )
  type(dose_factor_t), parameter :: Am_239  = dose_factor_t("Am-239", 11.9, "h", 0.005, 2.6E-09, 0.0005, &
    & [ 1.7E-09, 8.4E-10, 5.1E-10, 3.0E-10], 0.0005, 2.4E-10 )
  type(dose_factor_t), parameter :: Am_240  = dose_factor_t("Am-240", 50.8, "h", 0.005, 4.7E-09, 0.0005, &
    & [ 3.3E-09, 1.8E-09, 1.2E-09, 7.3E-10], 0.0005, 5.8E-10 )
  type(dose_factor_t), parameter :: Am_241  = dose_factor_t("Am-241", 432.2, "y", 0.005, 3.7E-06, 0.0005, &
    & [ 3.7E-07, 2.7E-07, 2.2E-07, 2.0E-07], 0.0005, 2.0E-07 )
  type(dose_factor_t), parameter :: Am_242  = dose_factor_t("Am-242", 16.02, "h", 0.005, 5.0E-09, 0.0005, &
    & [ 2.2E-09, 1.1E-09, 6.4E-10, 3.7E-10], 0.0005, 3.0E-10 )
  type(dose_factor_t), parameter :: Am_242m = dose_factor_t("Am-242m", 152, "y", 0.005, 3.1E-06, 0.0005, &
    & [ 3.0E-07, 2.3E-07, 2.0E-07, 1.9E-07], 0.0005, 1.9E-07 )
  type(dose_factor_t), parameter :: Am_243  = dose_factor_t("Am-243", 7380, "y", 0.005, 3.6E-06, 0.0005, &
    & [ 3.7E-07, 2.7E-07, 2.2E-07, 2.0E-07], 0.0005, 2.0E-07 )
  type(dose_factor_t), parameter :: Am_244  = dose_factor_t("Am-244", 10.1, "h", 0.005, 4.9E-09, 0.0005, &
    & [ 3.1E-09, 1.6E-09, 9.6E-10, 5.8E-10], 0.0005, 4.6E-10 )
  type(dose_factor_t), parameter :: Am_244m = dose_factor_t("Am-244m", 26, "m", 0.005, 3.7E-10, 0.0005, &
    & [ 2.0E-10, 9.6E-11, 5.5E-11, 3.7E-11], 0.0005, 2.9E-11 )
  type(dose_factor_t), parameter :: Am_245  = dose_factor_t("Am-245", 2.05, "h", 0.005, 6.8E-10, 0.0005, &
    & [ 4.5E-10, 2.2E-10, 1.3E-10, 7.9E-11], 0.0005, 6.2E-11 )
  type(dose_factor_t), parameter :: Am_246  = dose_factor_t("Am-246", 39, "m", 0.005, 6.7E-10, 0.0005, &
    & [ 3.8E-10, 1.9E-10, 1.1E-10, 7.3E-11], 0.0005, 5.8E-11 )
  type(dose_factor_t), parameter :: Am_246m = dose_factor_t("Am-246m", 25.0, "m", 0.005, 3.9E-10, 0.0005, &
    & [ 2.2E-10, 1.1E-10, 6.4E-11, 4.4E-11], 0.0005, 3.4E-11 )
  type(dose_factor_t), parameter :: Cm_238  = dose_factor_t("Cm-238", 2.4, "h", 0.005, 7.8E-10, 0.0005, &
    & [ 4.9E-10, 2.6E-10, 1.6E-10, 1.0E-10], 0.0005, 8.0E-11 )
  type(dose_factor_t), parameter :: Cm_240  = dose_factor_t("Cm-240", 27, "d", 0.005, 2.2E-07, 0.0005, &
    & [ 4.8E-08, 2.5E-08, 1.5E-08, 9.2E-09], 0.0005, 7.6E-09 )
  type(dose_factor_t), parameter :: Cm_241  = dose_factor_t("Cm-241", 32.8, "d", 0.005, 1.1E-08, 0.0005, &
    & [ 5.7E-09, 3.0E-09, 1.9E-09, 1.1E-09], 0.0005, 9.1E-10 )
  type(dose_factor_t), parameter :: Cm_242  = dose_factor_t("Cm-242", 162.8, "d", 0.005, 5.9E-07, 0.0005, &
    & [ 7.6E-08, 3.9E-08, 2.4E-08, 1.5E-08], 0.0005, 1.2E-08 )
  type(dose_factor_t), parameter :: Cm_243  = dose_factor_t("Cm-243", 28.5, "y", 0.005, 3.2E-06, 0.0005, &
    & [ 3.3E-07, 2.2E-07, 1.6E-07, 1.4E-07], 0.0005, 1.5E-07 )
  type(dose_factor_t), parameter :: Cm_244  = dose_factor_t("Cm-244", 18.11, "y", 0.005, 2.9E-06, 0.0005, &
    & [ 2.9E-07, 1.9E-07, 1.4E-07, 1.2E-07], 0.0005, 1.2E-07 )
  type(dose_factor_t), parameter :: Cm_245  = dose_factor_t("Cm-245", 8500, "y", 0.005, 3.7E-06, 0.0005, &
    & [ 3.7E-07, 2.8E-07, 2.3E-07, 2.1E-07], 0.0005, 2.1E-07 )
  type(dose_factor_t), parameter :: Cm_246  = dose_factor_t("Cm-246", 4730, "y", 0.005, 3.7E-06, 0.0005, &
    & [ 3.7E-07, 2.8E-07, 2.2E-07, 2.1E-07], 0.0005, 2.1E-07 )
  type(dose_factor_t), parameter :: Cm_247  = dose_factor_t("Cm-247", 1.56E7, "y", 0.005, 3.4E-06, 0.0005, &
    & [ 3.5E-07, 2.6E-07, 2.1E-07, 1.9E-07], 0.0005, 1.9E-07 )
  type(dose_factor_t), parameter :: Cm_248  = dose_factor_t("Cm-248", 3.39E5, "y", 0.005, 1.4E-05, 0.0005, &
    & [ 1.4E-06, 1.0E-06, 8.4E-07, 7.7E-07], 0.0005, 7.7E-07 )
  type(dose_factor_t), parameter :: Cm_249  = dose_factor_t("Cm-249", 64.15, "m", 0.005, 3.9E-10, 0.0005, &
    & [ 2.2E-10, 1.1E-10, 6.1E-11, 4.0E-11], 0.0005, 3.1E-11 )
  type(dose_factor_t), parameter :: Cm_250  = dose_factor_t("Cm-250", 6900, "y", 0.005, 7.8E-05, 0.0005, &
    & [ 8.2E-06, 6.0E-06, 4.9E-06, 4.4E-06], 0.0005, 4.4E-06 )
  type(dose_factor_t), parameter :: Bk_245  = dose_factor_t("Bk-245", 4.94, "d", 0.005, 6.1E-09, 0.0005, &
    & [ 3.9E-09, 2.0E-09, 1.2E-09, 7.2E-10], 0.0005, 5.7E-10 )
  type(dose_factor_t), parameter :: Bk_246  = dose_factor_t("Bk-246", 1.83, "d", 0.005, 3.7E-09, 0.0005, &
    & [ 2.6E-09, 1.4E-09, 9.4E-10, 6.0E-10], 0.0005, 4.8E-10 )
  type(dose_factor_t), parameter :: Bk_247  = dose_factor_t("Bk-247", 1380, "y", 0.005, 8.9E-06, 0.0005, &
    & [ 8.6E-07, 6.3E-07, 4.6E-07, 3.8E-07], 0.0005, 3.5E-07 )
  type(dose_factor_t), parameter :: Bk_249  = dose_factor_t("Bk-249", 320, "d", 0.005, 2.2E-08, 0.0005, &
    & [ 2.9E-09, 1.9E-09, 1.4E-09, 1.1E-09], 0.0005, 9.7E-10 )
  type(dose_factor_t), parameter :: Bk_250  = dose_factor_t("Bk-250", 3.222, "h", 0.005, 1.5E-09, 0.0005, &
    & [ 8.5E-10, 4.4E-10, 2.7E-10, 1.7E-10], 0.0005, 1.4E-10 )
  type(dose_factor_t), parameter :: Cf_244  = dose_factor_t("Cf-244", 19.4, "m", 0.005, 9.8E-10, 0.0005, &
    & [ 4.8E-10, 2.4E-10, 1.3E-10, 8.9E-11], 0.0005, 7.0E-11 )
  type(dose_factor_t), parameter :: Cf_246  = dose_factor_t("Cf-246", 35.7, "h", 0.005, 5.0E-08, 0.0005, &
    & [ 2.4E-08, 1.2E-08, 7.3E-09, 4.1E-09], 0.0005, 3.3E-09 )
  type(dose_factor_t), parameter :: Cf_248  = dose_factor_t("Cf-248", 333.5, "d", 0.005, 1.5E-06, 0.0005, &
    & [ 1.6E-07, 9.9E-08, 6.0E-08, 3.3E-08], 0.0005, 2.8E-08 )
  type(dose_factor_t), parameter :: Cf_249  = dose_factor_t("Cf-249", 350.6, "y", 0.005, 9.0E-06, 0.0005, &
    & [ 8.7E-07, 6.4E-07, 4.7E-07, 3.8E-07], 0.0005, 3.5E-07 )
  type(dose_factor_t), parameter :: Cf_250  = dose_factor_t("Cf-250", 13.08, "y", 0.005, 5.7E-06, 0.0005, &
    & [ 5.5E-07, 3.7E-07, 2.3E-07, 1.7E-07], 0.0005, 1.6E-07 )
  type(dose_factor_t), parameter :: Cf_251  = dose_factor_t("Cf-251", 898, "y", 0.005, 9.1E-06, 0.0005, &
    & [ 8.8E-07, 6.5E-07, 4.7E-07, 3.9E-07], 0.0005, 3.6E-07 )
  type(dose_factor_t), parameter :: Cf_252  = dose_factor_t("Cf-252", 2.638, "y", 0.005, 5.0E-06, 0.0005, &
    & [ 5.1E-07, 3.2E-07, 1.9E-07, 1.0E-07], 0.0005, 9.0E-08 )
  type(dose_factor_t), parameter :: Cf_253  = dose_factor_t("Cf-253", 17.81, "d", 0.005, 1.0E-07, 0.0005, &
    & [ 1.1E-08, 6.0E-09, 3.7E-09, 1.8E-09], 0.0005, 1.4E-09 )
  type(dose_factor_t), parameter :: Cf_254  = dose_factor_t("Cf-254", 60.5, "d", 0.005, 1.1E-05, 0.0005, &
    & [ 2.6E-06, 1.4E-06, 8.4E-07, 5.0E-07], 0.0005, 4.0E-07 )
  type(dose_factor_t), parameter :: Es_250m = dose_factor_t("Es-250m", 2.1, "h", 0.005, 2.3E-10, 0.0005, &
    & [ 9.9E-11, 5.7E-11, 3.7E-11, 2.6E-11], 0.0005, 2.1E-11 )
  type(dose_factor_t), parameter :: Es_251  = dose_factor_t("Es-251", 33, "h", 0.005, 1.9E-09, 0.0005, &
    & [ 1.2E-09, 6.1E-10, 3.7E-10, 2.2E-10], 0.0005, 1.7E-10 )
  type(dose_factor_t), parameter :: Es_253  = dose_factor_t("Es-253", 20.47, "d", 0.005, 1.7E-07, 0.0005, &
    & [ 4.5E-08, 2.3E-08, 1.4E-08, 7.6E-09], 0.0005, 6.1E-09 )
  type(dose_factor_t), parameter :: Es_254  = dose_factor_t("Es-254", 275.7, "d", 0.005, 1.4E-06, 0.0005, &
    & [ 1.6E-07, 9.8E-08, 6.0E-08, 3.3E-08], 0.0005, 2.8E-08 )
  type(dose_factor_t), parameter :: Es_254m = dose_factor_t("Es-254m", 39.3, "h", 0.005, 5.7E-08, 0.0005, &
    & [ 3.0E-08, 1.5E-08, 9.1E-09, 5.2E-09], 0.0005, 4.2E-09 )
  type(dose_factor_t), parameter :: Fm_252  = dose_factor_t("Fm-252", 22.7, "h", 0.005, 3.8E-08, 0.0005, &
    & [ 2.0E-08, 9.9E-09, 5.9E-09, 3.3E-09], 0.0005, 2.7E-09 )
  type(dose_factor_t), parameter :: Fm_253  = dose_factor_t("Fm-253", 3.00, "d", 0.005, 2.5E-08, 0.0005, &
    & [ 6.7E-09, 3.4E-09, 2.1E-09, 1.1E-09], 0.0005, 9.1E-10 )
  type(dose_factor_t), parameter :: Fm_254  = dose_factor_t("Fm-254", 3.240, "h", 0.005, 5.6E-09, 0.0005, &
    & [ 3.2E-09, 1.6E-09, 9.3E-10, 5.6E-10], 0.0005, 4.4E-10 )
  type(dose_factor_t), parameter :: Fm_255  = dose_factor_t("Fm-255", 20.07, "h", 0.005, 3.3E-08, 0.0005, &
    & [ 1.9E-08, 9.5E-09, 5.6E-09, 3.2E-09], 0.0005, 2.5E-09 )
  type(dose_factor_t), parameter :: Fm_257  = dose_factor_t("Fm-257", 100.5, "d", 0.005, 9.8E-07, 0.0005, &
    & [ 1.1E-07, 6.5E-08, 4.0E-08, 1.9E-08], 0.0005, 1.5E-08 )
  type(dose_factor_t), parameter :: Md_257  = dose_factor_t("Md-257", 5.2, "h", 0.005, 3.1E-09, 0.0005, &
    & [ 8.8E-10, 4.5E-10, 2.7E-10, 1.5E-10], 0.0005, 1.2E-10 )
  type(dose_factor_t), parameter :: Md_258  = dose_factor_t("Md-258", 55, "d", 0.005, 6.3E-07, 0.0005, &
    & [ 8.9E-08, 5.0E-08, 3.0E-08, 1.6E-08], 0.0005, 1.3E-08 )

  type(dose_factor_t) :: Dose_Factors(737)

  data dose_factors(  1: 10) / H_3    , Be_7   , Be_10  , C_11   , C_14   , F_18   , Na_22  , Na_24  , Mg_28  , Al_26   /
  data dose_factors( 11: 20) / Si_31  , Si_32  , P_32   , P_33   , S_35   , Cl_36  , Cl_38  , Cl_39  , K_40   , K_42    /
  data dose_factors( 21: 30) / K_43   , K_44   , K_45   , Ca_41  , Ca_45  , Ca_47  , Sc_43  , Sc_44  , Sc_44m , Sc_46   /
  data dose_factors( 31: 40) / Sc_47  , Sc_48  , Sc_49  , Ti_44  , Ti_45  , V_47   , V_48   , V_49   , Cr_48  , Cr_49   /
  data dose_factors( 41: 50) / Cr_51  , Mn_51  , Mn_52  , Mn_52m , Mn_53  , Mn_54  , Mn_56  , Fe_52  , Fe_55  , Fe_59   /
  data dose_factors( 51: 60) / Fe_60  , Co_55  , Co_56  , Co_57  , Co_58  , Co_58m , Co_60  , Co_60m , Co_61  , Co_62m  /
  data dose_factors( 61: 70) / Ni_56  , Ni_57  , Ni_59  , Ni_63  , Ni_65  , Ni_66  , Cu_60  , Cu_61  , Cu_64  , Cu_67   /
  data dose_factors( 71: 80) / Zn_62  , Zn_63  , Zn_65  , Zn_69  , Zn_69m , Zn_71m , Zn_72  , Ga_65  , Ga_66  , Ga_67   /
  data dose_factors( 81: 90) / Ga_68  , Ga_70  , Ga_72  , Ga_73  , Ge_66  , Ge_67  , Ge_68  , Ge_69  , Ge_71  , Ge_75   /
  data dose_factors( 91:100) / Ge_77  , Ge_78  , As_69  , As_70  , As_71  , As_72  , As_73  , As_74  , As_76  , As_77   /
  data dose_factors(101:110) / As_78  , Se_70  , Se_73  , Se_73m , Se_75  , Se_79  , Se_81  , Se_81m , Se_83  , Br_74   /
  data dose_factors(111:120) / Br_74m , Br_75  , Br_76  , Br_77  , Br_80  , Br_80m , Br_82  , Br_83  , Br_84  , Rb_79   /
  data dose_factors(121:130) / Rb_81  , Rb_81m , Rb_82m , Rb_83  , Rb_84  , Rb_86  , Rb_87  , Rb_88  , Rb_89  , Sr_80   /
  data dose_factors(131:140) / Sr_81  , Sr_82  , Sr_83  , Sr_85  , Sr_85m , Sr_87m , Sr_89  , Sr_90  , Sr_91  , Sr_92   /
  data dose_factors(141:150) / Y_86   , Y_86m  , Y_87   , Y_88   , Y_90   , Y_90m  , Y_91   , Y_91m  , Y_92   , Y_93    /
  data dose_factors(151:160) / Y_94   , Y_95   , Zr_86  , Zr_88  , Zr_89  , Zr_93  , Zr_95  , Zr_97  , Nb_88  , Nb_89   /
  data dose_factors(161:170) / Nb_89m , Nb_90  , Nb_93m , Nb_94  , Nb_95  , Nb_95m , Nb_96  , Nb_97  , Nb_98m , Mo_90   /
  data dose_factors(171:180) / Mo_93  , Mo_93m , Mo_99  , Mo_101 , Tc_93  , Tc_93m , Tc_94  , Tc_94m , Tc_95  , Tc_95m  /
  data dose_factors(181:190) / Tc_96  , Tc_96m , Tc_97  , Tc_97m , Tc_98  , Tc_99  , Tc_99m , Tc_101 , Tc_104 , Ru_94   /
  data dose_factors(191:200) / Ru_97  , Ru_103 , Ru_105 , Ru_106 , Rh_99  , Rh_99m , Rh_100 , Rh_101 , Rh_101m, Rh_102m /
  data dose_factors(201:210) / Rh_102 , Rh_103m, Rh_105 , Rh_106m, Rh_107 , Pd_100 , Pd_101 , Pd_103 , Pd_107 , Pd_109  /
  data dose_factors(211:220) / Ag_102 , Ag_103 , Ag_104 , Ag_104m, Ag_105 , Ag_106 , Ag_106m, Ag_108m, Ag_110m, Ag_111  /
  data dose_factors(221:230) / Ag_112 , Ag_115 , Cd_104 , Cd_107 , Cd_109 , Cd_113 , Cd_113m, Cd_115 , Cd_115m, Cd_117  /
  data dose_factors(231:240) / Cd_117m, In_109 , In_110 , In_110m, In_111 , In_112 , In_113m, In_114m, In_115 , In_115m /
  data dose_factors(241:250) / In_116m, In_117 , In_117m, In_119m, Sn_110 , Sn_111 , Sn_113 , Sn_117m, Sn_119m, Sn_121  /
  data dose_factors(251:260) / Sn_121m, Sn_123 , Sn_123m, Sn_125 , Sn_126 , Sn_127 , Sn_128 , Sb_115 , Sb_116 , Sb_116m /
  data dose_factors(261:270) / Sb_117 , Sb_118m, Sb_119 , Sb_120m, Sb_120 , Sb_122 , Sb_124 , Sb_124n, Sb_125 , Sb_126  /
  data dose_factors(271:280) / Sb_126m, Sb_127 , Sb_128 , Sb_128m, Sb_129 , Sb_130 , Sb_131 , Te_116 , Te_121 , Te_121m /
  data dose_factors(281:290) / Te_123 , Te_123m, Te_125m, Te_127 , Te_127m, Te_129 , Te_129m, Te_131 , Te_131m, Te_132  /
  data dose_factors(291:300) / Te_133 , Te_133m, Te_134 , I_120  , I_120m , I_121  , I_123  , I_124  , I_125  , I_126   /
  data dose_factors(301:310) / I_128  , I_129  , I_130  , I_131  , I_132  , I_132m , I_133  , I_134  , I_135  , Cs_125  /
  data dose_factors(311:320) / Cs_127 , Cs_129 , Cs_130 , Cs_131 , Cs_132 , Cs_134 , Cs_134m, Cs_135 , Cs_135m, Cs_136  /
  data dose_factors(321:330) / Cs_137 , Cs_138 , Ba_126 , Ba_128 , Ba_131 , Ba_131m, Ba_133 , Ba_133m, Ba_135m, Ba_139  /
  data dose_factors(331:340) / Ba_140 , Ba_141 , Ba_142 , La_131 , La_132 , La_135 , La_137 , La_138 , La_140 , La_141  /
  data dose_factors(341:350) / La_142 , La_143 , Ce_134 , Ce_135 , Ce_137 , Ce_137m, Ce_139 , Ce_141 , Ce_143 , Ce_144  /
  data dose_factors(351:360) / Pr_136 , Pr_137 , Pr_138m, Pr_139 , Pr_142 , Pr_142m, Pr_143 , Pr_144 , Pr_145 , Pr_147  /
  data dose_factors(361:370) / Nd_136 , Nd_138 , Nd_139 , Nd_139m, Nd_141 , Nd_147 , Nd_149 , Nd_151 , Pm_141 , Pm_143  /
  data dose_factors(371:380) / Pm_144 , Pm_145 , Pm_146 , Pm_147 , Pm_148 , Pm_148m, Pm_149 , Pm_150 , Pm_151 , Sm_141  /
  data dose_factors(381:390) / Sm_141m, Sm_142 , Sm_145 , Sm_146 , Sm_147 , Sm_151 , Sm_153 , Sm_155 , Sm_156 , Eu_145  /
  data dose_factors(391:400) / Eu_146 , Eu_147 , Eu_148 , Eu_149 , Eu_150 , Eu_150m, Eu_152 , Eu_152m, Eu_154 , Eu_155  /
  data dose_factors(401:410) / Eu_156 , Eu_157 , Eu_158 , Gd_145 , Gd_146 , Gd_147 , Gd_148 , Gd_149 , Gd_151 , Gd_152  /
  data dose_factors(411:420) / Gd_153 , Gd_159 , Tb_147 , Tb_149 , Tb_150 , Tb_151 , Tb_153 , Tb_154 , Tb_155 , Tb_156  /
  data dose_factors(421:430) / Tb_156n, Tb_156m, Tb_157 , Tb_158 , Tb_160 , Tb_161 , Dy_155 , Dy_157 , Dy_159 , Dy_165  /
  data dose_factors(431:440) / Dy_166 , Ho_155 , Ho_157 , Ho_159 , Ho_161 , Ho_162 , Ho_162m, Ho_164 , Ho_164m, Ho_166  /
  data dose_factors(441:450) / Ho_166m, Ho_167 , Er_161 , Er_165 , Er_169 , Er_171 , Er_172 , Tm_162 , Tm_166 , Tm_167  /
  data dose_factors(451:460) / Tm_170 , Tm_171 , Tm_172 , Tm_173 , Tm_175 , Yb_162 , Yb_166 , Yb_167 , Yb_169 , Yb_175  /
  data dose_factors(461:470) / Yb_177 , Yb_178 , Lu_169 , Lu_170 , Lu_171 , Lu_172 , Lu_173 , Lu_174 , Lu_174m, Lu_176  /
  data dose_factors(471:480) / Lu_176m, Lu_177 , Lu_177m, Lu_178 , Lu_178m, Lu_179 , Hf_170 , Hf_172 , Hf_173 , Hf_175  /
  data dose_factors(481:490) / Hf_177m, Hf_178m, Hf_179m, Hf_180m, Hf_181 , Hf_182 , Hf_182m, Hf_183 , Hf_184 , Ta_172  /
  data dose_factors(491:500) / Ta_173 , Ta_174 , Ta_175 , Ta_176 , Ta_177 , Ta_178m, Ta_179 , Ta_180 , Ta_182 , Ta_182m /
  data dose_factors(501:510) / Ta_183 , Ta_184 , Ta_185 , Ta_186 , W_176  , W_177  , W_178  , W_179  , W_181  , W_185   /
  data dose_factors(511:520) / W_187  , W_188  , Re_177 , Re_178 , Re_181 , Re_182 , Re_182m, Re_184 , Re_184m, Re_186  /
  data dose_factors(521:530) / Re_186m, Re_187 , Re_188 , Re_188m, Re_189 , Os_180 , Os_181 , Os_182 , Os_185 , Os_189m /
  data dose_factors(531:540) / Os_191 , Os_191m, Os_193 , Os_194 , Ir_182 , Ir_184 , Ir_185 , Ir_186 , Ir_186m, Ir_187  /
  data dose_factors(541:550) / Ir_188 , Ir_189 , Ir_190 , Ir_190n, Ir_190m, Ir_192 , Ir_192n, Ir_193m, Ir_194 , Ir_194m /
  data dose_factors(551:560) / Ir_195 , Ir_195m, Pt_186 , Pt_188 , Pt_189 , Pt_191 , Pt_193 , Pt_193m, Pt_195m, Pt_197  /
  data dose_factors(561:570) / Pt_197m, Pt_199 , Pt_200 , Au_193 , Au_194 , Au_195 , Au_198 , Au_198m, Au_199 , Au_200  /
  data dose_factors(571:580) / Au_200m, Au_201 , Hg_193 , Hg_193m, Hg_194 , Hg_195 , Hg_195m, Hg_197 , Hg_197m, Hg_199m /
  data dose_factors(581:590) / Hg_203 , Tl_194 , Tl_194m, Tl_195 , Tl_197 , Tl_198 , Tl_198m, Tl_199 , Tl_200 , Tl_201  /
  data dose_factors(591:600) / Tl_202 , Tl_204 , Pb_195m, Pb_198 , Pb_199 , Pb_200 , Pb_201 , Pb_202 , Pb_202m, Pb_203  /
  data dose_factors(601:610) / Pb_205 , Pb_209 , Pb_210 , Pb_211 , Pb_212 , Pb_214 , Bi_200 , Bi_201 , Bi_202 , Bi_203  /
  data dose_factors(611:620) / Bi_205 , Bi_206 , Bi_207 , Bi_210 , Bi_210m, Bi_212 , Bi_213 , Bi_214 , Po_203 , Po_205  /
  data dose_factors(621:630) / Po_207 , Po_210 , At_207 , At_211 , Fr_222 , Fr_223 , Ra_223 , Ra_224 , Ra_225 , Ra_226  /
  data dose_factors(631:640) / Ra_227 , Ra_228 , Ac_224 , Ac_225 , Ac_226 , Ac_227 , Ac_228 , Th_226 , Th_227 , Th_228  /
  data dose_factors(641:650) / Th_229 , Th_230 , Th_231 , Th_232 , Th_234 , Pa_227 , Pa_228 , Pa_230 , Pa_231 , Pa_232  /
  data dose_factors(651:660) / Pa_233 , Pa_234 , U_230  , U_231  , U_232  , U_233  , U_234  , U_235  , U_236  , U_237   /
  data dose_factors(661:670) / U_238  , U_239  , U_240  , Np_232 , Np_233 , Np_234 , Np_235 , Np_236 , Np_236m, Np_237  /
  data dose_factors(671:680) / Np_238 , Np_239 , Np_240 , Pu_234 , Pu_235 , Pu_236 , Pu_237 , Pu_238 , Pu_239 , Pu_240  /
  data dose_factors(681:690) / Pu_241 , Pu_242 , Pu_243 , Pu_244 , Pu_245 , Pu_246 , Am_237 , Am_238 , Am_239 , Am_240  /
  data dose_factors(691:700) / Am_241 , Am_242 , Am_242m, Am_243 , Am_244 , Am_244m, Am_245 , Am_246 , Am_246m, Cm_238  /
  data dose_factors(701:710) / Cm_240 , Cm_241 , Cm_242 , Cm_243 , Cm_244 , Cm_245 , Cm_246 , Cm_247 , Cm_248 , Cm_249  /
  data dose_factors(711:720) / Cm_250 , Bk_245 , Bk_246 , Bk_247 , Bk_249 , Bk_250 , Cf_244 , Cf_246 , Cf_248 , Cf_249  /
  data dose_factors(721:730) / Cf_250 , Cf_251 , Cf_252 , Cf_253 , Cf_254 , Es_250m, Es_251 , Es_253 , Es_254 , Es_254m /
  data dose_factors(731:737) / Fm_252 , Fm_253 , Fm_254 , Fm_255 , Fm_257 , Md_257 , Md_258  /

end module Dose_Factors_m
