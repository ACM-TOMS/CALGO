module ENDF_m

  ! Isotope data from an ENDF file such as JEFF 3.11

  implicit NONE
  public

  type :: ENDF_t
    character(2) :: Sym ! Atomic symbol
    integer :: Matn     ! ENDF material number
    integer :: A        ! Atomic mass index = Z + Neutrons
    integer :: Z        ! Atomic number
    character :: M      ! 'm' or blank, for metastable
    real :: AMU         ! Mass in AMU = grams/mole
    real :: Energy      ! Total decay energy (eV) for heat calculation purposes
    real :: RFS         ! Excitation state of daughter
    real :: RTYP        ! ENDF decay mode, including multiple modes, for decay
                        ! with largest branching ratio
                        ! 0 => Gamma (not in MT 457) (for stable isotopes here)
                        ! 1 => Beta- decay
                        ! 2 => e.c. Beta +, electron capture and/or positron emission
                        ! 3 => Isomeric transition
                        ! 4 => Alpha
                        ! 5 => Neutrons (not delayed neutron decay)
                        ! 6 => Spontaneous fission
                        ! 7 => Proton emission
                        ! Multiple particle decay indicated by more digits, e.g.
                        ! 1.4 => Beta followed by alpha, as in N-16 decay.
    real :: STYP        ! Decay radiation of first RTYP decay mode
                        ! 0 => Gamma
                        ! 1 => Beta -
                        ! 2 => e.c., Beta +
                        ! 4 => Alpha
                        ! 5 => Neutrons
                        ! 6 => Spontaneous fission
                        ! 7 => Protons
                        ! 8 => Discrete electrons
                        ! 9 => X-rays and annihilation radiation (photons not
                        !      arising as transitions between nuclear states)
    real :: T_Half      ! Half life in seconds, zero for stable isotopes
  end type ENDF_t


!            Sym   Matn  A   Z   M    AMU       Energy       RFS  RTYP    STYP   T_Half
  type(ENDF_t), parameter :: Cu_81 = &   ! AMU from JENDL 2000
    & ENDF_t("Cu",9999, 81, 29, " ",  80.9650 , 0.0        , 0.0, 1.0000,  1.0, 632.0000e-09)
  type(ENDF_t), parameter :: Se_85m = &  ! T_Half from JAERI 37002270.pdf
    & ENDF_t("Se",9998, 85, 34, "m",  84.1927 , 0.0        , 0.0, 3.0000,  3.0, 19.00000e-00)
  type(ENDF_t), parameter :: Br_86m = &  ! T_Half from JAERI 37002270.pdf
    & ENDF_t("Br",9997, 86, 35, "m",  85.1807 , 0.0        , 0.0, 3.0000,  3.0, 4.500000e-00)
  type(ENDF_t), parameter :: Rh_109m = & ! T_Half from JAERI 37002270.pdf
    & ENDF_t("Rh",9996,109, 45, "m",  107.973 , 0.0        , 0.0, 3.0000,  3.0, 50.00000e-00)
  type(ENDF_t), parameter :: Rh_123 = &  ! AMU from JENDL 2000
    & ENDF_t("Rh",9995,123, 45, " ",  121.944 , 0.0        , 0.0, 1.0000,  1.0, 42.00000e-03)
  type(ENDF_t), parameter :: Pd_126 = &  ! AMU from JENDL 2000
    & ENDF_t("Pd",9994,126, 46, " ",  125.940 , 0.0        , 0.0, 1.0000,  1.0, 48.60000e-03)
  ! Eu_151 half life (1.46e+26) seconds, Energy, and Energy_D (1.905 MeV)
  ! RTYP (4), and STYP (4) added manually.
  ! From P. Belli at al, "Search for \alpha decay of natural Europium,"
  ! Nuclear Physics 789 (2007) pp 15-29.
  type(ENDF_t), parameter :: Eu_151  = &
    & ENDF_t("Eu",2184,151, 63, " ",  149.623 , 1.905E+6   , 0.0, 4.0000,  4.0, 5.364678e+25)
  type(ENDF_t), parameter :: Tb_162m = & ! T_Half from JAERI 37002270.pdf
    & ENDF_t("Tb",9993,162, 65, "m",  160.538 , 0.0        , 0.0, 3.0000,  3.0, 8.028000e+03)
  type(ENDF_t), parameter :: Tm_170m = & ! AMU extrapolated
    & ENDF_t("Tm",9992,170, 69, "m",  168.476 , 0.0        , 0.0, 3.0000,  3.0, 4.100000e-06)

  include 'ENDF_data.f9h'

end module ENDF_m
