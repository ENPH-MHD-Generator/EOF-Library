!------------------------------------------------------------------------------
! MHDPlasma – Fortran interface to the mhd_plasma library (plasma/, Rust)
!
! The seeded-plasma physics (Saha ionization of the seed, electron transport
! from the case's electron_collisions.dat, elastic and inelastic losses, the
! local electron temperature) lives in plasma/src. This module only declares
! it for Fortran: two interoperable types, the constants MHDSolve uses, and
! procedures with the names and arguments MHDSolve has always called, each a
! direct call into the library. Nothing here computes anything; new physics
! goes in the library, with its declaration here and in
! plasma/include/mhd_plasma.h.
!------------------------------------------------------------------------------
MODULE MHDPlasma
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: kBoltz, eMass, eCharge, KelvinPerEV, Loschmidt
  PUBLIC :: SeedPlasma_t, PlasmaTransport_t
  PUBLIC :: ReadCollisionData, SahaDensity, SeedDensityDerivative, PlasmaTransport, &
      ElasticLossCoefficients, InelasticLoss, ElectronTemperature

  ! The library's constants (plasma/src/constants.rs), for MHDSolve's unit conversions
  REAL(KIND=C_DOUBLE), PARAMETER :: kBoltz = 1.380649d-23
  REAL(KIND=C_DOUBLE), PARAMETER :: eMass = 9.1093837015d-31
  REAL(KIND=C_DOUBLE), PARAMETER :: eCharge = 1.602176634d-19
  REAL(KIND=C_DOUBLE), PARAMETER :: KelvinPerEV = eCharge / kBoltz
  REAL(KIND=C_DOUBLE), PARAMETER :: Loschmidt = 2.6867811d25

  !> The seed and the model options (SeedPlasma in plasma/src/physics.rs)
  TYPE, BIND(C) :: SeedPlasma_t
    REAL(KIND=C_DOUBLE) :: SeedFrac      ! seed atoms per heavy particle
    REAL(KIND=C_DOUBLE) :: ChiJ          ! seed ionization energy [J]
    REAL(KIND=C_DOUBLE) :: WeightRatio   ! g_ion / g_neutral of the seed
    REAL(KIND=C_DOUBLE) :: IonReducedMobility = 0.0_C_DOUBLE  ! m^2/(V s) at Loschmidt density
    LOGICAL(KIND=C_BOOL) :: Lorentz = .FALSE.        ! Lorentz electron transport (else drifting)
    LOGICAL(KIND=C_BOOL) :: ExcitedAtGas = .FALSE.   ! K(4p) population at the gas temperature
  END TYPE SeedPlasma_t

  !> Local transport state at one point (PlasmaTransport in plasma/src/physics.rs)
  TYPE, BIND(C) :: PlasmaTransport_t
    REAL(KIND=C_DOUBLE) :: ne = 0.0_C_DOUBLE           ! electron density, m^-3
    REAL(KIND=C_DOUBLE) :: LnLambda = 0.0_C_DOUBLE     ! Coulomb logarithm
    REAL(KIND=C_DOUBLE) :: Mu0 = 0.0_C_DOUBLE, MuP = 0.0_C_DOUBLE, MuH = 0.0_C_DOUBLE   ! electron mobilities
    REAL(KIND=C_DOUBLE) :: MuI = 0.0_C_DOUBLE, MuIP = 0.0_C_DOUBLE, MuIH = 0.0_C_DOUBLE ! ion mobilities
    REAL(KIND=C_DOUBLE) :: CollisionFrequency = 0.0_C_DOUBLE  ! e / (me Mu0), 1/s
  END TYPE PlasmaTransport_t

  INTERFACE
    FUNCTION mhd_plasma_read_collision_data(file_name, lorentz, message, message_len) &
        BIND(C, NAME='mhd_plasma_read_collision_data') RESULT(status)
      IMPORT :: C_CHAR, C_BOOL, C_SIZE_T, C_INT32_T
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: file_name(*)
      LOGICAL(KIND=C_BOOL), VALUE :: lorentz
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(KIND=C_SIZE_T), VALUE :: message_len
      INTEGER(KIND=C_INT32_T) :: status
    END FUNCTION

    FUNCTION mhd_plasma_saha_density(pl, t, ns) BIND(C, NAME='mhd_plasma_saha_density') RESULT(ne)
      IMPORT :: SeedPlasma_t, C_DOUBLE
      TYPE(SeedPlasma_t), INTENT(IN) :: pl
      REAL(KIND=C_DOUBLE), VALUE :: t, ns
      REAL(KIND=C_DOUBLE) :: ne
    END FUNCTION

    FUNCTION mhd_plasma_seed_density_derivative(pl, t, ns, ne) &
        BIND(C, NAME='mhd_plasma_seed_density_derivative') RESULT(dne)
      IMPORT :: SeedPlasma_t, C_DOUBLE
      TYPE(SeedPlasma_t), INTENT(IN) :: pl
      REAL(KIND=C_DOUBLE), VALUE :: t, ns, ne
      REAL(KIND=C_DOUBLE) :: dne
    END FUNCTION

    SUBROUTINE mhd_plasma_transport(pl, te, n, ne, bmag, tr) BIND(C, NAME='mhd_plasma_transport')
      IMPORT :: SeedPlasma_t, PlasmaTransport_t, C_DOUBLE
      TYPE(SeedPlasma_t), INTENT(IN) :: pl
      REAL(KIND=C_DOUBLE), VALUE :: te, n, ne, bmag
      TYPE(PlasmaTransport_t), INTENT(OUT) :: tr
    END SUBROUTINE

    SUBROUTINE mhd_plasma_elastic_loss_coefficients(te, n, tr, cn, cei) &
        BIND(C, NAME='mhd_plasma_elastic_loss_coefficients')
      IMPORT :: PlasmaTransport_t, C_DOUBLE
      REAL(KIND=C_DOUBLE), VALUE :: te, n
      TYPE(PlasmaTransport_t), INTENT(IN) :: tr
      REAL(KIND=C_DOUBLE), INTENT(OUT) :: cn, cei
    END SUBROUTINE

    FUNCTION mhd_plasma_inelastic_loss(pl, te, tg, ne, nk) &
        BIND(C, NAME='mhd_plasma_inelastic_loss') RESULT(loss)
      IMPORT :: SeedPlasma_t, C_DOUBLE
      TYPE(SeedPlasma_t), INTENT(IN) :: pl
      REAL(KIND=C_DOUBLE), VALUE :: te, tg, ne, nk
      REAL(KIND=C_DOUBLE) :: loss
    END FUNCTION

    FUNCTION mhd_plasma_electron_temperature(pl, tg, n, ns, e_par2, e_perp2, bmag, te_max, at_max) &
        BIND(C, NAME='mhd_plasma_electron_temperature') RESULT(te)
      IMPORT :: SeedPlasma_t, C_DOUBLE, C_BOOL
      TYPE(SeedPlasma_t), INTENT(IN) :: pl
      REAL(KIND=C_DOUBLE), VALUE :: tg, n, ns, e_par2, e_perp2, bmag, te_max
      LOGICAL(KIND=C_BOOL), INTENT(OUT) :: at_max
      REAL(KIND=C_DOUBLE) :: te
    END FUNCTION
  END INTERFACE

CONTAINS

  !> Read the collision tables once; a later call with the same file does nothing
  SUBROUTINE ReadCollisionData( FileName, Lorentz, ErrorMessage )
    CHARACTER(LEN=*), INTENT(IN) :: FileName
    LOGICAL(KIND=C_BOOL), INTENT(IN) :: Lorentz
    CHARACTER(LEN=*), INTENT(OUT) :: ErrorMessage
    CHARACTER(KIND=C_CHAR) :: Message(1024)
    INTEGER :: i

    ErrorMessage = ''
    IF (mhd_plasma_read_collision_data(TRIM(FileName)//C_NULL_CHAR, Lorentz, Message, &
        INT(SIZE(Message), C_SIZE_T)) == 0) RETURN
    DO i = 1, MIN(SIZE(Message), LEN(ErrorMessage))
      IF (Message(i) == C_NULL_CHAR) EXIT
      ErrorMessage(i:i) = Message(i)
    END DO
  END SUBROUTINE ReadCollisionData

  !> Saha electron density at electron temperature T and seed density nS
  FUNCTION SahaDensity( Pl, T, nS ) RESULT(ne)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=C_DOUBLE), INTENT(IN) :: T, nS
    REAL(KIND=C_DOUBLE) :: ne
    ne = mhd_plasma_saha_density(Pl, T, nS)
  END FUNCTION SahaDensity

  !> d(ne)/dT of the Saha electron density at fixed seed density nS
  FUNCTION SeedDensityDerivative( Pl, T, nS, ne ) RESULT(dne)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=C_DOUBLE), INTENT(IN) :: T, nS, ne
    REAL(KIND=C_DOUBLE) :: dne
    dne = mhd_plasma_seed_density_derivative(Pl, T, nS, ne)
  END FUNCTION SeedDensityDerivative

  !> Electron and ion mobilities at Te, heavy-particle density N, ne and |B|
  SUBROUTINE PlasmaTransport( Pl, Te, N, ne, Bmag, Tr )
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=C_DOUBLE), INTENT(IN) :: Te, N, ne, Bmag
    TYPE(PlasmaTransport_t), INTENT(OUT) :: Tr
    CALL mhd_plasma_transport(Pl, Te, N, ne, Bmag, Tr)
  END SUBROUTINE PlasmaTransport

  !> Elastic loss coefficients: the loss is (Cn + Cei) (Te - Tg), W/m^3
  SUBROUTINE ElasticLossCoefficients( Te, N, Tr, Cn, Cei )
    REAL(KIND=C_DOUBLE), INTENT(IN) :: Te, N
    TYPE(PlasmaTransport_t), INTENT(IN) :: Tr
    REAL(KIND=C_DOUBLE), INTENT(OUT) :: Cn, Cei
    CALL mhd_plasma_elastic_loss_coefficients(Te, N, Tr, Cn, Cei)
  END SUBROUTINE ElasticLossCoefficients

  !> Net K 4s-4p inelastic energy loss of the electrons, W/m^3
  FUNCTION InelasticLoss( Pl, Te, Tg, ne, nK ) RESULT(L)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=C_DOUBLE), INTENT(IN) :: Te, Tg, ne, nK
    REAL(KIND=C_DOUBLE) :: L
    L = mhd_plasma_inelastic_loss(Pl, Te, Tg, ne, nK)
  END FUNCTION InelasticLoss

  !> Local electron temperature; AtMax is set when it is capped at TeMax
  FUNCTION ElectronTemperature( Pl, Tg, N, nS, Epar2, Eperp2, Bmag, TeMax, AtMax ) RESULT(Te)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=C_DOUBLE), INTENT(IN) :: Tg, N, nS, Epar2, Eperp2, Bmag, TeMax
    LOGICAL, INTENT(OUT) :: AtMax
    REAL(KIND=C_DOUBLE) :: Te
    LOGICAL(KIND=C_BOOL) :: Capped
    Te = mhd_plasma_electron_temperature(Pl, Tg, N, nS, Epar2, Eperp2, Bmag, TeMax, Capped)
    AtMax = Capped
  END FUNCTION ElectronTemperature

END MODULE MHDPlasma
