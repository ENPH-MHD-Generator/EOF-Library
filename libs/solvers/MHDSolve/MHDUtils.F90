!/******************************************************************************
! * MHDUtils.F90 – Utility and logging for MHDSolve
! * Split from MHDSolve.F90 to reduce solver bloat.
! *
! * Organization:
! *   MHDUtils     – Pure math (Invert3x3), no Elmer types.
! *   MHDLog       – Simple logging (LogLine, LogSection, LogSectionEnd).
! *   MHDDiagnostics – Electrode circuit/potential/current diagnostics and
! *                    helpers (GetElementNormal, ComputeElectrodeAreas, etc.).
! *****************************************************************************/

MODULE MHDUtils
  IMPLICIT NONE
  ! Local kind to avoid conflicting with Elmer Types%dp
  INTEGER, PARAMETER :: r8 = KIND(1.0d0)

CONTAINS

  SUBROUTINE Invert3x3(A, Ainv, ok)
    REAL(KIND=r8), INTENT(IN)  :: A(3,3)
    REAL(KIND=r8), INTENT(OUT) :: Ainv(3,3)
    LOGICAL, INTENT(OUT) :: ok
    REAL(KIND=r8) :: det

    det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

    IF (ABS(det) < 1.0d-30) THEN
      ok = .FALSE.
      Ainv = 0.0_r8
      RETURN
    END IF
    ok = .TRUE.

    Ainv(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
    Ainv(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
    Ainv(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det

    Ainv(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
    Ainv(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
    Ainv(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det

    Ainv(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
    Ainv(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
    Ainv(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
  END SUBROUTINE Invert3x3

END MODULE MHDUtils


!------------------------------------------------------------------------------
! MHDPlasma – seeded-plasma physics (no Elmer types)
!
! Argon seeded with potassium. Only the seed ionizes (two-temperature Saha at
! the electron temperature); heavy particles are at the gas temperature.
! Electron collision data come from plasma_collisions (a Maxwellian EEDF at Te),
! tabulated by case preparation into electron_collisions.dat: electron-neutral
! momentum transfer, recoil energy losses and K(4p) excitation. Electron-ion
! (Coulomb) collisions use the NRL Coulomb logarithm. Ions carry current with a
! polarization mobility (ion slip).
!------------------------------------------------------------------------------
MODULE MHDPlasma
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: r8 = KIND(1.0d0)

  REAL(KIND=r8), PARAMETER :: kBoltz = 1.380649d-23
  REAL(KIND=r8), PARAMETER :: eMass = 9.1093837015d-31
  REAL(KIND=r8), PARAMETER :: eCharge = 1.602176634d-19
  REAL(KIND=r8), PARAMETER, PRIVATE :: hPlanck = 6.62607015d-34
  REAL(KIND=r8), PARAMETER, PRIVATE :: Pi = 3.14159265358979323846d0
  REAL(KIND=r8), PARAMETER :: KelvinPerEV = eCharge / kBoltz
  !> Loschmidt number, m^-3 (273.15 K, 101325 Pa): reference for reduced mobility
  REAL(KIND=r8), PARAMETER :: Loschmidt = 2.6867811d25
  !> NRL electron-ion collision frequency coefficient, m^3 eV^(3/2)/s
  REAL(KIND=r8), PARAMETER, PRIVATE :: NrlCoefficient = 2.91d-12

  TYPE SeedPlasma_t
    REAL(KIND=r8) :: SeedFrac      ! seed atoms per heavy particle
    REAL(KIND=r8) :: ChiJ          ! seed ionization energy [J]
    REAL(KIND=r8) :: WeightRatio   ! g_ion / g_neutral of the seed
    REAL(KIND=r8) :: IonReducedMobility = 0.0_r8  ! m^2/(V s) at Loschmidt density
    LOGICAL :: Lorentz = .FALSE.   ! Lorentz electron transport (else drifting)
    LOGICAL :: ExcitedAtGas = .FALSE.  ! K(4p) population at the gas temperature
  END TYPE SeedPlasma_t

  !> Tables from electron_collisions.dat, natural logarithms of the values
  TYPE CollisionData_t
    LOGICAL :: Loaded = .FALSE.
    CHARACTER(LEN=512) :: FileName = ''
    INTEGER :: nT = 0, nA = 0, nW = 0
    REAL(KIND=r8), ALLOCATABLE :: LnT(:), LnA(:), LnW(:)
    REAL(KIND=r8), ALLOCATABLE :: LnNuN(:), LnElastic(:), LnEi(:), LnKexc(:), LnKsup(:)
    REAL(KIND=r8), ALLOCATABLE :: LnMu(:,:,:), LnMuP(:,:,:), LnMuH(:,:,:)
    REAL(KIND=r8) :: ExcEnergy = 0.0_r8, ExcWeight = 1.0_r8
  END TYPE CollisionData_t

  TYPE(CollisionData_t), SAVE :: CD

  !> Local transport state at one point
  TYPE PlasmaTransport_t
    REAL(KIND=r8) :: ne = 0.0_r8           ! electron density, m^-3
    REAL(KIND=r8) :: LnLambda = 0.0_r8     ! Coulomb logarithm
    REAL(KIND=r8) :: Mu0 = 0.0_r8, MuP = 0.0_r8, MuH = 0.0_r8   ! electron mobilities, m^2/(V s)
    REAL(KIND=r8) :: MuI = 0.0_r8, MuIP = 0.0_r8, MuIH = 0.0_r8 ! ion mobilities (Hall opposite)
    REAL(KIND=r8) :: CollisionFrequency = 0.0_r8  ! e / (me Mu0), 1/s
  END TYPE PlasmaTransport_t

CONTAINS

  !> Read the collision tables once; a later call with the same file does nothing
  SUBROUTINE ReadCollisionData( FileName, Lorentz, ErrorMessage )
    CHARACTER(LEN=*), INTENT(IN) :: FileName
    LOGICAL, INTENT(IN) :: Lorentz
    CHARACTER(LEN=*), INTENT(OUT) :: ErrorMessage
    CHARACTER(LEN=4096) :: Line
    CHARACTER(LEN=64) :: BlockName
    INTEGER :: Unit, Count, ios
    REAL(KIND=r8), ALLOCATABLE :: Values(:)
    REAL(KIND=r8) :: Scalar

    ErrorMessage = ''
    IF (CD % Loaded .AND. TRIM(CD % FileName) == TRIM(FileName)) RETURN
    CALL FreeCollisionData()

    OPEN(NEWUNIT=Unit, FILE=TRIM(FileName), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
      ErrorMessage = 'Cannot open electron collision data '//TRIM(FileName)//'; rerun mhd prepare'
      RETURN
    END IF
    DO
      READ(Unit, '(A)', IOSTAT=ios) Line
      IF (ios /= 0) EXIT
      Line = ADJUSTL(Line)
      IF (LEN_TRIM(Line) == 0 .OR. Line(1:1) == '!') CYCLE
      READ(Line, *, IOSTAT=ios) BlockName, Count
      IF (ios /= 0) THEN
        ErrorMessage = 'Malformed block header in '//TRIM(FileName)//': '//TRIM(Line)
        EXIT
      END IF
      IF (TRIM(BlockName) == 'end') EXIT
      IF (ALLOCATED(Values)) DEALLOCATE(Values)
      ALLOCATE(Values(Count))
      READ(Unit, *, IOSTAT=ios) Values
      IF (ios /= 0) THEN
        ErrorMessage = 'Short block '//TRIM(BlockName)//' in '//TRIM(FileName)
        EXIT
      END IF
      Scalar = Values(1)
      SELECT CASE (TRIM(BlockName))
      CASE ('format')
        IF (NINT(Scalar) /= 1) ErrorMessage = 'Unsupported collision data format in '//TRIM(FileName)
      CASE ('theta')
        CD % nT = Count
        CD % LnT = LOG(Values)
      CASE ('momentum_frequency_n')
        CD % LnNuN = SafeLog(Values)
      CASE ('elastic_recoil_n')
        CD % LnElastic = SafeLog(Values)
      CASE ('ei_recoil_per_coulomb')
        CD % LnEi = SafeLog(Values)
      CASE ('k_excitation_K')
        CD % LnKexc = SafeLog(Values)
      CASE ('k_superelastic_K')
        CD % LnKsup = SafeLog(Values)
      CASE ('excitation_energy_K')
        CD % ExcEnergy = Scalar
      CASE ('excitation_weight_ratio_K')
        CD % ExcWeight = Scalar
      CASE ('coulomb_parameter')
        CD % nA = Count
        CD % LnA = LOG(Values)
      CASE ('magnetic_parameter')
        CD % nW = Count
        CD % LnW = LOG(Values)
      CASE ('mobility_n')
        CD % LnMu = Cube(Values)
      CASE ('pedersen_mobility_n')
        CD % LnMuP = Cube(Values)
      CASE ('hall_mobility_n')
        CD % LnMuH = Cube(Values)
      END SELECT
      IF (LEN_TRIM(ErrorMessage) > 0) EXIT
    END DO
    CLOSE(Unit)
    IF (LEN_TRIM(ErrorMessage) > 0) RETURN

    IF (CD % nT < 2 .OR. .NOT. ALLOCATED(CD % LnNuN) .OR. .NOT. ALLOCATED(CD % LnElastic) &
        .OR. .NOT. ALLOCATED(CD % LnEi) .OR. .NOT. ALLOCATED(CD % LnKexc) &
        .OR. .NOT. ALLOCATED(CD % LnKsup)) THEN
      ErrorMessage = 'Incomplete electron collision data in '//TRIM(FileName)
      RETURN
    END IF
    IF (Lorentz .AND. .NOT. (ALLOCATED(CD % LnMu) .AND. ALLOCATED(CD % LnMuP) &
        .AND. ALLOCATED(CD % LnMuH))) THEN
      ErrorMessage = 'Lorentz transport needs the mobility tables; regenerate '//TRIM(FileName)// &
          ' with electron_transport_model: lorentz'
      RETURN
    END IF
    CD % FileName = FileName
    CD % Loaded = .TRUE.

  CONTAINS

    FUNCTION SafeLog(V) RESULT(L)
      REAL(KIND=r8), INTENT(IN) :: V(:)
      REAL(KIND=r8) :: L(SIZE(V))
      L = LOG(MAX(V, 1.0d-300))
    END FUNCTION SafeLog

    FUNCTION Cube(V) RESULT(L)
      REAL(KIND=r8), INTENT(IN) :: V(:)
      REAL(KIND=r8), ALLOCATABLE :: L(:,:,:)
      IF (CD % nT * CD % nA * CD % nW /= SIZE(V)) THEN
        ErrorMessage = 'Mobility table size does not match its grids in '//TRIM(FileName)
        ALLOCATE(L(1,1,1)); L = 0.0_r8
        RETURN
      END IF
      L = RESHAPE(SafeLog(V), (/ CD % nT, CD % nA, CD % nW /))
    END FUNCTION Cube

  END SUBROUTINE ReadCollisionData


  SUBROUTINE FreeCollisionData()
    IF (ALLOCATED(CD % LnT)) DEALLOCATE(CD % LnT)
    IF (ALLOCATED(CD % LnA)) DEALLOCATE(CD % LnA)
    IF (ALLOCATED(CD % LnW)) DEALLOCATE(CD % LnW)
    IF (ALLOCATED(CD % LnNuN)) DEALLOCATE(CD % LnNuN)
    IF (ALLOCATED(CD % LnElastic)) DEALLOCATE(CD % LnElastic)
    IF (ALLOCATED(CD % LnEi)) DEALLOCATE(CD % LnEi)
    IF (ALLOCATED(CD % LnKexc)) DEALLOCATE(CD % LnKexc)
    IF (ALLOCATED(CD % LnKsup)) DEALLOCATE(CD % LnKsup)
    IF (ALLOCATED(CD % LnMu)) DEALLOCATE(CD % LnMu)
    IF (ALLOCATED(CD % LnMuP)) DEALLOCATE(CD % LnMuP)
    IF (ALLOCATED(CD % LnMuH)) DEALLOCATE(CD % LnMuH)
    CD % nT = 0; CD % nA = 0; CD % nW = 0
    CD % Loaded = .FALSE.
  END SUBROUTINE FreeCollisionData


  !> Bracketing index and weight of x on an increasing grid, clamped to its ends
  PURE SUBROUTINE Locate( Grid, x, i, w )
    REAL(KIND=r8), INTENT(IN) :: Grid(:), x
    INTEGER, INTENT(OUT) :: i
    REAL(KIND=r8), INTENT(OUT) :: w
    INTEGER :: lo, hi, mid, n
    n = SIZE(Grid)
    IF (x <= Grid(1)) THEN
      i = 1; w = 0.0_r8; RETURN
    ELSE IF (x >= Grid(n)) THEN
      i = n - 1; w = 1.0_r8; RETURN
    END IF
    lo = 1; hi = n
    DO WHILE (hi - lo > 1)
      mid = (lo + hi) / 2
      IF (Grid(mid) <= x) THEN
        lo = mid
      ELSE
        hi = mid
      END IF
    END DO
    i = lo
    w = (x - Grid(lo)) / (Grid(lo+1) - Grid(lo))
  END SUBROUTINE Locate

  !> exp of a log-tabulated column at electron temperature Te [K] (log-log)
  FUNCTION Table1( LnValues, Te ) RESULT(V)
    REAL(KIND=r8), INTENT(IN) :: LnValues(:), Te
    REAL(KIND=r8) :: V, w
    INTEGER :: i
    CALL Locate( CD % LnT, LOG(Te / KelvinPerEV), i, w )
    V = EXP( (1.0_r8 - w) * LnValues(i) + w * LnValues(i+1) )
  END FUNCTION Table1

  !> exp of a log-tabulated cube at (theta, Coulomb parameter, magnetic parameter)
  FUNCTION Table3( LnValues, Te, A, W ) RESULT(V)
    REAL(KIND=r8), INTENT(IN) :: LnValues(:,:,:), Te, A, W
    REAL(KIND=r8) :: V, wt, wa, ww
    INTEGER :: it, ia, iw
    CALL Locate( CD % LnT, LOG(Te / KelvinPerEV), it, wt )
    CALL Locate( CD % LnA, LOG(MAX(A, 1.0d-300)), ia, wa )
    CALL Locate( CD % LnW, LOG(MAX(W, 1.0d-300)), iw, ww )
    V = EXP( &
        (1-ww) * ( (1-wa) * ((1-wt) * LnValues(it,ia,iw)   + wt * LnValues(it+1,ia,iw)) &
                 +    wa  * ((1-wt) * LnValues(it,ia+1,iw) + wt * LnValues(it+1,ia+1,iw)) ) &
      +    ww  * ( (1-wa) * ((1-wt) * LnValues(it,ia,iw+1)   + wt * LnValues(it+1,ia,iw+1)) &
                 +    wa  * ((1-wt) * LnValues(it,ia+1,iw+1) + wt * LnValues(it+1,ia+1,iw+1)) ) )
  END FUNCTION Table3


  !> Right-hand side S(T) of the Saha equation ne^2 / (nS - ne) = S [m^-3]
  FUNCTION SahaFactor( Pl, T ) RESULT(SahaS)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=r8), INTENT(IN) :: T
    REAL(KIND=r8) :: SahaS, Expo

    Expo = Pl % ChiJ / (kBoltz * T)
    IF (Expo > 700.0_r8) THEN
      SahaS = 0.0_r8
    ELSE
      SahaS = 2.0_r8 * Pl % WeightRatio * &
          ((2.0_r8*Pi*eMass*kBoltz*T) / (hPlanck*hPlanck))**1.5_r8 * EXP(-Expo)
    END IF
  END FUNCTION SahaFactor

  !> Saha electron density at electron temperature T and seed density nS
  FUNCTION SahaDensity( Pl, T, nS ) RESULT(ne)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=r8), INTENT(IN) :: T, nS
    REAL(KIND=r8) :: ne, SahaS
    ! Positive root of ne^2 + S ne - S nS = 0, in a form that avoids
    ! cancellation when S >> nS
    SahaS = SahaFactor( Pl, T )
    IF (SahaS > 0.0_r8) THEN
      ne = 2.0_r8 * SahaS * nS / (SahaS + SQRT(SahaS*SahaS + 4.0_r8*SahaS*nS))
    ELSE
      ne = 0.0_r8
    END IF
  END FUNCTION SahaDensity

  !> Derivative d(ne)/dT of the Saha electron density ne at electron
  !> temperature T and fixed seed density nS. From ne^2 + S ne - S nS = 0,
  !>   d(ne)/dS = (nS - ne) / (2 ne + S),  dS/dT = S (3/(2T) + chi/(kB T^2))
  FUNCTION SeedDensityDerivative( Pl, T, nS, ne ) RESULT(dne)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=r8), INTENT(IN) :: T, nS, ne
    REAL(KIND=r8) :: dne, SahaS

    dne = 0.0_r8
    SahaS = SahaFactor( Pl, T )
    IF (SahaS <= 0.0_r8 .OR. 2.0_r8*ne + SahaS <= 0.0_r8) RETURN
    dne = SahaS * (1.5_r8/T + Pl % ChiJ/(kBoltz*T*T)) * MAX(nS - ne, 0.0_r8) / (2.0_r8*ne + SahaS)
  END FUNCTION SeedDensityDerivative


  !> NRL electron-ion Coulomb logarithm (Z = 1, Te below 10 eV), floored at 1
  FUNCTION CoulombLogarithm( ne, Te ) RESULT(LnL)
    REAL(KIND=r8), INTENT(IN) :: ne, Te
    REAL(KIND=r8) :: LnL, theta
    theta = Te / KelvinPerEV
    IF (ne <= 0.0_r8) THEN
      LnL = 0.0_r8
      RETURN
    END IF
    LnL = MAX(23.0_r8 - LOG(SQRT(ne * 1.0d-6) * theta**(-1.5_r8)), 1.0_r8)
  END FUNCTION CoulombLogarithm


  !> Electron and ion mobilities at electron temperature Te, heavy-particle
  !> density N, electron density ne and |B|. Parallel (DC), Pedersen and Hall
  !> mobilities; for ions the Hall current runs the opposite way.
  !>
  !> Drifting: nu = N nu_en/N(theta) + 2.91e-12 ne lnL theta^-3/2 (NRL), and
  !>   the tensor of a single collision frequency, beta = mu |B|.
  !> Lorentz: collision frequencies add inside the velocity average, so the
  !>   Pedersen and Hall mobilities are tabulated on (theta, ne lnL / N,
  !>   omega_ce / N) and are not those of any single frequency.
  !> Ions: mu_i = mu_i0 N0 / N (polarization collisions, independent of speed).
  SUBROUTINE PlasmaTransport( Pl, Te, N, ne, Bmag, Tr )
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=r8), INTENT(IN) :: Te, N, ne, Bmag
    TYPE(PlasmaTransport_t), INTENT(OUT) :: Tr
    REAL(KIND=r8) :: theta, nu, beta, betaI, A, W

    theta = Te / KelvinPerEV
    Tr % ne = ne
    Tr % LnLambda = CoulombLogarithm( ne, Te )
    IF (Pl % Lorentz) THEN
      A = ne * Tr % LnLambda / N
      W = eCharge * Bmag / (eMass * N)
      Tr % Mu0 = Table3( CD % LnMu, Te, A, W ) / N
      Tr % MuP = Table3( CD % LnMuP, Te, A, W ) / N
      Tr % MuH = Table3( CD % LnMuH, Te, A, W ) / N
      IF (Bmag <= 0.0_r8) THEN
        Tr % MuP = Tr % Mu0
        Tr % MuH = 0.0_r8
      END IF
    ELSE
      nu = N * Table1( CD % LnNuN, Te ) + NrlCoefficient * ne * Tr % LnLambda * theta**(-1.5_r8)
      Tr % Mu0 = eCharge / (eMass * nu)
      beta = Tr % Mu0 * Bmag
      Tr % MuP = Tr % Mu0 / (1.0_r8 + beta*beta)
      Tr % MuH = Tr % Mu0 * beta / (1.0_r8 + beta*beta)
    END IF
    Tr % CollisionFrequency = eCharge / (eMass * Tr % Mu0)

    Tr % MuI = Pl % IonReducedMobility * Loschmidt / N
    betaI = Tr % MuI * Bmag
    Tr % MuIP = Tr % MuI / (1.0_r8 + betaI*betaI)
    Tr % MuIH = Tr % MuI * betaI / (1.0_r8 + betaI*betaI)
  END SUBROUTINE PlasmaTransport


  !> Elastic (recoil) energy loss coefficients of the electrons, W/(m^3 K), so
  !> that the loss is (Cn + Cei) (Te - Tg): to neutrals, e ne N G(theta) / Te,
  !> and to ions, e ne^2 lnL Gei(theta) / Te (the tables hold the
  !> electron-temperature integrals; (1 - Tg/Te) = (Te - Tg)/Te).
  SUBROUTINE ElasticLossCoefficients( Te, N, Tr, Cn, Cei )
    REAL(KIND=r8), INTENT(IN) :: Te, N
    TYPE(PlasmaTransport_t), INTENT(IN) :: Tr
    REAL(KIND=r8), INTENT(OUT) :: Cn, Cei
    Cn  = eCharge * Tr % ne * N * Table1( CD % LnElastic, Te ) / Te
    Cei = eCharge * Tr % ne * Tr % ne * Tr % LnLambda * Table1( CD % LnEi, Te ) / Te
  END SUBROUTINE ElasticLossCoefficients


  !> Net inelastic (K 4s-4p) energy loss of the electrons, W/m^3, for K(4p)
  !> populations at the gas temperature Tg (ExcitedAtGas). Zero when the
  !> excited states follow the electrons: excitation and de-excitation then
  !> balance by detailed balance.
  FUNCTION InelasticLoss( Pl, Te, Tg, ne, nK ) RESULT(L)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=r8), INTENT(IN) :: Te, Tg, ne, nK
    REAL(KIND=r8) :: L, r, yLow, yUp
    L = 0.0_r8
    IF (.NOT. Pl % ExcitedAtGas .OR. ne <= 0.0_r8 .OR. nK <= 0.0_r8) RETURN
    r = CD % ExcWeight * EXP(-MIN(CD % ExcEnergy * KelvinPerEV / Tg, 700.0_r8))
    yLow = 1.0_r8 / (1.0_r8 + r)
    yUp  = r / (1.0_r8 + r)
    L = eCharge * ne * nK * CD % ExcEnergy * &
        (yLow * Table1( CD % LnKexc, Te ) - yUp * Table1( CD % LnKsup, Te ))
  END FUNCTION InelasticLoss


  !> Electron temperature from the local balance of Joule heating of the
  !> electrons, e ne (mu_0 E'_par^2 + mu_P E'_perp^2), against elastic and
  !> inelastic losses, bracketed and bisected on [Tg, TeMax]. ne follows Te
  !> through the Saha equation. AtMax is set when the balance lies above TeMax
  !> and Te is capped there.
  FUNCTION ElectronTemperature( Pl, Tg, N, nS, Epar2, Eperp2, Bmag, TeMax, AtMax ) RESULT(Te)
    TYPE(SeedPlasma_t), INTENT(IN) :: Pl
    REAL(KIND=r8), INTENT(IN) :: Tg, N, nS, Epar2, Eperp2, Bmag, TeMax
    LOGICAL, INTENT(OUT) :: AtMax
    REAL(KIND=r8) :: Te, TeLo, TeHi
    INTEGER :: k

    AtMax = .FALSE.
    Te = Tg
    IF (Epar2 + Eperp2 <= 0.0_r8 .OR. TeMax <= Tg) RETURN

    ! Residual (heating - loss per electron) is positive at Tg; grow the upper
    ! bracket until it turns negative or reaches the cap
    TeLo = Tg
    TeHi = MIN(1.1_r8 * Tg, TeMax)
    DO k = 1, 60
      IF (Residual(TeHi) <= 0.0_r8) EXIT
      IF (TeHi >= TeMax) THEN
        AtMax = .TRUE.
        Te = TeMax
        RETURN
      END IF
      TeLo = TeHi
      TeHi = MIN(Tg + 2.0_r8 * (TeHi - Tg), TeMax)
    END DO

    DO k = 1, 100
      Te = 0.5_r8 * (TeLo + TeHi)
      IF (TeHi - TeLo < 1.0e-3_r8) EXIT
      IF (Residual(Te) > 0.0_r8) THEN
        TeLo = Te
      ELSE
        TeHi = Te
      END IF
    END DO

  CONTAINS

    FUNCTION Residual(T) RESULT(R)
      REAL(KIND=r8), INTENT(IN) :: T
      REAL(KIND=r8) :: R, ne, neEff, Cn, Cei
      TYPE(PlasmaTransport_t) :: Tr
      ne = SahaDensity( Pl, T, nS )
      ! A tiny floor keeps the per-electron balance defined in cold gas
      neEff = MAX(ne, 1.0d-12 * nS)
      CALL PlasmaTransport( Pl, T, N, neEff, Bmag, Tr )
      CALL ElasticLossCoefficients( T, N, Tr, Cn, Cei )
      R = eCharge * (Tr % Mu0 * Epar2 + Tr % MuP * Eperp2) &
          - ((Cn + Cei) * (T - Tg) + InelasticLoss( Pl, T, Tg, neEff, MAX(nS - ne, 0.0_r8) )) / neEff
    END FUNCTION Residual

  END FUNCTION ElectronTemperature

END MODULE MHDPlasma


MODULE MHDLog
  IMPLICIT NONE

CONTAINS

  SUBROUTINE LogLine(str)
    CHARACTER(LEN=*), INTENT(IN) :: str
    WRITE(*,*) TRIM(str)
  END SUBROUTINE LogLine

  SUBROUTINE LogSection(title)
    CHARACTER(LEN=*), INTENT(IN) :: title
    WRITE(*,*) '========================================'
    WRITE(*,*) TRIM(title)
    WRITE(*,*) '========================================'
  END SUBROUTINE LogSection

  SUBROUTINE LogSectionEnd()
    WRITE(*,*) '========================================'
  END SUBROUTINE LogSectionEnd

END MODULE MHDLog


!------------------------------------------------------------------------------
! MHDDiagnostics – electrode/logging subroutines that use Elmer types
!------------------------------------------------------------------------------
MODULE MHDDiagnostics
  USE DefUtils
  USE SolverUtils
  USE Types
  IMPLICIT NONE

CONTAINS

  SUBROUTINE LogElectrodeCktSolution(Solver, Potential, PotentialPerm, NPhi, NumPairs)
    TYPE(Solver_t), TARGET :: Solver
    REAL(dp), INTENT(IN)   :: Potential(:)
    INTEGER, INTENT(IN)    :: PotentialPerm(:)
    INTEGER, INTENT(IN)    :: NPhi, NumPairs

    INTEGER :: ep, cidxVp, cidxVm, cidxI
    INTEGER :: testNode, testRow
    REAL(dp) :: phiVar
    LOGICAL :: ok
    INTEGER :: multiplierSize
    TYPE(Variable_t), POINTER :: MultVar
    REAL(dp), POINTER :: MultiplierValues(:)

    IF (ParEnv % MyPE /= 0) RETURN

    ok = .FALSE.
    DO testNode = 1, SIZE(PotentialPerm)
      testRow = PotentialPerm(testNode)
      IF (testRow > 0 .AND. testRow <= NPhi) THEN
        phiVar = Potential(testRow)
        ok = .TRUE.
        EXIT
      END IF
    END DO

    IF (ok) THEN
      WRITE(*,'(A,ES18.10)') '[DBG] Sample potential value: phi=', phiVar
    ELSE
      WRITE(*,*) '[DBG] could not find a valid potential row to sanity-check'
    END IF

    MultVar => VariableGet(Solver % Mesh % Variables, 'Electrode Circuit Values')
    IF (.NOT. ASSOCIATED(MultVar)) THEN
      WRITE(*,*) '========================================='
      WRITE(*,*) '[Electrode Circuit Solution]'
      WRITE(*,'(A)') ' WARNING: Electrode Circuit Values variable not found!'
      WRITE(*,'(A)') '   Lagrange multipliers were not exported'
      WRITE(*,*) '========================================='
      RETURN
    END IF

    MultiplierValues => MultVar % Values
    multiplierSize = SIZE(MultiplierValues)

    WRITE(*,*) '========================================='
    WRITE(*,*) '[Electrode Circuit Solution]'
    WRITE(*,'(A,I6,A,I3,A,I6)') ' NPhi=', NPhi, '  Pairs=', NumPairs, &
        '  Multiplier size=', multiplierSize

    IF (multiplierSize < 3*NumPairs) THEN
      WRITE(*,'(A)') ' WARNING: Multiplier vector too small for constraint DOFs!'
      WRITE(*,'(A,I6,A,I6)') '   Expected at least ', 3*NumPairs, &
          ' but got ', multiplierSize
      WRITE(*,*) '========================================='
      RETURN
    END IF

    DO ep = 1, NumPairs
      cidxVp = 3*(ep-1) + 1
      cidxVm = 3*(ep-1) + 2
      cidxI  = 3*(ep-1) + 3
      WRITE(*,'(A,I3,A,3ES18.10)') ' EP=', ep, ' [Vp Vm I]=', &
        MultiplierValues(cidxVp), &
        MultiplierValues(cidxVm), &
        MultiplierValues(cidxI)
      WRITE(*,'(A,I3,A,ES18.10)') ' EP=', ep, ' Potential difference (Vp-Vm)=', &
        MultiplierValues(cidxVp) - MultiplierValues(cidxVm)
      WRITE(*,'(A,I3,A,ES18.10,A)') ' EP=', ep, ' Current magnitude: |I|=', &
        ABS(MultiplierValues(cidxI)), ' Amperes'
    END DO

    WRITE(*,*) '========================================'
  END SUBROUTINE LogElectrodeCktSolution


  SUBROUTINE DiagnoseElectrodeCurrents(Model, Solver, VolCurrent, CurrentPerm, Dim)
    TYPE(Model_t), INTENT(IN) :: Model
    TYPE(Solver_t), INTENT(IN) :: Solver
    REAL(dp), INTENT(IN) :: VolCurrent(:)
    INTEGER, INTENT(IN) :: CurrentPerm(:)
    INTEGER, INTENT(IN) :: Dim

    INTEGER :: bc, be, i, n, inode, jRow
    TYPE(Element_t), POINTER :: Elem
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(dp) :: Jx, Jy, Jz, Jmag
    REAL(dp) :: sumJx, sumJy, sumJz, sumJmag
    REAL(dp) :: minJmag, maxJmag, avgJmag
    INTEGER :: nNodes
    LOGICAL :: gotIt
    REAL(dp) :: currentDensityBC

    IF (ParEnv % MyPE /= 0) RETURN

    WRITE(*,*) '========================================'
    WRITE(*,*) '[Electrode Boundary Current Diagnostics]'
    WRITE(*,*) 'Nodal current density (averaged from bulk elements):'
    WRITE(*,*) '========================================'

    DO bc = 1, Model % NumberOfBCs
      i = ListGetInteger( Model % BCs(bc) % Values, 'Electrode Pair', gotIt )
      IF (.NOT. gotIt) THEN
        currentDensityBC = ListGetConstReal( Model % BCs(bc) % Values, &
          'Current Density', gotIt )
        IF (.NOT. gotIt) CYCLE
      END IF

      sumJx = 0.0_dp
      sumJy = 0.0_dp
      sumJz = 0.0_dp
      sumJmag = 0.0_dp
      minJmag = HUGE(minJmag)
      maxJmag = 0.0_dp
      nNodes = 0
      i = 0

      DO be = 1, Solver % Mesh % NumberOfBoundaryElements
        Elem => Solver % Mesh % Elements(Solver % Mesh % NumberOfBulkElements + be)
        IF (Elem % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag) CYCLE
        i = i + 1
        NodeIndexes => Elem % NodeIndexes
        n = Elem % TYPE % NumberOfNodes
        DO inode = 1, n
          jRow = CurrentPerm(NodeIndexes(inode))
          IF (jRow <= 0) CYCLE
          Jx = 0.0_dp
          Jy = 0.0_dp
          Jz = 0.0_dp
          IF (Dim >= 1) Jx = VolCurrent((jRow-1)*Dim + 1)
          IF (Dim >= 2) Jy = VolCurrent((jRow-1)*Dim + 2)
          IF (Dim >= 3) Jz = VolCurrent((jRow-1)*Dim + 3)
          Jmag = SQRT(Jx**2 + Jy**2 + Jz**2)
          sumJx = sumJx + Jx
          sumJy = sumJy + Jy
          sumJz = sumJz + Jz
          sumJmag = sumJmag + Jmag
          minJmag = MIN(minJmag, Jmag)
          maxJmag = MAX(maxJmag, Jmag)
          nNodes = nNodes + 1
        END DO
      END DO

      IF (nNodes > 0) THEN
        avgJmag = sumJmag / REAL(nNodes, dp)
        WRITE(*,'(A,I0,A,I5,A)') ' BC ', bc, ': ', nNodes, ' nodes'
        WRITE(*,'(A,3ES11.3,A)') '   Avg J = [', sumJx/REAL(nNodes,dp), &
          sumJy/REAL(nNodes,dp), sumJz/REAL(nNodes,dp), '] A/m²'
        WRITE(*,'(A,ES11.3,A,ES11.3,A)') '   |J| range: ', minJmag, ' to ', maxJmag, ' A/m²'
      END IF
    END DO
    WRITE(*,*) '========================================'
  END SUBROUTINE DiagnoseElectrodeCurrents


  ! --- Electrode geometry ---
  SUBROUTINE ComputeElectrodeAreas(Model, Solver, ElectrodePairOfBC, ElectrodeSignOfBC, &
      NumPairs, AreaPlus, AreaMinus)
    TYPE(Model_t), INTENT(IN) :: Model
    TYPE(Solver_t), INTENT(IN) :: Solver
    INTEGER, INTENT(IN) :: ElectrodePairOfBC(:), ElectrodeSignOfBC(:), NumPairs
    REAL(dp), INTENT(OUT) :: AreaPlus(:), AreaMinus(:)

    INTEGER :: ep, be, i, n, gp
    TYPE(Element_t), POINTER :: Elem
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Nodes_t) :: EN
    TYPE(GaussIntegrationPoints_t) :: Integ
    REAL(dp) :: Basis(MAX_ELEMENT_NODES), dBasisdx(MAX_ELEMENT_NODES,3)
    REAL(dp) :: SqrtElementMetric, s
    LOGICAL :: Stat

    ALLOCATE(EN % x(MAX_ELEMENT_NODES), EN % y(MAX_ELEMENT_NODES), EN % z(MAX_ELEMENT_NODES))
    AreaPlus = 0.0_dp
    AreaMinus = 0.0_dp

    DO be = 1, Solver % Mesh % NumberOfBoundaryElements
      Elem => Solver % Mesh % Elements(Solver % Mesh % NumberOfBulkElements + be)
      NodeIndexes => Elem % NodeIndexes
      DO i = 1, Model % NumberOfBCs
        IF (Elem % BoundaryInfo % Constraint /= Model % BCs(i) % Tag) CYCLE
        DO ep = 1, NumPairs
          IF (ElectrodePairOfBC(i) /= ep) CYCLE
          n = Elem % TYPE % NumberOfNodes
          EN % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes(1:n))
          EN % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes(1:n))
          EN % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes(1:n))
          Integ = GaussPoints(Elem)
          DO gp = 1, Integ % n
            Stat = ElementInfo(Elem, EN, Integ % u(gp), Integ % v(gp), Integ % w(gp), &
              SqrtElementMetric, Basis, dBasisdx)
            s = SqrtElementMetric * Integ % s(gp)
            IF (ElectrodeSignOfBC(i) == +1) THEN
              AreaPlus(ep) = AreaPlus(ep) + s
            ELSE
              AreaMinus(ep) = AreaMinus(ep) + s
            END IF
          END DO
        END DO
      END DO
    END DO

    IF (ParEnv % PEs > 1) THEN
      DO ep = 1, NumPairs
        AreaPlus(ep) = ParallelReduction(AreaPlus(ep))
        AreaMinus(ep) = ParallelReduction(AreaMinus(ep))
      END DO
    END IF

    IF (ParEnv % MyPE == 0) THEN
      WRITE(*,*) '[ComputeElectrodeAreas] Electrode areas computed:'
      DO ep = 1, NumPairs
        WRITE(*,'(A,I2,A,ES12.4,A,ES12.4,A)') '  Pair ', ep, ': A+ = ', AreaPlus(ep), &
          ' m², A- = ', AreaMinus(ep), ' m²'
      END DO
    END IF

    DEALLOCATE(EN % x, EN % y, EN % z)
  END SUBROUTINE ComputeElectrodeAreas


  SUBROUTINE UpdateLaggedCurrent(Solver, CurrentLagged, NumPairs, Iteration, &
      MaxChange, Converged, DampingFactor)
    TYPE(Solver_t), INTENT(IN) :: Solver
    REAL(dp), INTENT(INOUT) :: CurrentLagged(:)
    INTEGER, INTENT(IN) :: NumPairs, Iteration
    REAL(dp), INTENT(OUT) :: MaxChange
    LOGICAL, INTENT(OUT) :: Converged
    REAL(dp), INTENT(IN) :: DampingFactor

    INTEGER :: ep, cidxI, ierr, ConvergedInt
    REAL(dp) :: NewCurrent, NewCurrentDamped, Change
    TYPE(Variable_t), POINTER :: MultVar
    REAL(dp), POINTER :: MultiplierValues(:)
    REAL(KIND=dp) :: tolerance

    MaxChange = 0.0_dp
    Converged = .FALSE.
    ConvergedInt = 0

    IF (ParEnv % MyPE == 0) THEN
      MultVar => VariableGet(Solver % Mesh % Variables, 'Electrode Circuit Values')
      IF (.NOT. ASSOCIATED(MultVar)) THEN
        WRITE(*,*) 'WARNING: Cannot access electrode circuit values'
        Converged = .TRUE.
        ConvergedInt = 1
      ELSE
        MultiplierValues => MultVar % Values
        IF (SIZE(MultiplierValues) < 3*NumPairs) THEN
          WRITE(*,*) 'WARNING: Multiplier vector too small'
          Converged = .TRUE.
          ConvergedInt = 1
        ELSE
          WRITE(*,*) '========================================'
          WRITE(*,'(A,I0)') '[UpdateLaggedCurrent] Iteration ', Iteration
          WRITE(*,'(A,F6.3)') ' Damping factor: ', DampingFactor
          WRITE(*,*) 'Updating electrode currents for next iteration:'
          DO ep = 1, NumPairs
            cidxI = 3*(ep-1) + 3
            NewCurrent = MultiplierValues(cidxI)
            
            ! Apply damping: I_new_damped = I_old + damping * (I_new - I_old)
            NewCurrentDamped = CurrentLagged(ep) + DampingFactor * (NewCurrent - CurrentLagged(ep))
            
            Change = ABS(NewCurrentDamped - CurrentLagged(ep))
            MaxChange = MAX(MaxChange, Change)
            WRITE(*,'(A,I2,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A)') '  Pair ', ep, &
              ': I_old = ', CurrentLagged(ep), &
              ' → I_raw = ', NewCurrent, &
              ' → I_damped = ', NewCurrentDamped, &
              ' (ΔI = ', Change, ' A)'
            CurrentLagged(ep) = NewCurrentDamped
          END DO
          tolerance = 1.0_dp
          IF (MAXVAL(ABS(CurrentLagged)) > 10.0_dp) THEN
            tolerance = MAX(1.0_dp, 0.001_dp * MAXVAL(ABS(CurrentLagged)))
          END IF
          WRITE(*,'(A,ES12.4,A)') ' Max current change: ', MaxChange, ' A'
          WRITE(*,'(A,ES12.4,A)') ' Convergence tolerance: ', tolerance, ' A'
          IF (MaxChange < tolerance) THEN
            Converged = .TRUE.
            ConvergedInt = 1
            WRITE(*,*) ' ✓ Electrode current CONVERGED!'
          ELSE
            Converged = .FALSE.
            ConvergedInt = 0
            WRITE(*,*) ' ⚠ Continue iterating...'
          END IF
          WRITE(*,*) '========================================'
        END IF
      END IF
    END IF

    IF (ParEnv % PEs > 1) THEN
      CALL MPI_BCAST(CurrentLagged, NumPairs, MPI_DOUBLE_PRECISION, 0, &
          ELMER_COMM_WORLD, ierr)
      CALL MPI_BCAST(MaxChange, 1, MPI_DOUBLE_PRECISION, 0, &
          ELMER_COMM_WORLD, ierr)
      CALL MPI_BCAST(ConvergedInt, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
      Converged = (ConvergedInt /= 0)
    END IF
  END SUBROUTINE UpdateLaggedCurrent


  SUBROUTINE CheckBoundaryFlux(Model, Solver, Potential, PotentialPerm)
    TYPE(Model_t), INTENT(IN) :: Model
    TYPE(Solver_t), INTENT(IN) :: Solver
    REAL(dp), INTENT(IN) :: Potential(:)
    INTEGER, INTENT(IN) :: PotentialPerm(:)

    INTEGER :: bc, be, n, i, j, tg, N_Integ, matId
    TYPE(Element_t), POINTER :: Elem, Parent
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    REAL(dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
    REAL(dp) :: Basis(Model % MaxElementNodes)
    REAL(dp) :: dBasisdx(Model % MaxElementNodes, 3)
    REAL(dp) :: Normal(3), SqrtElementMetric, s, u, v, w
    REAL(dp) :: ElementPot(Model % MaxElementNodes)
    REAL(dp) :: GradPhi(3), Jn, FluxIntegral, BoundaryArea
    REAL(dp) :: sigma, Jgp(3)
    LOGICAL :: Stat, gotIt

    IF (ParEnv % MyPE /= 0) RETURN

    WRITE(*,*) '========================================'
    WRITE(*,*) '[HACKY DEBUG: Boundary Flux Check]'
    WRITE(*,*) 'Computing actual ∫J·n dS on boundaries'
    WRITE(*,*) 'NOTE: Simplified version - assumes isotropic conductivity'
    WRITE(*,*) '========================================'

    ALLOCATE(Nodes % x(Model % MaxElementNodes))
    ALLOCATE(Nodes % y(Model % MaxElementNodes))
    ALLOCATE(Nodes % z(Model % MaxElementNodes))

    DO bc = 1, Model % NumberOfBCs
      FluxIntegral = 0.0_dp
      BoundaryArea = 0.0_dp
      i = 0
      DO be = 1, Solver % Mesh % NumberOfBoundaryElements
        Elem => Solver % Mesh % Elements(Solver % Mesh % NumberOfBulkElements + be)
        IF (Elem % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag) CYCLE
        i = i + 1
        NodeIndexes => Elem % NodeIndexes
        n = Elem % TYPE % NumberOfNodes
        Parent => Elem % BoundaryInfo % Left
        IF (.NOT. ASSOCIATED(Parent)) Parent => Elem % BoundaryInfo % Right
        IF (.NOT. ASSOCIATED(Parent)) CYCLE
        matId = ListGetInteger(Model % Bodies(Parent % BodyId) % Values, &
          'Material', minv=1, maxv=Model % NumberOfMaterials)
        sigma = ListGetConstReal(Model % Materials(matId) % Values, &
          'Electric Conductivity', gotIt)
        IF (.NOT. gotIt) sigma = 1.0_dp
        ElementPot = 0.0_dp
        DO j = 1, n
          IF (PotentialPerm(NodeIndexes(j)) > 0) THEN
            ElementPot(j) = Potential(PotentialPerm(NodeIndexes(j)))
          END IF
        END DO
        Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        Nodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
        IntegStuff = GaussPoints(Elem)
        U_Integ => IntegStuff % u
        V_Integ => IntegStuff % v
        W_Integ => IntegStuff % w
        S_Integ => IntegStuff % s
        N_Integ = IntegStuff % n
        DO tg = 1, N_Integ
          u = U_Integ(tg)
          v = V_Integ(tg)
          w = W_Integ(tg)
          Stat = ElementInfo(Elem, Nodes, u, v, w, &
            SqrtElementMetric, Basis, dBasisdx)
          Normal = 0.0_dp
          CALL GetElementNormal(Elem, Nodes, u, v, Normal)
          s = SqrtElementMetric * S_Integ(tg)
          GradPhi = 0.0_dp
          DO j = 1, 3
            GradPhi(j) = SUM(dBasisdx(1:n, j) * ElementPot(1:n))
          END DO
          Jgp = -sigma * GradPhi
          Jn = DOT_PRODUCT(Jgp, Normal)
          FluxIntegral = FluxIntegral + Jn * s
          BoundaryArea = BoundaryArea + s
        END DO
      END DO
      IF (i > 0) THEN
        WRITE(*,'(A,I3,A,I5,A)') ' BC #', bc, ' (', i, ' elements)'
        WRITE(*,'(A,ES12.4,A)') '   Total flux ∫J·n dS = ', FluxIntegral, ' Amperes'
        WRITE(*,'(A,ES12.4,A)') '   Boundary area = ', BoundaryArea, ' m²'
        IF (BoundaryArea > 0) THEN
          WRITE(*,'(A,ES12.4,A)') '   Avg flux density J·n = ', FluxIntegral/BoundaryArea, ' A/m²'
        END IF
        IF (ABS(FluxIntegral) < 1.0_dp) THEN
          WRITE(*,*) '   CHECK: Boundary appears insulating (flux < 1 A)'
        ELSE
          WRITE(*,'(A,ES12.4,A)') '   CHECK: Non-zero flux detected: ', ABS(FluxIntegral), ' A'
        END IF
      END IF
    END DO
    WRITE(*,*) '========================================'
    DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
  END SUBROUTINE CheckBoundaryFlux


  SUBROUTINE GetElementNormal(Element, Nodes, u, v, Normal)
    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: u, v, Normal(3)
    REAL(KIND=dp) :: dBasisdx(Element % TYPE % NumberOfNodes, 3)
    REAL(KIND=dp) :: Basis(Element % TYPE % NumberOfNodes)
    REAL(KIND=dp) :: DetJ, Tangent1(3), Tangent2(3), NLen
    INTEGER :: n, i
    LOGICAL :: Stat

    n = Element % TYPE % NumberOfNodes
    Stat = ElementInfo(Element, Nodes, u, v, 0.0_dp, DetJ, Basis, dBasisdx)
    Tangent1 = 0.0_dp
    Tangent2 = 0.0_dp
    DO i = 1, n
      Tangent1(1) = Tangent1(1) + Nodes % x(i) * dBasisdx(i,1)
      Tangent1(2) = Tangent1(2) + Nodes % y(i) * dBasisdx(i,1)
      Tangent1(3) = Tangent1(3) + Nodes % z(i) * dBasisdx(i,1)
      Tangent2(1) = Tangent2(1) + Nodes % x(i) * dBasisdx(i,2)
      Tangent2(2) = Tangent2(2) + Nodes % y(i) * dBasisdx(i,2)
      Tangent2(3) = Tangent2(3) + Nodes % z(i) * dBasisdx(i,2)
    END DO
    Normal(1) = Tangent1(2)*Tangent2(3) - Tangent1(3)*Tangent2(2)
    Normal(2) = Tangent1(3)*Tangent2(1) - Tangent1(1)*Tangent2(3)
    Normal(3) = Tangent1(1)*Tangent2(2) - Tangent1(2)*Tangent2(1)
    NLen = SQRT(DOT_PRODUCT(Normal, Normal))
    IF (NLen > 1.0d-20) Normal = Normal / NLen
  END SUBROUTINE GetElementNormal


  SUBROUTINE DiagnoseElectrodePotentials(Model, Solver, Potential, PotentialPerm, &
      ElectrodePairOfBC, ElectrodeSignOfBC, NumPairs)
    TYPE(Model_t), INTENT(IN) :: Model
    TYPE(Solver_t), INTENT(IN) :: Solver
    REAL(dp), INTENT(IN) :: Potential(:)
    INTEGER, INTENT(IN) :: PotentialPerm(:)
    INTEGER, INTENT(IN) :: ElectrodePairOfBC(:), ElectrodeSignOfBC(:)
    INTEGER, INTENT(IN) :: NumPairs
    INTEGER :: ep, be, i, n, inode, pRow, sgn
    TYPE(Element_t), POINTER :: Elem
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(dp) :: sumPhi, avgPhi, minPhi, maxPhi
    INTEGER :: nNodes

    IF (ParEnv % MyPE /= 0) RETURN
    WRITE(*,*) '========================================'
    WRITE(*,*) '[Electrode Boundary Potential Diagnostics]'
    DO ep = 1, NumPairs
      DO sgn = -1, +1, 2
        sumPhi = 0.0_dp
        minPhi = HUGE(minPhi)
        maxPhi = -HUGE(maxPhi)
        nNodes = 0
        DO be = 1, Solver % Mesh % NumberOfBoundaryElements
          Elem => Solver % Mesh % Elements(Solver % Mesh % NumberOfBulkElements + be)
          NodeIndexes => Elem % NodeIndexes
          DO i = 1, Model % NumberOfBCs
            IF (Elem % BoundaryInfo % Constraint /= Model % BCs(i) % Tag) CYCLE
            IF (ElectrodePairOfBC(i) /= ep) CYCLE
            IF (ElectrodeSignOfBC(i) /= sgn) CYCLE
            n = Elem % TYPE % NumberOfNodes
            DO inode = 1, n
              pRow = PotentialPerm(NodeIndexes(inode))
              IF (pRow <= 0) CYCLE
              sumPhi = sumPhi + Potential(pRow)
              minPhi = MIN(minPhi, Potential(pRow))
              maxPhi = MAX(maxPhi, Potential(pRow))
              nNodes = nNodes + 1
            END DO
          END DO
        END DO
        IF (nNodes > 0) THEN
          avgPhi = sumPhi / REAL(nNodes, dp)
          IF (sgn == +1) THEN
            WRITE(*,'(A,I3,A,I5,A,ES12.4,A,ES12.4,A,ES12.4)') '  EP=', ep, ' (+) nodes=', nNodes, &
              ' phi: min=', minPhi, ' avg=', avgPhi, ' max=', maxPhi
          ELSE
            WRITE(*,'(A,I3,A,I5,A,ES12.4,A,ES12.4,A,ES12.4)') '  EP=', ep, ' (-) nodes=', nNodes, &
              ' phi: min=', minPhi, ' avg=', avgPhi, ' max=', maxPhi
          END IF
        END IF
      END DO
    END DO
    WRITE(*,*) '========================================'
  END SUBROUTINE DiagnoseElectrodePotentials


  SUBROUTINE DiagnoseBulkVsBoundaryCurrents(Model, Solver, VolCurrent, PotentialPerm, Dim)
    TYPE(Model_t), INTENT(IN) :: Model
    TYPE(Solver_t), INTENT(IN) :: Solver
    REAL(dp), INTENT(IN) :: VolCurrent(:)
    INTEGER, INTENT(IN) :: PotentialPerm(:)
    INTEGER, INTENT(IN) :: Dim
    TYPE(Element_t), POINTER :: Element, BoundaryElement
    INTEGER, POINTER :: NodeIndexes(:), BoundaryNodeIndexes(:)
    INTEGER :: t, be, i, j, n, pRow, bc_tag
    REAL(dp) :: Jmag, Jx, Jy, Jz
    REAL(dp) :: xc, yc, zc, dist_from_center
    REAL(dp) :: center_x, center_y, center_z
    REAL(dp) :: min_x, max_x, min_y, max_y, min_z, max_z
    REAL(dp) :: bulk_center_sum, bulk_center_max, bulk_center_count
    REAL(dp) :: bulk_boundary_sum, bulk_boundary_max, bulk_boundary_count
    REAL(dp) :: boundary_threshold
    INTEGER, PARAMETER :: MAX_BCS = 50
    REAL(dp) :: bc_current_sum(MAX_BCS), bc_current_max(MAX_BCS)
    INTEGER :: bc_elem_count(MAX_BCS)
    CHARACTER(LEN=MAX_NAME_LEN) :: bc_name
    LOGICAL :: bc_found(MAX_BCS)

    IF (ParEnv % MyPE /= 0) RETURN
    bulk_center_sum = 0.0_dp
    bulk_center_max = 0.0_dp
    bulk_center_count = 0.0_dp
    bulk_boundary_sum = 0.0_dp
    bulk_boundary_max = 0.0_dp
    bulk_boundary_count = 0.0_dp
    bc_current_sum = 0.0_dp
    bc_current_max = 0.0_dp
    bc_elem_count = 0
    bc_found = .FALSE.
    min_x = HUGE(min_x)
    max_x = -HUGE(max_x)
    min_y = HUGE(min_y)
    max_y = -HUGE(max_y)
    min_z = HUGE(min_z)
    max_z = -HUGE(max_z)
    DO i = 1, Model % NumberOfNodes
      min_x = MIN(min_x, Model % Nodes % x(i))
      max_x = MAX(max_x, Model % Nodes % x(i))
      min_y = MIN(min_y, Model % Nodes % y(i))
      max_y = MAX(max_y, Model % Nodes % y(i))
      min_z = MIN(min_z, Model % Nodes % z(i))
      max_z = MAX(max_z, Model % Nodes % z(i))
    END DO
    center_x = (min_x + max_x) / 2.0_dp
    center_y = (min_y + max_y) / 2.0_dp
    center_z = (min_z + max_z) / 2.0_dp
    boundary_threshold = 0.2_dp * MAX(max_x - min_x, max_y - min_y, max_z - min_z)
    WRITE(*,*) '========================================'
    WRITE(*,*) '[Bulk vs Boundary Current Diagnostics]'
    WRITE(*,'(A,3ES12.4)') '  Domain center: ', center_x, center_y, center_z
    WRITE(*,'(A,ES12.4)') '  Boundary threshold: ', boundary_threshold
    WRITE(*,*) ''
    DO t = 1, Solver % NumberOfActiveElements
      Element => Solver % Mesh % Elements(Solver % ActiveElements(t))
      NodeIndexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      xc = 0.0_dp
      yc = 0.0_dp
      zc = 0.0_dp
      DO i = 1, n
        xc = xc + Model % Nodes % x(NodeIndexes(i))
        yc = yc + Model % Nodes % y(NodeIndexes(i))
        zc = zc + Model % Nodes % z(NodeIndexes(i))
      END DO
      xc = xc / REAL(n, dp)
      yc = yc / REAL(n, dp)
      zc = zc / REAL(n, dp)
      Jx = 0.0_dp
      Jy = 0.0_dp
      Jz = 0.0_dp
      DO i = 1, n
        pRow = PotentialPerm(NodeIndexes(i))
        IF (pRow > 0) THEN
          Jx = Jx + VolCurrent(Dim*(pRow-1)+1)
          IF (Dim >= 2) Jy = Jy + VolCurrent(Dim*(pRow-1)+2)
          IF (Dim >= 3) Jz = Jz + VolCurrent(Dim*(pRow-1)+3)
        END IF
      END DO
      Jx = Jx / REAL(n, dp)
      Jy = Jy / REAL(n, dp)
      Jz = Jz / REAL(n, dp)
      Jmag = SQRT(Jx**2 + Jy**2 + Jz**2)
      dist_from_center = SQRT((xc - center_x)**2 + (yc - center_y)**2 + (zc - center_z)**2)
      IF (dist_from_center < boundary_threshold) THEN
        bulk_center_sum = bulk_center_sum + Jmag
        bulk_center_max = MAX(bulk_center_max, Jmag)
        bulk_center_count = bulk_center_count + 1.0_dp
      ELSE
        bulk_boundary_sum = bulk_boundary_sum + Jmag
        bulk_boundary_max = MAX(bulk_boundary_max, Jmag)
        bulk_boundary_count = bulk_boundary_count + 1.0_dp
      END IF
    END DO
    DO be = 1, Solver % Mesh % NumberOfBoundaryElements
      BoundaryElement => Solver % Mesh % Elements(Solver % Mesh % NumberOfBulkElements + be)
      BoundaryNodeIndexes => BoundaryElement % NodeIndexes
      n = BoundaryElement % TYPE % NumberOfNodes
      bc_tag = BoundaryElement % BoundaryInfo % Constraint
      IF (bc_tag <= 0 .OR. bc_tag > MAX_BCS) CYCLE
      Jx = 0.0_dp
      Jy = 0.0_dp
      Jz = 0.0_dp
      DO i = 1, n
        pRow = PotentialPerm(BoundaryNodeIndexes(i))
        IF (pRow > 0) THEN
          Jx = Jx + VolCurrent(Dim*(pRow-1)+1)
          IF (Dim >= 2) Jy = Jy + VolCurrent(Dim*(pRow-1)+2)
          IF (Dim >= 3) Jz = Jz + VolCurrent(Dim*(pRow-1)+3)
        END IF
      END DO
      Jx = Jx / REAL(n, dp)
      Jy = Jy / REAL(n, dp)
      Jz = Jz / REAL(n, dp)
      Jmag = SQRT(Jx**2 + Jy**2 + Jz**2)
      bc_current_sum(bc_tag) = bc_current_sum(bc_tag) + Jmag
      bc_current_max(bc_tag) = MAX(bc_current_max(bc_tag), Jmag)
      bc_elem_count(bc_tag) = bc_elem_count(bc_tag) + 1
      bc_found(bc_tag) = .TRUE.
    END DO
    WRITE(*,*) '--- Bulk Element Statistics ---'
    IF (bulk_center_count > 0) THEN
      WRITE(*,'(A,ES12.4,A,ES12.4)') '  Center region:  avg |J| = ', &
        bulk_center_sum / bulk_center_count, '  max |J| = ', bulk_center_max
      WRITE(*,'(A,I0)') '                  elements = ', INT(bulk_center_count)
    ELSE
      WRITE(*,*) '  Center region: No elements'
    END IF
    IF (bulk_boundary_count > 0) THEN
      WRITE(*,'(A,ES12.4,A,ES12.4)') '  Near boundary:  avg |J| = ', &
        bulk_boundary_sum / bulk_boundary_count, '  max |J| = ', bulk_boundary_max
      WRITE(*,'(A,I0)') '                  elements = ', INT(bulk_boundary_count)
    ELSE
      WRITE(*,*) '  Near boundary: No elements'
    END IF
    IF (bulk_center_count > 0 .AND. bulk_boundary_count > 0) THEN
      WRITE(*,'(A,F8.4)') '  Ratio (center/boundary): ', &
        (bulk_center_sum / bulk_center_count) / (bulk_boundary_sum / bulk_boundary_count)
    END IF
    WRITE(*,*) ''
    WRITE(*,*) '--- Boundary Element Statistics (by BC) ---'
    DO i = 1, Model % NumberOfBCs
      bc_tag = Model % BCs(i) % Tag
      IF (bc_tag <= 0 .OR. bc_tag > MAX_BCS) CYCLE
      IF (.NOT. bc_found(bc_tag)) CYCLE
      bc_name = ListGetString(Model % BCs(i) % Values, 'Name', bc_found(bc_tag))
      IF (.NOT. bc_found(bc_tag)) THEN
        WRITE(bc_name, '(A,I0)') 'BC_', bc_tag
      END IF
      IF (bc_elem_count(bc_tag) > 0) THEN
        WRITE(*,'(A,A,A)') '  BC: "', TRIM(bc_name), '"'
        WRITE(*,'(A,I0,A,I0)') '    Tag = ', bc_tag, '  Elements = ', bc_elem_count(bc_tag)
        WRITE(*,'(A,ES12.4,A,ES12.4)') '    avg |J| = ', &
          bc_current_sum(bc_tag) / REAL(bc_elem_count(bc_tag), dp), &
          '  max |J| = ', bc_current_max(bc_tag)
        IF (bulk_center_count > 0) THEN
          WRITE(*,'(A,F8.4)') '    Ratio (BC/center): ', &
            (bc_current_sum(bc_tag) / REAL(bc_elem_count(bc_tag), dp)) / &
            (bulk_center_sum / bulk_center_count)
        END IF
      END IF
    END DO
    WRITE(*,*) '========================================'
  END SUBROUTINE DiagnoseBulkVsBoundaryCurrents

END MODULE MHDDiagnostics
