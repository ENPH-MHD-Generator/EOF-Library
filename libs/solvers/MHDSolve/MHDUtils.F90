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
! MHDParams – central knobs for MHD solver behaviour
!
! This module collects tunable scalar parameters that should be shared by
! multiple MHD-related solvers. Keeping them here guarantees there is a
! single, well-documented “source of truth” instead of copy-pasted literals.
!------------------------------------------------------------------------------
MODULE MHDParams
  USE DefUtils, ONLY: dp
  IMPLICIT NONE

  ! REAL(KIND=dp), PARAMETER :: HallCoeffAlphaDefault = 0.070266_dp
  ! REAL(KIND=dp), PARAMETER :: HallCoeffAlphaDefault = 0.013333_dp

  ! https://pure.tue.nl/ws/files/4331881/7207091.pdf
  ! REAL(KIND=dp), PARAMETER :: HallCoeffAlphaDefault = -0.0403361_dp
  REAL(KIND=dp), PARAMETER :: HallCoeffAlphaDefault = 0.403361_dp


  ! Off
  ! REAL(KIND=dp), PARAMETER :: HallCoeffAlphaDefault = 0.0_dp

END MODULE MHDParams


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
