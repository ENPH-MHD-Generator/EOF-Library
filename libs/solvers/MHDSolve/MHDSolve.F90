!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Antti Pursula
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Aug 2002
! *
! *****************************************************************************/

!/******************************************************************************
! *
! *  Modified By: Felix Toft
! *  Email: felixtoft09@gmail.com
! *  Adapted from StatCurrent Solve
! *  https://github.com/ElmerCSC/elmerfem/blob/devel/fem/src/modules/StatCurrentSolve.F90
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. StatCurrentSolver.
!> ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver_Init( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE MHDUtils
  USE MHDLog
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)            :: Model
  TYPE(Solver_t), TARGET  :: Solver
  LOGICAL                 :: TransientSimulation
  REAL(KIND=dp)           :: dt
!------------------------------------------------------------------------------
  LOGICAL                 :: Found, Calculate
  TYPE(ValueList_t), POINTER :: Params
  INTEGER                 :: Dim, i
!------------------------------------------------------------------------------

  Params => GetSolverParams() 
  Dim    = CoordinateSystemDimension()

  !------------------------------------------------------------
  ! Exported variables
  !------------------------------------------------------------
  IF ( ListGetLogical( Params, 'Calculate Joule Heating', Found ) ) THEN
    CALL ListAddString( Params, &
         NextFreeKeyword('Exported Variable ', Params), &
         'Joule Heating' )
    ! Element-wise copy: nodal averaging smears the current and heating peaks at
    ! the electrode edges, and interpolating those onto larger OpenFOAM cells
    ! inflates their integrals. OpenFOAM gets these per element instead.
    CALL ListAddString( Params, &
         NextFreeKeyword('Exported Variable ', Params), &
         '-elem Element Joule Heating' )
  END IF

  IF ( ListGetLogical( Params, 'Calculate Nodal Heating', Found ) ) THEN
    CALL ListAddString( Params, &
         NextFreeKeyword('Exported Variable ', Params), &
         'Nodal Joule Heating' )
  END IF

  Calculate = ListGetLogical( Params, 'Calculate Volume Current', Found )
  IF ( Calculate ) THEN
    IF ( Dim == 2 ) THEN
      CALL ListAddString( Params, &
           NextFreeKeyword('Exported Variable ', Params), &
           'Volume Current[Volume Current:2]' )
    ELSE
      CALL ListAddString( Params, &
           NextFreeKeyword('Exported Variable ', Params), &
           'Volume Current[Volume Current:3]' )

    END IF

    DO i = 1, Dim
      CALL ListAddString( Params, &
           NextFreeKeyword('Exported Variable ', Params), &
           '-elem Element Volume Current '//TRIM(I2S(i)) )
    END DO
  END IF

  ! Seeded-plasma state computed from the Saha equation each call
  CALL ListAddString( Params, &
       NextFreeKeyword('Exported Variable ', Params), 'Ionization Fraction' )
  CALL ListAddString( Params, &
       NextFreeKeyword('Exported Variable ', Params), 'Electron Temperature' )
  CALL ListAddString( Params, &
       NextFreeKeyword('Exported Variable ', Params), 'Electron Density' )
  CALL ListAddString( Params, &
       NextFreeKeyword('Exported Variable ', Params), 'Electron Mobility' )
  CALL ListAddString( Params, &
       NextFreeKeyword('Exported Variable ', Params), &
       'Effective Electric Field[Effective Electric Field:3]' )

  ! Enable export of Lagrange multipliers (constraint DOF values)
  IF (.NOT. ListCheckPresent(Solver % Values, 'Export Lagrange Multiplier')) THEN
    CALL ListAddLogical(Solver % Values, 'Export Lagrange Multiplier', .TRUE.)
    CALL ListAddString(Solver % Values, 'Lagrange Multiplier Name', 'Electrode Circuit Values')
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE StatCurrentSolver_Init


    
!------------------------------------------------------------------------------
!>  Solve the Poisson equation for the electric potential and compute the 
!>  volume current and Joule heating
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ListMatrix
  USE MHDUtils
  USE MHDLog
  USE MHDDiagnostics
  USE MHDPlasma

  IMPLICIT NONE
!------------------------------------------------------------------------------ 
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET:: Solver
  REAL (KIND=DP) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Nodes_t) :: ElementNodes

  REAL (KIND=DP), POINTER :: ForceVector(:), Potential(:)
  REAL (KIND=DP), POINTER :: ElField(:), VolCurrent(:)
  REAL (KIND=DP), POINTER :: Heating(:), NodalHeating(:)
  REAL (KIND=DP), POINTER :: Cwrk(:,:,:)
  REAL (KIND=DP), ALLOCATABLE ::  Conductivity(:,:,:), &
    LocalStiffMatrix(:,:), Load(:), LocalForce(:)

  REAL (KIND=DP) :: Norm, HeatingTot, VolTot, CurrentTot, ControlTarget, ControlScaling = 1.0
  REAL (KIND=DP) :: PowerEmf, PowerLoad, UxAvg, CurrYAvg
  REAL (KIND=DP) :: Resistance, PotDiff
  REAL (KIND=DP) :: at, st, at0
#ifndef USE_ISO_C_BINDINGS
  REAL (KIND=DP) :: CPUTime, RealTime
#endif

  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: PotentialPerm(:)
  INTEGER :: i, j, k, n, t, istat, bf_id, LocalNodes, Dim, &
      iter, NonlinearIter

  LOGICAL :: AllocationsDone = .FALSE., gotIt, FluxBC
  LOGICAL :: CalculateField = .FALSE., ConstantWeights
  LOGICAL :: CalculateCurrent, CalculateHeating, CalculateNodalHeating
  LOGICAL :: ControlPower, ControlCurrent, Control

  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var

  CHARACTER(LEN=MAX_NAME_LEN) :: EquationName
  CHARACTER(LEN=256) :: LogMsg

  LOGICAL :: GetCondAtIp
  ! Gauss points for linear tetrahedra in assembly and current evaluation
  INTEGER :: TetraPoints
  TYPE(ValueHandle_t) :: CondAtIp_h
  REAL(KIND=dp) :: CondAtIp

  ! Velocity and Magnetic field values
  REAL(KIND=dp), POINTER :: UxVals(:), UyVals(:), UzVals(:)
  REAL(KIND=dp), POINTER :: BxVals(:), ByVals(:), BzVals(:)
  INTEGER, POINTER :: UxPerm(:), UyPerm(:), UzPerm(:)
  INTEGER, POINTER :: BxPerm(:), ByPerm(:), BzPerm(:)
  TYPE(Variable_t), POINTER :: UxVar, UyVar, UzVar, BxVar, ByVar, BzVar

  ! Pressure, gas temperature and conductivity (Saha ionization)
  REAL(KIND=dp), POINTER :: PVals(:), TgasVals(:), SigVals(:)
  INTEGER, POINTER :: PPerm(:), TgasPerm(:), SigPerm(:)
  TYPE(Variable_t), POINTER :: PVar, TgasVar, SigVar

  ! Seeded-plasma outputs (exported variables of this solver, one shared perm)
  REAL(KIND=dp), POINTER :: IonFrac(:), ElecTemp(:), ElecDens(:), ElecMob(:)
  REAL(KIND=dp), POINTER :: EffField(:)   ! E' = -grad(phi) + U x B, like VolCurrent

  ! Element-wise current and heating sent to OpenFOAM (shared Perm over elements)
  REAL(KIND=dp), POINTER :: ElemCurr1(:), ElemCurr2(:), ElemCurr3(:), ElemHeating(:)
  INTEGER, POINTER :: ElemPerm(:)
  INTEGER, POINTER :: PlasmaPerm(:)
  REAL(KIND=dp) :: TeChange, TeTol
  ! Per-node electron temperature relaxation and last step, for damping
  ! oscillating nodes (indexed like ElecTemp)
  REAL(KIND=dp), ALLOCATABLE :: TeNodeRelax(:), TeLastStep(:)
  ! Lumped nodal volumes (integral of each basis function), for the
  ! volume-weighted electron temperature convergence measure
  REAL(KIND=dp), ALLOCATABLE :: TeNodeWeight(:)

  ! Electrode Unknowns
  INTEGER, ALLOCATABLE :: ElectrodePairOfBC(:)
  INTEGER, ALLOCATABLE :: ElectrodeSignOfBC(:)
  CHARACTER(len=32) :: SignStr
  INTEGER :: NumElectrodePairs
  INTEGER :: sign
  
  ! Lagged iteration removed: current injection via circuit DOFs only

  ! Auxiliary matrix for electrode constraints
  TYPE(Matrix_t), POINTER, SAVE :: AuxMatrix => NULL()
  INTEGER :: NPhi, PermMax, GlobalNPhi

  ! Diagnostics for solution change (variables declared below in main declarations)
  REAL(dp), ALLOCATABLE :: OldPotential(:)

  ! Resistances
  REAL(dp), ALLOCATABLE :: ElectrodeResistance(:)
  LOGICAL :: gotItR
  REAL(dp) :: Rbc

  SAVE LocalStiffMatrix, Load, LocalForce, &
  ElementNodes, CalculateCurrent, CalculateHeating, &
  AllocationsDone, VolCurrent, Heating, Conductivity, &
  CalculateField, ConstantWeights, &
  Cwrk, ControlScaling, CalculateNodalHeating, &
  UxVar, UyVar, UzVar, BxVar, ByVar, BzVar, &
  UxVals, UyVals, UzVals, &
  BxVals, ByVals, BzVals, &
  UxPerm, UyPerm, UzPerm, &
  BxPerm, ByPerm, BzPerm, OldPotential, &
  PVar, TgasVar, SigVar, &
  PVals, TgasVals, SigVals, &
  PPerm, TgasPerm, SigPerm, &
  IonFrac, ElecTemp, ElecDens, ElecMob, PlasmaPerm, EffField, &
  TeNodeRelax, TeLastStep, TeNodeWeight, &
  ElemCurr1, ElemCurr2, ElemCurr3, ElemHeating, ElemPerm

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF(.NOT.ASSOCIATED(Solver % Matrix)) RETURN

  NumElectrodePairs = 0
  DO i = 1, Model % NumberOfBCs
    k = ListGetInteger( Model % BCs(i) % Values, 'Electrode Pair', gotIt )
    IF ( gotIt ) THEN
      NumElectrodePairs = MAX( NumElectrodePairs, k )
    END IF
  END DO

  ! Alocate electrode constraint memory
  IF (ALLOCATED(ElectrodeResistance)) DEALLOCATE(ElectrodeResistance)
  ALLOCATE( ElectrodeResistance(NumElectrodePairs) )
  ElectrodeResistance = -1.0_dp   ! sentinel = unset

  DO i = 1, Model % NumberOfBCs
    k = ListGetInteger(Model % BCs(i) % Values, 'Electrode Pair', gotIt)
    IF (.NOT. gotIt) CYCLE

    Rbc = GetCReal(Model % BCs(i) % Values, 'Electrode Resistance', gotItR)
    IF (.NOT. gotItR) CYCLE

    IF (ElectrodeResistance(k) < 0.0_dp) THEN
      ElectrodeResistance(k) = Rbc
    ELSE
      IF (ABS(ElectrodeResistance(k) - Rbc) > 1.0e-14_dp) THEN
        CALL Fatal('StatCurrentSolver','Electrode Resistance mismatch for pair index')
      END IF
    END IF
  END DO

  Potential     => Solver % Variable % Values
  PotentialPerm => Solver % Variable % Perm
  Params => GetSolverParams()

  LocalNodes = Model % NumberOfNodes
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  Norm = Solver % Variable % Norm
  DIM = CoordinateSystemDimension()

  ! We don't support 2 dimensions for MHD
  IF (Dim /= 3) THEN
    CALL Fatal( &
      'StatCurrentSolver', &
      'This solver requires a fully 3D coordinate system. ' // &
      'CoordinateSystemDimension() != 3. Aborting.' )
  END IF

  ControlTarget = GetCReal( Params,'Power Control',ControlPower)
  IF(ControlPower) THEN
    ControlCurrent = .FALSE.
  ELSE
    ControlTarget = GetCReal( Params,'Current Control',ControlCurrent)
  END IF
  Control = ControlPower .OR. ControlCurrent

  ! To obtain convergence rescale the potential to the original BCs
  IF( Control ) THEN
    Potential = Potential / ControlScaling
    Solver % Variable % Norm = Solver % Variable % Norm / ControlScaling
  END IF

  NonlinearIter = ListGetInteger( Params, &
      'Nonlinear System Max Iterations', GotIt )
  IF ( .NOT. GotIt ) NonlinearIter = 1

  GetCondAtIp = ListGetLogical( Params,'Conductivity At Ip',GotIt )

  ! Linear tetrahedra have constant basis gradients and our fields vary
  ! linearly, so the one-point centroid rule is the standard choice; Elmer's
  ! default of 4 points quadruples the element-loop cost. Assembly and current
  ! evaluation use the same rule, so the discrete energy balance still closes.
  TetraPoints = ListGetInteger( Params, 'Tetra Integration Points', GotIt )
  IF ( .NOT. GotIt ) TetraPoints = 1
  IF ( ALL(TetraPoints /= (/ 1, 4, 11 /)) ) &
      CALL Fatal('StatCurrentSolver', 'Tetra Integration Points must be 1, 4 or 11')

  !------------------------------------------------------------
  ! Electrode allocation and assignment
  !------------------------------------------------------------
  ALLOCATE( ElectrodePairOfBC( Model % NumberOfBCs ) )
  ALLOCATE( ElectrodeSignOfBC( Model % NumberOfBCs ) )

  ElectrodePairOfBC = 0
  ElectrodeSignOfBC = 0

  DO i = 1, Model % NumberOfBCs

    ! Is this BC an electrode?
    ElectrodePairOfBC(i) = ListGetInteger( &
        Model % BCs(i) % Values, 'Electrode Pair', gotIt )

    IF (.NOT. gotIt) CYCLE

    ! Get sign ONCE
    SignStr = ListGetString( Model % BCs(i) % Values, &
                            'Electrode Sign', gotIt )

    IF (.NOT. gotIt) THEN
      CALL Fatal( 'StatCurrentSolver', &
        'Electrode BC missing Electrode Sign (use "plus" or "minus")' )
    END IF

    SignStr = TRIM( SignStr )

    IF ( SignStr == 'plus' ) THEN
      ElectrodeSignOfBC(i) = +1
    ELSE IF ( SignStr == 'minus' ) THEN
      ElectrodeSignOfBC(i) = -1
    ELSE
      CALL Fatal( 'StatCurrentSolver', &
        'Electrode Sign must be "plus" or "minus" (lowercase)' )
    END IF
  END DO

     
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
    N = Model % MaxElementNodes

    IF(AllocationsDone) THEN
      DEALLOCATE( ElementNodes % x, &
                ElementNodes % y,   &
                ElementNodes % z,   &
                Conductivity,       &
                LocalForce,         &
                LocalStiffMatrix,   &
                Load )
    END IF

    ALLOCATE( ElementNodes % x(N),   &
              ElementNodes % y(N),   &
              ElementNodes % z(N),   &
              Conductivity(3,3,N),   &
              LocalForce(N),         &
              LocalStiffMatrix(N,N), &
              Load(N),               &
              STAT=istat )

    IF ( istat /= 0 ) THEN
      CALL Fatal( 'StatCurrentSolve', 'Memory allocation error.' )
    END IF

    NULLIFY( Cwrk )

    UxVar => VariableGet( Solver % Mesh % Variables, 'Ux' )
    UyVar => VariableGet( Solver % Mesh % Variables, 'Uy' )
    UzVar => VariableGet( Solver % Mesh % Variables, 'Uz' )

    BxVar => VariableGet( Solver % Mesh % Variables, 'Bx' )
    ByVar => VariableGet( Solver % Mesh % Variables, 'By' )
    BzVar => VariableGet( Solver % Mesh % Variables, 'Bz' )

    IF (.NOT.ASSOCIATED(UxVar)) CALL Fatal('StatCurrentSolver','Ux not found')
    IF (.NOT.ASSOCIATED(UyVar)) CALL Fatal('StatCurrentSolver','Uy not found')
    IF (.NOT.ASSOCIATED(UzVar)) CALL Fatal('StatCurrentSolver','Uz not found')

    IF (.NOT.ASSOCIATED(BxVar)) CALL Fatal('StatCurrentSolver','Bx not found')
    IF (.NOT.ASSOCIATED(ByVar)) CALL Fatal('StatCurrentSolver','By not found')
    IF (.NOT.ASSOCIATED(BzVar)) CALL Fatal('StatCurrentSolver','Bz not found')

    UxVals => UxVar % Values ; UxPerm => UxVar % Perm
    UyVals => UyVar % Values ; UyPerm => UyVar % Perm
    UzVals => UzVar % Values ; UzPerm => UzVar % Perm

    BxVals => BxVar % Values ; BxPerm => BxVar % Perm
    ByVals => ByVar % Values ; ByPerm => ByVar % Perm
    BzVals => BzVar % Values ; BzPerm => BzVar % Perm

    PVar     => VariableGet( Solver % Mesh % Variables, 'Pressure' )
    TgasVar  => VariableGet( Solver % Mesh % Variables, 'Gas Temperature' )
    SigVar   => VariableGet( Solver % Mesh % Variables, 'Electric Conductivity' )

    IF (.NOT.ASSOCIATED(PVar))     CALL Fatal('StatCurrentSolver','Pressure not found')
    IF (.NOT.ASSOCIATED(TgasVar))  CALL Fatal('StatCurrentSolver','Gas Temperature not found')
    IF (.NOT.ASSOCIATED(SigVar))   CALL Fatal('StatCurrentSolver','Electric Conductivity not found')

    PVals     => PVar    % Values ; PPerm     => PVar    % Perm
    TgasVals  => TgasVar % Values ; TgasPerm  => TgasVar % Perm
    SigVals   => SigVar  % Values ; SigPerm   => SigVar  % Perm

    Var => VariableGet( Solver % Mesh % Variables, 'Ionization Fraction' )
    IF (.NOT.ASSOCIATED(Var)) CALL Fatal('StatCurrentSolver','Ionization Fraction not found')
    IonFrac => Var % Values ; PlasmaPerm => Var % Perm

    Var => VariableGet( Solver % Mesh % Variables, 'Electron Temperature' )
    IF (.NOT.ASSOCIATED(Var)) CALL Fatal('StatCurrentSolver','Electron Temperature not found')
    ElecTemp => Var % Values

    Var => VariableGet( Solver % Mesh % Variables, 'Electron Density' )
    IF (.NOT.ASSOCIATED(Var)) CALL Fatal('StatCurrentSolver','Electron Density not found')
    ElecDens => Var % Values

    Var => VariableGet( Solver % Mesh % Variables, 'Electron Mobility' )
    IF (.NOT.ASSOCIATED(Var)) CALL Fatal('StatCurrentSolver','Electron Mobility not found')
    ElecMob => Var % Values

    Var => VariableGet( Solver % Mesh % Variables, 'Effective Electric Field' )
    IF (.NOT.ASSOCIATED(Var)) CALL Fatal('StatCurrentSolver','Effective Electric Field not found')
    EffField => Var % Values

    CalculateCurrent = ListGetLogical( Params, &
        'Calculate Volume Current', GotIt )
    IF ( CalculateCurrent ) THEN
      Var => VariableGet( Solver % Mesh % Variables,'Volume Current')
      IF( ASSOCIATED( Var) ) THEN
        VolCurrent => Var % Values
      ELSE
        CALL Fatal('StatCurrentSolver','Volume Current does not exist')
      END IF
    END IF

    CalculateHeating = ListGetLogicalAnyEquation( &
        Model,'Calculate Joule heating')
    IF ( .NOT. CalculateHeating )  &
        CalculateHeating = ListGetLogical( Params, &
        'Calculate Joule Heating', GotIt )
    IF ( CalculateHeating ) THEN
      Var => VariableGet( Solver % Mesh % Variables,'Joule Heating')
      IF( ASSOCIATED( Var) ) THEN
        Heating => Var % Values
      ELSE
        CALL Fatal('StatCurrentSolver','Joule Heating does not exist')
      END IF
    END IF

    NULLIFY( ElemCurr1, ElemCurr2, ElemCurr3, ElemHeating, ElemPerm )
    IF ( CalculateCurrent ) THEN
      Var => VariableGet( Solver % Mesh % Variables,'Element Volume Current 1')
      IF( .NOT. ASSOCIATED(Var) ) CALL Fatal('StatCurrentSolver','Element Volume Current 1 does not exist')
      ElemCurr1 => Var % Values ; ElemPerm => Var % Perm
      Var => VariableGet( Solver % Mesh % Variables,'Element Volume Current 2')
      IF( .NOT. ASSOCIATED(Var) ) CALL Fatal('StatCurrentSolver','Element Volume Current 2 does not exist')
      ElemCurr2 => Var % Values
      IF( DIM == 3 ) THEN
        Var => VariableGet( Solver % Mesh % Variables,'Element Volume Current 3')
        IF( .NOT. ASSOCIATED(Var) ) CALL Fatal('StatCurrentSolver','Element Volume Current 3 does not exist')
        ElemCurr3 => Var % Values
      END IF
    END IF

    IF ( CalculateHeating ) THEN
      Var => VariableGet( Solver % Mesh % Variables,'Element Joule Heating')
      IF( .NOT. ASSOCIATED(Var) ) CALL Fatal('StatCurrentSolver','Element Joule Heating does not exist')
      ElemHeating => Var % Values
      IF( .NOT. ASSOCIATED(ElemPerm) ) ElemPerm => Var % Perm
    END IF

    CalculateNodalHeating = ListGetLogical( Params, &
        'Calculate Nodal Heating', GotIt )
    IF ( CalculateNodalHeating ) THEN
      Var => VariableGet( Solver % Mesh % Variables,'Nodal Joule Heating')
      IF( ASSOCIATED( Var) ) THEN
        NodalHeating => Var % Values
      ELSE
        CALL Fatal('StatCurrentSolver','Nodal Joule Heating does not exist')
      END IF
    END IF

    ConstantWeights = ListGetLogical( Params, &
        'Constant Weights', GotIt )

!------------------------------------------------------------------------------

    IF ( .NOT.ASSOCIATED( StiffMatrix % MassValues ) ) THEN
      ALLOCATE( StiffMatrix % Massvalues( LocalNodes ) )
      StiffMatrix % MassValues = 0.0d0
    END IF

    ! Add electric field to variable list (disabled)
    IF ( CalculateField ) THEN
      CALL Info('StatCurrentSolver_bulk', '*** ABOUT TO ADD VARIABLE ***', Level=1)
      CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, &
            Solver, 'Electric Field', dim, ElField, PotentialPerm)
    END IF

    AllocationsDone = .TRUE.
  END IF
  

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

  EquationName = ListGetString( Params, 'Equation' )

  CALL Info( 'StatCurrentSolve', '-------------------------------------',Level=4 )
  CALL Info( 'StatCurrentSolve', 'STAT CURRENT SOLVER:  ', Level=4 )
  CALL Info( 'StatCurrentSolve', '-------------------------------------',Level=4 )

  CALL DefaultStart()

  ! The electron temperature depends on the current, so the plasma state is
  ! re-evaluated every nonlinear iteration from the previous iteration's
  ! current (or the previous coupling step's, on the first iteration).
  TeTol = ListGetCReal( Params, 'Nonlinear System Convergence Tolerance', GotIt )
  IF ( .NOT. GotIt ) TeTol = 1.0e-6_dp

  DO iter = 1, NonlinearIter
    at  = CPUTime()
    at0 = RealTime()

    CALL UpdateSeededPlasma( TeChange, iter == 1 )

    IF ( NonlinearIter > 1 ) THEN
      WRITE( Message, '(a,I0)' ) 'Static current iteration: ', iter
      CALL Info( 'StatCurrentSolve', Message, LEVEL=4 )
    END IF
    CALL Info( 'StatElecSolve', 'Starting Assembly...', Level=6 )

    CALL DefaultInitialize()
    
    !------------------------------------------------------------
    !    Do the assembly
    !------------------------------------------------------------

    IF( GetCondAtIp ) THEN
      CALL ListInitElementKeyword( CondAtIp_h,'Material','Electric Conductivity')
    END IF
      
    DO t = 1, Solver % NumberOfActiveElements

      IF ( RealTime() - at0 > 1.0 ) THEN
        WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
            (1.0*Solver % NumberOfActiveElements)), ' % done'

        CALL Info( 'StatCurrentSolve', Message, Level=5 )

        at0 = RealTime()
      END IF

      !------------------------------------------------------------------------------
      !        Check if this element belongs to a body where potential
      !        should be calculated
      !------------------------------------------------------------------------------
      CurrentElement => GetActiveElement(t)
      NodeIndexes => CurrentElement % NodeIndexes

      n = GetElementNOFNodes()

      ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

      bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
          Values, 'Body Force', gotIt, minv=1, maxv=Model % NumberOfBodyForces )

      Load  = 0.0d0
      IF ( gotIt ) THEN
        Load(1:n) = ListGetReal( Model % BodyForces(bf_id) % Values, &
            'Current Source',n,NodeIndexes, Gotit )
      END IF

      IF( .NOT. GetCondAtIp ) THEN

        CALL NodalConductivity( Conductivity, NodeIndexes, n )
      END IF

      !------------------------------------------------------------------------------
      !      Get element local matrix, and rhs vector
      !------------------------------------------------------------------------------
      CALL StatCurrentCompose( LocalStiffMatrix,LocalForce, &
          Conductivity,Load,CurrentElement,n,ElementNodes )
      !------------------------------------------------------------------------------
      !      Update global matrix and rhs vector from local matrix & vector
      !------------------------------------------------------------------------------

      CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

    END DO
    
    !-----------------------------------------------------------------------------
    !     Neumann boundary conditions
    !------------------------------------------------------------------------------
    DO t = Solver % Mesh % NumberOfBulkElements + 1, &
        Solver % Mesh % NumberOfBulkElements + &
        Solver % Mesh % NumberOfBoundaryElements

      CurrentElement => Solver % Mesh % Elements(t)

      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
          Model % BCs(i) % Tag ) THEN

          Model % CurrentElement => CurrentElement
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes
          IF ( ANY( PotentialPerm(NodeIndexes) <= 0 ) ) CYCLE

          !------------------------------------------------------------
          ! Electrode Boundary Neumann BCs
          !------------------------------------------------------------
          k = ListGetInteger(Model % BCs(i) % Values, 'Electrode Pair', gotIt)
          IF (gotIt) THEN
            ! This is an electrode boundary - skip explicit Neumann BC application
            ! Current is injected through circuit DOFs only
            CYCLE  ! Skip all Neumann handling for electrodes
          END IF

          FluxBC = ListGetLogical(Model % BCs(i) % Values, &
              'Current Density BC',gotIt)
          IF(GotIt .AND. .NOT. FluxBC) CYCLE

          ! BC: cond dPhi/dn = g
          Load = 0.0d0
          Load(1:n) = ListGetReal( Model % BCs(i) % Values,'Current Density', &
              n,NodeIndexes,gotIt )
          IF(.NOT. GotIt) CYCLE

          ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

          CALL StatCurrentBoundary( LocalStiffMatrix, LocalForce,  &
              Load, CurrentElement, n, ElementNodes )
          CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
        END IF ! of currentelement bc == bcs(i)
      END DO ! of i=1,model bcs
    END DO   ! Neumann BCs
    
    CALL DefaultFinishBulkAssembly()

    CALL DefaultFinishAssembly()

    CALL DefaultDirichletBCs()

    NPhi = Solver % Matrix % NumberOfRows
    IF (NPhi < 1) THEN
      CALL Fatal('StatCurrentSolver', 'Matrix NumberOfRows <= 0!')
    END IF
    
    IF (ParEnv % PEs > 1) THEN
      IF (ASSOCIATED(Solver % Matrix % ParallelInfo)) THEN
        IF (ASSOCIATED(Solver % Matrix % ParallelInfo % NeighbourList)) THEN
          IF (SIZE(Solver % Matrix % ParallelInfo % NeighbourList) /= NPhi) THEN
            IF (ParEnv % MyPE == 0) THEN
              WRITE(*,'(A,I0,A,I0,A)') '[StatCurrentSolver] ParallelInfo size ', &
                SIZE(Solver % Matrix % ParallelInfo % NeighbourList), ' /= ', NPhi, ', reinitializing'
            END IF
            CALL ParallelInitMatrix(Solver, Solver % Matrix)
          END IF
        END IF
      END IF
    END IF

    PermMax = MAXVAL(PotentialPerm, MASK=(PotentialPerm > 0))
    IF (ParEnv % PEs > 1) THEN
      GlobalNPhi = NINT(ParallelReduction(REAL(NPhi, dp), 2))  ! MPI_MAX
    ELSE
      GlobalNPhi = NPhi
    END IF

    ! Disconnect the old AddMatrix pointer before creating a new one.
    ! The old matrix is left for Elmer's memory management system to handle.
    ! Do NOT attempt to free it manually - Elmer has modified its internal
    ! structure and freeing it causes "invalid pointer" crashes.
    IF (ASSOCIATED(Solver % Matrix % AddMatrix)) THEN
      Solver % Matrix % AddMatrix => NULL()
    END IF

    IF (NumElectrodePairs > 0) THEN
      ! Build a fresh electrode constraint matrix for this iteration.
      ! BuildElectrodeAddMatrix will allocate a new matrix structure.
      CALL BuildElectrodeAddMatrix( Model, Solver, AuxMatrix, &
          PotentialPerm, ElectrodePairOfBC, ElectrodeSignOfBC, &
          ElectrodeResistance, NumElectrodePairs, NPhi )

      Solver % Matrix % AddMatrix => AuxMatrix
      
      ! Enable export of Lagrange multipliers (constraint DOF values)
      IF (.NOT. ListCheckPresent(Solver % Values, 'Export Lagrange Multiplier')) THEN
        CALL ListAddLogical(Solver % Values, 'Export Lagrange Multiplier', .TRUE.)
        CALL ListAddString(Solver % Values, 'Lagrange Multiplier Name', 'Electrode Circuit Values')
      END IF
    ELSE
      AuxMatrix => NULL()
      Solver % Matrix % AddMatrix => NULL()
      IF (ParEnv % MyPE == 0) THEN
        WRITE(*,'(A)') ' [StatCurrentSolver] NumElectrodePairs=0: skipping electrode AddMatrix assembly'
      END IF
    END IF

    at = CPUTime() - at
    WRITE( Message, * ) 'Assembly (s)          :',at
    CALL Info( 'StatCurrentSolve', Message, Level=5 )
    !------------------------------------------------------------------------------
    !    Solve the system and we are done.
    !------------------------------------------------------------------------------
    st = CPUTime()
    
    ! Store old potential for comparison
    IF (.NOT. ALLOCATED(OldPotential)) THEN
      ALLOCATE(OldPotential(SIZE(Potential)))
      OldPotential = 0.0_dp
    END IF
    IF (SIZE(OldPotential) /= SIZE(Potential)) THEN
      DEALLOCATE(OldPotential)
      ALLOCATE(OldPotential(SIZE(Potential)))
      OldPotential = 0.0_dp
    END IF
    
    ! Save old solution
    OldPotential = Potential
    
    Norm = DefaultSolve()

    st = CPUTime() - st
    WRITE( Message, * ) 'Solve (s)             :',st
    CALL Info( 'StatCurrentSolve', Message, Level=5 )
    
    ! Log electrode circuit solution
    CALL LogElectrodeCktSolution(Solver, Potential, PotentialPerm, NPhi, NumElectrodePairs)

!------------------------------------------------------------------------------
!    Compute the electric field from the potential: E = -grad Phi
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Compute the volume current from generalized Ohm law
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Compute the Joule heating from the Hall/EMF current model
!------------------------------------------------------------------------------
    IF ( Control .OR. CalculateCurrent .OR. CalculateHeating .OR. &
        CalculateNodalHeating ) THEN
      CALL GeneralCurrent( Model, Potential, PotentialPerm )
      
      ! Check current at electrode boundaries (only on last iteration)
      IF (CalculateCurrent .AND. iter == NonlinearIter) THEN
        CALL DiagnoseElectrodeCurrents(Model, Solver, VolCurrent, PotentialPerm, DIM)
        CALL DiagnoseBulkVsBoundaryCurrents(Model, Solver, VolCurrent, PotentialPerm, DIM)
      END IF

      WRITE( Message, * ) 'Total Heating Power   :', Heatingtot
      CALL Info( 'StatCurrentSolve', Message, Level=4 )

      ! Power balance: int J.(U x B) dV = Joule heating + power into the loads.
      ! Load power is sum (Vp - Vm) I = sum R I^2 from the circuit unknowns,
      ! which are owned by rank 0.
      PowerLoad = 0.0_dp
      Var => VariableGet( Solver % Mesh % Variables, 'Electrode Circuit Values' )
      IF ( ASSOCIATED(Var) ) THEN
        IF ( SIZE(Var % Values) >= 3*NumElectrodePairs ) THEN
          DO k = 1, NumElectrodePairs
            PowerLoad = PowerLoad + ( Var % Values(3*(k-1)+1) - Var % Values(3*(k-1)+2) ) &
                * Var % Values(3*(k-1)+3)
          END DO
        END IF
      END IF
      WRITE( Message, '(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A)' ) &
          'Power balance: P_emf ', PowerEmf, '  P_joule ', HeatingTot, &
          '  P_load ', PowerLoad, '  residual ', PowerEmf - HeatingTot - PowerLoad, ' W'
      CALL Info( 'StatCurrentSolve', Message, Level=4 )
      IF( VolTot > 0.0_dp ) THEN
        WRITE( Message, '(A,ES11.3,A,ES11.3,A,ES11.3)' ) &
            'Volume ', VolTot, ' m^3  mean Ux ', UxAvg / VolTot, &
            ' m/s  mean Jy ', CurrYAvg / VolTot
        CALL Info( 'StatCurrentSolve', Message, Level=4 )
      END IF
      CALL ListAddConstReal( Model % Simulation, &
          'RES: Total Joule Heating', Heatingtot )

      PotDiff = DirichletDofsRange( Solver )

      IF( PotDiff > 0 ) THEN
        Resistance = PotDiff**2 / HeatingTot
        WRITE( Message, * ) 'Effective Resistance  :', Resistance
        CALL Info( 'StatCurrentSolve', Message, Level=4 )
        CALL ListAddConstReal( Model % Simulation, &
            'RES: Effective Resistance', Resistance )
      END IF
    END IF

    IF(Control ) THEN
      WRITE( Message, * ) 'Total Volume          :', VolTot
      CALL Info( 'StatCurrentSolve', Message, Level=4 )

      ControlScaling = 1.0_dp
      IF( ControlPower ) THEN
        ControlScaling = SQRT( ControlTarget / HeatingTot )
      ELSE IF( ControlCurrent ) THEN
        IF( PotDiff > 0.0d0 ) THEN
          CurrentTot = HeatingTot / PotDiff
          ControlScaling = ControlTarget / CurrentTot
          WRITE( Message, * ) 'Total Current         :', CurrentTot
          CALL Info( 'StatCurrentSolve', Message, Level=4 )
          CALL ListAddConstReal( Model % Simulation, &
              'RES: TotalCurrent', CurrentTot )
        ELSE
          CALL Warn('StatCurrentSolver','Current cannot be determined without pot. difference')
        END IF
      END IF

      WRITE( Message, * ) 'Control Scaling       :', ControlScaling
      CALL Info( 'StatCurrentSolve', Message, Level=4 )
      CALL ListAddConstReal( Model % Simulation, &
          'RES: CurrentSolver Scaling', ControlScaling )
      Potential = ControlScaling * Potential
      ! Solver % Variable % Norm = ControlScaling * Solver % Variable % Norm

      IF ( CalculateHeating )     Heating = ControlScaling**2 * Heating
      IF ( CalculateNodalHeating) NodalHeating = ControlScaling**2 * NodalHeating
      IF ( CalculateCurrent )     VolCurrent = ControlScaling * VolCurrent
    END IF


    ! Converged when both the potential and the electron temperature settle
    IF( Solver % Variable % NonlinConverged > 0 .AND. TeChange < TeTol ) EXIT
  END DO

  CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Potential')

  IF ( CalculateCurrent ) THEN
    CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Volume Current')
  END IF

  IF ( CalculateHeating ) THEN
    CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Joule Heating')
  END IF

  IF ( CalculateNodalHeating ) THEN
    CALL InvalidateVariable( Model % Meshes, Solver % Mesh, &
        'Nodal Joule Heating')
  END IF

  ! Deallocate electrode stuff
  IF (ALLOCATED(ElectrodePairOfBC)) DEALLOCATE(ElectrodePairOfBC)
  IF (ALLOCATED(ElectrodeSignOfBC)) DEALLOCATE(ElectrodeSignOfBC)
  IF (ALLOCATED(ElectrodeResistance)) DEALLOCATE(ElectrodeResistance)

  CALL DefaultFinish()

  CONTAINS

!------------------------------------------------------------------------------
!> Evaluate the seeded-plasma state at every node. Only the alkali seed
!> ionizes; the carrier gas (e.g. argon) is treated as fully neutral.
!> Heavy particles are at the gas temperature Tg, electrons at Te.
!>
!>   n     = (p + p_ref) / (kB Tg)                  heavy-particle density
!>   n_s   = x_s n                                  seed density
!>   ne^2 / (n_s - ne) = S(Te)                      two-temperature Saha
!>   S     = 2 (g_i/g_n) (2 pi me kB Te / h^2)^(3/2) exp(-chi / (kB Te))
!>   nu_c  = vth (n - n_s) Q_c,  nu_s = vth (n_s - ne) Q_s,  nu = nu_c + nu_s
!>   mu_e  = e / (me nu),   sigma = e ne mu_e,   vth = sqrt(8 kB Te / (pi me))
!>
!> Two-temperature mode (Kerrebrock): Joule heating of the electrons balances
!> their elastic collisional losses to heavy particles,
!>   J^2/sigma = 3 delta ne me kB (Te - Tg) sum_s nu_s / M_s,
!> which with E = |J|/sigma becomes
!>   Te - Tg = e^2 E^2 / (3 delta kB me^2 nu sum_s nu_s / M_s).
!> E is taken from the previous current solution, the root is found by
!> bisection, and the update is under-relaxed. Equilibrium mode sets Te = Tg.
!>
!> Fills Electric Conductivity (clamped to [Sigma Min, Sigma Max]), Ionization
!> Fraction (ne/n), Electron Temperature, Electron Density and Electron
!> Mobility. The mobility sets the Hall coefficient, see HallCoefficient.
!> MaxRelChange returns the largest relative change of Te over all nodes.
!------------------------------------------------------------------------------
SUBROUTINE UpdateSeededPlasma( RmsRelChange, FirstIteration )
  !> Volume-weighted RMS of the relative electron temperature change. This,
  !> not the largest nodal change, decides convergence: the singular current
  !> concentration at electrode edges keeps a handful of cells changing, and
  !> letting them decide stalls every update at the iteration limit on fine
  !> meshes while the solution elsewhere has long converged.
  REAL(KIND=dp), INTENT(OUT) :: RmsRelChange
  LOGICAL, INTENT(IN) :: FirstIteration
  REAL(KIND=dp) :: MaxRelChange, RelChange, SumW, SumWC2

  TYPE(ValueList_t), POINTER :: Mat
  TYPE(SeedPlasma_t) :: Pl
  INTEGER :: i, ipT, ipP, ipS, ipX, ipJ, matId
  INTEGER :: nOwned, nSigClamped, nTeClamped
  REAL(KIND=dp) :: Pref, TeRelax, TeMax, SigmaMin, SigmaMax
  REAL(KIND=dp) :: Tg, Te, TeOld, TeStep, Pabs, nHeavy, nSeed
  REAL(KIND=dp) :: WorstChange, WorstPos(3)
  INTEGER :: nDamped
  REAL(KIND=dp) :: ne, nu, NuOverMass, mu, sigma
  REAL(KIND=dp) :: Ep(3), Bn(3), Bmag, Epar2, Eperp2
  REAL(KIND=dp) :: IonMin, IonMax, SigMin, SigMax, TeMinSeen, TeMaxSeen, DTeMax
  LOGICAL :: Found, TwoTemperature, Owned, AtMax

  REAL(KIND=dp), PARAMETER :: NA = 6.02214076d23

  ! Single conducting body: use the material of Body 1
  matId = ListGetInteger(Model % Bodies(1) % Values, 'Material', Found, &
      minv=1, maxv=Model % NumberOfMaterials)
  IF (.NOT. Found) CALL Fatal('UpdateSeededPlasma','Could not get Material id from Body 1')
  Mat => Model % Materials(matId) % Values

  Pl % SeedFrac    = RequiredMaterialReal(Mat, 'Seed Mole Fraction')
  Pl % ChiJ        = RequiredMaterialReal(Mat, 'Seed Ionization Energy') * eCharge   ! eV -> J
  Pl % WeightRatio = RequiredMaterialReal(Mat, 'Seed Statistical Weight Ratio')
  Pl % Qseed       = RequiredMaterialReal(Mat, 'Seed Electron Neutral Cross Section')
  Pl % Qcarrier    = RequiredMaterialReal(Mat, 'Carrier Electron Neutral Cross Section')
  Pl % Mseed       = 1.0_dp
  Pl % Mcarrier    = 1.0_dp
  Pl % LossFactor  = 1.0_dp

  IF (Pl % SeedFrac <= 0.0_dp .OR. Pl % SeedFrac >= 1.0_dp) &
      CALL Fatal('UpdateSeededPlasma','Seed Mole Fraction must be in (0, 1)')
  IF (Pl % WeightRatio <= 0.0_dp) &
      CALL Fatal('UpdateSeededPlasma','Seed Statistical Weight Ratio must be positive')
  IF (Pl % Qseed < 0.0_dp .OR. Pl % Qcarrier < 0.0_dp .OR. Pl % Qseed + Pl % Qcarrier <= 0.0_dp) &
      CALL Fatal('UpdateSeededPlasma','Electron-neutral cross sections must be non-negative and not both zero')

  TwoTemperature = ListGetLogical(Mat, 'Two Temperature', Found)
  IF (TwoTemperature) THEN
    IF (.NOT. CalculateCurrent) CALL Fatal('UpdateSeededPlasma', &
        'Two Temperature requires Calculate Volume Current = True')
    Pl % Mseed      = RequiredMaterialReal(Mat, 'Seed Molar Mass') * 1.0e-3_dp / NA      ! g/mol -> kg
    Pl % Mcarrier   = RequiredMaterialReal(Mat, 'Carrier Molar Mass') * 1.0e-3_dp / NA  ! g/mol -> kg
    Pl % LossFactor = RequiredMaterialReal(Mat, 'Electron Energy Loss Factor')
    TeMax   = RequiredMaterialReal(Mat, 'Electron Temperature Max')
    TeRelax = RequiredMaterialReal(Mat, 'Electron Temperature Relaxation')
    IF (Pl % Mseed <= 0.0_dp .OR. Pl % Mcarrier <= 0.0_dp .OR. Pl % LossFactor <= 0.0_dp .OR. &
        TeMax <= 0.0_dp .OR. TeRelax <= 0.0_dp .OR. TeRelax > 1.0_dp) &
        CALL Fatal('UpdateSeededPlasma', 'Invalid two-temperature parameters')
  END IF

  ! OpenFOAM pressure is gauge; this is the absolute pressure it is relative to
  Pref = ListGetCReal(Mat, 'Reference Pressure', Found)
  IF (.NOT. Found) Pref = 101325.0_dp
  SigmaMin = ListGetCReal(Mat, 'Sigma Min', Found)
  IF (.NOT. Found) SigmaMin = 1.0d-2
  SigmaMax = ListGetCReal(Mat, 'Sigma Max', Found)
  IF (.NOT. Found) SigmaMax = 1.0d6

  IonMin = HUGE(1.0_dp); IonMax = 0.0_dp
  SigMin = HUGE(1.0_dp); SigMax = 0.0_dp
  TeMinSeen = HUGE(1.0_dp); TeMaxSeen = 0.0_dp
  DTeMax = 0.0_dp
  MaxRelChange = 0.0_dp
  RmsRelChange = 0.0_dp
  SumW = 0.0_dp; SumWC2 = 0.0_dp
  nOwned = 0; nSigClamped = 0; nTeClamped = 0
  nDamped = 0
  WorstChange = -1.0_dp; WorstPos = 0.0_dp

  IF (TwoTemperature) THEN
    IF (.NOT. ALLOCATED(TeNodeRelax)) THEN
      ALLOCATE( TeNodeRelax(SIZE(ElecTemp)), TeLastStep(SIZE(ElecTemp)) )
    ELSE IF (SIZE(TeNodeRelax) /= SIZE(ElecTemp)) THEN
      DEALLOCATE( TeNodeRelax, TeLastStep )
      ALLOCATE( TeNodeRelax(SIZE(ElecTemp)), TeLastStep(SIZE(ElecTemp)) )
    END IF
    IF (FirstIteration) THEN
      TeNodeRelax = TeRelax
      TeLastStep = 0.0_dp
    END IF
    IF (.NOT. ALLOCATED(TeNodeWeight)) THEN
      CALL ComputeNodeWeights()
    ELSE IF (SIZE(TeNodeWeight) /= SIZE(ElecTemp)) THEN
      DEALLOCATE( TeNodeWeight )
      CALL ComputeNodeWeights()
    END IF
  END IF

  DO i = 1, Solver % Mesh % NumberOfNodes
    ipT = TgasPerm(i)
    ipP = PPerm(i)
    ipS = SigPerm(i)
    ipX = PlasmaPerm(i)
    IF (ipT <= 0 .OR. ipP <= 0 .OR. ipS <= 0 .OR. ipX <= 0) CYCLE

    ! Count each node once across partitions (the first neighbour owns it)
    Owned = .TRUE.
    IF (ParEnv % PEs > 1) Owned = &
        Solver % Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE

    Tg   = TgasVals(ipT)
    Pabs = PVals(ipP) + Pref

    ! Non-physical input (T <= 0, NaN, or negative absolute pressure):
    ! non-conducting, no Hall effect
    IF (.NOT. (Tg > 0.0_dp) .OR. .NOT. (Pabs > 0.0_dp)) THEN
      Te = MAX(Tg, 0.0_dp)
      ne = 0.0_dp
      nHeavy = 1.0_dp
      mu = 0.0_dp
      sigma = SigmaMin
    ELSE
      nHeavy = Pabs / (kBoltz * Tg)
      nSeed  = Pl % SeedFrac * nHeavy
      Te = Tg

      IF (TwoTemperature) THEN
        ! Field seen by the electrons, E' = -grad(phi) + U x B, from the
        ! previous potential solution, split along and across B
        Ep = 0.0_dp
        ipJ = PotentialPerm(i)
        IF (ipJ > 0) Ep(1:DIM) = EffField(DIM*(ipJ-1)+1 : DIM*(ipJ-1)+DIM)
        Bn = 0.0_dp
        IF (BxPerm(i) > 0) Bn(1) = BxVals(BxPerm(i))
        IF (ByPerm(i) > 0) Bn(2) = ByVals(ByPerm(i))
        IF (BzPerm(i) > 0) Bn(3) = BzVals(BzPerm(i))
        Bmag = SQRT(SUM(Bn**2))
        Epar2 = 0.0_dp
        IF (Bmag > 0.0_dp) Epar2 = (SUM(Ep*Bn) / Bmag)**2
        Eperp2 = MAX(SUM(Ep**2) - Epar2, 0.0_dp)

        Te = ElectronTemperature( Pl, Tg, nHeavy, nSeed, Epar2, Eperp2, Bmag, TeMax, AtMax )
        IF (AtMax .AND. Owned) nTeClamped = nTeClamped + 1

        ! Under-relax against the previous iterate (Tg on the first call).
        ! A node whose step reverses direction is oscillating around its fixed
        ! point (strong Te-sigma-current feedback); halve its relaxation so the
        ! oscillation decays instead of stalling the whole nonlinear solve.
        TeOld = ElecTemp(ipX)
        IF (TeOld <= 0.0_dp) TeOld = Tg
        TeStep = Te - TeOld
        ! Only reversals larger than the convergence tolerance count: converged
        ! nodes flip sign at noise level and must not be slowed down
        IF (TeStep * TeLastStep(ipX) < 0.0_dp .AND. ABS(TeStep) > TeTol * TeOld) &
            TeNodeRelax(ipX) = MAX(0.5_dp * TeNodeRelax(ipX), 0.1_dp * TeRelax)
        TeLastStep(ipX) = TeStep
        IF (TeNodeRelax(ipX) < TeRelax .AND. Owned) nDamped = nDamped + 1

        ! Convergence is judged on the applied (relaxed) change; the damping
        ! floor of 0.1 x the base relaxation bounds how much it can understate
        ! the distance to the energy-balance solution
        IF (Owned) THEN
          RelChange = TeNodeRelax(ipX) * ABS(TeStep) / TeOld
          SumW = SumW + TeNodeWeight(ipX)
          SumWC2 = SumWC2 + TeNodeWeight(ipX) * RelChange**2
          IF (RelChange > MaxRelChange) THEN
            MaxRelChange = RelChange
            IF (MaxRelChange > WorstChange) THEN
              WorstChange = MaxRelChange
              WorstPos = (/ Solver % Mesh % Nodes % x(i), Solver % Mesh % Nodes % y(i), &
                            Solver % Mesh % Nodes % z(i) /)
            END IF
          END IF
        END IF
        Te = MAX(Tg, TeOld + TeNodeRelax(ipX) * TeStep)
      END IF

      CALL SeedPlasmaState( Pl, Te, nHeavy, nSeed, ne, nu, NuOverMass )
      IF (nu > 0.0_dp) THEN
        mu    = eCharge / (eMass * nu)
        sigma = eCharge * ne * mu
      ELSE
        mu    = 0.0_dp
        sigma = SigmaMin
      END IF
    END IF

    IF (Owned) THEN
      nOwned = nOwned + 1
      IF (sigma < SigmaMin .OR. sigma > SigmaMax) nSigClamped = nSigClamped + 1
    END IF
    sigma = MIN(MAX(sigma, SigmaMin), SigmaMax)

    SigVals(ipS)  = sigma
    IonFrac(ipX)  = ne / nHeavy
    ElecTemp(ipX) = Te
    ElecDens(ipX) = ne
    ElecMob(ipX)  = mu

    IonMin = MIN(IonMin, IonFrac(ipX)); IonMax = MAX(IonMax, IonFrac(ipX))
    SigMin = MIN(SigMin, sigma);        SigMax = MAX(SigMax, sigma)
    TeMinSeen = MIN(TeMinSeen, Te);     TeMaxSeen = MAX(TeMaxSeen, Te)
    DTeMax = MAX(DTeMax, Te - Tg)
  END DO

  IF (ParEnv % PEs > 1) THEN
    IonMin = ParallelReduction(IonMin, 1)
    IonMax = ParallelReduction(IonMax, 2)
    SigMin = ParallelReduction(SigMin, 1)
    SigMax = ParallelReduction(SigMax, 2)
    TeMinSeen = ParallelReduction(TeMinSeen, 1)
    TeMaxSeen = ParallelReduction(TeMaxSeen, 2)
    DTeMax = ParallelReduction(DTeMax, 2)
    MaxRelChange = ParallelReduction(MaxRelChange, 2)
    SumW = ParallelReduction(SumW)
    SumWC2 = ParallelReduction(SumWC2)
    nOwned = NINT(ParallelReduction(REAL(nOwned, dp)))
    nSigClamped = NINT(ParallelReduction(REAL(nSigClamped, dp)))
    nTeClamped = NINT(ParallelReduction(REAL(nTeClamped, dp)))
    nDamped = NINT(ParallelReduction(REAL(nDamped, dp)))
  END IF

  IF (SumW > 0.0_dp) RmsRelChange = SQRT(SumWC2 / SumW)

  WRITE(Message,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A,I0,A)') &
      'Ionization fraction [', IonMin, ', ', IonMax, &
      ']  sigma [', SigMin, ', ', SigMax, '] S/m, clamped at ', &
      nSigClamped, ' of ', nOwned, ' nodes'
  CALL Info('UpdateSeededPlasma', Message, Level=4)
  IF (TwoTemperature) THEN
    WRITE(Message,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A,I0,A)') &
        'Te [', TeMinSeen, ', ', TeMaxSeen, '] K  max(Te-Tg) ', DTeMax, &
        ' K  max rel. change ', MaxRelChange, '  capped at Te max: ', nTeClamped, &
        ' nodes  damped: ', nDamped, ' nodes'
    CALL Info('UpdateSeededPlasma', Message, Level=4)
    WRITE(Message,'(A,ES10.3,A,ES10.3,A)') &
        'Te convergence: volume-weighted RMS change ', RmsRelChange, &
        ' (tolerance ', TeTol, ')'
    CALL Info('UpdateSeededPlasma', Message, Level=4)
    IF (WorstChange > 0.0_dp .AND. MaxRelChange > 0.05_dp) THEN
      WRITE(Message,'(A,ES10.3,A,3F8.2,A)') '  largest change on this rank ', WorstChange, &
          ' at (', 1.0e3_dp * WorstPos, ') mm'
      CALL Info('UpdateSeededPlasma', Message, Level=4)
    END IF
  END IF

END SUBROUTINE UpdateSeededPlasma

!> Lumped nodal volumes, integral of each nodal basis function over the
!> local active elements, stored with the plasma variable permutation
SUBROUTINE ComputeNodeWeights()
  TYPE(Element_t), POINTER :: Elem
  TYPE(Nodes_t) :: EN
  TYPE(GaussIntegrationPoints_t) :: GIP
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3), DetJ
  INTEGER :: t, g, k, nn, ipW
  LOGICAL :: stat

  ALLOCATE( TeNodeWeight(SIZE(ElecTemp)) )
  TeNodeWeight = 0.0_dp
  ALLOCATE( EN % x(Model % MaxElementNodes), EN % y(Model % MaxElementNodes), &
            EN % z(Model % MaxElementNodes) )

  DO t = 1, Solver % NumberOfActiveElements
    Elem => Solver % Mesh % Elements( Solver % ActiveElements(t) )
    nn = Elem % TYPE % NumberOfNodes
    EN % x(1:nn) = Solver % Mesh % Nodes % x(Elem % NodeIndexes(1:nn))
    EN % y(1:nn) = Solver % Mesh % Nodes % y(Elem % NodeIndexes(1:nn))
    EN % z(1:nn) = Solver % Mesh % Nodes % z(Elem % NodeIndexes(1:nn))
    GIP = GaussPoints( Elem )
    DO g = 1, GIP % n
      stat = ElementInfo( Elem, EN, GIP % u(g), GIP % v(g), GIP % w(g), DetJ, Basis, dBasisdx )
      DO k = 1, nn
        ipW = PlasmaPerm(Elem % NodeIndexes(k))
        IF (ipW > 0) TeNodeWeight(ipW) = TeNodeWeight(ipW) + DetJ * GIP % s(g) * Basis(k)
      END DO
    END DO
  END DO

  DEALLOCATE( EN % x, EN % y, EN % z )
END SUBROUTINE ComputeNodeWeights


!------------------------------------------------------------------------------
!> Isotropic nodal conductivity of an element, read straight from the
!> Electric Conductivity variable that UpdateSeededPlasma fills. Going through
!> the Material keyword instead evaluates its MATC expression for every node of
!> every element, which dominated the assembly and current-evaluation time.
!------------------------------------------------------------------------------
SUBROUTINE NodalConductivity( Cond, NodeIndexes, n )
  REAL(KIND=dp) :: Cond(:,:,:)
  INTEGER :: NodeIndexes(:), n
  INTEGER :: p, d, ipS

  Cond = 0.0_dp
  DO p = 1, n
    ipS = SigPerm(NodeIndexes(p))
    IF (ipS <= 0) CYCLE
    DO d = 1, 3
      Cond(d,d,p) = SigVals(ipS)
    END DO
  END DO
END SUBROUTINE NodalConductivity


!------------------------------------------------------------------------------
!> Integration points for the potential assembly and current evaluation:
!> TetraPoints for linear tetrahedra, Elmer's default rule otherwise (a single
!> point would admit hourglass modes on hexahedra).
!------------------------------------------------------------------------------
FUNCTION ElementGaussPoints( Elem ) RESULT( IP )
  TYPE(Element_t) :: Elem
  TYPE(GaussIntegrationPoints_t) :: IP

  IF ( Elem % TYPE % ElementCode == 504 ) THEN
    IP = GaussPoints( Elem, TetraPoints )
  ELSE
    IP = GaussPoints( Elem )
  END IF
END FUNCTION ElementGaussPoints


FUNCTION RequiredMaterialReal(Mat, Name) RESULT(Val)
  TYPE(ValueList_t), POINTER :: Mat
  CHARACTER(LEN=*) :: Name
  REAL(KIND=dp) :: Val
  LOGICAL :: Found

  Val = ListGetCReal(Mat, Name, Found)
  IF (.NOT. Found) CALL Fatal('UpdateSeededPlasma', 'Missing Material keyword: '//TRIM(Name))
END FUNCTION RequiredMaterialReal


!------------------------------------------------------------------------------
!> Hall coefficient 1/(ne e) of the generalized Ohm's law
!>   E + U x B = eta J + (1/(ne e)) J x B
!> written as mu_e * eta, which equals 1/(ne e) when sigma = e ne mu_e and
!> keeps the Hall parameter beta = mu_e |B| physical where sigma is clamped.
!------------------------------------------------------------------------------
FUNCTION HallCoefficient(Basis, NodeIndexes, n, Eta) RESULT(Alpha)
  REAL(KIND=dp) :: Basis(:), Eta
  INTEGER :: NodeIndexes(:), n
  REAL(KIND=dp) :: Alpha, MobGp
  INTEGER :: i, ip

  MobGp = 0.0_dp
  DO i = 1, n
    ip = PlasmaPerm(NodeIndexes(i))
    IF (ip > 0) MobGp = MobGp + Basis(i) * ElecMob(ip)
  END DO
  Alpha = MobGp * Eta
END FUNCTION HallCoefficient

!------------------------------------------------------------------------------
!> Compute the Current and Joule Heating at model nodes.
!------------------------------------------------------------------------------
  SUBROUTINE GeneralCurrent( Model, Potential, Reorder )
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: Potential(:)
    INTEGER :: Reorder(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes 
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
    REAL(KIND=dp), ALLOCATABLE :: SumOfWeights(:), tmp(:)
    REAL(KIND=dp) :: Conductivity(3,3,Model % MaxElementNodes)
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes)
    REAL(KIND=dp) :: dBasisdx(Model % MaxElementNodes,3)
    REAL(KIND=DP) :: SqrtElementMetric, ElemVol
    REAL(KIND=dp) :: ElementPot(Model % MaxElementNodes)
    REAL(KIND=dp) :: Current(3), EStar(3)
    REAL(KIND=dp) :: s, ug, vg, wg, Grad(3)
    REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
    REAL(KIND=dp) :: HeatingDensity, x, y, z
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: N_Integ, t, tg, i, j, k
    LOGICAL :: Stat

    REAL(KIND=dp) :: Ugp(3), Bgp(3), UxBgp(3)
    INTEGER :: ip
    REAL(KIND=dp) :: Cgp(3,3)

    REAL(KIND=dp) :: HallCoeffAlpha, EtaGP, SigmaIso, JouleGp
    REAL(KIND=dp) :: RHS(3), Jgp(3)
    REAL(KIND=dp) :: M(3,3), Minv(3,3)

    ALLOCATE( Nodes % x( Model % MaxElementNodes ) )
    ALLOCATE( Nodes % y( Model % MaxElementNodes ) )
    ALLOCATE( Nodes % z( Model % MaxElementNodes ) )

    IF( CalculateHeating .OR. CalculateCurrent ) THEN
      ALLOCATE( SumOfWeights( Model % NumberOfNodes ) )
      SumOfWeights = 0.0d0
    END IF

    HeatingTot = 0.0d0
    PowerEmf = 0.0d0
    UxAvg = 0.0d0
    CurrYAvg = 0.0d0
    VolTot = 0.0d0
    IF ( CalculateHeating )  Heating = 0.0d0
    IF ( CalculateNodalHeating)  NodalHeating = 0.0d0
    IF ( CalculateCurrent )  VolCurrent = 0.0d0
    IF ( CalculateCurrent )  EffField = 0.0d0

    IF (CalculateCurrent) THEN
      IF (.NOT. ASSOCIATED(VolCurrent)) THEN
        CALL Fatal('GeneralCurrent','DBG: VolCurrent NA')
      END IF
      IF (MOD(SIZE(VolCurrent), DIM) /= 0) THEN
        CALL Fatal('GeneralCurrent','DBG: VC bad size')
      END IF
    END IF


    IF( GetCondAtIp ) THEN
      CALL ListInitElementKeyword( CondAtIp_h,'Material','Electric Conductivity')
    END IF
     
!------------------------------------------------------------------------------
!   Go through model elements, we will compute on average of elementwise
!   fluxes to nodes of the model
!------------------------------------------------------------------------------
    DO t = 1,Solver % NumberOfActiveElements
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where electrostatics
!        should be calculated
!------------------------------------------------------------------------------
       Element => Solver % Mesh % Elements( Solver % ActiveElements( t ) )
       Model % CurrentElement => Element
       NodeIndexes => Element % NodeIndexes

       IF ( Element % PartIndex /= ParEnv % MyPE ) CYCLE

       n = Element % TYPE % NumberOfNodes

       IF ( ANY(Reorder(NodeIndexes) == 0) ) CYCLE

       ElementPot(1:n) = Potential( Reorder( NodeIndexes(1:n) ) )
       
       Nodes % x(1:n) = Model % Nodes % x( NodeIndexes )
       Nodes % y(1:n) = Model % Nodes % y( NodeIndexes )
       Nodes % z(1:n) = Model % Nodes % z( NodeIndexes )

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
       IntegStuff = ElementGaussPoints( Element )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------

       IF( .NOT. GetCondAtIp ) THEN
         CALL NodalConductivity( Conductivity, NodeIndexes, n )
       END IF
         
!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------

       HeatingDensity = 0.0d0
       Current = 0.0d0
       EStar = 0.0d0
       ElemVol = 0.0d0


       DO tg=1,N_Integ

          ug = U_Integ(tg)
          vg = V_Integ(tg)
          wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element, Nodes,ug,vg,wg, &
               SqrtElementMetric,Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
          s = SqrtElementMetric * S_Integ(tg)

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
            x = SUM( Nodes % x(1:n)*Basis(1:n) )
            y = SUM( Nodes % y(1:n)*Basis(1:n) )
            z = SUM( Nodes % z(1:n)*Basis(1:n) )
            
            CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
            s = s * SqrtMetric * 2 * PI
          END IF

!------------------------------------------------------------------------------

          DO j = 1, DIM
            Grad(j) = SUM( dBasisdx(1:n,j) * ElementPot(1:n) )
          END DO

          Ugp = 0.0_dp
          Bgp = 0.0_dp

          DO i = 1, n
            ! Velocity components
            ip = UxPerm(NodeIndexes(i))
            IF (ip > 0) Ugp(1) = Ugp(1) + Basis(i) * UxVals(ip)

            ip = UyPerm(NodeIndexes(i))
            IF (ip > 0) Ugp(2) = Ugp(2) + Basis(i) * UyVals(ip)

            IF (DIM == 3) THEN
              ip = UzPerm(NodeIndexes(i))
              IF (ip > 0) Ugp(3) = Ugp(3) + Basis(i) * UzVals(ip)
            END IF

            ! Magnetic field components
            ip = BxPerm(NodeIndexes(i))
            IF (ip > 0) Bgp(1) = Bgp(1) + Basis(i) * BxVals(ip)

            ip = ByPerm(NodeIndexes(i))
            IF (ip > 0) Bgp(2) = Bgp(2) + Basis(i) * ByVals(ip)

            IF (DIM == 3) THEN
              ip = BzPerm(NodeIndexes(i))
              IF (ip > 0) Bgp(3) = Bgp(3) + Basis(i) * BzVals(ip)
            END IF
          END DO

          ! Cross product U x B
          UxBgp(1) = Ugp(2)*Bgp(3) - Ugp(3)*Bgp(2)
          UxBgp(2) = Ugp(3)*Bgp(1) - Ugp(1)*Bgp(3)
          UxBgp(3) = Ugp(1)*Bgp(2) - Ugp(2)*Bgp(1)


          Cgp = 0.0_dp
          IF( GetCondAtIp ) THEN
            CondAtIp = ListGetElementReal( CondAtIp_h, Basis, Element, Stat, GaussPoint = tg )
            DO i = 1, dim
              Cgp(i,i) = CondAtIp
            END DO
          ELSE
            DO i = 1, dim
              DO j = 1, dim
                Cgp(i,j) = SUM( Conductivity(i,j,1:n) * Basis(1:n) )
              END DO
            END DO
          END IF
        
          ! Caluclate resistivity from conductivity
          SigmaIso = (Cgp(1,1) + Cgp(2,2) + Cgp(3,3)) / REAL(dim,dp)
          IF (SigmaIso > 0.0_dp) THEN
            EtaGP = 1.0_dp / SigmaIso
          ELSE
            EtaGP = 0.0_dp
          END IF
          HallCoeffAlpha = HallCoefficient( Basis, NodeIndexes, n, EtaGP )

          ! Build the M matrix
          M = 0.0_dp
          DO i=1,dim
            M(i,i) = EtaGP
          END DO

          M(1,2) = M(1,2) - HallCoeffAlpha * (-Bgp(3))
          M(1,3) = M(1,3) - HallCoeffAlpha * ( Bgp(2))
          M(2,1) = M(2,1) - HallCoeffAlpha * ( Bgp(3))
          M(2,3) = M(2,3) - HallCoeffAlpha * (-Bgp(1))
          M(3,1) = M(3,1) - HallCoeffAlpha * (-Bgp(2))
          M(3,2) = M(3,2) - HallCoeffAlpha * ( Bgp(1))

          RHS = 0.0_dp
          DO j=1,dim
            RHS(j) = -Grad(j) + UxBgp(j)
          END DO


          VolTot = VolTot + s

          IF( Control .OR. CalculateHeating .OR. CalculateCurrent .OR. CalculateNodalHeating ) THEN
            ! Invert the hall matrix
            CALL Invert3x3(M, Minv, Stat)
              
            IF (.NOT. Stat) THEN
              WRITE(*,*) 'Hall matrix inversion failed at Gauss point'
              CALL Fatal( &
                'GeneralCurrent / Hall MHD', &
                'Hall conductivity matrix is singular or ill-conditioned at Gauss point.' )
            END IF
            Jgp = 0.0_dp
            DO i=1,3
              DO j=1,3
                Jgp(i) = Jgp(i) + Minv(i,j) * RHS(j)
              END DO
            END DO

            ! Resistive Joule heating for generalized Ohm law:
            ! J·(E + UxB) = eta*|J|^2 since Hall term is non-dissipative.
            JouleGp = EtaGP * SUM( Jgp(1:DIM) * Jgp(1:DIM) )
            HeatingTot = HeatingTot + s * JouleGp
            PowerEmf = PowerEmf + s * SUM( Jgp(1:DIM) * UxBgp(1:DIM) )
            UxAvg = UxAvg + s * Ugp(1)
            CurrYAvg = CurrYAvg + s * Jgp(2)
            HeatingDensity = HeatingDensity + s * JouleGp

            DO j=1,dim
              Current(j) = Current(j) + Jgp(j) * s
              EStar(j) = EStar(j) + RHS(j) * s
            END DO

            ElemVol = ElemVol + s
          END IF


       END DO! of the Gauss integration points

!------------------------------------------------------------------------------
!   Element-wise values for OpenFOAM, before the nodal averaging below
!------------------------------------------------------------------------------
       IF ( ASSOCIATED(ElemPerm) .AND. ElemVol > 0.0d0 ) THEN
         i = Solver % ActiveElements(t)
         IF ( i <= SIZE(ElemPerm) ) THEN
           ip = ElemPerm(i)
           IF ( ip > 0 ) THEN
             IF ( CalculateCurrent ) THEN
               ElemCurr1(ip) = Current(1) / ElemVol
               ElemCurr2(ip) = Current(2) / ElemVol
               IF ( DIM == 3 ) ElemCurr3(ip) = Current(3) / ElemVol
             END IF
             IF ( CalculateHeating ) ElemHeating(ip) = HeatingDensity / ElemVol
           END IF
         END IF
       END IF

!------------------------------------------------------------------------------
!   Weight with element area if required
!------------------------------------------------------------------------------

       IF( CalculateHeating .OR. CalculateCurrent ) THEN
         IF ( ConstantWeights ) THEN
           HeatingDensity = HeatingDensity / ElemVol
           Current(1:Dim) = Current(1:Dim) / ElemVol
           EStar(1:Dim) = EStar(1:Dim) / ElemVol
           SumOfWeights( Reorder( NodeIndexes(1:n) ) ) = &
               SumOfWeights( Reorder( NodeIndexes(1:n) ) ) + 1
         ELSE
           SumOfWeights( Reorder( NodeIndexes(1:n) ) ) = &
               SumOfWeights( Reorder( NodeIndexes(1:n) ) ) + ElemVol
         END IF
       END IF
         
       IF ( CalculateHeating ) THEN
         Heating( Reorder(NodeIndexes(1:n)) ) = &
             Heating( Reorder(NodeIndexes(1:n)) ) + HeatingDensity
       END IF
       
       IF ( CalculateNodalHeating ) THEN
         NodalHeating( Reorder(NodeIndexes(1:n)) ) = &
             NodalHeating( Reorder(NodeIndexes(1:n)) ) + HeatingDensity
       END IF
         
       IF ( CalculateCurrent ) THEN
         DO j=1,DIM 
           VolCurrent(DIM*(Reorder(NodeIndexes(1:n))-1)+j) = &
               VolCurrent(DIM*(Reorder(NodeIndexes(1:n))-1)+j) + &
               Current(j)
           EffField(DIM*(Reorder(NodeIndexes(1:n))-1)+j) = &
               EffField(DIM*(Reorder(NodeIndexes(1:n))-1)+j) + &
               EStar(j)
         END DO
       END IF

    END DO! of the bulk elements

    IF ( CalculateHeating .OR. CalculateCurrent) THEN
      IF ( ParEnv % PEs > 1) THEN
        VolTot     = ParallelReduction(VolTot)
        HeatingTot = ParallelReduction(HeatingTot)
        PowerEmf   = ParallelReduction(PowerEmf)
        UxAvg      = ParallelReduction(UxAvg)
        CurrYAvg   = ParallelReduction(CurrYAvg)
        
        IF ( CalculateCurrent) THEN
          ALLOCATE(tmp(SIZE(VolCurrent)/dim))
          DO i=1,dim
            tmp = VolCurrent(i::dim)
            CALL ParallelSumVector(Solver % Matrix, tmp)
            Volcurrent(i::dim) = tmp
            tmp = EffField(i::dim)
            CALL ParallelSumVector(Solver % Matrix, tmp)
            EffField(i::dim) = tmp
          END DO
        END IF
        IF (CalculateHeating ) CALL ParallelSumVector(Solver % Matrix, Heating)
        CALL ParallelSumVector(Solver % Matrix, SumOfWeights)
      END IF
      
!------------------------------------------------------------------------------
!   Finally, compute average of the fluxes at nodes
!------------------------------------------------------------------------------
      DO i = 1, Model % NumberOfNodes
        IF ( ABS( SumOfWeights(i) ) > 0.0D0 ) THEN
          IF ( CalculateHeating )  Heating(i) = Heating(i) / SumOfWeights(i)
          DO j = 1, DIM
            IF ( CalculateCurrent )  VolCurrent(DIM*(i-1)+j) = &
                VolCurrent(DIM*(i-1)+j) /  SumOfWeights(i)
            IF ( CalculateCurrent )  EffField(DIM*(i-1)+j) = &
                EffField(DIM*(i-1)+j) /  SumOfWeights(i)
          END DO
        END IF
      END DO
      DEALLOCATE( SumOfWeights ) 
    END IF
      
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

!------------------------------------------------------------------------------
   END SUBROUTINE GeneralCurrent
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
    SUBROUTINE StatCurrentCompose( StiffMatrix,Force,Conductivity, &
                            Load,Element,n,Nodes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: StiffMatrix(:,:),Force(:),Load(:), Conductivity(:,:,:)
      INTEGER :: n
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
      REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L,C(3,3),x,y,z
      LOGICAL :: Stat

      INTEGER :: i,j,p,q,t,DIM
 
      TYPE(GaussIntegrationPoints_t) :: IntegStuff

      REAL(KIND=dp) :: Ugp(3), Bgp(3), UxBgp(3)
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER :: ip

      REAL(KIND=dp) :: HallCoeffAlpha, EtaGP, SigmaIso
      REAL(KIND=dp) :: M(3,3), Minv(3,3)

      ! Guard against element is not part of this rank
      IF ( Element % PartIndex /= ParEnv % MyPE ) THEN
        Force = 0.0_dp
        StiffMatrix = 0.0_dp
        RETURN
      END IF

!------------------------------------------------------------------------------
      DIM = CoordinateSystemDimension()

      Force = 0.0d0
      StiffMatrix = 0.0d0
!------------------------------------------------------------------------------

      NodeIndexes => Element % NodeIndexes
 
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
      IntegStuff = ElementGaussPoints( Element )

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                  Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
          y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
          z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
        END IF

        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )

        S = S * SqrtElementMetric * SqrtMetric

        L = SUM( Load(1:n) * Basis )

        IF( GetCondAtIp ) THEN
          CondAtIp = ListGetElementReal( CondAtIp_h, Basis, Element, Stat, GaussPoint = t )
          C(1:dim,1:dim) = 0.0_dp
          DO i=1,dim
            C(i,i) = CondAtIp
          END DO
        ELSE
          DO i=1,DIM
            DO j=1,DIM
              C(i,j) = SUM( Conductivity(i,j,1:n) * Basis(1:n) )
            END DO
          END DO
        END IF

        ! Caluclate resistivity from conductivity
        SigmaIso = (C(1,1) + C(2,2) + C(3,3)) / REAL(dim,dp)
        IF (SigmaIso > 0.0_dp) THEN
          EtaGP = 1.0_dp / SigmaIso
        ELSE
          EtaGP = 0.0_dp
        END IF


        ! Reset Gauss-point velocity and magnetic field
        Ugp = 0.0_dp
        Bgp = 0.0_dp

        DO i = 1, n
          ! Velocity components: use each variable's own perm
          ip = UxPerm(NodeIndexes(i))
          IF (ip > 0) Ugp(1) = Ugp(1) + Basis(i) * UxVals(ip)

          ip = UyPerm(NodeIndexes(i))
          IF (ip > 0) Ugp(2) = Ugp(2) + Basis(i) * UyVals(ip)

          IF (DIM == 3) THEN
            ip = UzPerm(NodeIndexes(i))
            IF (ip > 0) Ugp(3) = Ugp(3) + Basis(i) * UzVals(ip)
          END IF

          ! Magnetic field components: use each variable's own perm
          ip = BxPerm(NodeIndexes(i))
          IF (ip > 0) Bgp(1) = Bgp(1) + Basis(i) * BxVals(ip)

          ip = ByPerm(NodeIndexes(i))
          IF (ip > 0) Bgp(2) = Bgp(2) + Basis(i) * ByVals(ip)

          IF (DIM == 3) THEN
            ip = BzPerm(NodeIndexes(i))
            IF (ip > 0) Bgp(3) = Bgp(3) + Basis(i) * BzVals(ip)
          END IF
        END DO

        HallCoeffAlpha = HallCoefficient( Basis, NodeIndexes, n, EtaGP )

        ! Build the M matrix
        M = 0.0_dp
        DO i=1,dim
          M(i,i) = EtaGP
        END DO

        M(1,2) = M(1,2) - HallCoeffAlpha * (-Bgp(3))
        M(1,3) = M(1,3) - HallCoeffAlpha * ( Bgp(2))
        M(2,1) = M(2,1) - HallCoeffAlpha * ( Bgp(3))
        M(2,3) = M(2,3) - HallCoeffAlpha * (-Bgp(1))
        M(3,1) = M(3,1) - HallCoeffAlpha * (-Bgp(2))
        M(3,2) = M(3,2) - HallCoeffAlpha * ( Bgp(1))

        ! Invert the hall matrix
        CALL Invert3x3(M, Minv, Stat)


        IF (.NOT. Stat) THEN
          WRITE(*,*) 'Hall matrix inversion failed at Gauss point'
          CALL Fatal('GeneralCurrent / Hall MHD', &
                    'Hall conductivity matrix is singular or ill-conditioned at Gauss point.')
        END IF

        ! Cross product UxB = U x B
        UxBgp(1) = Ugp(2)*Bgp(3) - Ugp(3)*Bgp(2)
        UxBgp(2) = Ugp(3)*Bgp(1) - Ugp(1)*Bgp(3)
        UxBgp(3) = Ugp(1)*Bgp(2) - Ugp(2)*Bgp(1)

!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------
        DO p=1,n
          DO q=1,n
            A = 0.d0
            DO i=1,DIM
              DO J=1,DIM
                A = A + dBasisdx(p,i) * Minv(i,j) * dBasisdx(q,j)
              END DO
            END DO
            StiffMatrix(p,q) = StiffMatrix(p,q) + S*A
          END DO
          Force(p) = Force(p) + S*L*Basis(p)
          ! Add forcing terms from the hall and motional EMF
          DO i=1,dim
            DO j=1, dim
              Force(p) = Force(p) + S * dBasisdx(p,i) * Minv(i,j) * UxBgp(j)
            END DO
          END DO
        END DO
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE StatCurrentCompose
!------------------------------------------------------------------------------

  LOGICAL FUNCTION OwnerCircuitRow(gid)
    USE DefUtils
    INTEGER, INTENT(IN) :: gid
    INTEGER :: loc

    ! Circuit rows are [NPhi+1 ... NPhi+NX]
    loc = gid - NPhi
    IF (loc <= 0) THEN
      OwnerCircuitRow = .FALSE.
      RETURN
    END IF

    OwnerCircuitRow = (MOD(loc-1, ParEnv % PEs) == ParEnv % MyPE)
  END FUNCTION OwnerCircuitRow


!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for boundary conditions
!>  of the electrostatic equation. 
!------------------------------------------------------------------------------
  SUBROUTINE StatCurrentBoundary( BoundaryMatrix, BoundaryVector, &
        LoadVector, Element, n, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), LoadVector(:)
    TYPE(Nodes_t)   :: Nodes
    TYPE(Element_t) :: Element
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n)
    REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

    REAL(KIND=dp) :: u,v,w,s,x,y,z
    REAL(KIND=dp) :: Force
    REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

    INTEGER :: t,q,N_Integ

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    LOGICAL :: stat
!------------------------------------------------------------------------------

    BoundaryVector = 0.0d0
    BoundaryMatrix = 0.0d0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
    DO t=1,N_Integ
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivates at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
        y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
        z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
      END IF

      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )

      s = S_Integ(t) * SqrtElementMetric * SqrtMetric
!------------------------------------------------------------------------------
      Force = SUM( LoadVector(1:n)*Basis )

      DO q=1,N
        BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
      END DO
    END DO
  END SUBROUTINE StatCurrentBoundary
!------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !> Electrode circuits as extra unknowns added to the potential system. Pair ep
  !> has electrode voltages Vp, Vm and the load current I through its resistor:
  !>   rows NPhi+3(ep-1)+1 = Vp,  +2 = Vm,  +3 = I.
  !> Each electrode surface is tied to its voltage V by a stiff contact
  !> conductance (Robin condition J.n = Gc (phi - V)), assembled per boundary
  !> element with lumped weights kappa = Gc int psi_p dS:
  !>   phi rows:  + kappa (phi_p - V)                     current into the electrode
  !>   V rows:    sum kappa (V - phi_p) + I = 0  (plus)   Kirchhoff at the electrode
  !>              sum kappa (V - phi_p) - I = 0  (minus)
  !>   I row:     R I - (Vp - Vm) = 0                     Ohm's law of the load
  !> Gc = Electrode Penalty Factor * sigma_max / h, with sigma_max the largest
  !> conductivity in the domain, so the contact is at least that factor stiffer
  !> than any plasma element and each electrode is equipotential to roughly
  !> 1/factor of one element's potential drop. (Scaling by the local wall
  !> conductivity instead makes the contact a real resistor when the gas at the
  !> wall is cold, dissipating power that belongs to the load.) Unlike per-node Lagrange
  !> multipliers this keeps a positive diagonal in every row (no saddle point),
  !> which ILU-preconditioned Krylov solvers handle well, and adds only three
  !> rows per pair.
  !------------------------------------------------------------------------------
  SUBROUTINE BuildElectrodeAddMatrix( Model, Solver, AuxMatrix, &
      PotentialPerm, ElectrodePairOfBC, ElectrodeSignOfBC, &
      ElectrodeResistance, NumElectrodePairs, NPhi )

    USE DefUtils
    USE SolverUtils
    USE ListMatrix
    IMPLICIT NONE

    TYPE(Model_t),  INTENT(IN)    :: Model
    TYPE(Solver_t), INTENT(IN)    :: Solver
    TYPE(Matrix_t), POINTER       :: AuxMatrix
    INTEGER,        INTENT(IN)    :: PotentialPerm(:)
    INTEGER,        INTENT(IN)    :: ElectrodePairOfBC(:)
    INTEGER,        INTENT(IN)    :: ElectrodeSignOfBC(:)
    REAL(dp),       INTENT(IN)    :: ElectrodeResistance(:)
    INTEGER,        INTENT(IN)    :: NumElectrodePairs
    INTEGER,        INTENT(IN)    :: NPhi       ! Local Solver % Matrix % NumberOfRows

    INTEGER :: NX, ep, be, i, k, n, inode, pRow, gp, gidV, gidVp, gidVm, gidI, sideIdx, maxN
    TYPE(Element_t), POINTER :: Elem
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Nodes_t) :: EN
    TYPE(GaussIntegrationPoints_t) :: Integ
    REAL(dp) :: Basis(MAX_ELEMENT_NODES), dBasisdx(MAX_ELEMENT_NODES,3)
    REAL(dp) :: SqrtElementMetric, s, ElemArea, SigRef, Gc, kappa, PenaltyFactor
    REAL(dp), ALLOCATABLE :: SideArea(:), SideKappa(:)
    LOGICAL :: stat, Found

    IF (NumElectrodePairs <= 0) THEN
      AuxMatrix => NULL()
      RETURN
    END IF

    PenaltyFactor = ListGetCReal(Solver % Values, 'Electrode Penalty Factor', Found)
    IF (.NOT. Found) PenaltyFactor = 1.0e3_dp

    SigRef = 0.0_dp
    IF (SIZE(SigVals) > 0) SigRef = MAXVAL(SigVals)
    IF (ParEnv % PEs > 1) SigRef = ParallelReduction(SigRef, 2)
    IF (SigRef <= 0.0_dp) SigRef = 1.0_dp

    NX = 3 * NumElectrodePairs

    AuxMatrix => AllocateMatrix()
    AuxMatrix % FORMAT = MATRIX_LIST
    AuxMatrix % Symmetric = .FALSE.
    AuxMatrix % NumberOfRows = NPhi + NX
    ALLOCATE(AuxMatrix % RHS(AuxMatrix % NumberOfRows))
    AuxMatrix % RHS = 0.0_dp

    ! Extra (circuit) rows are owned by rank 0; ParallelInitMatrix builds the
    ! neighbour lists from RowOwner.
    IF (ParEnv % PEs > 1) THEN
      ALLOCATE(AuxMatrix % RowOwner(NPhi + NX))
      AuxMatrix % RowOwner = 0
    END IF

    maxN = Model % MaxElementNodes
    ALLOCATE( EN % x(maxN), EN % y(maxN), EN % z(maxN) )
    ALLOCATE( SideArea(2*NumElectrodePairs), SideKappa(2*NumElectrodePairs) )
    SideArea = 0.0_dp
    SideKappa = 0.0_dp

    !------------------------------------------------------------
    ! Contact conductance between each electrode surface and its voltage
    !------------------------------------------------------------
    DO be = 1, Solver % Mesh % NumberOfBoundaryElements
      Elem => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements + be )
      NodeIndexes => Elem % NodeIndexes
      n = Elem % TYPE % NumberOfNodes

      DO i = 1, Model % NumberOfBCs
        IF (Elem % BoundaryInfo % Constraint /= Model % BCs(i) % Tag) CYCLE
        ep = ElectrodePairOfBC(i)
        IF (ep <= 0) CYCLE

        IF (ElectrodeSignOfBC(i) == +1) THEN
          gidV = NPhi + 3*(ep-1) + 1
          sideIdx = 2*(ep-1) + 2
        ELSE
          gidV = NPhi + 3*(ep-1) + 2
          sideIdx = 2*(ep-1) + 1
        END IF

        EN % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes(1:n))
        EN % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes(1:n))
        EN % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes(1:n))

        Integ = GaussPoints(Elem)
        ElemArea = 0.0_dp
        DO gp = 1, Integ % n
          stat = ElementInfo(Elem, EN, Integ % u(gp), Integ % v(gp), Integ % w(gp), &
                            SqrtElementMetric, Basis, dBasisdx)
          ElemArea = ElemArea + SqrtElementMetric * Integ % s(gp)
        END DO
        IF (ElemArea <= 0.0_dp) CYCLE

        Gc = PenaltyFactor * SigRef / SQRT(ElemArea)

        DO gp = 1, Integ % n
          stat = ElementInfo(Elem, EN, Integ % u(gp), Integ % v(gp), Integ % w(gp), &
                            SqrtElementMetric, Basis, dBasisdx)
          s = SqrtElementMetric * Integ % s(gp)
          SideArea(sideIdx) = SideArea(sideIdx) + s

          DO inode = 1, n
            pRow = PotentialPerm(NodeIndexes(inode))
            IF (pRow <= 0) CYCLE
            kappa = Gc * s * Basis(inode)
            CALL AddToMatrixElement(AuxMatrix, pRow, pRow,  kappa)
            CALL AddToMatrixElement(AuxMatrix, pRow, gidV, -kappa)
            CALL AddToMatrixElement(AuxMatrix, gidV, pRow, -kappa)
            CALL AddToMatrixElement(AuxMatrix, gidV, gidV,  kappa)
            SideKappa(sideIdx) = SideKappa(sideIdx) + kappa
          END DO
        END DO
      END DO
    END DO

    IF (ParEnv % PEs > 1) THEN
      DO sideIdx = 1, 2*NumElectrodePairs
        SideArea(sideIdx)  = ParallelReduction(SideArea(sideIdx))
        SideKappa(sideIdx) = ParallelReduction(SideKappa(sideIdx))
      END DO
    END IF
    DO sideIdx = 1, 2*NumElectrodePairs
      IF (SideArea(sideIdx) <= 0.0_dp) THEN
        CALL Fatal('BuildElectrodeAddMatrix', 'Electrode side has zero area. Check BC tagging/mesh.')
      END IF
    END DO

    !------------------------------------------------------------
    ! Kirchhoff at the electrodes and Ohm's law of each load (rank 0)
    !------------------------------------------------------------
    IF (ParEnv % MyPE == 0) THEN
      DO ep = 1, NumElectrodePairs
        gidVp = NPhi + 3*(ep-1) + 1
        gidVm = NPhi + 3*(ep-1) + 2
        gidI  = NPhi + 3*(ep-1) + 3
        CALL AddToMatrixElement(AuxMatrix, gidVp, gidI,  1.0_dp)
        CALL AddToMatrixElement(AuxMatrix, gidVm, gidI, -1.0_dp)
        CALL AddToMatrixElement(AuxMatrix, gidI, gidVp, -1.0_dp)
        CALL AddToMatrixElement(AuxMatrix, gidI, gidVm,  1.0_dp)
        CALL AddToMatrixElement(AuxMatrix, gidI, gidI,  ElectrodeResistance(ep))
      END DO

      ! Gauge: ground the minus electrode of pair 1 through a conductance of
      ! the same size as its contact. The total current into the plasma is
      ! zero, so no current flows to ground and Vm(1) = 0 exactly.
      gidVm = NPhi + 2
      CALL AddToMatrixElement(AuxMatrix, gidVm, gidVm, SideKappa(1))
    END IF

    ! Zero diagonals on all ranks keep the extra row structures alive through
    ! the LIST -> CRS conversion; real entries are added on top.
    DO k = NPhi + 1, NPhi + NX
      CALL AddToMatrixElement(AuxMatrix, k, k, 0.0_dp)
    END DO

    CALL List_ToCRSMatrix(AuxMatrix)

    DEALLOCATE( SideArea, SideKappa )
    DEALLOCATE( EN % x, EN % y, EN % z )
  END SUBROUTINE BuildElectrodeAddMatrix

END SUBROUTINE StatCurrentSolver


!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver_post( Model, Solver, dt, Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  ! No-op: postprocessing is done inside StatCurrentSolver (GeneralCurrent).
  RETURN
END SUBROUTINE StatCurrentSolver_post

