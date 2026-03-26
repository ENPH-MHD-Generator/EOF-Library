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
  INTEGER                 :: Dim
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
  END IF
  
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
  USE MHDParams, ONLY: HallCoeffAlphaDefault

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
  TYPE(ValueHandle_t) :: CondAtIp_h
  REAL(KIND=dp) :: CondAtIp

  ! Velocity and Magnetic field values
  REAL(KIND=dp), POINTER :: UxVals(:), UyVals(:), UzVals(:)
  REAL(KIND=dp), POINTER :: BxVals(:), ByVals(:), BzVals(:)
  INTEGER, POINTER :: UxPerm(:), UyPerm(:), UzPerm(:)
  INTEGER, POINTER :: BxPerm(:), ByPerm(:), BzPerm(:)
  TYPE(Variable_t), POINTER :: UxVar, UyVar, UzVar, BxVar, ByVar, BzVar

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
  BxPerm, ByPerm, BzPerm, OldPotential

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
  
  DO iter = 1, NonlinearIter
    at  = CPUTime()
    at0 = RealTime()

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

        k = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
            Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

        !------------------------------------------------------------------------------
        !      Read conductivity values (might be a tensor)
        !------------------------------------------------------------------------------
        
        CALL ListGetRealArray( Model % Materials(k) % Values, &
            'Electric Conductivity', Cwrk, n, NodeIndexes )

        Conductivity = 0.0d0
        IF ( SIZE(Cwrk,1) == 1 ) THEN
          DO i=1,3
            Conductivity( i,i,1:n ) = Cwrk( 1,1,1:n )
          END DO
        ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk,1))
            Conductivity(i,i,1:n) = Cwrk(i,1,1:n)
          END DO
        ELSE
          DO i=1,MIN(3,SIZE(Cwrk,1))
            DO j=1,MIN(3,SIZE(Cwrk,2))
              Conductivity( i,j,1:n ) = Cwrk(i,j,1:n)
            END DO
          END DO
        END IF
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


    ! Standard convergence check
    IF( Solver % Variable % NonlinConverged > 0 ) EXIT
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
    REAL(KIND=dp) :: Current(3)
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

    HallCoeffAlpha = HallCoeffAlphaDefault

    HeatingTot = 0.0d0
    VolTot = 0.0d0
    IF ( CalculateHeating )  Heating = 0.0d0
    IF ( CalculateNodalHeating)  NodalHeating = 0.0d0
    IF ( CalculateCurrent )  VolCurrent = 0.0d0

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
       IntegStuff = GaussPoints( Element )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------

       IF( .NOT. GetCondAtIp ) THEN
         k = ListGetInteger( Model % Bodies( Element % BodyId ) % &
             Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

         CALL ListGetRealArray( Model % Materials(k) % Values, &
             'Electric Conductivity', Cwrk, n, NodeIndexes, gotIt )

         Conductivity = 0.0d0
         IF ( SIZE(Cwrk,1) == 1 ) THEN
           DO i=1,3
             Conductivity( i,i,1:n ) = Cwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Cwrk,1))
             Conductivity(i,i,1:n) = Cwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Cwrk,1))
             DO j=1,MIN(3,SIZE(Cwrk,2))
               Conductivity( i,j,1:n ) = Cwrk(i,j,1:n)
             END DO
           END DO
         END IF
       END IF
         
!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------

       HeatingDensity = 0.0d0
       Current = 0.0d0
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
            HeatingDensity = HeatingDensity + s * JouleGp

            DO j=1,dim
              Current(j) = Current(j) + Jgp(j) * s
            END DO

            ElemVol = ElemVol + s
          END IF


       END DO! of the Gauss integration points

!------------------------------------------------------------------------------
!   Weight with element area if required
!------------------------------------------------------------------------------

       IF( CalculateHeating .OR. CalculateCurrent ) THEN
         IF ( ConstantWeights ) THEN
           HeatingDensity = HeatingDensity / ElemVol
           Current(1:Dim) = Current(1:Dim) / ElemVol
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
         END DO
       END IF

    END DO! of the bulk elements

    IF ( CalculateHeating .OR. CalculateCurrent) THEN
      IF ( ParEnv % PEs > 1) THEN
        VolTot     = ParallelReduction(VolTot)
        HeatingTot = ParallelReduction(HeatingTot)
        
        IF ( CalculateCurrent) THEN
          ALLOCATE(tmp(SIZE(VolCurrent)/dim))
          DO i=1,dim
            tmp = VolCurrent(i::dim)
            CALL ParallelSumVector(Solver % Matrix, tmp)
            Volcurrent(i::dim) = tmp
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
      USE MHDParams, ONLY: HallCoeffAlphaDefault

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

      HallCoeffAlpha = HallCoeffAlphaDefault


!------------------------------------------------------------------------------
      DIM = CoordinateSystemDimension()

      Force = 0.0d0
      StiffMatrix = 0.0d0
!------------------------------------------------------------------------------

      NodeIndexes => Element % NodeIndexes
 
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Element )

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

    INTEGER :: NX, ep, sgn, be, i, n, inode, pRow, gidV, gidVp, gidVm, gidI, p, gPhi
    INTEGER :: nSides, sideIdx, gidC, NumNodeConstraints, globalSideCount, globalNPhi
    INTEGER :: gaugeRow, gaugeNodeId
    LOGICAL :: UseNodalConstraints, gotIt
    TYPE(Element_t), POINTER :: Elem
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: maxN
    INTEGER :: gp

    TYPE(Nodes_t) :: EN
    TYPE(GaussIntegrationPoints_t) :: Integ
    REAL(dp) :: Basis(MAX_ELEMENT_NODES)
    REAL(dp) :: dBasisdx(MAX_ELEMENT_NODES,3)
    REAL(dp) :: SqrtElementMetric, s
    LOGICAL :: stat

    REAL(dp), ALLOCATABLE :: SideRowBuf(:)
    INTEGER, ALLOCATABLE :: SideConstraintCount(:), SideIsElectrodeRow(:,:), RowMembershipCount(:)
    
    ! Variables for boundary current coupling
    REAL(dp) :: AreaPlus, AreaMinus, areaCoeff
    INTEGER :: nCoupled

    ! Defensive no-op path: no electrode circuit DOFs requested.
    IF (NumElectrodePairs <= 0) THEN
      AuxMatrix => NULL()
      RETURN
    END IF

    maxN = Model % MaxElementNodes
    ALLOCATE( EN % x(maxN), EN % y(maxN), EN % z(maxN) )

    ! Keep legacy electrode circuit indices first:
    !   Vp(ep) = NPhi + 3*(ep-1) + 1
    !   Vm(ep) = NPhi + 3*(ep-1) + 2
    !   I(ep)  = NPhi + 3*(ep-1) + 3
    !
    ! Append nodal equipotential constraints after these 3*NumElectrodePairs rows.
    nSides = 2 * NumElectrodePairs
    UseNodalConstraints = .TRUE.
    NumNodeConstraints = 0
    IF (UseNodalConstraints) THEN
      ALLOCATE(SideConstraintCount(nSides))
      SideConstraintCount = 0
      ALLOCATE(SideIsElectrodeRow(nSides, NPhi))
      SideIsElectrodeRow = 0
      ALLOCATE(RowMembershipCount(NPhi))
      RowMembershipCount = 0
      ALLOCATE(SideRowBuf(NPhi))

      DO ep = 1, NumElectrodePairs
        DO sgn = -1, +1, 2
          SideRowBuf = 0.0_dp

          DO be = 1, Solver % Mesh % NumberOfBoundaryElements
            Elem => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements + be )
            NodeIndexes => Elem % NodeIndexes

            DO i = 1, Model % NumberOfBCs
              IF (Elem % BoundaryInfo % Constraint /= Model % BCs(i) % Tag) CYCLE
              IF (ElectrodePairOfBC(i) /= ep) CYCLE
              IF (ElectrodeSignOfBC(i) /= sgn) CYCLE

              n = Elem % TYPE % NumberOfNodes
              DO inode = 1, n
                pRow = PotentialPerm(NodeIndexes(inode))
              IF (pRow <= 0) CYCLE
              IF (ParEnv % PEs > 1 .AND. ALLOCATED(Solver % Matrix % RowOwner)) THEN
                IF (Solver % Matrix % RowOwner(pRow) /= ParEnv % MyPE) CYCLE
              END IF
              SideRowBuf(pRow) = 1.0_dp
              END DO
            END DO
          END DO

          IF (sgn == -1) THEN
            sideIdx = 2*(ep-1) + 1
          ELSE
            sideIdx = 2*(ep-1) + 2
          END IF

          DO pRow = 1, NPhi
            IF (SideRowBuf(pRow) > 0.5_dp) THEN
              NumNodeConstraints = NumNodeConstraints + 1
              SideConstraintCount(sideIdx) = SideConstraintCount(sideIdx) + 1
              SideIsElectrodeRow(sideIdx, pRow) = 1
            END IF
          END DO
        END DO
      END DO

      ! A potential DOF must belong to at most one electrode side.
      ! Multiple memberships indicate broken BC tagging/meshing, not recoverable here.
      DO pRow = 1, NPhi
        RowMembershipCount(pRow) = SUM(SideIsElectrodeRow(:, pRow))
        IF (RowMembershipCount(pRow) > 1) THEN
          CALL Fatal('BuildElectrodeAddMatrix', &
            'Potential DOF belongs to multiple electrode sides. Check BC tagging/mesh.')
        END IF
      END DO

      ! Sanity check (GLOBAL in MPI): each electrode side must contribute
      ! at least one constrained potential DOF across all ranks.
      DO sideIdx = 1, nSides
        IF (ParEnv % PEs > 1) THEN
          globalSideCount = NINT(ParallelReduction(REAL(SideConstraintCount(sideIdx), dp)))
        ELSE
          globalSideCount = SideConstraintCount(sideIdx)
        END IF

        IF (ParEnv % PEs > 1 .AND. SideConstraintCount(sideIdx) <= 0) THEN
          WRITE(*,'(A,I4,A,I4,A,I8)') ' [BuildElectrodeAddMatrix][rank ', ParEnv % MyPE, &
            '] local side ', sideIdx, ' has zero constrained DOFs (global=', globalSideCount
        END IF

        IF (ParEnv % MyPE == 0) THEN
          WRITE(*,'(A,I4,A,I8)') ' [BuildElectrodeAddMatrix] side ', sideIdx, &
            ' global constrained DOFs = ', globalSideCount
        END IF

        IF (globalSideCount <= 0) THEN
          CALL Fatal('BuildElectrodeAddMatrix', 'Electrode side has zero constrained potential DOFs')
        END IF
      END DO
    END IF


    IF (ParEnv % PEs > 1) THEN
      globalNPhi = NINT(ParallelReduction(REAL(NPhi, dp), 2))  ! MPI_MAX
    ELSE
      globalNPhi = NPhi
    END IF

    ! Append one nodal row per (sideIdx, global potential DOF).
    IF (UseNodalConstraints) THEN
      NX = 3 * NumElectrodePairs + nSides * globalNPhi
    ELSE
      NX = 3 * NumElectrodePairs
    END IF

    IF (ParEnv % PEs > 1) THEN
      WRITE(*,'(A,I4,A,I8,A,I8,A,I8)') ' [BuildElectrodeAddMatrix][rank ', ParEnv % MyPE, &
        '] active local nodal constraints=', NumNodeConstraints, ' local NPhi=', NPhi, &
        ' globalNPhi=', globalNPhi
      WRITE(*,'(A,I4,A,I8)') ' [BuildElectrodeAddMatrix][rank ', ParEnv % MyPE, &
        '] local NX=', NX
    END IF


    AuxMatrix => AllocateMatrix()
    AuxMatrix % FORMAT = MATRIX_LIST
    AuxMatrix % Symmetric = .FALSE.
    AuxMatrix % NumberOfRows = NPhi + NX

    ALLOCATE(AuxMatrix % RHS(AuxMatrix % NumberOfRows))
    AuxMatrix % RHS = 0.0_dp

    ! Initialize parallel info for AddMatrix following CircuitUtils pattern
    IF (ParEnv % PEs > 1) THEN
      ! Allocate RowOwner and initialize to -1 (like CircuitUtils.F90:1742)
      ALLOCATE(AuxMatrix % RowOwner(NPhi + NX))
      AuxMatrix % RowOwner = -1
      
      ! Set base circuit row ownership (Vp/Vm/I): owned by rank 0
      DO i = 1, NX
        IF (i <= 3*NumElectrodePairs) THEN
          AuxMatrix % RowOwner(NPhi + i) = 0
        ELSE
          AuxMatrix % RowOwner(NPhi + i) = 0
        END IF
      END DO

      ! Re-own nodal constraint rows by owner of each corresponding potential DOF.
      ! This avoids rank-local row index aliasing in MPI.
      IF (UseNodalConstraints .AND. ASSOCIATED(Solver % Matrix % ParallelInfo) .AND. &
          ASSOCIATED(Solver % Matrix % ParallelInfo % GlobalDOFs) .AND. &
          ALLOCATED(Solver % Matrix % RowOwner)) THEN
        DO ep = 1, NumElectrodePairs
          DO sgn = -1, +1, 2
            IF (sgn == -1) THEN
              sideIdx = 2*(ep-1) + 1
            ELSE
              sideIdx = 2*(ep-1) + 2
            END IF
            DO pRow = 1, NPhi
              IF (SideIsElectrodeRow(sideIdx, pRow) /= 1) CYCLE
              IF (Solver % Matrix % RowOwner(pRow) /= ParEnv % MyPE) CYCLE
              gPhi = Solver % Matrix % ParallelInfo % GlobalDOFs(pRow)
              IF (gPhi < 1 .OR. gPhi > globalNPhi) CYCLE
              gidC = NPhi + 3*NumElectrodePairs + (sideIdx-1)*globalNPhi + gPhi
              AuxMatrix % RowOwner(gidC) = ParEnv % MyPE
            END DO
          END DO
        END DO
      END IF
      
      ! Allocate ParallelInfo structure
      ALLOCATE(AuxMatrix % ParallelInfo)
      ALLOCATE(AuxMatrix % ParallelInfo % NeighbourList(NPhi + NX))
      
      ! Initialize phi rows (1:NPhi) - set to NULL since we use Solver%Matrix for phi
      DO i = 1, NPhi
        AuxMatrix % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
      END DO
      
      ! Initialize constraint DOF neighbor lists (NPhi+1:NPhi+NX)
      DO i = NPhi+1, NPhi+NX
        ALLOCATE(AuxMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1))
        IF (ALLOCATED(AuxMatrix % RowOwner) .AND. AuxMatrix % RowOwner(i) >= 0) THEN
          AuxMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1) = AuxMatrix % RowOwner(i)
        ELSE
          AuxMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1) = 0
        END IF
      END DO
    END IF

    !------------------------------------------------------------
    ! Exact nodal equipotential constraints for each electrode boundary DOF:
    !   phi(pRow) - Vside = 0
    !
    ! For each constraint row gidC:
    !   A(gidC, pRow) = +1
    !   A(gidC, gidV) = -1
    !
    ! And the transpose coupling for saddle-point symmetry:
    !   A(pRow, gidC) = +1
    !------------------------------------------------------------
    IF (UseNodalConstraints) THEN
      DO ep = 1, NumElectrodePairs
        DO sgn = -1, +1, 2

          IF (sgn == +1) THEN
            gidV = NPhi + 3*(ep-1) + 1   ! Vp (use local NPhi for local matrix indexing)
            sideIdx = 2*(ep-1) + 2
          ELSE
            gidV = NPhi + 3*(ep-1) + 2   ! Vm (use local NPhi for local matrix indexing)
            sideIdx = 2*(ep-1) + 1
          END IF

          DO pRow = 1, NPhi
            IF (ParEnv % PEs > 1 .AND. ASSOCIATED(Solver % Matrix % ParallelInfo) .AND. &
                ASSOCIATED(Solver % Matrix % ParallelInfo % GlobalDOFs)) THEN
              gPhi = Solver % Matrix % ParallelInfo % GlobalDOFs(pRow)
            ELSE
              gPhi = pRow
            END IF
            IF (gPhi < 1 .OR. gPhi > globalNPhi) CYCLE
            gidC = NPhi + 3*NumElectrodePairs + (sideIdx-1)*globalNPhi + gPhi

            IF (SideIsElectrodeRow(sideIdx, pRow) == 1) THEN
              IF (ParEnv % PEs > 1 .AND. ALLOCATED(Solver % Matrix % RowOwner)) THEN
                IF (Solver % Matrix % RowOwner(pRow) == ParEnv % MyPE) THEN
                  CALL AddToMatrixElement(AuxMatrix, pRow, gidC, 1.0_dp)
                  CALL AddToMatrixElement(AuxMatrix, gidC, pRow, 1.0_dp)
                  CALL AddToMatrixElement(AuxMatrix, gidC, gidV, -1.0_dp)
                END IF
              ELSE
                CALL AddToMatrixElement(AuxMatrix, pRow, gidC, 1.0_dp)
                CALL AddToMatrixElement(AuxMatrix, gidC, pRow, 1.0_dp)
                CALL AddToMatrixElement(AuxMatrix, gidC, gidV, -1.0_dp)
              END IF

              IF (ParEnv % MyPE == 0) THEN
                CALL AddToMatrixElement(AuxMatrix, gidV, gidC, -1.0_dp)
              END IF
            ELSE
              ! This reserved row is not part of the current electrode side;
              ! keep it harmless and non-singular.
              IF (ParEnv % MyPE == 0) THEN
                CALL AddToMatrixElement(AuxMatrix, gidC, gidC, 1.0_dp)
              END IF
            END IF
          END DO
        END DO
      END DO
    END IF

    !    Ohm law per pair: (Vp - Vm) - R * I = 0   (row gidI)
    !    Keep these rows on rank 0 (simple and stable).
    DO ep = 1, NumElectrodePairs
      gidVp = NPhi + 3*(ep-1) + 1  ! Use local NPhi for local matrix indexing
      gidVm = NPhi + 3*(ep-1) + 2  
      gidI  = NPhi + 3*(ep-1) + 3  ! Use local NPhi for local matrix indexing

      IF (ParEnv % MyPE == 0) THEN
        CALL AddToMatrixElement(AuxMatrix, gidI, gidVp,  1.0_dp)
        CALL AddToMatrixElement(AuxMatrix, gidI, gidVm, -1.0_dp)
        CALL AddToMatrixElement(AuxMatrix, gidI, gidI,  -ElectrodeResistance(ep))
      END IF
    END DO

    ! Pin potential=0 at one mesh node on the first cathode surface.
    ! One gauge point removes the 1D null space for any number of pairs.
    gaugeRow = 0
    IF (UseNodalConstraints) THEN
      DO pRow = 1, NPhi
        IF (SideIsElectrodeRow(1, pRow) == 1) THEN
          gaugeRow = pRow
          EXIT
        END IF
      END DO
    END IF
    IF (gaugeRow > 0) THEN
      CALL AddToMatrixElement(AuxMatrix, gaugeRow, gaugeRow, 1.0e8_dp)
      WRITE(*,'(A,I6)') ' [BuildElectrodeAddMatrix] Gauge: pinning Phi=0 at mesh DOF ', gaugeRow
    ELSE
      gidVm = NPhi + 2
      IF (ParEnv % MyPE == 0) THEN
        CALL AddToMatrixElement(AuxMatrix, gidVm, gidVm, 1.0e8_dp)
        WRITE(*,'(A,I6)') ' [BuildElectrodeAddMatrix] Gauge: fallback pinning Vm=0 at row ', gidVm
      END IF
    END IF

    ! Add zero diagonals for all constraint rows on non-owner ranks
    ! to keep row structures alive through LIST->CRS conversion.
    IF (ParEnv % MyPE /= 0) THEN
      DO p = NPhi + 1, NPhi + NX
        CALL AddToMatrixElement(AuxMatrix, p, p, 0.0_dp)
      END DO
    END IF
    
    IF (ParEnv % MyPE == 0) THEN
      WRITE(*,'(A)') ' [BuildElectrodeAddMatrix] Starting boundary current coupling'
      WRITE(*,'(A,I3)') '   NumElectrodePairs = ', NumElectrodePairs
      WRITE(*,'(A,I6)') '   NPhi = ', NPhi
      WRITE(*,'(A,I8)') '   Active nodal electrode constraints = ', NumNodeConstraints
    END IF
    
    nCoupled = 0
    DO ep = 1, NumElectrodePairs
      gidI = NPhi + 3*(ep-1) + 3
      
      IF (ParEnv % MyPE == 0) THEN
        WRITE(*,'(A,I3,A,I6)') '   Processing electrode pair ', ep, ', gidI = ', gidI
      END IF
      
      ! Compute electrode areas (same as earlier area computation)
      AreaPlus = 0.0_dp
      AreaMinus = 0.0_dp
      
      DO be = 1, Solver % Mesh % NumberOfBoundaryElements
        Elem => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements + be )
        NodeIndexes => Elem % NodeIndexes

        DO i = 1, Model % NumberOfBCs
          IF (Elem % BoundaryInfo % Constraint /= Model % BCs(i) % Tag) CYCLE
          IF (ElectrodePairOfBC(i) /= ep) CYCLE

          n = Elem % TYPE % NumberOfNodes
          EN % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes(1:n))
          EN % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes(1:n))
          EN % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes(1:n))

          Integ = GaussPoints(Elem)
          DO gp = 1, Integ % n
            stat = ElementInfo(Elem, EN, Integ % u(gp), Integ % v(gp), Integ % w(gp), &
                              SqrtElementMetric, Basis, dBasisdx)
            s = SqrtElementMetric * Integ % s(gp)
            IF (ElectrodeSignOfBC(i) == +1) THEN
              AreaPlus = AreaPlus + s
            ELSE
              AreaMinus = AreaMinus + s
            END IF
          END DO
        END DO
      END DO
      
      ! Parallel reduction for areas
      IF (ParEnv % PEs > 1) THEN
        AreaPlus = ParallelReduction(AreaPlus)
        AreaMinus = ParallelReduction(AreaMinus)
      END IF
      
      IF (ParEnv % MyPE == 0) THEN
        WRITE(*,'(A,I3,A,ES12.4,A,ES12.4)') ' [BoundaryCoupling] Pair ', ep, &
          ': AreaPlus=', AreaPlus, ' AreaMinus=', AreaMinus
        WRITE(*,'(A,I6)') '   gidI (current DOF index) = ', gidI
      END IF
      

      DO be = 1, Solver % Mesh % NumberOfBoundaryElements
        Elem => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements + be )
        NodeIndexes => Elem % NodeIndexes

        DO i = 1, Model % NumberOfBCs
          IF (Elem % BoundaryInfo % Constraint /= Model % BCs(i) % Tag) CYCLE
          IF (ElectrodePairOfBC(i) /= ep) CYCLE
          
          n = Elem % TYPE % NumberOfNodes
          EN % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes(1:n))
          EN % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes(1:n))
          EN % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes(1:n))

          ! Determine area coefficient
          IF (ElectrodeSignOfBC(i) == +1 .AND. AreaPlus > 1.0e-20) THEN
            areaCoeff = 1.0_dp / AreaPlus  ! Positive: I/A
          ELSE IF (ElectrodeSignOfBC(i) == -1 .AND. AreaMinus > 1.0e-20) THEN
            areaCoeff = -1.0_dp / AreaMinus  ! Negative: -I/A
          ELSE
            CYCLE
          END IF

          Integ = GaussPoints(Elem)
          DO gp = 1, Integ % n
            stat = ElementInfo(Elem, EN, Integ % u(gp), Integ % v(gp), Integ % w(gp), &
                              SqrtElementMetric, Basis, dBasisdx)
            s = SqrtElementMetric * Integ % s(gp)
            
            DO inode = 1, n
              pRow = PotentialPerm(NodeIndexes(inode))
              IF (pRow <= 0) CYCLE
              
              ! FORWARD coupling: phi equation gets I contribution
              ! Weak form: K·φ = ∫(I/A)ψ dS  -->  K·φ - B·I = 0
              ! So: A(pRow, gidI) = -∫ψ/A dS  (NEGATIVE!)
              IF (ParEnv % PEs > 1 .AND. ALLOCATED(Solver % Matrix % RowOwner)) THEN
                IF (Solver % Matrix % RowOwner(pRow) == ParEnv % MyPE) THEN
                  CALL AddToMatrixElement(AuxMatrix, pRow, gidI, -s * Basis(inode) * areaCoeff)
                  IF (nCoupled < 5 .AND. ParEnv % MyPE == 0) THEN
                    WRITE(*,'(A,I3,A,I6,A,I6,A,ES12.4,A,ES12.4)') &
                      '   [Forward] EP=', ep, ' A(', pRow, ',', gidI, ')=', &
                      -s * Basis(inode) * areaCoeff, ' (areaCoeff=', areaCoeff, ')'
                  END IF
                  nCoupled = nCoupled + 1
                END IF
              ELSE
                CALL AddToMatrixElement(AuxMatrix, pRow, gidI, -s * Basis(inode) * areaCoeff)
                IF (nCoupled < 5) THEN
                  WRITE(*,'(A,I3,A,I6,A,I6,A,ES12.4,A,ES12.4)') &
                    '   [Forward] EP=', ep, ' A(', pRow, ',', gidI, ')=', &
                    -s * Basis(inode) * areaCoeff, ' (areaCoeff=', areaCoeff, ')'
                END IF
                nCoupled = nCoupled + 1
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
    
    IF (ParEnv % MyPE == 0) THEN
      WRITE(*,'(A,I6)') ' [BuildElectrodeAddMatrix] Added boundary current coupling entries: ', nCoupled
    END IF

    ! Convert AddMatrix to CRS like Circuits does
    CALL List_ToCRSMatrix(AuxMatrix)

    ! Debug output
    IF (ParEnv % MyPE == 0) THEN
      WRITE(*,'(A,I6,A,I8)') ' [BuildElectrodeAddMatrix] AuxMatrix % NumberOfRows = ', &
          AuxMatrix % NumberOfRows, '  (NPhi + NX, NX=', NX
      IF (ALLOCATED(AuxMatrix % RowOwner)) THEN
        WRITE(*,'(A)') ' [BuildElectrodeAddMatrix] AuxMatrix % RowOwner allocated successfully'
      ELSE
        WRITE(*,'(A)') ' [BuildElectrodeAddMatrix] WARNING: AuxMatrix % RowOwner NOT allocated'
      END IF
      
      ! Log some matrix statistics
      WRITE(*,'(A,I6)') ' [BuildElectrodeAddMatrix] Total coupling entries added: ', nCoupled
      WRITE(*,'(A)') ' [BuildElectrodeAddMatrix] Matrix structure summary:'
      WRITE(*,'(A,I8)') '   Total non-zeros in AddMatrix: ', &
        COUNT(AuxMatrix % Values /= 0.0_dp)
      
      ! Check a sample coupling entry by scanning constraint row entries
      DO ep = 1, MIN(1, NumElectrodePairs)
        gidI = NPhi + 3*(ep-1) + 3
        WRITE(*,'(A,I3,A,I6,A)') '   Checking constraint row gidI=', gidI, &
          ' (index ', gidI, ') for non-zero entries...'
        ! Note: After CRS conversion, we can't easily inspect individual entries
        ! without walking the CRS structure, so this is just informational
      END DO      
    END IF

    IF (UseNodalConstraints) THEN
      DEALLOCATE(RowMembershipCount)
      DEALLOCATE(SideIsElectrodeRow)
      DEALLOCATE(SideConstraintCount)
      DEALLOCATE(SideRowBuf)
    END IF
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

