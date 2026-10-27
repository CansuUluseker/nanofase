!> Module containing the LakeSegment and ReservoirSegment types.
!!
!! Lakes and reservoirs are represented as *lake segments*, contained within grid
!! cells and thus with a maximum surface area equal to the grid cell's area. Lakes
!! or reservoirs spanning multiple grid cells are split into multiple lake segments,
!! connected to each other by the normal reach routing (advection transfers water,
!! SPM and contaminants to/from neighbouring segments). Each lake segment has a bed
!! sediment.
!!
!! Processes transferring material to/from each segment:
!!  - Advection to/from neighbouring segments and reaches (inflow/outflow pointers)
!!  - Surface runoff (water) and the soil erosion it causes (sediment, contaminants)
!!  - Point sources and atmospheric deposition (via the diffuse source machinery)
!!  - Sediment deposition and resuspension to/from the bed sediment
!!  - Settling from the epilimnion to the hypolimnion when thermally stratified
!!  - Mixing of water and associated SPM/contaminants between the two layers
!! Contaminant transformations (heteroaggregation to SPM, dissolution, chemical
!! alteration, ...) are handled by the segment's Reactor.
!!
!! Thermal stratification into an epilimnion and hypolimnion is diagnosed each
!! timestep from meteorological data (wind speed and solar radiation, when
!! provided; a water-temperature proxy otherwise) using a simple empirical scheme
!! whose parameters live in the optional &lake group of the constants file.
!! Rather than duplicating the mass stores per layer, the canonical stores
!! (me%m_spm, me%m_contaminant) always hold whole-segment totals - so the reactor,
!! outputs and checkpointing see the segment exactly like any reach - and the
!! epilimnion's share of each store is tracked as a fraction (f_epi_*). When the
!! segment is unstratified the shares are 1.
module LakeSegmentModule
    use GlobalsModule
    use ReachModule
    use UtilModule
    use ResultModule
    use BedSedimentModule
    use LoggerModule, only: LOGR
    use DataInputModule, only: DATASET
    use ReactorModule
    use BiotaWaterModule
    use ContaminantModule
    use PFASConstantsModule, only: PFAS_AQ, PFAS_SPM, PFAS_SOL
    implicit none

    ! Types for reservoir segments
    integer, parameter, public :: SPILLWAY_UNCONTROLLED = 0
    integer, parameter, public :: SPILLWAY_CONTROLLED = 1
    integer, parameter, public :: SPILLWAY_COMBINATION = 2

    ! Constants for the wind resuspension parameterisation
    real(dp), parameter :: rho_air = 1.225_dp       ! Density of air [kg/m3]
    real(dp), parameter :: C_drag_wind = 1.3e-3_dp  ! Wind drag coefficient at 10 m [-]

    !> A LakeSegment is the part of a lake contained within one grid cell
    type, public, extends(Reach) :: LakeSegment
        real(dp) :: maxDepth = 0.0_dp           !! Depth at the outlet/spillway level [m]
        real(dp) :: areaFraction = 1.0_dp       !! Fraction of the grid cell covered by this segment [-]
        integer  :: lakeId = 0                  !! Identifier of the lake this segment belongs to (segments of a multi-cell lake share an id; 0 = unassigned)
        ! Thermal stratification state
        logical  :: isStratified = .false.      !! Is the segment currently thermally stratified?
        real(dp) :: epilimnionDepth = 0.0_dp    !! Depth of the epilimnion (whole water column when unstratified) [m]
        ! Epilimnion share of the canonical (whole-segment) mass stores
        real(dp), allocatable :: f_epi_spm(:)   !! Epilimnion share of each SPM size class [-]
        real(dp), allocatable :: f_epi_cont(:,:,:) !! Epilimnion share of each contaminant (species, form, phase) element [-]
        real(dp) :: f_epi_diss = 1.0_dp         !! Epilimnion share of the dissolved contaminant store [-]
      contains
        ! Create
        procedure :: create => createLakeSegment
        procedure :: allocateAndInitialise => allocateAndInitialiseLakeSegment
        ! Simulators
        procedure :: update => updateLakeSegment
        procedure :: updateStratification => updateStratificationLakeSegment
        procedure :: setWindResuspensionRate => setWindResuspensionRateLakeSegment
        procedure :: calculateOutflowVolume => calculateOutflowVolumeLakeSegment
        ! Data handlers
        procedure :: parseInputData => parseInputDataLakeSegment
        ! Helpers
        procedure :: initialiseLakeState
        procedure :: setLakeDimensions
        procedure :: windSpeedAt
        procedure :: finalise => finaliseLakeSegment
    end type

    !> A ReservoirSegment is a LakeSegment whose outflow is controlled by a dam
    !! and spillway
    type, public, extends(LakeSegment) :: ReservoirSegment
        integer  :: spillwayType = SPILLWAY_UNCONTROLLED    !! Uncontrolled, controlled or a combination
        real(dp) :: releaseRate = 0.0_dp                    !! Controlled spillway release rate [m3/s]
      contains
        procedure :: create => createReservoirSegment
        procedure :: calculateOutflowVolume => calculateOutflowVolumeReservoirSegment
    end type

  contains

    !> Create this LakeSegment with the provided grid cell and waterbody indices
    !! (x, y, w) and sediment size class distribution
    function createLakeSegment(me, x, y, w, distributionSediment) result(rslt)
        class(LakeSegment), intent(inout) :: me
        integer, intent(in) :: x, y, w
        real(dp), intent(in) :: distributionSediment(C%nSizeClassesSpm)
        type(Result) :: rslt
        integer :: i, s
        real(dp) :: T0, rho_s

        ! Create base waterbody and set ref (ReservoirSegment overrides the ref after this)
        call rslt%addErrors(.errors. me%WaterBody%create(x, y, w, distributionSediment))
        me%ref = trim(ref("LakeSegment", x, y, w))

        ! Parse grid-based inputs (geometry, resuspension/sediment transport params, topology)
        call rslt%addErrors(.errors. me%parseInputData())

        ! Create the contaminant object for WATER (sizes, rates, etc.)
        call rslt%addErrors(.errors. me%m_contaminant%create_from_data( &
            compartment='water', &
            contaminantDensity=DATASET%contaminantDensity, &
            soilAttachmentEfficiency=DATASET%soilAttachmentEfficiencyConstant, &
            riverAttachmentEfficiency=DATASET%riverAttachmentEfficiency, &
            estuaryAttachmentEfficiency=DATASET%estuaryAttachmentEfficiency, &
            k_diss_pristine=DATASET%contaminant_k_diss_pristine, &
            k_diss_transformed=DATASET%contaminant_k_diss_transformed, &
            k_transform_pristine=DATASET%contaminant_k_transform_pristine, &
            waterTemperature=DATASET%waterTemperature(C%startDate%yearday()), x=me%x, y=me%y ))

        ! The initial water state is applied in WaterBody%finaliseCreate

        ! Settling velocities for SPM size classes
        if (allocated(me%W_settle_spm)) deallocate(me%W_settle_spm)
        allocate(me%W_settle_spm(C%nSizeClassesSpm))
        T0 = DATASET%waterTemperature(C%startDate%yearday())
        do s = 1, C%nSizeClassesSpm
            if (allocated(DATASET%spmDensityBySizeClass)) then
                rho_s = DATASET%spmDensityBySizeClass(min(s, size(DATASET%spmDensityBySizeClass)))
            else
                rho_s = C%sedimentParticleDensities(min(s, size(C%sedimentParticleDensities)))
            end if
            me%W_settle_spm(s) = me%calculateSettlingVelocity( &
                C%d_spm(s), rho_s, T0, &
                DATASET%depositionAlphaConstant, DATASET%depositionBetaConstant)
        end do

        ! Ensure an SPM vector exists for the reactor
        if (allocated(me%C_spm)) deallocate(me%C_spm)
        allocate(me%C_spm(C%nSizeClassesSpm))
        me%C_spm = 0.0_dp

        ! Lake state (dimensions, stratification shares); starts full to the outlet level
        call me%initialiseLakeState()

        ! Create the bed sediment and reactor. Standing water: zero velocity for the reactor.
        allocate(BedSediment :: me%bedSediment)
        allocate(Reactor     :: me%reactor)
        call rslt%addErrors([ &
            .errors. me%bedSediment%create(me%x, me%y, me%w), &
            .errors. me%reactor%create( &
                me%x, me%y, 'water', &
                me%m_contaminant, me%volume, &
                DATASET%waterTemperature(C%startDate%yearday()), &
                C_spm=me%C_spm, W_settle_spm=me%W_settle_spm, &
                G=DATASET%shearRate, velocity=0.0_dp) ])

        ! Water biota
        allocate(me%biotaIndices(0))
        if (DATASET%hasBiota) then
            do i = 1, DATASET%nBiota
                if (trim(DATASET%biotaCompartment(i)) == 'water') then
                    me%nBiota = me%nBiota + 1
                    me%biotaIndices = [me%biotaIndices, i]
                end if
            end do
        end if
        allocate(me%biota(me%nBiota))
        do i = 1, me%nBiota
            call rslt%addErrors(.errors. me%biota(i)%create(me%biotaIndices(i)))
        end do

        call rslt%addToTrace('Creating ' // trim(me%ref))
        if (me%lakeId > 0) then
            call LOGR%toFile("Creating " // trim(me%ref) // " (lake id " &
                // trim(str(me%lakeId)) // "): success")
        else
            call LOGR%toFile("Creating " // trim(me%ref) // ": success")
        end if
    end function

    !> Create this ReservoirSegment: a LakeSegment plus the dam/spillway parameters
    function createReservoirSegment(me, x, y, w, distributionSediment) result(rslt)
        class(ReservoirSegment), intent(inout) :: me
        integer, intent(in) :: x, y, w
        real(dp), intent(in) :: distributionSediment(C%nSizeClassesSpm)
        type(Result) :: rslt

        rslt = me%LakeSegment%create(x, y, w, distributionSediment)
        me%ref = trim(ref("ReservoirSegment", x, y, w))

        ! Spillway operating parameters, from data if provided, else defaults
        me%spillwayType = SPILLWAY_UNCONTROLLED
        if (allocated(DATASET%reservoirSpillwayType)) then
            if (me%x <= size(DATASET%reservoirSpillwayType, 1) &
                .and. me%y <= size(DATASET%reservoirSpillwayType, 2)) then
                me%spillwayType = DATASET%reservoirSpillwayType(me%x, me%y)
            end if
        end if
        me%releaseRate = DATASET%reservoirDefaultReleaseRate
        if (allocated(DATASET%reservoirReleaseRate)) then
            if (me%x <= size(DATASET%reservoirReleaseRate, 1) &
                .and. me%y <= size(DATASET%reservoirReleaseRate, 2)) then
                me%releaseRate = real(DATASET%reservoirReleaseRate(me%x, me%y), dp)
            end if
        end if
    end function

    !> Allocate memory and set initial values. Called from the base class's
    !! finaliseCreate, which zeroes the volume - so the lake state (full to the
    !! outlet level) is re-imposed afterwards.
    subroutine allocateAndInitialiseLakeSegment(me)
        class(LakeSegment), intent(inout) :: me
        type(Result) :: r
        ! Mirror Reach%allocateAndInitialise (the parent binding can't be called
        ! directly as Reach is abstract): base waterbody allocation, roughness
        ! default and the WATER contaminant state
        call me%WaterBody%allocateAndInitialise()
        me%n = C%n_river
        call r%addErrors(.errors. me%m_contaminant%create_from_data('water', DATASET%contaminantDensity, &
            DATASET%soilAttachmentEfficiencyConstant, DATASET%riverAttachmentEfficiency, &
            DATASET%estuaryAttachmentEfficiency, DATASET%contaminant_k_diss_pristine, &
            DATASET%contaminant_k_diss_transformed, DATASET%contaminant_k_transform_pristine, &
            DATASET%waterTemperature(C%startDate%yearday()), x=me%x, y=me%y))
        ! Re-impose the lake state (WaterBody%allocateAndInitialise zeroes the volume)
        call me%initialiseLakeState()
    end subroutine

    !> Set the lake segment's geometry and layer state to its initial (full,
    !! unstratified) condition
    subroutine initialiseLakeState(me)
        class(LakeSegment), intent(inout) :: me
        ! Start full to the outlet/spillway level, so through-flow is immediate
        ! rather than spending years filling an empty basin
        me%depth = me%maxDepth
        call me%setLakeDimensions()
        me%isStratified = .false.
        me%epilimnionDepth = me%depth
        ! Everything in the (single) surface layer
        if (allocated(me%f_epi_spm)) deallocate(me%f_epi_spm)
        allocate(me%f_epi_spm(C%nSizeClassesSpm))
        me%f_epi_spm = 1.0_dp
        if (allocated(me%f_epi_cont)) deallocate(me%f_epi_cont)
        allocate(me%f_epi_cont(C%contaminantDim(1), C%contaminantDim(2), &
            max(PFAS_NPHASES, C%contaminantDim(3))))
        me%f_epi_cont = 1.0_dp
        me%f_epi_diss = 1.0_dp
    end subroutine

    !> Set the segment's dimensions from its current depth. The surface area is
    !! fixed (vertical-sided basin covering areaFraction of the cell); depth and
    !! volume vary with storage.
    subroutine setLakeDimensions(me)
        class(LakeSegment), intent(inout) :: me
        me%surfaceArea = me%areaFraction * DATASET%gridRes(1) * DATASET%gridRes(2)
        me%bedArea = me%surfaceArea
        me%volume = me%depth * me%surfaceArea
        ! Effective plan dimensions, used for erosion delivery scaling and outputs
        me%length = sqrt(me%surfaceArea)
        me%width = sqrt(me%surfaceArea)
        me%xsArea = me%width * me%depth
        me%velocity = 0.0_dp
    end subroutine

    !> Parse data from the input file for this lake segment
    function parseInputDataLakeSegment(me) result(rslt)
        class(LakeSegment), intent(inout) :: me
        type(Result) :: rslt
        integer :: nx, ny
        logical :: ok

        me%f_m = 1.0_dp                             ! No meandering for a basin
        me%alpha_hetero = DATASET%riverAttachmentEfficiency
        me%n = C%n_river

        ! Which lake does this segment belong to? Segments of a multi-cell lake
        ! share an id (0 if no 'lake_id' input variable was provided)
        if (allocated(DATASET%lakeId)) then
            if (me%x <= size(DATASET%lakeId, 1) .and. me%y <= size(DATASET%lakeId, 2)) then
                me%lakeId = DATASET%lakeId(me%x, me%y)
            end if
        end if

        ! Depth at the outlet/spillway level, from data if provided
        me%maxDepth = DATASET%lakeDefaultDepth
        if (allocated(DATASET%lakeDepth)) then
            if (me%x <= size(DATASET%lakeDepth, 1) .and. me%y <= size(DATASET%lakeDepth, 2)) then
                if (DATASET%lakeDepth(me%x, me%y) > 0.0) then
                    me%maxDepth = real(DATASET%lakeDepth(me%x, me%y), dp)
                end if
            end if
        end if

        ! Fraction of the cell covered. A lake segment cannot be larger than its cell.
        me%areaFraction = min(max(DATASET%lakeDefaultAreaFraction, 0.0_dp), 1.0_dp)
        if (allocated(DATASET%lakeAreaFraction)) then
            if (me%x <= size(DATASET%lakeAreaFraction, 1) .and. me%y <= size(DATASET%lakeAreaFraction, 2)) then
                if (DATASET%lakeAreaFraction(me%x, me%y) > 0.0) then
                    me%areaFraction = min(real(DATASET%lakeAreaFraction(me%x, me%y), dp), 1.0_dp)
                end if
            end if
        end if

        ! Resuspension calibration params with grid -> scalar fallback (river values)
        ok = .false.
        if (allocated(DATASET%resuspensionAlpha)) then
            nx = size(DATASET%resuspensionAlpha, 1); ny = size(DATASET%resuspensionAlpha, 2)
            if (nx > 0 .and. ny > 0 .and. me%x <= nx .and. me%y <= ny) then
                me%alpha_resus = DATASET%resuspensionAlpha(me%x, me%y)
                ok = .true.
            end if
        end if
        if (.not. ok) me%alpha_resus = DATASET%waterResuspensionAlpha
        ok = .false.
        if (allocated(DATASET%resuspensionBeta)) then
            nx = size(DATASET%resuspensionBeta, 1); ny = size(DATASET%resuspensionBeta, 2)
            if (nx > 0 .and. ny > 0 .and. me%x <= nx .and. me%y <= ny) then
                me%beta_resus = DATASET%resuspensionBeta(me%x, me%y)
                ok = .true.
            end if
        end if
        if (.not. ok) me%beta_resus = DATASET%waterResuspensionBeta

        ! Sediment transport capacity params with grid -> constant fallback
        ok = .false.
        if (allocated(DATASET%sedimentTransport_a)) then
            nx = size(DATASET%sedimentTransport_a, 1); ny = size(DATASET%sedimentTransport_a, 2)
            if (nx > 0 .and. ny > 0 .and. me%x <= nx .and. me%y <= ny) then
                me%a_stc = DATASET%sedimentTransport_a(me%x, me%y)
                ok = .true.
            end if
        end if
        if (.not. ok) me%a_stc = DATASET%sedimentTransport_aConstant
        ok = .false.
        if (allocated(DATASET%sedimentTransport_b)) then
            nx = size(DATASET%sedimentTransport_b, 1); ny = size(DATASET%sedimentTransport_b, 2)
            if (nx > 0 .and. ny > 0 .and. me%x <= nx .and. me%y <= ny) then
                me%b_stc = DATASET%sedimentTransport_b(me%x, me%y)
                ok = .true.
            end if
        end if
        if (.not. ok) me%b_stc = DATASET%sedimentTransport_bConstant
        ok = .false.
        if (allocated(DATASET%sedimentTransport_c)) then
            nx = size(DATASET%sedimentTransport_c, 1); ny = size(DATASET%sedimentTransport_c, 2)
            if (nx > 0 .and. ny > 0 .and. me%x <= nx .and. me%y <= ny) then
                me%c_stc = DATASET%sedimentTransport_c(me%x, me%y)
                ok = .true.
            end if
        end if
        if (.not. ok) me%c_stc = DATASET%sedimentTransport_cConstant

        ! Water temperature (vector over day-of-year)
        me%T_water = DATASET%waterTemperature

        ! Inflow/outflow topology, exactly like a reach - this is what connects
        ! neighbouring lake segments (and up/downstream reaches) by advection
        call rslt%addErrors(.errors. me%parseInflowsAndOutflow())
        ! Standing water: slope is meaningless but keep it defined for safety
        me%slope = C%minStreamSlope

        call rslt%addToTrace('Parsing input data (LakeSegment)')
    end function

    !> The wind speed forcing for this segment on timestep t [m/s]: the gridded
    !! time series if provided, else the scalar constant
    function windSpeedAt(me, t) result(u10)
        class(LakeSegment), intent(in) :: me
        integer, intent(in) :: t
        real(dp) :: u10
        if (allocated(DATASET%windSpeedTimeSeries)) then
            u10 = real(DATASET%windSpeedTimeSeries(me%x, me%y, min(t, size(DATASET%windSpeedTimeSeries, 3))), dp)
        else
            u10 = DATASET%windSpeed
        end if
        u10 = max(0.0_dp, u10)
    end function

    !> Diagnose thermal stratification for this timestep from meteorological data.
    !! Stratification forms when there is enough surface heating (solar radiation
    !! above a threshold, or a water temperature proxy if no radiation data) and
    !! not too much wind, and the basin is deep enough to hold two layers. The
    !! epilimnion deepens with wind speed. On the transition mixed -> stratified,
    !! each store's epilimnion share is initialised to the epilimnion volume
    !! fraction; on stratified -> mixed (turnover), the layers are fully merged.
    subroutine updateStratificationLakeSegment(me, t, T_water_t)
        class(LakeSegment), intent(inout) :: me
        integer, intent(in) :: t                    !! Current timestep
        real(dp), intent(in) :: T_water_t           !! Water temperature today [deg C]
        real(dp) :: u10, z_epi_pot, f0
        logical :: heating, stratifiedNow

        u10 = me%windSpeedAt(t)

        ! Is there enough surface heating to sustain a thermocline?
        if (allocated(DATASET%solarRadiation)) then
            heating = DATASET%solarRadiation(me%x, me%y, min(t, size(DATASET%solarRadiation, 3))) &
                > DATASET%lakeStratificationRadiationThreshold
        else
            heating = T_water_t > DATASET%lakeStratificationTemperatureThreshold
        end if

        ! Potential (wind-mixed) epilimnion depth
        z_epi_pot = max(DATASET%lakeEpilimnionMinDepth, &
            DATASET%lakeEpilimnionBaseDepth + DATASET%lakeEpilimnionWindCoefficient * u10)

        stratifiedNow = heating &
            .and. u10 < DATASET%lakeStratificationWindThreshold &
            .and. z_epi_pot < 0.75_dp * me%depth

        if (stratifiedNow) then
            me%epilimnionDepth = min(z_epi_pot, 0.75_dp * me%depth)
            if (.not. me%isStratified) then
                ! Onset of stratification: split all stores by volume fraction
                f0 = min(max(divideCheckZero(me%epilimnionDepth, me%depth), 0.0_dp), 1.0_dp)
                me%f_epi_spm = f0
                me%f_epi_cont = f0
                me%f_epi_diss = f0
            end if
        else
            if (me%isStratified) then
                ! Turnover: the whole column mixes into a single layer
                me%f_epi_spm = 1.0_dp
                me%f_epi_cont = 1.0_dp
                me%f_epi_diss = 1.0_dp
            end if
            me%epilimnionDepth = me%depth
        end if
        me%isStratified = stratifiedNow
    end subroutine

    !> Set the resuspension rate from wind-driven turbulence: the wind power per
    !! unit surface area (rho_air * C_d * u10^3) is attenuated exponentially with
    !! depth before reaching the bed, and scaled by the (calibratable) beta_resus.
    !! When the segment is stratified, the thermocline shields the bed from
    !! wind-driven turbulence and resuspension is switched off.
    subroutine setWindResuspensionRateLakeSegment(me, t)
        class(LakeSegment), intent(inout) :: me
        integer, intent(in) :: t
        real(dp) :: u10, omega_bed

        if (me%isStratified .or. me%depth <= C%epsilon) then
            me%k_resus = 0.0_dp
            return
        end if
        u10 = me%windSpeedAt(t)
        ! Wind power per unit area reaching the bed [W/m2]
        omega_bed = rho_air * C_drag_wind * u10**3 &
            * exp(-me%depth / max(C%epsilon, DATASET%lakeWindResuspensionAttenuationDepth))
        ! beta_resus [s2/kg] is reused as the calibration factor, as for reaches
        me%k_resus = me%beta_resus * omega_bed
    end subroutine

    !> The volume of water leaving this lake segment this timestep, given the
    !! available storage V_avail [m3]. A lake behaves as an uncontrolled spillway
    !! at its natural outlet: any water above the outlet level leaves immediately,
    !! keeping the lake at (or below) maxDepth.
    function calculateOutflowVolumeLakeSegment(me, V_avail, dt) result(V_out)
        class(LakeSegment), intent(in) :: me
        real(dp), intent(in) :: V_avail             !! Storage + inflows this timestep [m3]
        real(dp), intent(in) :: dt                  !! Timestep length [s]
        real(dp) :: V_out                           !! Outflow volume [m3/timestep]
        real(dp) :: V_cap
        V_cap = me%surfaceArea * me%maxDepth
        V_out = max(0.0_dp, V_avail - V_cap)
    end function

    !> The outflow volume for a reservoir segment, controlled by the dam and
    !! spillway. See the module description for the three spillway types.
    function calculateOutflowVolumeReservoirSegment(me, V_avail, dt) result(V_out)
        class(ReservoirSegment), intent(in) :: me
        real(dp), intent(in) :: V_avail             !! Storage + inflows this timestep [m3]
        real(dp), intent(in) :: dt                  !! Timestep length [s]
        real(dp) :: V_out                           !! Outflow volume [m3/timestep]
        real(dp) :: V_cap, V_release
        V_cap = me%surfaceArea * me%maxDepth
        V_release = max(0.0_dp, me%releaseRate) * dt
        select case (me%spillwayType)
            case (SPILLWAY_CONTROLLED)
                ! Regulated outflow only; the reservoir may surcharge above the
                ! spillway level if inflows exceed the operating release
                V_out = min(V_avail, V_release)
            case (SPILLWAY_COMBINATION)
                ! Controlled release plus uncontrolled spill above the spillway
                V_out = min(V_avail, V_release)
                V_out = V_out + max(0.0_dp, V_avail - V_out - V_cap)
            case default     ! SPILLWAY_UNCONTROLLED
                V_out = max(0.0_dp, V_avail - V_cap)
        end select
    end function

    !> Run the lake segment simulation for this timestep
    subroutine updateLakeSegment(me, t, q_runoff, q_overland, j_spm_runoff, j_contaminant_runoff, &
                                 contributingArea, isWarmUp)
        class(LakeSegment), intent(inout) :: me                    !! This LakeSegment instance
        integer, intent(in) :: t                                   !! The current timestep
        real(dp), intent(in) :: q_runoff                           !! Runoff from the hydrological model [m3/m2/timestep]
        real(dp), intent(in) :: q_overland                         !! Overland runoff [m3/m2/timestep]
        real(dp), intent(in) :: j_spm_runoff(:)                    !! Eroded sediment runoff to this segment [kg/timestep]
        type(Contaminant), intent(in) :: j_contaminant_runoff      !! Contaminant runoff to this segment [kg/timestep]
        real(dp), intent(in) :: contributingArea                   !! Area contributing to this segment [m2]
        logical, intent(in) :: isWarmUp                            !! Are we in a warm up period?
        type(Result) :: rslt
        type(datetime) :: currentDate
        real(dp) :: T_water_t, dt
        integer :: i
        real(dp) :: V_avail, V_out, V_epi, V_hyp, z_epi, z_hyp, k_out
        real(dp) :: frac_settle_cont, frac_dep_cont, fd, fu, V_bot
        real(dp) :: m_spm_old(C%nSizeClassesSpm), bot_share_spm(C%nSizeClassesSpm)
        real(dp) :: dj_spm_in(C%nSizeClassesSpm), dj_spm_deposit(C%nSizeClassesSpm)
        real(dp) :: dj_spm_out(C%nSizeClassesSpm), frac_settle(C%nSizeClassesSpm)
        real(dp) :: frac_dep(C%nSizeClassesSpm)
        real(dp) :: dj_spm_resus(C%nSizeClassesSpm)
        real(dp) :: dj_spm_resus_perArea(C%nSizeClassesSpm), dj_spm_resus_perArea_(C%nSizeClassesSpm)
        real(dp) :: dj_spm_deposit_perArea(C%nSizeClassesSpm)
        type(Contaminant) :: dj_cont_in, dj_cont_dep, dj_cont_out
        type(Contaminant) :: cont_resus_spm, j_contam_dep_perArea, m_bottom, dj_exchange
        type(Contaminant) :: c_env_contaminant
        real(dp), allocatable :: c_old(:,:,:), bot_share_cont(:,:,:)
        real(dp) :: diss_old
        integer, parameter :: sp = kind(1.0)

        dt = real(C%timeStep, dp)

        ! Reset all flows to zero and hard-reset contaminant flux holders
        call me%emptyFlows()
        call rslt%addErrors(.errors. me%j_contaminant_inflow%create())
        call rslt%addErrors(.errors. me%j_contaminant_runoff%create())
        call rslt%addErrors(.errors. me%j_contaminant_transfers%create())
        call rslt%addErrors(.errors. me%j_contaminant_deposition%create())
        call rslt%addErrors(.errors. me%j_contaminant_resuspension%create())
        call rslt%addErrors(.errors. me%j_contaminant_outflow%create())
        call rslt%addErrors(.errors. dj_cont_in%create())
        call rslt%addErrors(.errors. dj_cont_dep%create())
        call rslt%addErrors(.errors. dj_cont_out%create())
        call rslt%addErrors(.errors. cont_resus_spm%create())
        call rslt%addErrors(.errors. j_contam_dep_perArea%create())

        currentDate = C%startDate + timedelta(t-1)
        T_water_t = me%T_water(currentDate%yearday())

        ! --------------------------------------------------------------------
        ! 1. Advection in from neighbouring lake segments and upstream reaches
        ! --------------------------------------------------------------------
        do i = 1, me%nInflows
            me%Q%inflow     = me%Q%inflow     - me%inflows(i)%item%Q%outflow
            me%j_spm%inflow = me%j_spm%inflow - me%inflows(i)%item%j_spm%outflow
            call me%j_contaminant_inflow%add(me%inflows(i)%item%get_j_contaminant_outflow())
        end do

        ! Surface runoff provides water to the segment
        me%Q%runoff = q_runoff * contributingArea
        me%Q_in_total = me%Q%inflow + me%Q%runoff

        ! --------------------------------------------------------------------
        ! 2. Soil erosion delivers sediment and contaminants to the segment,
        !    limited by the sediment transport capacity of the overland flow.
        ! --------------------------------------------------------------------
        call me%setErosionYields(j_spm_runoff, q_overland, contributingArea, j_contaminant_runoff)
        me%j_spm%bankErosion = 0.0_dp

        ! --------------------------------------------------------------------
        ! 3. Point sources and atmospheric deposition (diffuse source 2)
        ! --------------------------------------------------------------------
        if (.not. C%ignoreContaminant .and. .not. isWarmUp) then
            call me%updateSources(t)
        end if

        ! --------------------------------------------------------------------
        ! 4. Water balance: level-pool storage with outlet/spillway outflow
        ! --------------------------------------------------------------------
        V_avail = max(0.0_dp, me%volume + me%Q_in_total)
        V_out = me%calculateOutflowVolume(V_avail, dt)
        V_out = min(max(V_out, 0.0_dp), V_avail)
        me%depth = divideCheckZero(V_avail - V_out, me%surfaceArea)
        call me%setLakeDimensions()

        ! --------------------------------------------------------------------
        ! 5. Thermal stratification (epilimnion/hypolimnion), from wind speed
        !    and solar radiation
        ! --------------------------------------------------------------------
        call me%updateStratification(t, T_water_t)
        z_epi = min(me%epilimnionDepth, me%depth)
        z_hyp = max(me%depth - z_epi, 0.0_dp)
        V_epi = z_epi * me%surfaceArea
        V_hyp = z_hyp * me%surfaceArea

        ! Settling velocities/rates for today's temperature, and wind resuspension
        call me%setSettlingRate(T_water_t)
        call me%setWindResuspensionRate(t)

        ! --------------------------------------------------------------------
        ! 6. Mass bookkeeping. The canonical stores hold whole-segment totals;
        !    f_epi_* track the epilimnion's share of each store.
        ! --------------------------------------------------------------------
        if (me%volume > C%epsilon) then

            ! --- Inflows enter the epilimnion (surface layer) ---
            dj_spm_in = me%j_spm%inflow + me%j_spm%soilErosion
            call dj_cont_in%add(me%j_contaminant_inflow)
            call dj_cont_in%add(me%j_contaminant_soilErosion)
            call dj_cont_in%add(me%j_contaminant_pointSources)
            call dj_cont_in%add(me%j_contaminant_diffuseSources)
            call dj_cont_in%add(me%j_contaminant_dissolvedRunoff)

            m_spm_old = me%m_spm
            me%m_spm = flushToZero(max(me%m_spm + dj_spm_in, 0.0_dp))
            me%f_epi_spm = shareAfterSurfaceGain(me%f_epi_spm, m_spm_old, dj_spm_in)
            if (.not. C%ignoreContaminant .and. allocated(me%m_contaminant%c) &
                .and. allocated(dj_cont_in%c)) then
                c_old = me%m_contaminant%c
                diss_old = me%m_contaminant%m_dissolved
                call me%m_contaminant%add(dj_cont_in)
                me%f_epi_cont = shareAfterSurfaceGain(me%f_epi_cont, c_old, dj_cont_in%c)
                me%f_epi_diss = shareAfterSurfaceGain(me%f_epi_diss, diss_old, dj_cont_in%m_dissolved)
            end if

            ! --- Settling: epilimnion -> hypolimnion (internal transfer) ---
            if (me%isStratified .and. z_epi > C%epsilon) then
                frac_settle = min(max(me%W_settle_spm * dt / z_epi, 0.0_dp), 1.0_dp)
                ! Only SPM-sorbed contaminant settles with the particles, in proportion
                ! to the epilimnion SPM mass that settles (summed over size classes)
                frac_settle_cont = divideCheckZero(sum(me%m_spm * me%f_epi_spm * frac_settle), &
                                                   sum(me%m_spm * me%f_epi_spm))
                me%f_epi_spm = me%f_epi_spm * (1.0_dp - frac_settle)
                me%f_epi_cont(:,:,PFAS_SPM) = me%f_epi_cont(:,:,PFAS_SPM) * (1.0_dp - frac_settle_cont)
            end if

            ! --- Deposition: bottom layer -> bed sediment ---
            if (me%isStratified) then
                z_hyp = max(z_hyp, C%epsilon)
                frac_dep = min(max(me%W_settle_spm * dt / z_hyp, 0.0_dp), 1.0_dp)
                bot_share_spm = 1.0_dp - me%f_epi_spm
                allocate(bot_share_cont, mold=me%f_epi_cont)
                bot_share_cont = 1.0_dp - me%f_epi_cont
            else
                frac_dep = min(max(me%W_settle_spm * dt / max(me%depth, C%epsilon), 0.0_dp), 1.0_dp)
                bot_share_spm = 1.0_dp
                allocate(bot_share_cont, mold=me%f_epi_cont)
                bot_share_cont = 1.0_dp
            end if
            dj_spm_deposit = flushToZero(min(me%m_spm * bot_share_spm * frac_dep, me%m_spm))
            ! SPM-sorbed contaminant deposits in proportion to the bottom SPM mass deposited
            frac_dep_cont = divideCheckZero(sum(dj_spm_deposit), sum(me%m_spm * bot_share_spm))

            ! Resuspension demand as an area flux; the bed returns what it can't supply
            dj_spm_resus_perArea = flushToZero(me%k_resus * me%bedSediment%Mf_bed_by_size() * dt)
            dj_spm_resus_perArea_ = dj_spm_resus_perArea
            call rslt%addErrors(.errors. me%bedSediment%resuspend(dj_spm_resus_perArea_))
            dj_spm_resus_perArea = dj_spm_resus_perArea - dj_spm_resus_perArea_
            dj_spm_resus = dj_spm_resus_perArea * me%bedArea

            ! Move deposited SPM into the bed sediment
            call rslt%addErrors(.errors. me%depositToBed(dj_spm_deposit))
            ! depositToBed adjusts depth for pore water uptake; keep volume consistent
            me%volume = me%depth * me%surfaceArea

            ! SPM store update (deposition out of the bottom layer; the epilimnion
            ! mass is unchanged, so its share grows)
            m_spm_old = me%m_spm
            me%m_spm = flushToZero(max(me%m_spm - dj_spm_deposit, 0.0_dp))
            me%f_epi_spm = shareAfterBottomLoss(me%f_epi_spm, m_spm_old, dj_spm_deposit)

            if (.not. C%ignoreContaminant .and. allocated(me%m_contaminant%c)) then
                ! SPM-sorbed contaminant deposits with the settling particles. The
                ! package keeps the mass in the PFAS_SPM slot: that is the phase it
                ! leaves the water from, and the slot transferContaminant reads.
                dj_cont_dep%c(:,:,PFAS_SPM) = me%m_contaminant%c(:,:,PFAS_SPM) &
                    * bot_share_cont(:,:,PFAS_SPM) * frac_dep_cont
                c_old = me%m_contaminant%c
                me%m_contaminant%c(:,:,PFAS_SPM) = me%m_contaminant%c(:,:,PFAS_SPM) &
                    - dj_cont_dep%c(:,:,PFAS_SPM)
                me%f_epi_cont(:,:,PFAS_SPM) = shareAfterBottomLoss(me%f_epi_cont(:,:,PFAS_SPM), &
                    c_old(:,:,PFAS_SPM), dj_cont_dep%c(:,:,PFAS_SPM))

                ! Contaminant resuspension is handled solely by resuspend_spm below,
                ! which removes mass from the bed pools; the old scalar proxy on the
                ! bed total created mass and double-counted the resuspend_spm flux.

                ! Bed coupling, one matrix per step (as in the reaches): transfer the
                ! deposit into the bed and the resuspended layers' PFAS into the
                ! resuspension pool, then export that pool to the bottom layer
                if (isZero(me%bedArea)) then
                    dj_spm_deposit_perArea = 0.0_dp
                    dj_spm_resus_perArea = 0.0_dp
                else
                    dj_spm_deposit_perArea = dj_spm_deposit / me%bedArea
                end if
                call me%bedSediment%getmatrix(dj_spm_deposit_perArea, dj_spm_resus_perArea)
                if (me%bedArea > C%epsilon) then
                    j_contam_dep_perArea = dj_cont_dep * (1.0_dp / me%bedArea)
                end if
                call rslt%addErrors(.errors. me%bedSediment%transferContaminant(j_contam_dep_perArea))
                call rslt%addErrors(.errors. me%bedSediment%resuspend_spm(me%bedArea, cont_resus_spm))
                if (allocated(cont_resus_spm%c)) then
                    c_old = me%m_contaminant%c
                    call me%m_contaminant%add(cont_resus_spm)
                    me%f_epi_cont = shareAfterBottomGain(me%f_epi_cont, c_old, cont_resus_spm%c)
                end if
                ! Flush FP noise on the canonical store
                call me%m_contaminant%clamp_negative()
            end if

            ! SPM resuspension into the bottom layer
            m_spm_old = me%m_spm
            me%m_spm = flushToZero(me%m_spm + dj_spm_resus)
            me%f_epi_spm = shareAfterBottomGain(me%f_epi_spm, m_spm_old, dj_spm_resus)

            ! --- Mixing between epilimnion and hypolimnion ---
            if (me%isStratified .and. V_epi > C%epsilon .and. V_hyp > C%epsilon) then
                fd = min(max(DATASET%lakeInterlayerMixingVelocity / 86400.0_dp &
                    * me%surfaceArea * dt / V_epi, 0.0_dp), 1.0_dp)
                fu = min(max(DATASET%lakeInterlayerMixingVelocity / 86400.0_dp &
                    * me%surfaceArea * dt / V_hyp, 0.0_dp), 1.0_dp)
                me%f_epi_spm = shareAfterMixing(me%f_epi_spm, fd, fu)
                if (allocated(me%f_epi_cont)) then
                    me%f_epi_cont = shareAfterMixing(me%f_epi_cont, fd, fu)
                end if
                me%f_epi_diss = shareAfterMixing(me%f_epi_diss, fd, fu)
            end if

            ! --- Outflow, drawn from the epilimnion (surface outlet/spillway) ---
            k_out = min(max(divideCheckZero(V_out, max(V_epi, C%epsilon)), 0.0_dp), 1.0_dp)
            dj_spm_out = flushToZero(min(me%m_spm * me%f_epi_spm * k_out, me%m_spm))
            m_spm_old = me%m_spm
            me%m_spm = flushToZero(max(me%m_spm - dj_spm_out, 0.0_dp))
            me%f_epi_spm = shareAfterSurfaceLoss(me%f_epi_spm, m_spm_old, dj_spm_out)
            if (.not. C%ignoreContaminant .and. allocated(me%m_contaminant%c)) then
                dj_cont_out%c = me%m_contaminant%c * me%f_epi_cont * k_out
                dj_cont_out%m_dissolved = me%m_contaminant%m_dissolved * me%f_epi_diss * k_out
                c_old = me%m_contaminant%c
                diss_old = me%m_contaminant%m_dissolved
                me%m_contaminant%c = max(me%m_contaminant%c - dj_cont_out%c, 0.0_dp)
                me%m_contaminant%m_dissolved = max(me%m_contaminant%m_dissolved &
                    - dj_cont_out%m_dissolved, 0.0_dp)
                me%f_epi_cont = shareAfterSurfaceLoss(me%f_epi_cont, c_old, dj_cont_out%c)
                me%f_epi_diss = shareAfterSurfaceLoss(me%f_epi_diss, diss_old, dj_cont_out%m_dissolved)
            end if

            ! --- Record flows, using the reach sign conventions ---
            me%Q%outflow = -V_out
            me%j_spm%deposition = me%j_spm%deposition - dj_spm_deposit
            me%j_spm%resuspension = me%j_spm%resuspension + dj_spm_resus
            me%j_spm%outflow = me%j_spm%outflow - dj_spm_out
            ! Deposition diagnostics keep the legacy -ve (loss) convention; the
            ! outflow holder is stored positive, since downstream reaches and grid
            ! cell getters accumulate it with a plain add
            call me%j_contaminant_deposition%add_scaled(dj_cont_dep, -1.0_dp)
            call me%j_contaminant_resuspension%add(cont_resus_spm)
            call me%j_contaminant_outflow%add(dj_cont_out)
        else
            ! Dry basin: zero the stores
            me%m_spm = 0.0_dp
            call me%m_contaminant%empty()
            me%f_epi_spm = 1.0_dp
            if (allocated(me%f_epi_cont)) me%f_epi_cont = 1.0_dp
            me%f_epi_diss = 1.0_dp
            me%Q%outflow = -V_out
        end if

        ! --------------------------------------------------------------------
        ! 7. Concentrations, transformations (reactor) and biota
        ! --------------------------------------------------------------------
        me%C_spm = divideCheckZero(me%m_spm, me%volume)

        if (.not. C%ignoreContaminant .and. .not. isZero(me%volume)) then
            ! The reactor holds a pointer to me%m_contaminant (whole-segment total),
            ! so heteroaggregation, dissolution and other transformations run on the
            ! total. The layer shares are left unchanged by transformations - i.e.
            ! both layers are assumed to react at the same rate.
            ! Bed porewater partitioning/transformation and the dissolved sediment-water
            ! exchange, against the water in contact with the bed: the hypolimnion
            ! when stratified, else the whole segment
            if (allocated(me%bedSediment) .and. allocated(me%m_contaminant%c)) then
                call rslt%addErrors(.errors. m_bottom%create())
                if (me%isStratified .and. V_hyp > C%epsilon .and. allocated(me%f_epi_cont)) then
                    m_bottom%c(:,:,PFAS_AQ) = me%m_contaminant%c(:,:,PFAS_AQ) &
                        * (1.0_dp - me%f_epi_cont(:,:,PFAS_AQ))
                    V_bot = V_hyp
                else
                    m_bottom%c(:,:,PFAS_AQ) = me%m_contaminant%c(:,:,PFAS_AQ)
                    V_bot = me%volume
                end if
                call me%bedSediment%update_pfas(dt, m_bottom, V_bot, me%bedArea, dj_exchange)
                c_old = me%m_contaminant%c
                call me%applyBedExchange(dj_exchange)
                ! The exchange is at the bottom, so it changes the epilimnion share
                if (allocated(me%f_epi_cont)) then
                    me%f_epi_cont(:,:,PFAS_AQ) = shareAfterBottomLoss(me%f_epi_cont(:,:,PFAS_AQ), &
                        c_old(:,:,PFAS_AQ), max(0.0_dp, dj_exchange%c(:,:,PFAS_AQ)))
                    c_old(:,:,PFAS_AQ) = c_old(:,:,PFAS_AQ) - max(0.0_dp, dj_exchange%c(:,:,PFAS_AQ))
                    me%f_epi_cont(:,:,PFAS_AQ) = shareAfterBottomGain(me%f_epi_cont(:,:,PFAS_AQ), &
                        c_old(:,:,PFAS_AQ), max(0.0_dp, -dj_exchange%c(:,:,PFAS_AQ)))
                end if
            end if

            ! Refresh the reactor's copies of this step's volume, SPM, settling,
            ! temperature and surface (air-water interface) area before it runs
            call me%reactor%setState(me%volume, T_water_t, me%C_spm, me%W_settle_spm, me%surfaceArea)
            call rslt%addErrors(.errors. me%reactor%update(dt=dt))
            call me%recordSurfaceLosses()
            if (allocated(me%m_contaminant%c)) then
                me%m_contaminant%m_dissolved = sum(me%m_contaminant%c(:,:,PFAS_AQ))
            end if
        end if

        do i = 1, me%nBiota
            c_env_contaminant = me%m_contaminant%divideCheckZero(me%volume)
            ! Dietary exposure from the prey's current burden (food web)
            call update_diets(me%biota)
            call rslt%addErrors(.errors. me%biota(i)%update(t, c_env_contaminant))
            call c_env_contaminant%finalise()
        end do

        ! Finalise temporaries
        call dj_cont_in%finalise()
        call dj_cont_dep%finalise()
        call dj_cont_out%finalise()
        call m_bottom%finalise()
        call dj_exchange%finalise()
        call cont_resus_spm%finalise()
        call j_contam_dep_perArea%finalise()

        call rslt%addToTrace("Updating " // trim(me%ref) // " on timestep #" // trim(str(t)))
        call LOGR%toFile(errors = .errors. rslt)
        call ERROR_HANDLER%trigger(errors = .errors. rslt)
    end subroutine updateLakeSegment

    subroutine finaliseLakeSegment(me)
        class(LakeSegment), intent(inout) :: me
        ! Mirror Reach%finalise (the parent binding can't be called directly as
        ! Reach is abstract): base waterbody finalisation plus the Reach members
        call me%WaterBody%finalise()
        if (allocated(me%inflowsArr)) deallocate(me%inflowsArr)
        if (allocated(me%inflows)) deallocate(me%inflows)
        if (allocated(me%domainOutflow)) deallocate(me%domainOutflow)
        if (allocated(me%f_epi_spm)) deallocate(me%f_epi_spm)
        if (allocated(me%f_epi_cont)) deallocate(me%f_epi_cont)
    end subroutine

!---------------------------------!
!-- LAYER SHARE HELPER FUNCTIONS --!
!---------------------------------!
! f is the epilimnion's share of a store's total mass m. These update f for mass
! entering or leaving one specific layer, keeping 0 <= f <= 1 and defaulting to a
! fully-mixed state (f = 1) when the store empties.

    !> Share after mass dj is added to the epilimnion (inflows, runoff, sources)
    pure elemental function shareAfterSurfaceGain(f, m_old, dj) result(fnew)
        real(dp), intent(in) :: f, m_old, dj
        real(dp) :: fnew
        if (m_old + dj > C%epsilon) then
            fnew = min(max((f * m_old + max(dj, 0.0_dp)) / (m_old + dj), 0.0_dp), 1.0_dp)
        else
            fnew = 1.0_dp
        end if
    end function

    !> Share after mass dj is added to the hypolimnion/bottom layer (resuspension)
    pure elemental function shareAfterBottomGain(f, m_old, dj) result(fnew)
        real(dp), intent(in) :: f, m_old, dj
        real(dp) :: fnew
        if (m_old + dj > C%epsilon) then
            fnew = min(max(f * m_old / (m_old + dj), 0.0_dp), 1.0_dp)
        else
            fnew = 1.0_dp
        end if
    end function

    !> Share after mass dj leaves the hypolimnion/bottom layer (deposition to bed)
    pure elemental function shareAfterBottomLoss(f, m_old, dj) result(fnew)
        real(dp), intent(in) :: f, m_old, dj
        real(dp) :: fnew
        if (m_old - dj > C%epsilon) then
            fnew = min(max(f * m_old / (m_old - dj), 0.0_dp), 1.0_dp)
        else
            fnew = 1.0_dp
        end if
    end function

    !> Share after mass dj leaves the epilimnion (outflow)
    pure elemental function shareAfterSurfaceLoss(f, m_old, dj) result(fnew)
        real(dp), intent(in) :: f, m_old, dj
        real(dp) :: fnew
        if (m_old - dj > C%epsilon) then
            fnew = min(max((f * m_old - dj) / (m_old - dj), 0.0_dp), 1.0_dp)
        else
            fnew = 1.0_dp
        end if
    end function

    !> Share after two-way mixing: fraction fd of the epilimnion's mass moves down
    !! and fraction fu of the hypolimnion's mass moves up. Total mass is unchanged.
    pure elemental function shareAfterMixing(f, fd, fu) result(fnew)
        real(dp), intent(in) :: f, fd, fu
        real(dp) :: fnew
        fnew = min(max(f * (1.0_dp - fd) + (1.0_dp - f) * fu, 0.0_dp), 1.0_dp)
    end function

end module
