module BiotaSoilModule
    use AbstractBiotaModule
    use ResultModule, only: Result
    use ErrorInstanceModule
    use GlobalsModule
    use DataInputModule
    use datetime_module
    use ContaminantModule
    use PFASConstantsModule, only: PFAS_AQ, PFAS_SOL, PFAS_AWI
    use UtilModule, only: str
    implicit none

    type, public, extends(AbstractBiota) :: BiotaSoil
    contains
        procedure :: create => createBiotaSoil
        procedure :: update => updateBiotaSoil
        procedure :: update_soil => updateSoilBiotaSoil
        procedure :: parseInputData => parseInputDataBiotaSoil
    end type

contains

    function createBiotaSoil(me, biotaIndex) result(rslt)
        class(BiotaSoil), intent(inout) :: me
        integer, intent(in)             :: biotaIndex
        type(Result)                    :: rslt

        me%ref = "BiotaSoil_" // trim(str(biotaIndex))
        call rslt%addErrors(.errors. createAbstractBiota(me, biotaIndex))

        ! Default: porewater exposure
        if (len_trim(me%uptakeFromPhase) == 0 .or. trim(me%uptakeFromPhase) == "free" .or. &
            trim(me%uptakeFromPhase) == "water_aq") then
            me%uptakeFromPhase = "soil_aq"
        end if

        call rslt%addToTrace('Creating BiotaSoil')
    end function createBiotaSoil


    !> Update from a concentration object, for compatibility with the AbstractBiota
    !! interface: the AQ slot is taken as the porewater concentration [kg m-3] and the
    !! SOL + AWI slots as concentrations per kg dry soil [kg kg-1]. Soil layers call
    !! update_soil instead, which derives these from the layer's masses.
    function updateBiotaSoil(me, t, C_env_contaminant) result(rslt)
        class(BiotaSoil), intent(inout) :: me
        integer, intent(in)             :: t
        type(Contaminant), intent(in)   :: C_env_contaminant
        type(Result)                    :: rslt
        real(dp), allocatable :: C_pw(:), C_bulk(:)
        integer :: i
        if (.not. allocated(C_env_contaminant%c)) then
            call rslt%addError(ErrorInstance(code=105, message="Contaminant array not allocated in BiotaSoil update"))
            return
        end if
        allocate(C_pw(size(C_env_contaminant%c, 1)), C_bulk(size(C_env_contaminant%c, 1)))
        do i = 1, size(C_pw)
            C_pw(i) = sum(C_env_contaminant%c(i,:,PFAS_AQ))
            C_bulk(i) = sum(C_env_contaminant%c(i,:,PFAS_SOL)) + sum(C_env_contaminant%c(i,:,PFAS_AWI))
        end do
        call rslt%addErrors(.errors. integrateBurden(me, t, soil_biota_exposure(C_pw, C_bulk, trim(me%uptakeFromPhase))))
    end function updateBiotaSoil

    !> Update from a soil layer's PFAS masses. Exposure is either the porewater
    !! concentration C_pw = M_aq / V_pw [kg m-3] ("soil_aq", the default: plant and
    !! soil organism uptake is mainly from porewater, and sorption limits it), or the
    !! total soil concentration C_soil = (M_aq + M_sol + M_awi) / M_soil [kg kg-1 dry
    !! soil] ("soil_total", for BSAF-type parameterisations). k_uptake must be in
    !! units matching the chosen exposure.
    function updateSoilBiotaSoil(me, t, m_contaminant, V_porewater, M_soil) result(rslt)
        class(BiotaSoil), intent(inout) :: me
        integer, intent(in)             :: t
        type(Contaminant), intent(in)   :: m_contaminant     !! Layer PFAS masses [kg]
        real(dp), intent(in)            :: V_porewater       !! Layer porewater volume [m3]
        real(dp), intent(in)            :: M_soil            !! Layer dry soil mass [kg]
        type(Result)                    :: rslt
        real(dp), allocatable :: C_pw(:), C_bulk(:)
        integer :: i
        if (.not. allocated(m_contaminant%c)) then
            call rslt%addError(ErrorInstance(code=105, message="Contaminant array not allocated in BiotaSoil update"))
            return
        end if
        allocate(C_pw(size(m_contaminant%c, 1)), C_bulk(size(m_contaminant%c, 1)))
        C_pw = 0.0_dp
        C_bulk = 0.0_dp
        do i = 1, size(C_pw)
            if (V_porewater > C%epsilon) C_pw(i) = sum(m_contaminant%c(i,:,PFAS_AQ)) / V_porewater
            if (M_soil > C%epsilon) C_bulk(i) = (sum(m_contaminant%c(i,:,PFAS_AQ)) &
                + sum(m_contaminant%c(i,:,PFAS_SOL)) + sum(m_contaminant%c(i,:,PFAS_AWI))) / M_soil
        end do
        call rslt%addErrors(.errors. integrateBurden(me, t, soil_biota_exposure(C_pw, C_bulk, trim(me%uptakeFromPhase))))
    end function updateSoilBiotaSoil

    !> First-order uptake/elimination/growth-dilution of the active and stored
    !! burdens over one time step, solved exactly for a constant exposure.
    function integrateBurden(me, t, exposure) result(rslt)
        class(BiotaSoil), intent(inout) :: me
        integer, intent(in)  :: t
        real(dp), intent(in) :: exposure(:)
        type(Result)         :: rslt
        real(dp) :: dt_days, lambda_active, lambda_stored, target_active, target_stored
        type(datetime) :: currentDate
        integer :: i, nSpecies

        nSpecies = size(exposure)
        if (.not. allocated(me%C_active)) then
            call rslt%addErrors(.errors. me%create(me%biotaIndex))
        else if (size(me%C_active) /= nSpecies) then
            call me%finalise()
            call rslt%addErrors(.errors. me%create(me%biotaIndex))
        end if
        if (.not. allocated(me%C_active)) return

        currentDate = C%startDate + timedelta(days=t-1)
        if (me%harvestInMonth > 0 .and. me%harvestInMonth == currentDate%getMonth()) then
            me%C_active = 0.0_dp
            me%C_stored = 0.0_dp
            return
        end if

        dt_days = real(C%timeStep, dp) / 86400.0_dp
        do i = 1, nSpecies
            lambda_active = max(0.0_dp, me%k_elim(i) + me%k_growth + me%k_death)
            if (lambda_active > C%epsilon) then
                target_active = me%uptake(i, exposure(i)) * (1.0_dp - me%storedFraction) / lambda_active
                me%C_active(i) = target_active + (me%C_active(i) - target_active) * exp(-lambda_active * dt_days)
            else
                me%C_active(i) = me%C_active(i) + &
                                 me%uptake(i, exposure(i)) * (1.0_dp - me%storedFraction) * dt_days
            end if
            lambda_stored = max(0.0_dp, me%k_growth + me%k_death)
            if (lambda_stored > C%epsilon) then
                target_stored = me%uptake(i, exposure(i)) * me%storedFraction / lambda_stored
                me%C_stored(i) = target_stored + (me%C_stored(i) - target_stored) * exp(-lambda_stored * dt_days)
            else
                me%C_stored(i) = me%C_stored(i) + me%uptake(i, exposure(i)) * me%storedFraction * dt_days
            end if
        end do
        ! Metabolism of precursors in the organism
        call me%transformBurden(dt_days)
    end function integrateBurden

    function parseInputDataBiotaSoil(me) result(rslt)
        class(BiotaSoil), intent(inout) :: me
        type(Result)                    :: rslt

        call rslt%addErrors(.errors. parseInputDataAbstractBiota(me))

        if (len_trim(me%uptakeFromPhase) == 0 .or. trim(me%uptakeFromPhase) == "water_aq" .or. &
            trim(me%uptakeFromPhase) == "free") then
            me%uptakeFromPhase = "soil_aq"
        end if
    end function parseInputDataBiotaSoil


    !> Exposure concentration by the configured phase selector.
    function soil_biota_exposure(C_pw, C_bulk, source) result(exposure)
        real(dp), intent(in)         :: C_pw(:)      !! Porewater concentration [kg m-3]
        real(dp), intent(in)         :: C_bulk(:)    !! Total soil concentration [kg kg-1 dry soil]
        character(len=*), intent(in) :: source
        real(dp) :: exposure(size(C_pw))
        select case(trim(source))
        case("soil_total", "soil_aq_sol", "soil_aq_sol_awi", "free_and_attached")
            exposure = C_bulk
        case default
            exposure = C_pw
        end select
    end function soil_biota_exposure
end module BiotaSoilModule
