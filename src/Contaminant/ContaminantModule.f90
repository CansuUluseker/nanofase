!> PFAS contaminant module.
!!
!! A `Contaminant` object holds PFAS MASS [kg] (or mass per unit area [kg/m2] in the
!! bed sediment pools, or a concentration when built with divideCheckZero) in the
!! array c(species, form, phase):
!!   - species: the PFAS compounds (e.g. PFOS, PFOA, a precursor), n_pfas_species;
!!   - form:    a legacy grouping axis kept from NanoFASE (usually 1);
!!   - phase:   where the mass sits, indexed by the PFASConstantsModule parameters
!!              PFAS_AQ   dissolved in water or porewater,
!!              PFAS_SOL  sorbed to soil or bed-sediment solids,
!!              PFAS_SPM  sorbed to suspended particulate matter (water column),
!!              PFAS_AWI  adsorbed at the air-water interface (soil pores, water surface),
!!              PFAS_FOAM exported to foam (only used in flux packages),
!!              PFAS_AIR  lost to the atmosphere (only used in flux packages).


module ContaminantModule
    use GlobalsModule, only: dp, C
    use ResultModule, only: Result, Result0D
    use ErrorInstanceModule
    use LoggerModule, only: LOGR
    use DataInputModule, only: DATASET
    use PFASConstantsModule
    use PFASProcessModule
    use PFASParameterisationModule, only: pfas_kd_soil, pfas_kd_spm, pfas_kd_sediment, pfas_kia_water
    implicit none

    private
    public :: Contaminant
    public :: operator(+), operator(-), operator(*)

    !> PFAS mass (or concentration) by species, form and phase, plus the per-species
    !! process parameters that apply to it. See the module description.
    type, public :: Contaminant
        real(dp), allocatable :: c(:,:,:)          !! PFAS mass c(species, form, phase) [kg] (bed pools [kg/m2])
        real(dp) :: m_dissolved = 0.0_dp           !! Legacy alias kept in sync with sum(c(:,:,PFAS_AQ)) [kg]

        ! Legacy NanoFASE scalar fields, retained only for backwards-compatible interfaces;
        ! none of them drives a PFAS process except k_transform_pristine (see create_from_data).
        real(dp) :: rho_contaminant = 0.0_dp
        real(dp) :: k_diss_pristine = 0.0_dp
        real(dp) :: k_diss_transformed = 0.0_dp
        real(dp) :: k_transform_pristine = 0.0_dp
        real(dp) :: alpha_hetero = 0.0_dp
        real(dp) :: alpha_att = 0.0_dp
        character(len=100) :: compartment = ''    !! 'soil', 'water', 'estuary', 'sediment' or 'atmospheric'

        ! PFAS species properties and process parameters, one entry per species.
        character(len=PFAS_MAX_NAME), allocatable :: species_name(:)  !! Species names (&pfas species_names, else PFAS_i)
        real(dp), allocatable :: molecular_weight(:), charge(:), pka(:)  !! Informational properties [g/mol], [-], [-]
        real(dp), allocatable :: kd_solid(:)       !! Solid-water partition coefficient for soil/bed solids, this cell [m3 kg-1]
        real(dp), allocatable :: kd_spm(:)         !! SPM-water partition coefficient, this cell (incl. salinity) [m3 kg-1]
        real(dp), allocatable :: kia_awi(:)        !! Air-water interfacial adsorption coefficient, Gamma = Kia*C_w [m]
        real(dp), allocatable :: k_mt_solid(:), k_mt_spm(:), k_mt_awi(:)  !! First-order mass-transfer rates towards
                                                   !! sorption equilibrium; 0 = instantaneous equilibrium [s-1]
        real(dp), allocatable :: k_react(:,:)      !! Transformation rate from species i to species j [s-1]
        real(dp), allocatable :: reaction_yield(:,:)  !! Yield of species j per unit mass of species i transformed [-]

        ! Diagnostic per-species fluxes of the last process call on this object
        ! [kg timestep-1]; reset by reset_fluxes / the update procedures.
        real(dp), allocatable :: j_to_groundwater(:)  !! Leached from a soil layer (set by leach)
        real(dp), allocatable :: j_to_downstream(:)   !! Not currently set by any procedure
        real(dp), allocatable :: j_to_sediment(:)     !! Set by sediment_exchange (legacy, unused)
        real(dp), allocatable :: j_to_atmosphere(:)   !! Volatilised + aerosolised (set by atmosphere_exchange)
        real(dp), allocatable :: j_to_foam(:)         !! Exported to foam (set by foam_exchange)
    contains
        ! Construction, destruction and housekeeping
        procedure :: create => contaminant_create
        procedure :: create_from_data => contaminant_create_from_data
        ! Arithmetic on the mass array
        procedure :: add => contaminant_add
        procedure :: add_scaled => contaminant_add_scaled
        procedure :: multiply_scalar => contaminant_multiply_scalar
        procedure :: divideCheckZero => contaminant_divideCheckZero
        procedure :: finalise => contaminant_finalise
        procedure :: empty => contaminant_empty
        procedure :: clamp_negative => contaminant_clamp_negative
        ! Compartment updates (partitioning + transformation over a time step)
        procedure :: update => contaminant_update
        procedure :: update_water => contaminant_update_water
        procedure :: update_sediment => contaminant_update_sediment
        procedure :: update_soil => contaminant_update_soil
        ! Individual processes
        procedure :: partition_phases => contaminant_partition_phases
        procedure :: move_phase => contaminant_move_phase
        procedure :: transformation => contaminant_transformation
        procedure :: transform_network => contaminant_transformation   ! Alias of transformation
        ! Flux packages exchanged between compartments
        procedure :: leach => contaminant_leach
        procedure :: erosion_export => contaminant_erosion_export
        procedure :: deposition => contaminant_deposition
        procedure :: sediment_exchange => contaminant_sediment_exchange
        procedure :: outflow_split => contaminant_outflow_split
        procedure :: foam_exchange => contaminant_foam_exchange
        procedure :: atmosphere_exchange => contaminant_atmosphere_exchange
        ! Getters
        procedure :: getConcentration => contaminant_getConcentration
        procedure :: get_free => contaminant_get_free
        procedure :: get_attached => contaminant_get_attached
        procedure :: get_phase => contaminant_get_phase
        procedure :: mass_total => contaminant_mass_total
        procedure :: reset_fluxes => contaminant_reset_fluxes

        ! Legacy NanoFASE nanoparticle procedures, kept so old call sites still compile.
        ! They do nothing (or return zero) for PFAS.
        procedure :: heteroaggregation => contaminant_noop_result
        procedure :: dissolution => contaminant_noop_result_dt
        procedure :: attachment => contaminant_noop_result_att
        procedure :: calculateCollisionRate => contaminant_calculateCollisionRate
        procedure :: calculateParticleConcentration => contaminant_calculateParticleConcentration
        procedure :: calculateAttachmentRate => contaminant_calculateAttachmentRate
        procedure :: calculateSettlingVelocity => contaminant_calculateSettlingVelocity
    end type

    !> Contaminant * scalar: a new object with every mass scaled (e.g. per-area conversion)
    interface operator(*)
        module procedure multiply_contaminant_scalar
    end interface

    !> Contaminant: a new object with every mass negated
    interface operator(-)
        module procedure negate_contaminant
    end interface

    !> Contaminant + Contaminant: a new object with the masses summed elementwise
    interface operator(+)
        module procedure add_contaminant
    end interface

contains

    !> A local copy of UtilModule's str, to avoid
    !! a UtilModule dependency (and possible circular module dependencies).
    !! Not currently called.
    function int_to_string(i) result(s)
        integer, intent(in) :: i
        character(len=20) :: s
        write(s, '(I0)') i
        s = trim(adjustl(s))
    end function


    !> Allocate an empty Contaminant: all masses, fluxes and process parameters zero.
    !! Array sizes come from the configuration: species and forms from
    !! C%contaminantDim(1:2), and at least PFAS_NPHASES phases. 
    function contaminant_create(this) result(r)
        class(Contaminant), intent(inout) :: this
        type(Result) :: r
        type(ErrorInstance) :: err
        integer :: ns, nf, np, stat

        ns = max(1, C%contaminantDim(1))
        nf = max(1, C%contaminantDim(2))
        np = max(PFAS_NPHASES, C%contaminantDim(3))

        call this%finalise()
        allocate(this%c(ns,nf,np), this%species_name(ns), this%molecular_weight(ns), &
                 this%charge(ns), this%pka(ns), this%kd_solid(ns), this%kd_spm(ns), &
                 this%kia_awi(ns), this%k_mt_solid(ns), this%k_mt_spm(ns), this%k_mt_awi(ns), &
                 this%k_react(ns,ns), this%reaction_yield(ns,ns), &
                 this%j_to_groundwater(ns), this%j_to_downstream(ns), this%j_to_sediment(ns), &
                 this%j_to_atmosphere(ns), this%j_to_foam(ns), stat=stat)
        if (stat /= 0) then
            err = ErrorInstance(code=901, message='P-FASE Contaminant allocation failed')
            call r%addError(err)
            call LOGR%toFile(errors=r%errors)
            return
        end if

        call this%empty()
        this%species_name = ''
        this%molecular_weight = 0.0_dp
        this%charge = 0.0_dp
        this%pka = 0.0_dp
        this%kd_solid = 0.0_dp
        this%kd_spm = 0.0_dp
        this%kia_awi = 0.0_dp
        this%k_mt_solid = 0.0_dp
        this%k_mt_spm = 0.0_dp
        this%k_mt_awi = 0.0_dp
        this%k_react = 0.0_dp
        this%reaction_yield = 0.0_dp
    end function


    !> Create a Contaminant store for a compartment, with its PFAS process parameters.
    !! 1. Allocates an empty object (create).
    !! 2. Copies the species properties and process rates from the &pfas constants group
    !!    (DATASET): Kd, Kia, mass-transfer rates, transformation rates and yields. A
    !!    'sediment' store takes the bed-specific rates reaction_rate_sediment.
    !! 3. If the grid cell (x, y) is given, replaces the constants with the spatially
    !!    explicit values for that cell (PFASParameterisationModule): soil Kd for 'soil',
    !!    bed Kd for 'sediment', SPM Kd and salinity-adjusted Kia for 'water'/'estuary'.
    !! 4. Names unnamed species PFAS_i and, if the legacy k_transform_pristine > 0, links
    !!    species i -> i+1 with that rate wherever reaction_rate gives no rate.
    function contaminant_create_from_data(this, compartment, contaminantDensity, soilAttachmentEfficiency, &
                                          riverAttachmentEfficiency, estuaryAttachmentEfficiency, &
                                          k_diss_pristine, k_diss_transformed, k_transform_pristine, &
                                          waterTemperature, x, y) result(r)
        class(Contaminant), intent(inout) :: this
        character(len=*), intent(in) :: compartment
        real(dp), intent(in), optional :: contaminantDensity, soilAttachmentEfficiency
        real(dp), intent(in), optional :: riverAttachmentEfficiency, estuaryAttachmentEfficiency
        real(dp), intent(in), optional :: k_diss_pristine, k_diss_transformed, k_transform_pristine
        real(dp), intent(in), optional :: waterTemperature
        integer, intent(in), optional :: x, y       !! Grid cell, for the spatially explicit Kd
        type(Result) :: r
        integer :: i

        r = this%create()
        if (r%hasCriticalError()) return

        this%compartment = trim(compartment)
        if (present(contaminantDensity)) this%rho_contaminant = contaminantDensity
        if (present(k_diss_pristine)) this%k_diss_pristine = k_diss_pristine
        if (present(k_diss_transformed)) this%k_diss_transformed = k_diss_transformed
        if (present(k_transform_pristine)) this%k_transform_pristine = k_transform_pristine

        select case(trim(compartment))
        case('soil','atmospheric')
            if (present(soilAttachmentEfficiency)) this%alpha_att = soilAttachmentEfficiency
        case('water')
            if (present(riverAttachmentEfficiency)) this%alpha_att = riverAttachmentEfficiency
        case('estuary','sediment')
            if (present(estuaryAttachmentEfficiency)) this%alpha_att = estuaryAttachmentEfficiency
        end select


        ! Species-resolved PFAS properties from the &pfas constants group
        call copy_species(DATASET%pfasMolecularWeight, this%molecular_weight)
        call copy_species(DATASET%pfasCharge, this%charge)
        call copy_species(DATASET%pfasPka, this%pka)
        call copy_species(DATASET%pfasKdSoil, this%kd_solid)
        call copy_species(DATASET%pfasKdSpm, this%kd_spm)
        call copy_species(DATASET%pfasKia, this%kia_awi)
        call copy_species(DATASET%pfasKMtSoil, this%k_mt_solid)
        call copy_species(DATASET%pfasKMtSpm, this%k_mt_spm)
        call copy_species(DATASET%pfasKMtAwi, this%k_mt_awi)
        if (trim(compartment) == 'sediment' .and. allocated(DATASET%pfasReactionRateSediment)) then
            ! Bed sediment is largely anoxic, so it may have its own (slower/different) rates
            i = min(size(this%k_react,1), size(DATASET%pfasReactionRateSediment,1))
            this%k_react(:i,:i) = DATASET%pfasReactionRateSediment(:i,:i)
        else if (allocated(DATASET%pfasReactionRate)) then
            i = min(size(this%k_react,1), size(DATASET%pfasReactionRate,1))
            this%k_react(:i,:i) = DATASET%pfasReactionRate(:i,:i)
        end if
        if (allocated(DATASET%pfasReactionYield)) then
            i = min(size(this%reaction_yield,1), size(DATASET%pfasReactionYield,1))
            this%reaction_yield(:i,:i) = DATASET%pfasReactionYield(:i,:i)
        end if

        ! Kd for this grid cell, from the soil/water property maps or the constants
        if (present(x) .and. present(y)) then
            select case(trim(compartment))
            case('soil')
                this%kd_solid = pfas_kd_soil(x, y)
            case('sediment')
                this%kd_solid = pfas_kd_sediment(x, y)
            case('water', 'estuary')
                this%kd_spm = pfas_kd_spm(x, y)
                this%kia_awi = pfas_kia_water(x, y)
            end select
        end if

        do i = 1, size(this%c,1)
            if (allocated(DATASET%pfasSpeciesNames)) then
                if (i <= size(DATASET%pfasSpeciesNames)) this%species_name(i) = DATASET%pfasSpeciesNames(i)
            end if
            if (len_trim(this%species_name(i)) == 0) write(this%species_name(i),'(A,I0)') 'PFAS_', i
            ! Legacy single-chain transformation, only where &pfas reaction_rate gives none
            if (i < size(this%c,1) .and. this%k_transform_pristine > 0.0_dp) then
                if (this%k_react(i,i+1) <= 0.0_dp) then
                    this%k_react(i,i+1) = this%k_transform_pristine
                    this%reaction_yield(i,i+1) = 1.0_dp
                end if
            end if
        end do

    contains

        !> Copy a species-resolved constants array into a parameter array, as far as
        !! both sizes allow; leaves dst unchanged if the constants are absent.
        subroutine copy_species(src, dst)
            real(dp), allocatable, intent(in) :: src(:)
            real(dp), intent(inout) :: dst(:)
            integer :: n
            if (.not. allocated(src)) return
            n = min(size(src), size(dst))
            dst(:n) = src(:n)
        end subroutine
    end function

    !> Deallocate all arrays and clear the scalar state. The object can be reused
    !! after a new create / create_from_data.
    subroutine contaminant_finalise(this)
        class(Contaminant), intent(inout) :: this
        if (allocated(this%c)) deallocate(this%c)
        if (allocated(this%species_name)) deallocate(this%species_name)
        if (allocated(this%molecular_weight)) deallocate(this%molecular_weight)
        if (allocated(this%charge)) deallocate(this%charge)
        if (allocated(this%pka)) deallocate(this%pka)
        if (allocated(this%kd_solid)) deallocate(this%kd_solid)
        if (allocated(this%kd_spm)) deallocate(this%kd_spm)
        if (allocated(this%kia_awi)) deallocate(this%kia_awi)
        if (allocated(this%k_mt_solid)) deallocate(this%k_mt_solid)
        if (allocated(this%k_mt_spm)) deallocate(this%k_mt_spm)
        if (allocated(this%k_mt_awi)) deallocate(this%k_mt_awi)
        if (allocated(this%k_react)) deallocate(this%k_react)
        if (allocated(this%reaction_yield)) deallocate(this%reaction_yield)
        if (allocated(this%j_to_groundwater)) deallocate(this%j_to_groundwater)
        if (allocated(this%j_to_downstream)) deallocate(this%j_to_downstream)
        if (allocated(this%j_to_sediment)) deallocate(this%j_to_sediment)
        if (allocated(this%j_to_atmosphere)) deallocate(this%j_to_atmosphere)
        if (allocated(this%j_to_foam)) deallocate(this%j_to_foam)
        this%m_dissolved = 0.0_dp
        this%compartment = ''
    end subroutine contaminant_finalise

    !> Set all masses and diagnostic fluxes to zero, keeping the allocation and the
    !! process parameters (e.g. to reset a flux holder at the start of a time step).
    subroutine contaminant_empty(this)
        class(Contaminant), intent(inout) :: this
        if (allocated(this%c)) this%c = 0.0_dp
        this%m_dissolved = 0.0_dp
        call this%reset_fluxes()
    end subroutine contaminant_empty

    !> Zero the per-species diagnostic fluxes (j_to_*), leaving the masses untouched.
    subroutine contaminant_reset_fluxes(this)
        class(Contaminant), intent(inout) :: this
        if (allocated(this%j_to_groundwater)) this%j_to_groundwater = 0.0_dp
        if (allocated(this%j_to_downstream))  this%j_to_downstream  = 0.0_dp
        if (allocated(this%j_to_sediment))    this%j_to_sediment    = 0.0_dp
        if (allocated(this%j_to_atmosphere))  this%j_to_atmosphere  = 0.0_dp
        if (allocated(this%j_to_foam))        this%j_to_foam        = 0.0_dp
    end subroutine contaminant_reset_fluxes

    !> Keep the legacy m_dissolved scalar in sync with the aqueous phase mass.
    !! Called at the end of every procedure that changes c.
    subroutine update_m_dissolved(this)
        class(Contaminant), intent(inout) :: this
        if (allocated(this%c)) then
            this%m_dissolved = sum(this%c(:,:,PFAS_AQ))
        else
            this%m_dissolved = 0.0_dp
        end if
    end subroutine update_m_dissolved


    subroutine contaminant_add(this, addition)
        class(Contaminant), intent(inout) :: this
        type(Contaminant), intent(in) :: addition
        if (.not. allocated(this%c) .or. .not. allocated(addition%c)) return
        this%c = this%c + addition%c
        call update_m_dissolved(this)
    end subroutine contaminant_add

    !> Apply a flux package with a sign (scale = -1 to remove it) or a fraction of it.
    subroutine contaminant_add_scaled(this, addition, scale)
        class(Contaminant), intent(inout) :: this
        type(Contaminant), intent(in) :: addition
        real(dp), intent(in) :: scale
        if (.not. allocated(this%c) .or. .not. allocated(addition%c)) return
        this%c = this%c + addition%c * scale
        call update_m_dissolved(this)
    end subroutine contaminant_add_scaled

    !> Flush negative masses (from floating-point noise in flux arithmetic) to
    !! zero. 
    subroutine contaminant_clamp_negative(this)
        class(Contaminant), intent(inout) :: this
        if (allocated(this%c)) this%c = max(0.0_dp, this%c)
        call update_m_dissolved(this)
    end subroutine contaminant_clamp_negative

    !> Implements operator(+): returns a new object holding this%c + other%c, with the
    !! process parameters of the first operand.
    function add_contaminant(this, other) result(out)
        type(Contaminant), intent(in) :: this, other
        type(Contaminant) :: out
        type(Result) :: r
        r = out%create()
        if (r%hasCriticalError()) return
        if (allocated(this%c) .and. allocated(other%c)) out%c = this%c + other%c
        call copy_params(this, out)
        call update_m_dissolved(out)
    end function add_contaminant

    !> Implements unary operator(-): returns a new object holding -this%c.
    function negate_contaminant(this) result(out)
        type(Contaminant), intent(in) :: this
        type(Contaminant) :: out
        type(Result) :: r
        r = out%create()
        if (r%hasCriticalError()) return
        if (allocated(this%c)) out%c = -this%c
        call copy_params(this, out)
        call update_m_dissolved(out)
    end function negate_contaminant

    !> Implements operator(*): returns a new object holding this%c * scalar, e.g. to
    !! convert an absolute package to per unit area (scalar = 1/bedArea).
    function multiply_contaminant_scalar(this, scalar) result(out)
        type(Contaminant), intent(in) :: this
        real(dp), intent(in) :: scalar
        type(Contaminant) :: out
        type(Result) :: r
        r = out%create()
        if (r%hasCriticalError()) return
        if (allocated(this%c)) out%c = this%c * scalar
        call copy_params(this, out)
        call update_m_dissolved(out)
    end function multiply_contaminant_scalar

    !> Re-create this object as source * scalar (masses scaled, parameters copied).
    !! Typical use: this reach's share of a cell-level package, scalar = length ratio.
    subroutine contaminant_multiply_scalar(this, source, scalar)
        class(Contaminant), intent(inout) :: this
        type(Contaminant), intent(in) :: source
        real(dp), intent(in) :: scalar
        type(Result) :: r
        r = this%create()
        if (r%hasCriticalError()) return
        if (allocated(source%c)) this%c = source%c * scalar
        call copy_params(source, this)
        call update_m_dissolved(this)
    end subroutine contaminant_multiply_scalar

    !> Returns a new object holding this%c / denominator, or all zeros when the
    !! denominator is (near) zero. 
    function contaminant_divideCheckZero(this, denominator) result(divided)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: denominator
        type(Contaminant) :: divided
        type(Result) :: r
        r = divided%create()
        if (r%hasCriticalError()) return
        if (allocated(this%c) .and. abs(denominator) > C%epsilon) divided%c = this%c / denominator
        call copy_params(this, divided)
        call update_m_dissolved(divided)
    end function contaminant_divideCheckZero

    !> Copy the process parameters and legacy scalars (not the masses) from src to dst,
    !! so that derived objects (sums, scaled copies, packages) keep the cell-specific
    !! Kd, Kia and rates of the store they came from.
    subroutine copy_params(src, dst)
        type(Contaminant), intent(in) :: src
        type(Contaminant), intent(inout) :: dst
        dst%compartment = src%compartment
        dst%rho_contaminant = src%rho_contaminant
        dst%k_diss_pristine = src%k_diss_pristine
        dst%k_diss_transformed = src%k_diss_transformed
        dst%k_transform_pristine = src%k_transform_pristine
        dst%alpha_hetero = src%alpha_hetero
        dst%alpha_att = src%alpha_att
        if (allocated(src%species_name) .and. allocated(dst%species_name)) dst%species_name = src%species_name
        if (allocated(src%molecular_weight) .and. allocated(dst%molecular_weight)) dst%molecular_weight = src%molecular_weight
        if (allocated(src%charge) .and. allocated(dst%charge)) dst%charge = src%charge
        if (allocated(src%pka) .and. allocated(dst%pka)) dst%pka = src%pka
        if (allocated(src%kd_solid) .and. allocated(dst%kd_solid)) dst%kd_solid = src%kd_solid
        if (allocated(src%kd_spm) .and. allocated(dst%kd_spm)) dst%kd_spm = src%kd_spm
        if (allocated(src%kia_awi) .and. allocated(dst%kia_awi)) dst%kia_awi = src%kia_awi
        if (allocated(src%k_mt_solid) .and. allocated(dst%k_mt_solid)) dst%k_mt_solid = src%k_mt_solid
        if (allocated(src%k_mt_spm) .and. allocated(dst%k_mt_spm)) dst%k_mt_spm = src%k_mt_spm
        if (allocated(src%k_mt_awi) .and. allocated(dst%k_mt_awi)) dst%k_mt_awi = src%k_mt_awi
        if (allocated(src%k_react) .and. allocated(dst%k_react)) dst%k_react = src%k_react
        if (allocated(src%reaction_yield) .and. allocated(dst%reaction_yield)) dst%reaction_yield = src%reaction_yield
    end subroutine copy_params

    !> Run the in-compartment PFAS processes over one time step dt [s]: partitioning
    !! between phases, then transformation. Dispatches on the compartment name to
    !! update_soil, update_water (also used for 'estuary') or update_sediment.
    !! This is what the Reactor calls.
    function contaminant_update(this, dt, T_water, C_spm, W_settle_spm, G, volume, compartment, &
                                k_att, alpha_att, mass_solid, awi_area) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt, T_water, C_spm(:), W_settle_spm(:), G, volume
        character(len=*), intent(in) :: compartment
        real(dp), intent(in), optional :: k_att(:), alpha_att
        real(dp), intent(in), optional :: mass_solid    !! Solid mass available for Kd partitioning [kg]
        real(dp), intent(in), optional :: awi_area      !! Air-water interfacial area (soil) or water surface area [m2]
        type(Result) :: r
        select case(trim(compartment))
        case('soil')
            call r%addErrors(.errors. this%update_soil(dt,T_water,C_spm,W_settle_spm,G,volume,k_att,alpha_att, &
                                                       mass_solid,awi_area))
        case('water','estuary')
            call r%addErrors(.errors. this%update_water(dt,T_water,C_spm,W_settle_spm,G,volume,awi_area))
        case('sediment')
            call r%addErrors(.errors. this%update_sediment(dt,T_water,C_spm,W_settle_spm,G,volume,mass_solid))
        end select
        call update_m_dissolved(this)
    end function contaminant_update

    !> Soil layer step: three-phase partitioning between porewater (AQ), soil solids
    !! (SOL, Kd = kd_solid) and the pore air-water interface (AWI, Kia), then the
    !! transformation network. SoilLayer%update calls partition_phases and
    !! transformation directly (with leaching in between); this wrapper serves the
    !! Reactor 'soil' path.
    function contaminant_update_soil(this, dt, T_water, C_spm, W_settle_spm, G, volume, k_att, alpha_att, &
                                     mass_solid, awi_area) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt, T_water, C_spm(:), W_settle_spm(:), G, volume
        real(dp), intent(in), optional :: k_att(:), alpha_att
        real(dp), intent(in), optional :: mass_solid    !! Soil solid mass, bulkDensity*area*depth [kg]
        real(dp), intent(in), optional :: awi_area      !! Air-water interfacial area [m2]
        type(Result) :: r
        real(dp) :: soil_mass, awi
        call this%reset_fluxes()
        ! volume is the porewater volume [m3]. The soil solid mass must come from the
        ! caller (bulk density x layer volume)
        if (present(mass_solid)) then
            soil_mass = max(C%epsilon, mass_solid)
        else
            soil_mass = max(C%epsilon, volume)
        end if
        ! AWI adsorption is disabled unless the caller supplies an interfacial area
        awi = 0.0_dp
        if (present(awi_area)) awi = max(0.0_dp, awi_area)
        call this%partition_phases(dt, volume, PFAS_SOL, soil_mass, awi)
        call r%addErrors(.errors. this%transformation(dt))
        call update_m_dissolved(this)
    end function contaminant_update_soil

    !> Water column step (rivers, lakes, estuaries): three-phase partitioning between
    !! dissolved (AQ), SPM-sorbed (SPM, Kd = kd_spm, SPM mass = sum(C_spm)*volume) and
    !! the water-surface microlayer (AWI, Kia, area = awi_area), then transformation.
    function contaminant_update_water(this, dt, T_water, C_spm, W_settle_spm, G, volume, awi_area) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt, T_water, C_spm(:), W_settle_spm(:), G, volume
        real(dp), intent(in), optional :: awi_area      !! Water surface area [m2]
        type(Result) :: r
        real(dp) :: spm_mass, awi
        call this%reset_fluxes()
        spm_mass = max(0.0_dp, sum(C_spm) * max(volume, 0.0_dp))
        ! The surface microlayer is the air-water interface of an open water body
        awi = 0.0_dp
        if (present(awi_area)) awi = max(0.0_dp, awi_area)
        call this%partition_phases(dt, volume, PFAS_SPM, spm_mass, awi)
        call r%addErrors(.errors. this%transformation(dt))
        call update_m_dissolved(this)
    end function contaminant_update_water

    !> Bed sediment step: two-phase partitioning between porewater (AQ) and bed solids
    !! (SOL, Kd = kd_solid), then transformation. The live bed path is
    !! AbstractBedSediment%update_pfas, which calls partition_phases and transformation
    !! per layer (per unit area).
    function contaminant_update_sediment(this, dt, T_water, C_spm, W_settle_spm, G, volume, mass_solid) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt, T_water, C_spm(:), W_settle_spm(:), G, volume
        real(dp), intent(in), optional :: mass_solid    !! Bed solid mass for Kd partitioning [kg]
        type(Result) :: r
        real(dp) :: solid
        call this%reset_fluxes()
        ! volume is the porewater volume [m3]; callers should pass the bed solid
        ! mass, otherwise it falls back to the (dimensionally wrong) legacy placeholder
        if (present(mass_solid)) then
            solid = max(C%epsilon, mass_solid)
        else
            solid = max(C%epsilon, volume)
        end if
        call this%partition_phases(dt, max(C%epsilon, volume), PFAS_SOL, solid, 0.0_dp)
        call r%addErrors(.errors. this%transformation(dt))
        call update_m_dissolved(this)
    end function contaminant_update_sediment

    !> Partition between the aqueous phase, one sorbed phase and the air-water
    !! interface, using this object's (cell-specific) Kd, Kia and mass-transfer rates.
    !! At equilibrium the total is shared in proportion to the capacities
    !! V_w : Kd*M_sorbent : Kia*A_aw; with a mass-transfer rate k > 0 each pool relaxes
    !! towards its equilibrium mass as exp(-k*dt) (rate-limited sorption). Mass is
    !! conserved per species and form. The sorbed phase is PFAS_SPM (uses kd_spm,
    !! k_mt_spm) or a solid phase (uses kd_solid, k_mt_solid). With dt = 0, equilibrium
    !! pools re-partition instantly and rate-limited pools are left unchanged.
    !! PFASProcessModule%pfas_partition_phases.
    subroutine contaminant_partition_phases(this, dt, volume_water, sorbed_phase, mass_sorbent, awi_area)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: volume_water    !! Water volume [m3]
        integer, intent(in)  :: sorbed_phase    !! PFAS_SOL (soil/bed) or PFAS_SPM (water column)
        real(dp), intent(in) :: mass_sorbent    !! Solid or SPM mass [kg]
        real(dp), intent(in) :: awi_area        !! Air-water interfacial area [m2]
        if (.not. allocated(this%c)) return
        if (sorbed_phase == PFAS_SPM) then
            call pfas_partition_phases(this%c, dt, volume_water, PFAS_SPM, mass_sorbent, this%kd_spm, &
                                       this%k_mt_spm, awi_area, this%kia_awi, this%k_mt_awi)
        else
            call pfas_partition_phases(this%c, dt, volume_water, sorbed_phase, mass_sorbent, this%kd_solid, &
                                       this%k_mt_solid, awi_area, this%kia_awi, this%k_mt_awi)
        end if
        call update_m_dissolved(this)
    end subroutine contaminant_partition_phases

    !> Move all mass from one phase slot to another, e.g. soil-sorbed PFAS on eroded
    !! soil becoming SPM-sorbed PFAS when it enters a water body.
    subroutine contaminant_move_phase(this, from_phase, to_phase)
        class(Contaminant), intent(inout) :: this
        integer, intent(in) :: from_phase, to_phase
        if (.not. allocated(this%c) .or. from_phase == to_phase) return
        this%c(:,:,to_phase) = this%c(:,:,to_phase) + this%c(:,:,from_phase)
        this%c(:,:,from_phase) = 0.0_dp
        call update_m_dissolved(this)
    end subroutine contaminant_move_phase

    !> Transform precursors to other species over dt [s] with the first-order network
    !! k_react(i,j) and yields reaction_yield(i,j): species i loses a fraction
    !! 1 - exp(-sum_j k_react(i,j)*dt), shared between products in proportion to the
    !! rates, and product j gains that loss times the yield.
    !! PFASProcessModule%pfas_transform_network.
    function contaminant_transformation(this, dt) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt
        type(Result) :: r
        if (allocated(this%c)) call pfas_transform_network(this%c, dt, this%k_react, this%reaction_yield)
        call update_m_dissolved(this)
    end function contaminant_transformation

    !> Leaching from a soil layer: removes the fraction water_fraction (percolating
    !! water / layer water, V_perc/V_w) of the DISSOLVED mass only and returns it in
    !! j_leached (AQ slot), to be passed to the layer below or, from the bottom layer,
    !! to the groundwater reservoir. Sorbed and interface PFAS stay behind, which is
    !! what retards PFAS leaching. Also sets j_to_groundwater.
    subroutine contaminant_leach(this, water_fraction, j_leached)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: water_fraction
        type(Contaminant), intent(inout) :: j_leached
        type(Result) :: rr
        integer :: f
        rr = j_leached%create()
        call copy_params(this, j_leached)
        if (.not. allocated(this%c)) return
        do f = 1, size(this%c,2)
            call pfas_remove_by_fraction(this%c(:,f,PFAS_AQ), water_fraction, j_leached%c(:,f,PFAS_AQ))
        end do
        this%j_to_groundwater = sum(j_leached%c(:,:,PFAS_AQ), dim=2)
        call update_m_dissolved(this)
        call update_m_dissolved(j_leached)
    end subroutine contaminant_leach

    !> Soil erosion: removes the fraction erosion_fraction (eroded soil mass / layer soil
    !! mass) of the soil-SORBED mass and returns it in j_eroded (SOL slot). Only PFAS
    !! sorbed to the soil leaves with the eroded particles. 
    subroutine contaminant_erosion_export(this, erosion_fraction, j_eroded)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: erosion_fraction
        type(Contaminant), intent(inout) :: j_eroded
        type(Result) :: rr
        integer :: f
        rr = j_eroded%create()
        call copy_params(this, j_eroded)
        if (.not. allocated(this%c)) return
        do f = 1, size(this%c,2)
            call pfas_remove_by_fraction(this%c(:,f,PFAS_SOL), erosion_fraction, j_eroded%c(:,f,PFAS_SOL))
        end do
        call update_m_dissolved(this)
        call update_m_dissolved(j_eroded)
    end subroutine contaminant_erosion_export

    !> Deposition from the water column to the bed
    subroutine contaminant_deposition(this, spm_deposit_fraction, dj_dep)
        !! SPM-sorbed PFAS that settles with the depositing SPM. PFAS is assumed to be
        !! spread over the SPM in proportion to its mass, so the fraction deposited is
        !! the fraction of the SPM mass deposited (sum over size classes of the deposit
        !! divided by the SPM mass it was taken from), not a mean settling velocity.
        !! The package keeps the mass in the PFAS_SPM slot so that (a) subtracting it
        !! from the water column hits the phase it came from and (b) the bed sediment's
        !! transferContaminant, which reads the PFAS_SPM slot of incoming deposits,
        !! receives it. The bed itself stores it as PFAS_SOL.
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: spm_deposit_fraction    !! Fraction of the SPM mass deposited [-]
        type(Contaminant), intent(inout) :: dj_dep
        type(Result) :: rr
        rr = dj_dep%create()
        call copy_params(this, dj_dep)
        if (.not. allocated(this%c)) return
        dj_dep%c(:,:,PFAS_SPM) = this%c(:,:,PFAS_SPM) * pfas_clamp01(spm_deposit_fraction)
        call update_m_dissolved(dj_dep)
    end subroutine contaminant_deposition

    !> Legacy deposition helper, not currently called: removes deposition_fraction of
    !! the SPM-sorbed mass into j_sediment (SOL slot). resuspension_fraction is ignored.
    !! Use deposition with the bed sediment's transferContaminant instead, which keeps
    !! the deposit in the SPM slot.
    subroutine contaminant_sediment_exchange(this, deposition_fraction, resuspension_fraction, j_sediment)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: deposition_fraction, resuspension_fraction
        type(Contaminant), intent(inout) :: j_sediment
        type(Result) :: rr
        integer :: f
        rr = j_sediment%create()
        call copy_params(this, j_sediment)
        if (.not. allocated(this%c)) return
        do f = 1, size(this%c,2)
            call pfas_remove_by_fraction(this%c(:,f,PFAS_SPM), deposition_fraction, j_sediment%c(:,f,PFAS_SOL))
        end do
        this%j_to_sediment = sum(j_sediment%c(:,:,PFAS_SOL), dim=2)
        call update_m_dissolved(this)
        call update_m_dissolved(j_sediment)
    end subroutine contaminant_sediment_exchange

    !> Outflow from a water body on one displacement: builds the package that leaves
    !! downstream, without changing this object. Dissolved and surface-film (AWI) PFAS
    !! leave with the water, fraction k_outflow (outflow volume / reach volume); SPM-sorbed
    !! PFAS leaves with the SPM, fraction sum(dj_spm_outflow)/sum(m_spm)
    subroutine contaminant_outflow_split(this, k_outflow, dj_spm_outflow, m_spm, dj_out)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: k_outflow, dj_spm_outflow(:), m_spm(:)
        type(Contaminant), intent(inout) :: dj_out
        type(Result) :: rr
        real(dp) :: spm_frac
        rr = dj_out%create()
        call copy_params(this, dj_out)
        if (.not. allocated(this%c)) return
        spm_frac = 0.0_dp
        if (sum(m_spm) > C%epsilon) spm_frac = pfas_clamp01(sum(dj_spm_outflow) / sum(m_spm))
        dj_out%c(:,:,PFAS_AQ)  = this%c(:,:,PFAS_AQ)  * pfas_clamp01(k_outflow)
        dj_out%c(:,:,PFAS_SPM) = this%c(:,:,PFAS_SPM) * spm_frac
        ! The surface film (air-water interface) is carried downstream with the water
        dj_out%c(:,:,PFAS_AWI) = this%c(:,:,PFAS_AWI) * pfas_clamp01(k_outflow)
        call update_m_dissolved(dj_out)
    end subroutine contaminant_outflow_split

    !> Foam export from the water surface: removes the fraction foam_fraction of the
    !! air-water interface (surface microlayer) mass and returns it in j_foam (FOAM
    !! slot). The Reactor computes the fraction as 1 - exp(-k_foam*dt); species
    !! differences in foam enrichment come from their interface loading (Kia).
    !! Also sets j_to_foam.
    subroutine contaminant_foam_exchange(this, foam_fraction, j_foam)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: foam_fraction
        type(Contaminant), intent(inout) :: j_foam
        type(Result) :: rr
        integer :: f
        rr = j_foam%create()
        call copy_params(this, j_foam)
        if (.not. allocated(this%c)) return
        do f = 1, size(this%c,2)
            call pfas_remove_by_fraction(this%c(:,f,PFAS_AWI), foam_fraction, j_foam%c(:,f,PFAS_FOAM))
        end do
        this%j_to_foam = sum(j_foam%c(:,:,PFAS_FOAM), dim=2)
        call update_m_dissolved(this)
        call update_m_dissolved(j_foam)
    end subroutine contaminant_foam_exchange

    !> Loss to the atmosphere from a water body: volatilisation of dissolved PFAS and
    !! sea-spray aerosol export of the surface-film (AWI) mass. Removed mass goes to the
    !! AIR slot of j_air; also sets j_to_atmosphere.
    subroutine contaminant_atmosphere_exchange(this, volatilisation_fraction, aerosol_fraction, j_air)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: volatilisation_fraction(:), aerosol_fraction
        type(Contaminant), intent(inout) :: j_air
        type(Result) :: rr
        integer :: i, f
        real(dp) :: afrac
        rr = j_air%create()
        call copy_params(this, j_air)
        if (.not. allocated(this%c)) return
        afrac = pfas_clamp01(aerosol_fraction)
        do f = 1, size(this%c,2)
            do i = 1, min(size(this%c,1), size(volatilisation_fraction))
                j_air%c(i,f,PFAS_AIR) = this%c(i,f,PFAS_AQ) * pfas_clamp01(volatilisation_fraction(i))
                this%c(i,f,PFAS_AQ) = this%c(i,f,PFAS_AQ) - j_air%c(i,f,PFAS_AIR)
            end do
            j_air%c(:,f,PFAS_AIR) = j_air%c(:,f,PFAS_AIR) + this%c(:,f,PFAS_AWI) * afrac
            this%c(:,f,PFAS_AWI) = this%c(:,f,PFAS_AWI) * (1.0_dp - afrac)
        end do
        this%j_to_atmosphere = sum(j_air%c(:,:,PFAS_AIR), dim=2)
        call update_m_dissolved(this)
        call update_m_dissolved(j_air)
    end subroutine contaminant_atmosphere_exchange

    !> Total concentration: the mass summed over ALL species, forms and phases divided
    !! by volume [kg m-3], returned as Result0D data; zero for a dry water body.
    function contaminant_getConcentration(this, volume) result(r)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: volume
        type(Result0D) :: r
        real(dp) :: C_total
        ! A waterbody with no water has zero concentration
        if (volume > C%epsilon .and. allocated(this%c)) then
            C_total = sum(this%c) / volume
        else
            C_total = 0.0_dp
        end if
        allocate(r%data, source=C_total)
    end function contaminant_getConcentration

    !> The dissolved (AQ) phase, by species and form. Despite the result name this is
    !! whatever c holds (mass, or concentration for a divideCheckZero result).
    function contaminant_get_free(this) result(C_free)
        class(Contaminant), intent(in) :: this
        real(dp) :: C_free(C%contaminantDim(1), C%contaminantDim(2))
        C_free = 0.0_dp
        if (allocated(this%c)) C_free = this%c(:,:,PFAS_AQ)
    end function contaminant_get_free

    !> The soil/bed-sorbed (SOL) phase, by species and form. Note: in the water column
    !! particle-bound PFAS is in the SPM phase, so use get_phase(PFAS_SPM) there.
    function contaminant_get_attached(this) result(C_attached)
        class(Contaminant), intent(in) :: this
        real(dp) :: C_attached(C%contaminantDim(1), C%contaminantDim(2))
        C_attached = 0.0_dp
        if (allocated(this%c)) C_attached = this%c(:,:,PFAS_SOL)
    end function contaminant_get_attached

    !> Any one phase (a PFAS_* index), by species and form; zeros for an invalid index.
    function contaminant_get_phase(this, phase) result(m)
        class(Contaminant), intent(in) :: this
        integer, intent(in) :: phase
        real(dp) :: m(C%contaminantDim(1), C%contaminantDim(2))
        m = 0.0_dp
        if (allocated(this%c)) then
            if (phase >= 1 .and. phase <= size(this%c,3)) m = this%c(:,:,phase)
        end if
    end function contaminant_get_phase

    !> Total mass over all species, forms and phases [kg] (zero if unallocated).
    function contaminant_mass_total(this) result(m)
        class(Contaminant), intent(in) :: this
        real(dp) :: m
        if (allocated(this%c)) then
            m = sum(this%c)
        else
            m = 0.0_dp
        end if
    end function contaminant_mass_total

    !> Legacy heteroaggregation (nanoparticle-SPM collision): no-op for PFAS, whose
    !! association with SPM is handled by partition_phases.
    function contaminant_noop_result(this, dt, T_water, C_spm, W_settle_spm, C_spm_particle) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt, T_water, C_spm(:), W_settle_spm(:), C_spm_particle(:)
        type(Result) :: r
    end function contaminant_noop_result

    !> Legacy nanoparticle dissolution: no-op for PFAS (it is already a solute).
    function contaminant_noop_result_dt(this, dt) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt
        type(Result) :: r
    end function contaminant_noop_result_dt

    !> Legacy nanoparticle attachment to soil grains: no-op for PFAS, whose sorption
    !! to soil is handled by partition_phases.
    function contaminant_noop_result_att(this, dt, k_att, alpha_att) result(r)
        class(Contaminant), intent(inout) :: this
        real(dp), intent(in) :: dt, k_att(:), alpha_att
        type(Result) :: r
    end function contaminant_noop_result_att

    !> Legacy nanoparticle-SPM collision rate: always zero for PFAS.
    function contaminant_calculateCollisionRate(this, T_water, G, W_settle_spm) result(k_coll)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: T_water, G, W_settle_spm(:)
        real(dp) :: k_coll(C%contaminantDim(1), C%nSizeClassesSpm)
        k_coll = 0.0_dp
    end function contaminant_calculateCollisionRate

    !> Legacy particle number concentration: always zero for PFAS.
    function contaminant_calculateParticleConcentration(this, C_mass, rho_particle, d) result(C_particle)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: C_mass, rho_particle, d
        real(dp) :: C_particle
        C_particle = 0.0_dp
    end function contaminant_calculateParticleConcentration

    !> Legacy nanoparticle attachment rate: always zero for PFAS.
    function contaminant_calculateAttachmentRate(this, T_water, porosity, d_grain, velocity) result(k_att)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: T_water, porosity, d_grain
        real(dp), intent(in), optional :: velocity
        real(dp) :: k_att(C%contaminantDim(1))
        k_att = 0.0_dp
    end function contaminant_calculateAttachmentRate

    !> Legacy free-particle settling velocity: always zero for PFAS, which settles
    !! only when sorbed to SPM.
    function contaminant_calculateSettlingVelocity(this, d, rho_particle, T_water) result(W_settle)
        class(Contaminant), intent(in) :: this
        real(dp), intent(in) :: d, rho_particle, T_water
        real(dp) :: W_settle
        W_settle = 0.0_dp
    end function contaminant_calculateSettlingVelocity
end module
