!> The DataInputModule contains the Database type, which is responsible for
!! data input to the model, as well as a variable (DATASET) of type(Database), which
!! can be imported into other modules, thus making the data parsed by the Database type
!! accessible throughout the model.
module DataInputModule
    use mo_netcdf
    use DefaultsModule
    use ConstantsDefaultsModule
    use GlobalsModule, only: dp, C
    use PFASConstantsModule, only: PFAS_AQ, PFAS_SOL, PFAS_SPM, PFAS_AWI, PFAS_FOAM, PFAS_AIR, PFAS_NPHASES
    use ResultModule, only: Result
    use ErrorInstanceModule, only: ErrorInstance
    use LoggerModule, only: LOGR
    use UtilModule
    implicit none
    ! P-FASE UPDATE:
    ! DataInput now treats Contaminant%c(:,:,phase) as PFAS phase pools.
    ! Legacy input variables are still accepted, but are mapped into PFAS_AQ by default.
    ! New PFAS-ready NetCDF/constants files should provide species x form x phase arrays.

    !> The Database type is responsible for data input to the model. It parses data
    !! from the NetCDF and constant namelist files.
    type, public :: Database
        type(NcDataset)     :: nc                               ! The NetCDF dataset

        ! CONSTANTS
        ! ---------
        ! Contaminant
        real(dp)          :: contaminantDensity                    ! Density of the contaminant [kg/m3]
        real(dp), allocatable :: contaminantSizeClasses(:)             ! Legacy field; for P-FASE use as PFAS species support array if still present
        real, allocatable :: defaultDistributionContaminant(:) ! Default distribution to split contaminant across size classes
        real, allocatable :: defaultContaminantFormDistribution(:)
        integer           :: nContaminantSizeClasses               ! Number of PFAS species / legacy contaminant size classes
        ! Sediment
        real, allocatable :: defaultSpmSizeDistribution(:)      ! Default distribution to split SPM across size classes
        real(dp), allocatable :: spmDensityBySizeClass(:)           ! Density of sediment in each size class [kg/m3]
        real(dp), allocatable :: spmSizeClasses(:)                  ! Diameter of each SPM size class [m]
        real, allocatable :: defaultMatrixEmbeddedDistributionToSpm(:)  ! Default distribution to proportion matrix-embedded releases to SPM size classes
        integer :: nSizeClassesSpm                              ! Number of SPM size classes
        real(dp) :: sedimentEnrichment_k                        ! Clay enrichment scaling factor
        real(dp) :: sedimentEnrichment_a                        ! Clay enrichment skew factor
        ! Soil
        real(dp) :: soilDarcyVelocity                           ! Darcy velocity in soil [m/s]
        real(dp) :: soilDefaultPorosity                         ! Default porosity [-]
        real(dp) :: soilHamakerConstant                         ! Hamaker constant for soil [J]
        real(dp) :: soilParticleDensity                         ! Particle density of soil [kg m-3]
        real(dp) :: soilErosivity_a1                            ! Erosivity a1 parameter [-]
        real(dp) :: soilErosivity_a2                            ! Erosivity a2 parameter [-]
        real(dp) :: soilErosivity_a3                            ! Erosivity a3 parameter [-]
        real(dp) :: soilErosivity_b                             ! Erosivity b parameter [-]
        real(dp) :: soilAttachmentEfficiencyConstant            ! Attachment efficiency to soil matrix [-]
        real(dp) :: sedimentTransport_aConstant                 ! Sediment transport capacity a parameter (scaling factor) [kg/m2/km2]
        real(dp) :: sedimentTransport_bConstant                 ! Sediment transport capacity b parameter (overland flow threshold) [m2/s]
        real(dp) :: sedimentTransport_cConstant                 ! Sediment transport capacity c parameter (non-linear coefficient) [-]
        ! Earthworms
        integer :: earthwormDensityArable                       ! Earthworm density arable [individuals/m2]
        integer :: earthwormDensityConiferous                   ! Earthworm density coniferous [individuals/m2]
        integer :: earthwormDensityDeciduous                    ! Earthworm density deciduous [individuals/m2]
        integer :: earthwormDensityGrassland                    ! Earthworm density grassland [individuals/m2]
        integer :: earthwormDensityHeathland                    ! Earthworm density heathland [individuals/m2]
        integer :: earthwormDensityUrbanCapped                  ! Earthworm density urban_capped [individuals/m2]
        integer :: earthwormDensityUrbanGardens                 ! Earthworm density urban_gardens [individuals/m2]
        integer :: earthwormDensityUrbanParks                   ! Earthworm density urban_parks [individuals/m2]
        real, allocatable :: earthwormVerticalDistribution(:)   ! Vertical distribution of earthworms
        ! Biota
        character(len=100), allocatable :: biotaName(:)
        character(len=100), allocatable :: biotaCompartment(:)
        real(dp), allocatable :: biotaInitial_C_org(:)
        real(dp), allocatable :: biota_k_growth(:)
        real(dp), allocatable :: biota_k_death(:)                ! Loss/death dilution [/day]; optional k_death in &biota
        integer, allocatable :: biotaPrey(:)                     ! Index of the biota this one eats (same compartment), 0 = none
        real(dp), allocatable :: biota_k_diet(:)                 ! Dietary uptake = feeding rate x assimilation efficiency [/day]
        real(dp), allocatable :: pfasBiotaTransformationRate(:,:)  ! In-organism precursor transformation, species i -> j [/day]
        real(dp), allocatable :: biota_k_uptake_contaminant(:,:) ! Uptake rates for contaminant forms [nBiota, nContaminantForms]
        real(dp), allocatable :: biota_k_elim_contaminant(:,:)   ! Elimination rates for contaminant forms [nBiota, nContaminantForms]
        real(dp), allocatable :: biota_k_uptake_dissolved(:)     ! Uptake rate for dissolved contaminant
        real(dp), allocatable :: biota_k_elim_dissolved(:)       ! Elimination rate for dissolved contaminant
        real, allocatable :: biotaStoredFraction(:)
        character(len=17), allocatable :: biotaUptakeFromForm(:)
        integer, allocatable :: biotaHarvestInMonth(:)
        logical :: hasBiota = .false.
        integer :: nBiota = 0
        ! Water
        real     :: riverMeanderingFactor           ! Meandering factor for rivers (not estuaries) [-]
        real(dp) :: waterResuspensionAlpha          ! Resuspension parameter alpha
        real(dp) :: waterResuspensionBeta           ! Resuspension parameter beta
        real(dp) :: waterResuspensionAlphaEstuary   ! Resuspension parameter alpha for estuary
        real(dp) :: waterResuspensionBetaEstuary    ! Resuspension parameter beta for estuary
        real(dp) :: depositionAlphaConstant         ! Deposition parameter alpha - constant if spatial variable not supplied
        real(dp) :: depositionBetaConstant          ! Deposition parameter beta - constant if spatial variable not supplied
        real(dp) :: bankErosionAlphaConstant        ! Bank erosion parameter alpha - constant if spatial variable not supplied
        real(dp) :: bankErosionBetaConstant         ! Bank erosion parameter beta - constant if spatial variable not supplied
        real(dp) :: contaminant_k_diss_pristine     ! Dissolution rate constant for pristine contaminant [/s]
        real(dp) :: contaminant_k_diss_transformed  ! Dissolution rate constant for transformed contaminant [/s]
        real(dp) :: contaminant_k_transform_pristine ! Transformation rate constant for pristine contaminant [/s]
        real(dp) :: shearRate                       ! Shear rate [/s]
        real(dp) :: waterTemperature(366)           ! Temporally varying water temperature [deg C]
        real(dp) :: riverAttachmentEfficiency
        ! Estuary
        real(dp) :: estuaryAttachmentEfficiency     ! Attachment efficiency for NM to SPM in estuaries [-]
        real :: estuaryTidalM2                      ! Estuary tidal harmonics parameter M2
        real :: estuaryTidalS2                      ! Estuary tidal harmonics parameter S2
        real :: estuaryMeanDepthExpA                ! Estuary mean depth exponential parameter A
        real :: estuaryMeanDepthExpB                ! Estuary mean depth exponential parameter B
        real :: estuaryWidthExpA                    ! Estuary width exponential parameter A
        real :: estuaryWidthExpB                    ! Estuary width exponential parameter B
        real :: estuaryMeanderingFactor             ! Estuary meandering factor, used to calculate distance to mouth [-]
        real :: estuaryMouthCoords(2)               ! Coordinates of the estuary mouth, used to calculate distance to mouth
        ! Sediment
        real, allocatable :: sedimentPorosity(:)    ! Porosity of the bed sediment layers [-]
        real(dp), allocatable :: sedimentInitialMass(:) ! Initial mass of each sediment size class [kg/m2]
        real, allocatable :: sedimentFractionalComposition(:) ! Distribution of sediment amongst fractional compositions [-]

        ! NETCDF VARIABLES
        ! ----------------
        ! Grid and coordinate variables
        integer, allocatable :: gridShape(:)        ! Number of grid cells along each grid axis [-]
        real, allocatable :: gridRes(:)             ! Resolution of grid cells [m]
        real, allocatable :: gridBounds(:)          ! Bounding box of grid, indexed as left, bottom, right, top [m]
        logical, allocatable :: gridMask(:,:)       ! Logical mask for extent of grid [-]
        character(len=2000) :: crsWKT               ! The well-known text of the CRS, to be saved to output data
        integer :: epsgCode                         ! The EPSG code of the CRS
        real, allocatable :: x(:)                   ! Centre of cell [m]
        real, allocatable :: x_l(:)                 ! Left side of cell [m]
        real, allocatable :: y(:)                   ! Centre of cell [m]
        real, allocatable :: y_u(:)                 ! Upper side of cell [m]
        integer, allocatable :: t(:)                ! Seconds since start date [s]
        integer :: nTimesteps                       ! Number of timesteps for the model run [-]
        integer, allocatable :: dem(:,:)            ! Digital elevation model [dm asl]
        ! Simulation mask
        logical, allocatable :: simulationMask(:,:) ! Boolean mask to apply to the simulation
        integer :: nNonMaskedCells                  ! Number of non-masked grid cells
        ! Routing variables
        integer, allocatable :: outflow(:,:,:)
        integer, allocatable :: inflows(:,:,:,:)
        logical, allocatable :: isHeadwater(:,:)
        integer, allocatable :: nWaterbodies(:,:)
        logical, allocatable :: isEstuary(:,:)
        integer :: maxNWaterbodies                  ! Maximum number of waterbodies per cell in the model domain
        ! Lakes and reservoirs. Lake segments are contained within grid cells (max size = cell
        ! area); waterbodies spanning multiple cells are split into one segment per cell.
        logical, allocatable :: isLake(:,:)                 ! Does this cell contain a lake segment?
        logical, allocatable :: isReservoir(:,:)            ! Does this cell contain a reservoir segment?
        real, allocatable :: lakeDepth(:,:)                 ! Lake segment depth at outlet/spillway level [m]
        real, allocatable :: lakeAreaFraction(:,:)          ! Fraction of cell area covered by the lake segment [-]
        real, allocatable :: reservoirReleaseRate(:,:)      ! Controlled spillway release rate [m3/s]
        integer, allocatable :: reservoirSpillwayType(:,:)  ! 0 = uncontrolled, 1 = controlled, 2 = combination
        integer, allocatable :: lakeId(:,:)                 ! Lake identifier; segments of the same (multi-cell) lake share an id, 0 = no lake
        logical :: hasLakeRouting = .false.                 ! Was a 'lake_outflow' routing network provided and merged in?
        ! Spatiotemporal variables
        real, allocatable :: runoff(:,:,:)
        real, allocatable :: quickflow(:,:,:)
        real, allocatable :: precip(:,:,:)
        real, allocatable :: evap(:,:,:)
        real, allocatable :: windSpeedTimeSeries(:,:,:)     ! Wind speed at 10 m [m/s], used for lake stratification/resuspension
        real, allocatable :: solarRadiation(:,:,:)          ! Downward solar radiation [W/m2], used for lake stratification
        ! Spatial variables
        real, allocatable :: soilBulkDensity(:,:)
        real, allocatable :: soilWaterContentFieldCapacity(:,:)
        real, allocatable :: soilWaterContentSaturation(:,:)
        real, allocatable :: soilHydraulicConductivity(:,:)
        real, allocatable :: soilTextureClayContent(:,:)
        real, allocatable :: soilTextureSandContent(:,:)
        real, allocatable :: soilTextureSiltContent(:,:)
        real, allocatable :: soilTextureCoarseFragContent(:,:)
        real, allocatable :: soilAttachmentEfficiency(:,:)
        real, allocatable :: soilAttachmentRate(:,:)
        real(dp), allocatable :: soilUsleCFactor(:,:)
        real(dp), allocatable :: soilUslePFactor(:,:)
        real(dp), allocatable :: soilUsleLSFactor(:,:)
        real(dp), allocatable :: resuspensionAlpha(:,:)
        real(dp), allocatable :: resuspensionBeta(:,:)
        real(dp), allocatable :: depositionAlpha(:,:)
        real(dp), allocatable :: depositionBeta(:,:)
        real(dp), allocatable :: bankErosionAlpha(:,:)
        real(dp), allocatable :: bankErosionBeta(:,:)
        real(dp), allocatable :: sedimentTransport_a(:,:)                       ! Sediment transport capacity a parameter (scaling factor) [kg/m2/km2]
        real(dp), allocatable :: sedimentTransport_b(:,:)                       ! Sediment transport capacity b parameter (overland flow threshold) [m2/s]
        real(dp), allocatable :: sedimentTransport_c(:,:)                       ! Sediment transport capacity c parameter (non-linear coefficient) [-]
        ! Initial concentrations
        real(dp), allocatable :: initialContaminantConcsSoil(:,:,:,:,:)
        real(dp), allocatable :: initialContaminantConcsWater(:,:,:,:,:)
        real(dp), allocatable :: initialContaminantConcsSediment(:,:,:,:,:)
        real(dp), allocatable :: initialDissolvedConcsSoil(:,:)
        real(dp), allocatable :: initialDissolvedConcsWater(:,:)
        real(dp), allocatable :: initialDissolvedConcsSediment(:,:)
        ! Emissions - areal
        real(dp), allocatable :: emissionsArealSoilContaminant(:,:,:,:,:)
        real(dp), allocatable :: emissionsArealWaterContaminant(:,:,:,:,:)
        real(dp), allocatable :: emissionsArealSoilDissolvedContaminant(:,:)
        real(dp), allocatable :: emissionsArealWaterDissolvedContaminant(:,:)
        ! Emissions - atmospheric depo
        real(dp), allocatable :: emissionsAtmosphericDryDepoContaminant(:,:,:,:,:,:)
        real(dp), allocatable :: emissionsAtmosphericWetDepoContaminant(:,:,:,:,:,:)
        real(dp), allocatable :: emissionsAtmosphericDryDepoDissolvedContaminant(:,:,:)
        real(dp), allocatable :: emissionsAtmosphericWetDepoDissolvedContaminant(:,:,:)
        ! Emissions - point
        real(dp), allocatable :: emissionsPointWaterCoords(:,:,:,:)
        real(dp), allocatable :: emissionsPointWaterContaminant(:,:,:,:,:,:,:)  ! (x, y, t, p, species, form, phase)
        real(dp), allocatable :: emissionsPointWaterDissolvedContaminant(:,:,:) 
        integer, allocatable :: nPointSources(:,:)
        integer :: maxPointSources                                              ! Maximum number of point sources in a cell in the whole environment
        ! Spatial 1D variables
        real, allocatable :: landUse(:,:,:)

        ! P-FASE PFAS species properties, read from the optional &pfas group of the
        ! constants file (all species-resolved, sized C%contaminantDim(1)). Partition
        ! coefficients are stored in SI [m3 kg-1]; the namelist takes L kg-1.
        character(len=64), allocatable :: pfasSpeciesNames(:)
        real(dp), allocatable :: pfasMolecularWeight(:), pfasCharge(:), pfasPka(:)
        real(dp) :: pfasPkaUnset = -900.0_dp            ! pKa at or below this = not given (PFSA pKa can be < 0)
        real(dp), allocatable :: pfasKdSoil(:)          ! Constant soil/bed Kd, fallback when no map and no Koc [m3 kg-1]
        real(dp), allocatable :: pfasKdSpm(:)           ! Constant SPM Kd, fallback when no map and no Koc [m3 kg-1]
        real(dp), allocatable :: pfasKoc(:)             ! Organic-carbon normalised Kd; > 0 enables Kd = Koc*foc [m3 kg-1 OC]
        real(dp), allocatable :: pfasKdPhSlope(:)       ! d(log10 Kd)/d(pH) for the Koc*foc route [-]
        real(dp) :: pfasKdPhReference = 7.0_dp          ! pH at which Koc applies unmodified [-]
        real(dp), allocatable :: pfasKia(:)             ! Air-water interfacial adsorption coefficient [m]
        real(dp), allocatable :: pfasKMtSoil(:)         ! Soil/bed sorption mass-transfer rate, 0 = equilibrium [s-1]
        real(dp), allocatable :: pfasKMtSpm(:)          ! SPM sorption mass-transfer rate, 0 = equilibrium [s-1]
        real(dp), allocatable :: pfasKMtAwi(:)          ! AWI adsorption mass-transfer rate, 0 = equilibrium [s-1]
        real(dp), allocatable :: pfasReactionRate(:,:), pfasReactionYield(:,:)
        real(dp), allocatable :: pfasKdSediment(:)      ! Constant bed-sediment Kd; negative = use the SPM Kd [m3 kg-1]
        real(dp), allocatable :: pfasReactionRateSediment(:,:)  ! Bed (anoxic) transformation rates; unset = pfasReactionRate [s-1]
        real(dp), allocatable :: pfasDiffusivity(:)     ! Free-water molecular diffusivity [m2 s-1]
        real(dp) :: sedimentWaterExchangeVelocity = -1.0_dp  ! Sediment-water mass-transfer velocity; negative = from diffusivity [m s-1]
        real(dp), allocatable :: pfasKdSalinityCoefficient(:)   ! d(log10 Kd)/d(salinity), salting-out [psu-1]
        real(dp), allocatable :: pfasKiaSalinityCoefficient(:)  ! d(log10 Kia)/d(salinity) [psu-1]
        real(dp) :: estuarySalinityConstant = 0.0_dp        ! Salinity of estuary cells where 'water_salinity' has no value [psu]
        real(dp) :: runoffMixingDepth = 0.01_dp             ! Soil depth that exchanges with overland flow; <= 0 = no dissolved runoff [m]
        real(dp) :: groundwaterResidenceTimeConstant = -1.0_dp  ! Groundwater reservoir residence time; <= 0 = no reservoir [days]
        ! Soil and surface water properties driving the PFAS partitioning. The maps are
        ! (x,y) [or (x,y,species)] and optional; cells with no map value (missing, fill
        ! or negative) fall back to the model-computed value, then to the constants.
        real(dp), allocatable :: pfasKdSoilMap(:,:,:)       ! 'pfas_kd_soil' (x,y,species) [m3 kg-1]
        real(dp), allocatable :: pfasKdSpmMap(:,:,:)        ! 'pfas_kd_spm' (x,y,species) [m3 kg-1]
        real(dp), allocatable :: pfasKdSedimentMap(:,:,:)   ! 'pfas_kd_sediment' (x,y,species) [m3 kg-1]
        real(dp), allocatable :: sedimentOrganicCarbon(:,:) ! 'sediment_organic_carbon' (x,y) [kg OC kg-1]
        real(dp), allocatable :: waterSalinity(:,:)         ! 'water_salinity' (x,y) [psu]
        real(dp), allocatable :: groundwaterResidenceTime(:,:)  ! 'groundwater_residence_time' (x,y) [days]
        real(dp), allocatable :: soilOrganicCarbon(:,:)     ! 'soil_organic_carbon' (x,y) [kg OC kg-1]
        real(dp), allocatable :: soilPh(:,:)                ! 'soil_ph' (x,y) [-]
        real(dp), allocatable :: spmOrganicCarbon(:,:)      ! 'spm_organic_carbon' (x,y) [kg OC kg-1]
        real(dp), allocatable :: waterPh(:,:)               ! 'water_ph' (x,y) [-]
        real(dp), allocatable :: soilAwiSpecificArea(:,:)   ! 'soil_awi_specific_area' (x,y), at zero saturation [m2 m-3]
        real(dp) :: soilOrganicCarbonConstant = -1.0_dp     ! Negative = not supplied
        real(dp) :: soilPhConstant = -1.0_dp
        real(dp) :: spmOrganicCarbonConstant = -1.0_dp
        real(dp) :: waterPhConstant = -1.0_dp
        real(dp) :: sedimentOrganicCarbonConstant = -1.0_dp
        real(dp) :: soilAwiSpecificAreaConstant = 0.0_dp    ! Used only when soilAwiGeometricFactor <= 0 [m2 m-3]
        real(dp) :: soilAwiGeometricFactor = 1.0_dp         ! Multiplier on the smooth-sphere solid area 6(1-n)/d [-]
        ! Scalar forcing defaults used by water/estuary/biota modules when species-resolved
        ! arrays are absent or when a compact parameterisation is preferred.
        real(dp) :: pfasFoamCoefficient = 0.0_dp
        real(dp) :: pfasVolatilisationRateScalar = 0.0_dp
        real(dp) :: pfasSeaSprayAerosolRate = 0.0_dp
        real(dp) :: windSpeed = 0.0_dp

        ! Lake/reservoir process constants. Defaults here; can be overridden via the
        ! optional &lake namelist group in the constants file.
        real(dp) :: lakeDefaultDepth = 5.0_dp                       ! Lake segment depth when 'lake_depth' isn't in the input data [m]
        real(dp) :: lakeDefaultAreaFraction = 1.0_dp                ! Fraction of cell covered when 'lake_area_fraction' isn't in the input data [-]
        real(dp) :: lakeStratificationRadiationThreshold = 120.0_dp ! Solar radiation above which stratification can form [W/m2]
        real(dp) :: lakeStratificationWindThreshold = 8.0_dp        ! Wind speed above which stratification breaks down [m/s]
        real(dp) :: lakeStratificationTemperatureThreshold = 12.0_dp ! Water temperature proxy threshold when no solar radiation data [deg C]
        real(dp) :: lakeEpilimnionBaseDepth = 2.0_dp                ! Epilimnion depth at zero wind [m]
        real(dp) :: lakeEpilimnionWindCoefficient = 0.5_dp          ! Epilimnion deepening per unit wind speed [m/(m/s)]
        real(dp) :: lakeEpilimnionMinDepth = 1.0_dp                 ! Minimum epilimnion depth [m]
        real(dp) :: lakeInterlayerMixingVelocity = 0.05_dp          ! Diffusive exchange velocity across the thermocline [m/day]
        real(dp) :: lakeWindResuspensionAttenuationDepth = 3.0_dp   ! e-folding depth for wind energy reaching the lake bed [m]
        real(dp) :: reservoirDefaultReleaseRate = 0.0_dp            ! Controlled spillway release when 'reservoir_release_rate' isn't in the input data [m3/s]
contains
        procedure, public   :: init => initDatabase
        procedure, public   :: update => updateDatabase
        procedure, public   :: readBatchVariables => readBatchVariablesDatabase
        procedure, private  :: parseHydrology => parseHydrologyDatabase
        procedure, private  :: mergeLakeRouting => mergeLakeRoutingDatabase
        procedure, private  :: parseConstants => parseConstantsDatabase
        procedure, private  :: mask => maskDatabase
        procedure, public   :: inModelDomain => inModelDomainDatabase
        procedure, private  :: calculateNPointSources => calculateNPointSourcesDatabase
        procedure, private  :: calculateMeanderingFactorFromCellSize => calculateMeanderingFactorFromCellSizeDatabase
        procedure, public   :: coordsToCellIndex => coordsToCellIndexDatabase
        procedure, public   :: coordsToFractionalCellIndex => coordsToFractionalCellIndexDatabase
        procedure, private  :: calculateWaterTemperatureTimeSeries => calculateWaterTemperatureTimeSeriesWaterBody
        ! Auditing
        procedure, private :: audit => auditDatabase
    end type

    type(Database) :: DATASET

  contains

    !> Initialise the database, by opening the NetCDF and constants files, and parsing
    !! their deposition
    subroutine initDatabase(me, inputFile, constantsFile)
        class(Database)     :: me
        type(NcDataset)     :: nc_simulationMask
        type(NcVariable)    :: var
        character(len=*)    :: inputFile, constantsFile
        type(Result)        :: rslt
        ! temps returned by mo_netcdf in Fortran order (reversed NetCDF dims)
        integer, allocatable :: outflow_dxy(:,:,:)
        integer, allocatable :: inflows_dwxy(:,:,:,:)
        integer, allocatable :: isHeadwaterInt_xy(:,:), isEstuaryInt_xy(:,:)
        integer, allocatable :: lakeInt_xy(:,:), lakeOutflow_dxy(:,:,:)
        logical :: hadIsLakeVar
        integer, allocatable :: nWaterbodies_xy(:,:)
        integer, allocatable :: simulationMask(:,:)
        integer :: nx, ny

        ! Open the dataset and parse constants
        me%nc = NcDataset(inputFile, 'r')
        call me%parseConstants(constantsFile)

        ! GRID / COORDS
        var = me%nc%getVariable('grid_shape');  call var%getData(me%gridShape)
        var = me%nc%getVariable('grid_res');    call var%getData(me%gridRes)
        var = me%nc%getVariable('grid_bounds'); call var%getData(me%gridBounds)
        var = me%nc%getVariable('x');           call var%getData(me%x)
        var = me%nc%getVariable('y');           call var%getData(me%y)
        allocate(me%x_l(size(me%x))); me%x_l = me%x - 0.5 * me%gridRes(1)
        allocate(me%y_u(size(me%y))); me%y_u = me%y + 0.5 * me%gridRes(2)
        var = me%nc%getVariable('crs'); call var%getAttribute('crs_wkt', me%crsWKT)

        nx = me%gridShape(1)
        ny = me%gridShape(2)

        ! ROUTING (getData already reversed dims to Fortran order)
        ! outflow: file (y,x,d) -> returned (d,x,y) => model (d,x,y)
        var = me%nc%getVariable('outflow');  call var%getData(outflow_dxy)
        if (allocated(me%outflow)) deallocate(me%outflow)
        allocate(me%outflow( size(outflow_dxy,1), size(outflow_dxy,2), size(outflow_dxy,3) ))
        me%outflow = outflow_dxy
        deallocate(outflow_dxy)

        ! inflows: file (y,x,w,d) -> returned (d,w,x,y) => model (d,w,x,y)
        var = me%nc%getVariable('inflows');  call var%getData(inflows_dwxy)
        if (allocated(me%inflows)) deallocate(me%inflows)
        allocate(me%inflows( size(inflows_dwxy,1), size(inflows_dwxy,2), &
                            size(inflows_dwxy,3), size(inflows_dwxy,4) ))
        me%inflows = inflows_dwxy
        deallocate(inflows_dwxy)

        ! headwater / n_waterbodies / estuary: file (y,x) -> returned (x,y) => model (x,y)
        var = me%nc%getVariable('is_headwater');  call var%getData(isHeadwaterInt_xy)
        me%isHeadwater = ulgcl(isHeadwaterInt_xy)
        deallocate(isHeadwaterInt_xy)

        var = me%nc%getVariable('n_waterbodies'); call var%getData(nWaterbodies_xy)
        me%nWaterbodies = nWaterbodies_xy
        deallocate(nWaterbodies_xy)
        me%maxNWaterbodies = maxval(me%nWaterbodies)

        if (C%includeEstuary) then
            var = me%nc%getVariable('is_estuary'); call var%getData(isEstuaryInt_xy)
            me%isEstuary = ulgcl(isEstuaryInt_xy)
            deallocate(isEstuaryInt_xy)
        else
            allocate(me%isEstuary(nx, ny)); me%isEstuary = .false.
        end if

        ! Lake and reservoir segment masks. Both are optional; cells not flagged
        ! (or with the masks absent, or include_lakes = .false.) are simulated as rivers.
        if (allocated(me%isLake)) deallocate(me%isLake)
        allocate(me%isLake(nx, ny)); me%isLake = .false.
        if (allocated(me%isReservoir)) deallocate(me%isReservoir)
        allocate(me%isReservoir(nx, ny)); me%isReservoir = .false.
        if (allocated(me%lakeId)) deallocate(me%lakeId)
        allocate(me%lakeId(nx, ny)); me%lakeId = 0
        me%hasLakeRouting = .false.
        if (C%includeLakes) then
            hadIsLakeVar = me%nc%hasVariable('is_lake')
            if (hadIsLakeVar) then
                var = me%nc%getVariable('is_lake'); call var%getData(lakeInt_xy)
                me%isLake = ulgcl(lakeInt_xy)
                deallocate(lakeInt_xy)
                call LOGR%add("DataInput: Read 'is_lake' mask (" // trim(str(count(me%isLake))) // " lake cells)")
            end if
            if (me%nc%hasVariable('is_reservoir')) then
                var = me%nc%getVariable('is_reservoir'); call var%getData(lakeInt_xy)
                me%isReservoir = ulgcl(lakeInt_xy)
                deallocate(lakeInt_xy)
                call LOGR%add("DataInput: Read 'is_reservoir' mask (" &
                    // trim(str(count(me%isReservoir))) // " reservoir cells)")
            end if

            ! Lake identifiers: segments of the same (multi-cell) lake or reservoir
            ! share an id. If no 'is_lake' mask was given, lake cells are derived
            ! from the ids instead, so 'lake_id' alone is a sufficient mask.
            if (me%nc%hasVariable('lake_id')) then
                var = me%nc%getVariable('lake_id'); call var%getData(lakeInt_xy)
                where (lakeInt_xy > 0) me%lakeId = lakeInt_xy
                deallocate(lakeInt_xy)
                if (.not. hadIsLakeVar) then
                    me%isLake = me%lakeId > 0
                    call LOGR%add("DataInput: Derived lake cells from 'lake_id' (" &
                        // trim(str(count(me%isLake))) // " cells)")
                end if
            end if
            ! A reservoir is a lake segment with a dam; don't double count
            where (me%isReservoir) me%isLake = .false.

            ! Lake routing network: an optional per-cell outflow override for lake and
            ! reservoir segments, connecting the segments of a multi-cell lake to each
            ! other and to the river network at the lake outlet. 
            if (me%nc%hasVariable('lake_outflow')) then
                ! Same layout as 'outflow': file (y,x,d) -> returned (d,x,y)
                var = me%nc%getVariable('lake_outflow'); call var%getData(lakeOutflow_dxy)
                call me%mergeLakeRouting(lakeOutflow_dxy)
                deallocate(lakeOutflow_dxy)
            end if
        end if

        ! Optional lake/reservoir spatial parameters, (x,y) indexed like the
        ! resuspension/sediment transport parameters (no transpose)
        if (me%nc%hasVariable('lake_depth')) then
            var = me%nc%getVariable('lake_depth')
            if (allocated(me%lakeDepth)) deallocate(me%lakeDepth)
            call var%getData(me%lakeDepth)
        end if
        if (me%nc%hasVariable('lake_area_fraction')) then
            var = me%nc%getVariable('lake_area_fraction')
            if (allocated(me%lakeAreaFraction)) deallocate(me%lakeAreaFraction)
            call var%getData(me%lakeAreaFraction)
        end if
        if (me%nc%hasVariable('reservoir_release_rate')) then
            var = me%nc%getVariable('reservoir_release_rate')
            if (allocated(me%reservoirReleaseRate)) deallocate(me%reservoirReleaseRate)
            call var%getData(me%reservoirReleaseRate)
        end if
        if (me%nc%hasVariable('reservoir_spillway_type')) then
            var = me%nc%getVariable('reservoir_spillway_type')
            if (allocated(me%reservoirSpillwayType)) deallocate(me%reservoirSpillwayType)
            call var%getData(me%reservoirSpillwayType)
        end if

        ! Digital elevation model [dm asl]
        ! Reach%parseInputData uses this to derive the stream slope from the elevation drop
        ! along the reach; if it isn't allocated, every reach silently falls back to
        ! defaultSlope
        if (me%nc%hasVariable('dem')) then
            var = me%nc%getVariable('dem'); call var%getData(me%dem)
        else
            call LOGR%toFile(errors=[ &
                ErrorInstance(message='Digital elevation model (dem) not found in input file. ' // &
                    'Default slope of ' // trim(str(defaultSlope)) // ' m/m will be used.', iscritical=.false.) &
            ])
        end if

        ! Grid mask from nWaterbodies
        allocate(me%gridMask(nx, ny))
        me%gridMask = me%mask(me%nWaterbodies)

        ! Derive meandering factors from grid size if not set in constants
        if (isZero(me%riverMeanderingFactor))   me%riverMeanderingFactor   = &
            me%calculateMeanderingFactorFromCellSize()
        if (isZero(me%estuaryMeanderingFactor)) me%estuaryMeanderingFactor = &
            me%calculateMeanderingFactorFromCellSize()

        ! Chunk-varying variables
        call me%readBatchVariables()

        ! Close input dataset
        call me%nc%close()

        ! Simulation mask (same reversal: file (y,x) -> returned (x,y))
        if (C%hasSimulationMask) then
            nc_simulationMask = NcDataset(C%simulationMaskPath, 'r')
            var = nc_simulationMask%getVariable('simulation_mask')
            call var%getData(simulationMask)
            me%simulationMask = ulgcl(simulationMask)
            me%nNonMaskedCells = count(me%simulationMask)
        else
            allocate(me%simulationMask(nx, ny))
            me%simulationMask = .true.
            me%nNonMaskedCells = count(.not. me%gridMask)
        end if

        ! Audit & log
        call rslt%addErrors(.errors. me%audit())
        call rslt%addToTrace('Initialising database')
        call ERROR_HANDLER%trigger(errors=.errors.rslt)
        call LOGR%toFile("Initialising database: success")
        call LOGR%toConsole("Initialising database: "//COLOR_GREEN//"success"//COLOR_RESET)
    end subroutine

    !> Update the database based on data for a new chunk (k), or for the only chunk if this
    !! isn't a batch run.
    subroutine updateDatabase(me, k)
        class(Database) :: me
        integer         :: k

        ! Get the config options for this chunk
        C%inputFile   = C%batchInputFiles(k)
        C%constantsFile = C%batchConstantFiles(k)
        C%nTimeSteps  = C%batchNTimesteps(k)
        C%startDate   = C%batchStartDates(k)

        call me%parseConstants(C%constantsFile)

        me%nc = NcDataset(C%inputFile, 'r')

        ! Deallocate previous-chunk vars
        if (allocated(me%t)) deallocate(me%t)
        if (allocated(me%soilAttachmentRate)) deallocate(me%soilAttachmentRate)
        if (allocated(me%soilAttachmentEfficiency)) deallocate(me%soilAttachmentEfficiency)
        if (allocated(me%emissionsArealSoilContaminant)) deallocate(me%emissionsArealSoilContaminant)
        if (allocated(me%emissionsArealWaterContaminant)) deallocate(me%emissionsArealWaterContaminant)
        if (allocated(me%emissionsAtmosphericDryDepoContaminant)) deallocate(me%emissionsAtmosphericDryDepoContaminant)
        if (allocated(me%emissionsAtmosphericWetDepoContaminant)) deallocate(me%emissionsAtmosphericWetDepoContaminant)
        if (allocated(me%emissionsPointWaterContaminant)) deallocate(me%emissionsPointWaterContaminant)
        if (allocated(me%emissionsPointWaterCoords)) deallocate(me%emissionsPointWaterCoords)
        if (allocated(me%resuspensionAlpha)) deallocate(me%resuspensionAlpha)
        if (allocated(me%resuspensionBeta)) deallocate(me%resuspensionBeta)
        if (allocated(me%depositionAlpha)) deallocate(me%depositionAlpha)
        if (allocated(me%depositionBeta)) deallocate(me%depositionBeta)
        if (allocated(me%bankErosionAlpha)) deallocate(me%bankErosionAlpha)
        if (allocated(me%bankErosionBeta)) deallocate(me%bankErosionBeta)
        if (allocated(me%sedimentTransport_a)) deallocate(me%sedimentTransport_a)
        if (allocated(me%sedimentTransport_b)) deallocate(me%sedimentTransport_b)
        if (allocated(me%sedimentTransport_c)) deallocate(me%sedimentTransport_c)
        if (allocated(me%initialContaminantConcsSoil)) deallocate(me%initialContaminantConcsSoil)
        if (allocated(me%initialContaminantConcsWater)) deallocate(me%initialContaminantConcsWater)
        if (allocated(me%initialContaminantConcsSediment)) deallocate(me%initialContaminantConcsSediment)
        if (allocated(me%initialDissolvedConcsSoil)) deallocate(me%initialDissolvedConcsSoil)
        if (allocated(me%initialDissolvedConcsWater)) deallocate(me%initialDissolvedConcsWater)
        if (allocated(me%initialDissolvedConcsSediment)) deallocate(me%initialDissolvedConcsSediment)
        if (allocated(me%emissionsArealSoilDissolvedContaminant)) deallocate(me%emissionsArealSoilDissolvedContaminant)
        if (allocated(me%emissionsArealWaterDissolvedContaminant)) deallocate(me%emissionsArealWaterDissolvedContaminant)
        if (allocated(me%emissionsAtmosphericDryDepoDissolvedContaminant)) &
            deallocate(me%emissionsAtmosphericDryDepoDissolvedContaminant)
        if (allocated(me%emissionsAtmosphericWetDepoDissolvedContaminant)) &
            deallocate(me%emissionsAtmosphericWetDepoDissolvedContaminant)
        if (allocated(me%emissionsPointWaterDissolvedContaminant)) &
            deallocate(me%emissionsPointWaterDissolvedContaminant)

        call me%readBatchVariables()
        call me%nc%close()
    end subroutine

    !> Read variables in for the new chunk as part of a batch run
    subroutine readBatchVariablesDatabase(me)
        class(Database)     :: me
        type(NcVariable)    :: var
        type(NcDimension)   :: p_dim
        logical :: haveCoordVar
        integer :: n
        integer :: alloc_stat
        integer :: nx, ny, nt, nforms, nsizes, np
        integer :: f_pris, f_mat, f_tra

        ! temp arrays with explicit ranks that match legacy file vars
        real(dp), allocatable :: A2(:,:)          ! (x,y)
        real(dp), allocatable :: A3(:,:,:)        ! (x,y,t)
        real(dp), allocatable :: COORD4(:,:,:,:)  ! (x,y,p,d)
        real(dp), allocatable :: A4(:,:,:,:)      ! (x,y,t,p)
        real(dp), allocatable :: T2(:,:)          ! temp for spatial (x,y) → transpose → (y,x)

        nx     = me%gridShape(1)
        ny     = me%gridShape(2)
        nt     = C%nTimeSteps
        nsizes = C%contaminantDim(1)
        nforms = C%contaminantDim(2)

        ! Legacy form indices (cap to available number of forms)
        f_pris = 1
        f_mat  = merge(2, 1, nforms >= 2)
        f_tra  = merge(3, 1, nforms >= 3)

        ! ----------------------
        ! SPATIAL SOIL VARIABLES
        ! ----------------------
        ! File vars are 2-D (y,x); the NetCDF helper returns (x,y) on getData.
        ! We read into T2(nx,ny) and store as (y,x) to match NetCDFOutput expectations.

        ! Soil bulk density [kg/m3]
        if (me%nc%hasVariable('soil_bulk_density')) then
            var = me%nc%getVariable('soil_bulk_density')
            if (allocated(me%soilBulkDensity)) deallocate(me%soilBulkDensity)
            allocate(T2(nx,ny)); call var%getData(T2)        ! (x,y)
            allocate(me%soilBulkDensity(ny,nx))              ! (y,x)
            me%soilBulkDensity = transpose(T2)               ! → (y,x)
            deallocate(T2)
        else
            if (allocated(me%soilBulkDensity)) deallocate(me%soilBulkDensity)
            allocate(me%soilBulkDensity(ny,nx))
            me%soilBulkDensity = nf90_fill_real
        end if

        ! Soil water content at field capacity [cm3/cm3]
        if (me%nc%hasVariable('soil_water_content_field_capacity')) then
            var = me%nc%getVariable('soil_water_content_field_capacity')
            if (allocated(me%soilWaterContentFieldCapacity)) deallocate(me%soilWaterContentFieldCapacity)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilWaterContentFieldCapacity(ny,nx))
            me%soilWaterContentFieldCapacity = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilWaterContentFieldCapacity)) deallocate(me%soilWaterContentFieldCapacity)
            allocate(me%soilWaterContentFieldCapacity(ny,nx))
            me%soilWaterContentFieldCapacity = nf90_fill_real
        end if

        ! Soil water content at saturation [cm3/cm3]
        if (me%nc%hasVariable('soil_water_content_saturation')) then
            var = me%nc%getVariable('soil_water_content_saturation')
            if (allocated(me%soilWaterContentSaturation)) deallocate(me%soilWaterContentSaturation)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilWaterContentSaturation(ny,nx))
            me%soilWaterContentSaturation = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilWaterContentSaturation)) deallocate(me%soilWaterContentSaturation)
            allocate(me%soilWaterContentSaturation(ny,nx))
            me%soilWaterContentSaturation = nf90_fill_real
        end if

        ! Soil hydraulic conductivity [m/s]
        if (me%nc%hasVariable('soil_hydraulic_conductivity')) then
            var = me%nc%getVariable('soil_hydraulic_conductivity')
            if (allocated(me%soilHydraulicConductivity)) deallocate(me%soilHydraulicConductivity)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilHydraulicConductivity(ny,nx))
            me%soilHydraulicConductivity = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilHydraulicConductivity)) deallocate(me%soilHydraulicConductivity)
            allocate(me%soilHydraulicConductivity(ny,nx))
            me%soilHydraulicConductivity = nf90_fill_real
        end if

        ! Soil texture [%] — clay
        if (me%nc%hasVariable('soil_texture_clay_content')) then
            var = me%nc%getVariable('soil_texture_clay_content')
            if (allocated(me%soilTextureClayContent)) deallocate(me%soilTextureClayContent)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilTextureClayContent(ny,nx))
            me%soilTextureClayContent = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilTextureClayContent)) deallocate(me%soilTextureClayContent)
            allocate(me%soilTextureClayContent(ny,nx))
            me%soilTextureClayContent = nf90_fill_real
        end if

        ! Soil texture [%] — sand
        if (me%nc%hasVariable('soil_texture_sand_content')) then
            var = me%nc%getVariable('soil_texture_sand_content')
            if (allocated(me%soilTextureSandContent)) deallocate(me%soilTextureSandContent)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilTextureSandContent(ny,nx))
            me%soilTextureSandContent = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilTextureSandContent)) deallocate(me%soilTextureSandContent)
            allocate(me%soilTextureSandContent(ny,nx))
            me%soilTextureSandContent = nf90_fill_real
        end if

        ! Soil texture [%] — silt
        if (me%nc%hasVariable('soil_texture_silt_content')) then
            var = me%nc%getVariable('soil_texture_silt_content')
            if (allocated(me%soilTextureSiltContent)) deallocate(me%soilTextureSiltContent)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilTextureSiltContent(ny,nx))
            me%soilTextureSiltContent = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilTextureSiltContent)) deallocate(me%soilTextureSiltContent)
            allocate(me%soilTextureSiltContent(ny,nx))
            me%soilTextureSiltContent = nf90_fill_real
        end if

        ! Soil texture [%] — coarse fragments
        if (me%nc%hasVariable('soil_texture_coarse_frag_content')) then
            var = me%nc%getVariable('soil_texture_coarse_frag_content')
            if (allocated(me%soilTextureCoarseFragContent)) deallocate(me%soilTextureCoarseFragContent)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilTextureCoarseFragContent(ny,nx))
            me%soilTextureCoarseFragContent = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilTextureCoarseFragContent)) deallocate(me%soilTextureCoarseFragContent)
            allocate(me%soilTextureCoarseFragContent(ny,nx))
            me%soilTextureCoarseFragContent = nf90_fill_real
        end if

        ! USLE factors [-] — C
        if (me%nc%hasVariable('soil_usle_c_factor')) then
            var = me%nc%getVariable('soil_usle_c_factor')
            if (allocated(me%soilUsleCFactor)) deallocate(me%soilUsleCFactor)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilUsleCFactor(ny,nx))
            me%soilUsleCFactor = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilUsleCFactor)) deallocate(me%soilUsleCFactor)
            allocate(me%soilUsleCFactor(ny,nx))
            me%soilUsleCFactor = nf90_fill_real
        end if

        ! USLE factors [-] — LS
        if (me%nc%hasVariable('soil_usle_ls_factor')) then
            var = me%nc%getVariable('soil_usle_ls_factor')
            if (allocated(me%soilUsleLSFactor)) deallocate(me%soilUsleLSFactor)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilUsleLSFactor(ny,nx))
            me%soilUsleLSFactor = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilUsleLSFactor)) deallocate(me%soilUsleLSFactor)
            allocate(me%soilUsleLSFactor(ny,nx))
            me%soilUsleLSFactor = nf90_fill_real
        end if

        ! USLE factors [-] — P
        if (me%nc%hasVariable('soil_usle_p_factor')) then
            var = me%nc%getVariable('soil_usle_p_factor')
            if (allocated(me%soilUslePFactor)) deallocate(me%soilUslePFactor)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilUslePFactor(ny,nx))
            me%soilUslePFactor = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilUslePFactor)) deallocate(me%soilUslePFactor)
            allocate(me%soilUslePFactor(ny,nx))
            me%soilUslePFactor = nf90_fill_real
        end if

        ! Soil attachment — prefer explicit rate; otherwise efficiency (default fallback elsewhere)
        if (me%nc%hasVariable('soil_attachment_rate')) then
            var = me%nc%getVariable('soil_attachment_rate')
            if (allocated(me%soilAttachmentRate)) deallocate(me%soilAttachmentRate)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilAttachmentRate(ny,nx))
            me%soilAttachmentRate = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilAttachmentRate)) deallocate(me%soilAttachmentRate)
            allocate(me%soilAttachmentRate(ny,nx))
            me%soilAttachmentRate = nf90_fill_real
        end if

        if (me%nc%hasVariable('soil_attachment_efficiency')) then
            var = me%nc%getVariable('soil_attachment_efficiency')
            if (allocated(me%soilAttachmentEfficiency)) deallocate(me%soilAttachmentEfficiency)
            allocate(T2(nx,ny)); call var%getData(T2)
            allocate(me%soilAttachmentEfficiency(ny,nx))
            me%soilAttachmentEfficiency = transpose(T2)
            deallocate(T2)
        else
            if (allocated(me%soilAttachmentEfficiency)) deallocate(me%soilAttachmentEfficiency)
            allocate(me%soilAttachmentEfficiency(ny,nx))
            me%soilAttachmentEfficiency = me%soilAttachmentEfficiencyConstant
        end if


        ! ----------------------
        ! SEDIMENT TRANSPORT
        ! ----------------------
        !
        ! If a spatial var is available in the NetCDF file, this is used. If not,
        ! the constant value is used for the whole grid. If a constant isn't available
        ! in the constant namelist, this will have already been set to the default
        ! from ConstantsDefaultsModule

        ! Deposition alpha
        if (me%nc%hasVariable('deposition_alpha')) then
            var = me%nc%getVariable('deposition_alpha')
            call var%getData(me%depositionAlpha)
        else
            allocate(me%depositionAlpha(me%gridShape(1), me%gridShape(2)))
            me%depositionAlpha = me%depositionAlphaConstant
        end if
        ! Deposition beta
        if (me%nc%hasVariable('deposition_beta')) then
            var = me%nc%getVariable('deposition_beta')
            call var%getData(me%depositionBeta)
        else
            allocate(me%depositionBeta(me%gridShape(1), me%gridShape(2)))
            me%depositionBeta = me%depositionBetaConstant
        end if

        ! Resuspension alpha
        ! TODO add check that at least one of the variables exists, or put a default in
        if (me%nc%hasVariable('resuspension_alpha')) then
            var = me%nc%getVariable('resuspension_alpha')
            call var%getData(me%resuspensionAlpha)
        else
            allocate(me%resuspensionAlpha(me%gridShape(1), me%gridShape(2)))
            ! Use the estuary mask to get a different resuspension alpha value in
            ! estuaries. If no specific alpha value is given for estuaries, it
            ! defaults to that for rivers
            where (me%isEstuary)
                me%resuspensionAlpha = me%waterResuspensionAlphaEstuary
            elsewhere
                me%resuspensionAlpha = me%waterResuspensionAlpha
            end where
        end if
        ! Resuspension beta
        if (me%nc%hasVariable('resuspension_beta')) then
            var = me%nc%getVariable('resuspension_beta')
            call var%getData(me%resuspensionBeta)
        else
            allocate(me%resuspensionBeta(me%gridShape(1), me%gridShape(2)))
            ! Same as for alpha, use different value for estuaries if available
            where (me%isEstuary)
                me%resuspensionBeta = me%waterResuspensionBetaEstuary
            elsewhere
                me%resuspensionBeta = me%waterResuspensionBeta
            end where
        end if

        ! Sediment transport param a 
        if (me%nc%hasVariable('sediment_transport_a')) then 
            var = me%nc%getVariable('sediment_transport_a') 
            call var%getData(me%sedimentTransport_a) 
        else 
            allocate(me%sedimentTransport_a(me%gridShape(1), me%gridShape(2))) 
            me%sedimentTransport_a = me%sedimentTransport_aConstant 
        end if 

        ! Sediment transport param b 
        if (me%nc%hasVariable('sediment_transport_b')) then 
            var = me%nc%getVariable('sediment_transport_b') 
            call var%getData(me%sedimentTransport_b) 
        else 
            allocate(me%sedimentTransport_b(me%gridShape(1), me%gridShape(2))) 
            me%sedimentTransport_b = me%sedimentTransport_bConstant 
        end if 

        ! Sediment transport param b 
        if (me%nc%hasVariable('sediment_transport_c')) then 
            var = me%nc%getVariable('sediment_transport_c') 
            call var%getData(me%sedimentTransport_c) 
        else 
            allocate(me%sedimentTransport_c(me%gridShape(1), me%gridShape(2))) 
            me%sedimentTransport_c = me%sedimentTransport_cConstant 
        end if 

        ! Bank erosion alpha
        if (me%nc%hasVariable('bank_erosion_alpha')) then
            var = me%nc%getVariable('bank_erosion_alpha')
            call var%getData(me%bankErosionAlpha)
        else
            allocate(me%bankErosionAlpha(me%gridShape(1), me%gridShape(2)))
            me%bankErosionAlpha = me%bankErosionAlphaConstant
        end if

        ! Bank erosion beta
        if (me%nc%hasVariable('bank_erosion_beta')) then
            var = me%nc%getVariable('bank_erosion_beta')
            call var%getData(me%bankErosionBeta)
        else
            allocate(me%bankErosionBeta(me%gridShape(1), me%gridShape(2)))
            me%bankErosionBeta = me%bankErosionBetaConstant
        end if
        

        !----------------------
        ! BASIC TIME SERIES
        !----------------------
        ! Total runoff (me%runoff) and surface runoff (me%quickflow) are derived
        ! from whichever runoff components the dataset provides
        call me%parseHydrology()

        if (me%nc%hasVariable('precip')) then
            var = me%nc%getVariable('precip')
            if (allocated(me%precip)) deallocate(me%precip)
            allocate(me%precip(nx,ny,nt))
            call var%getData(me%precip)
        else
            if (allocated(me%precip)) deallocate(me%precip)
            allocate(me%precip(nx,ny,nt))
            me%precip = 0.0_dp
        end if

        if (me%nc%hasVariable('evap')) then
            var = me%nc%getVariable('evap')
            if (allocated(me%evap)) deallocate(me%evap)
            allocate(me%evap(nx,ny,nt))
            call var%getData(me%evap)
        else
            if (allocated(me%evap)) deallocate(me%evap)
            allocate(me%evap(nx,ny,nt))
            me%evap = 0.0_dp
        end if

        ! Meteorological time series used by lake segments (stratification and wind
        ! resuspension). Both optional: wind speed falls back to the scalar windSpeed
        ! constant, solar radiation.
        if (me%nc%hasVariable('wind_speed')) then
            var = me%nc%getVariable('wind_speed')
            if (allocated(me%windSpeedTimeSeries)) deallocate(me%windSpeedTimeSeries)
            allocate(me%windSpeedTimeSeries(nx,ny,nt))
            call var%getData(me%windSpeedTimeSeries)
        else
            if (allocated(me%windSpeedTimeSeries)) deallocate(me%windSpeedTimeSeries)
        end if

        if (me%nc%hasVariable('solar_radiation')) then
            var = me%nc%getVariable('solar_radiation')
            if (allocated(me%solarRadiation)) deallocate(me%solarRadiation)
            allocate(me%solarRadiation(nx,ny,nt))
            call var%getData(me%solarRadiation)
        else
            if (allocated(me%solarRadiation)) deallocate(me%solarRadiation)
        end if

        !----------------------
        ! POINT-SOURCE DIM
        !----------------------
        if (me%nc%hasDimension('p')) then
            p_dim              = me%nc%getDimension('p')
            me%maxPointSources = p_dim%getLength()
        else
            me%maxPointSources = 0
        end if
        np = max(1, me%maxPointSources)

        !----------------------
        ! DEALLOC & ALLOC EMISSIONS
        !----------------------
        if (allocated(me%emissionsArealSoilContaminant)) &
            deallocate(me%emissionsArealSoilContaminant)
        if (allocated(me%emissionsArealWaterContaminant)) &
            deallocate(me%emissionsArealWaterContaminant)
        if (allocated(me%emissionsAtmosphericDryDepoContaminant)) &
            deallocate(me%emissionsAtmosphericDryDepoContaminant)
        if (allocated(me%emissionsAtmosphericWetDepoContaminant)) &
            deallocate(me%emissionsAtmosphericWetDepoContaminant)
        if (allocated(me%emissionsPointWaterContaminant)) &
            deallocate(me%emissionsPointWaterContaminant)
        if (allocated(me%emissionsArealSoilDissolvedContaminant)) &
            deallocate(me%emissionsArealSoilDissolvedContaminant)
        if (allocated(me%emissionsArealWaterDissolvedContaminant)) &
            deallocate(me%emissionsArealWaterDissolvedContaminant)
        if (allocated(me%emissionsAtmosphericDryDepoDissolvedContaminant)) &
            deallocate(me%emissionsAtmosphericDryDepoDissolvedContaminant)
        if (allocated(me%emissionsAtmosphericWetDepoDissolvedContaminant)) &
            deallocate(me%emissionsAtmosphericWetDepoDissolvedContaminant)
        if (allocated(me%emissionsPointWaterDissolvedContaminant)) &
            deallocate(me%emissionsPointWaterDissolvedContaminant)
        if (allocated(me%emissionsPointWaterCoords)) &
            deallocate(me%emissionsPointWaterCoords)

        allocate( &
            me%emissionsArealSoilContaminant( nx, ny, nsizes, nforms, C%contaminantDim(3) ), &
            me%emissionsArealWaterContaminant( nx, ny, nsizes, nforms, C%contaminantDim(3) ), &
            me%emissionsAtmosphericDryDepoContaminant( nx, ny, nt, nsizes, nforms, &
                                                    C%contaminantDim(3) ), &
            me%emissionsAtmosphericWetDepoContaminant( nx, ny, nt, nsizes, nforms, &
                                                    C%contaminantDim(3) ), &
            me%emissionsPointWaterContaminant( nx, ny, nt, np, nsizes, nforms, &
                                            C%contaminantDim(3) ), &
            me%emissionsArealSoilDissolvedContaminant( nx, ny ), &
            me%emissionsArealWaterDissolvedContaminant( nx, ny ), &
            me%emissionsAtmosphericDryDepoDissolvedContaminant( nx, ny, nt ), &
            me%emissionsAtmosphericWetDepoDissolvedContaminant( nx, ny, nt ), &
            me%emissionsPointWaterDissolvedContaminant( nx, ny, nt ), &
            stat=alloc_stat )

        if (alloc_stat /= 0) then
            call ERROR_HANDLER%trigger( &
                error=ErrorInstance(message='Emission allocation failed') )
            return
        end if

        me%emissionsArealSoilContaminant                   = 0.0_dp
        me%emissionsArealWaterContaminant                  = 0.0_dp
        me%emissionsAtmosphericDryDepoContaminant          = 0.0_dp
        me%emissionsAtmosphericWetDepoContaminant          = 0.0_dp
        me%emissionsPointWaterContaminant                  = 0.0_dp
        me%emissionsArealSoilDissolvedContaminant          = 0.0_dp
        me%emissionsArealWaterDissolvedContaminant         = 0.0_dp
        me%emissionsAtmosphericDryDepoDissolvedContaminant = 0.0_dp
        me%emissionsAtmosphericWetDepoDissolvedContaminant = 0.0_dp
        me%emissionsPointWaterDissolvedContaminant         = 0.0_dp

        !-----------------------------------
        ! AREAL EMISSIONS (2-D y,x -> (x,y))
        !-----------------------------------

        ! 1. SOIL EMISSIONS (PRISTINE) - Updated with Fallback
        ! ----------------------------------------------------
        if (me%nc%hasVariable('emissions_areal_soil_pristine')) then
            ! Try NEW name first
            var = me%nc%getVariable('emissions_areal_soil_pristine')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            do n = 1, nsizes
                me%emissionsArealSoilContaminant(:,:,n,f_pris,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
            call LOGR%add("DataInput: Read 'emissions_areal_soil_pristine'")

        else if (me%nc%hasVariable('emissions_areal_soil_nm')) then
            ! FALLBACK: Try OLD name (legacy support)
            var = me%nc%getVariable('emissions_areal_soil_nm')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            
            ! Map legacy 'nm' emissions to 'pristine' form
            do n = 1, nsizes
                me%emissionsArealSoilContaminant(:,:,n,f_pris,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
            call LOGR%add("DataInput: Read legacy 'emissions_areal_soil_nm' as pristine")
        end if

        ! 2. SOIL EMISSIONS (MATRIX EMBEDDED)
        ! -----------------------------------
        if (me%nc%hasVariable('emissions_areal_soil_matrixembedded')) then
            var = me%nc%getVariable('emissions_areal_soil_matrixembedded')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            do n = 1, nsizes
                me%emissionsArealSoilContaminant(:,:,n,f_mat,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
        end if

        ! 3. SOIL EMISSIONS (TRANSFORMED)
        ! -------------------------------
        if (me%nc%hasVariable('emissions_areal_soil_transformed')) then
            var = me%nc%getVariable('emissions_areal_soil_transformed')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            do n = 1, nsizes
                me%emissionsArealSoilContaminant(:,:,n,f_tra,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
        end if

        ! 4. WATER EMISSIONS (PRISTINE) - Updated with Fallback
        ! -----------------------------------------------------
        if (me%nc%hasVariable('emissions_areal_water_pristine')) then
            ! Try NEW name first
            var = me%nc%getVariable('emissions_areal_water_pristine')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            do n = 1, nsizes
                me%emissionsArealWaterContaminant(:,:,n,f_pris,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
            call LOGR%add("DataInput: Read 'emissions_areal_water_pristine'")

        else if (me%nc%hasVariable('emissions_areal_water_nm')) then
            ! FALLBACK: Try OLD name (legacy support)
            var = me%nc%getVariable('emissions_areal_water_nm')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            
            ! Map legacy 'nm' emissions to 'pristine' form
            do n = 1, nsizes
                me%emissionsArealWaterContaminant(:,:,n,f_pris,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
            call LOGR%add("DataInput: Read legacy 'emissions_areal_water_nm' as pristine")
        end if

        ! 5. WATER EMISSIONS (MATRIX EMBEDDED)
        ! ------------------------------------
        if (me%nc%hasVariable('emissions_areal_water_matrixembedded')) then
            var = me%nc%getVariable('emissions_areal_water_matrixembedded')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            do n = 1, nsizes
                me%emissionsArealWaterContaminant(:,:,n,f_mat,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
        end if

        ! 6. WATER EMISSIONS (TRANSFORMED)
        ! --------------------------------
        if (me%nc%hasVariable('emissions_areal_water_transformed')) then
            var = me%nc%getVariable('emissions_areal_water_transformed')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            do n = 1, nsizes
                me%emissionsArealWaterContaminant(:,:,n,f_tra,PFAS_AQ) = &
                    A2 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A2)
        end if

        ! 7. DISSOLVED EMISSIONS
        ! ----------------------
        if (me%nc%hasVariable('emissions_areal_soil_dissolved')) then
            var = me%nc%getVariable('emissions_areal_soil_dissolved')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            me%emissionsArealSoilDissolvedContaminant = A2
            deallocate(A2)
        end if

        if (me%nc%hasVariable('emissions_areal_water_dissolved')) then
            var = me%nc%getVariable('emissions_areal_water_dissolved')
            allocate(A2(nx,ny)); call var%getData(A2); A2 = zeroFillValue(A2)
            me%emissionsArealWaterDissolvedContaminant = A2
            deallocate(A2)
        end if

        !-----------------------------------------
        ! ATMOSPHERIC DEPOSITION (3-D t,y,x → (x,y,t))
        !-----------------------------------------
        if (me%nc%hasVariable('emissions_atmospheric_drydepo_pristine')) then
            var = me%nc%getVariable('emissions_atmospheric_drydepo_pristine')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            do n = 1, nsizes
                me%emissionsAtmosphericDryDepoContaminant(:,:,:,n,f_pris,PFAS_AQ) = &
                    A3 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_drydepo_matrixembedded')) then
            var = me%nc%getVariable('emissions_atmospheric_drydepo_matrixembedded')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            do n = 1, nsizes
                me%emissionsAtmosphericDryDepoContaminant(:,:,:,n,f_mat,PFAS_AQ) = &
                    A3 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_drydepo_transformed')) then
            var = me%nc%getVariable('emissions_atmospheric_drydepo_transformed')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            do n = 1, nsizes
                me%emissionsAtmosphericDryDepoContaminant(:,:,:,n,f_tra,PFAS_AQ) = &
                    A3 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_wetdepo_pristine')) then
            var = me%nc%getVariable('emissions_atmospheric_wetdepo_pristine')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            do n = 1, nsizes
                me%emissionsAtmosphericWetDepoContaminant(:,:,:,n,f_pris,PFAS_AQ) = &
                    A3 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_wetdepo_matrixembedded')) then
            var = me%nc%getVariable('emissions_atmospheric_wetdepo_matrixembedded')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            do n = 1, nsizes
                me%emissionsAtmosphericWetDepoContaminant(:,:,:,n,f_mat,PFAS_AQ) = &
                    A3 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_wetdepo_transformed')) then
            var = me%nc%getVariable('emissions_atmospheric_wetdepo_transformed')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            do n = 1, nsizes
                me%emissionsAtmosphericWetDepoContaminant(:,:,:,n,f_tra,PFAS_AQ) = &
                    A3 * me%defaultDistributionContaminant(n)
            end do
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_drydepo_dissolved')) then
            var = me%nc%getVariable('emissions_atmospheric_drydepo_dissolved')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            me%emissionsAtmosphericDryDepoDissolvedContaminant = A3
            deallocate(A3)
        end if

        if (me%nc%hasVariable('emissions_atmospheric_wetdepo_dissolved')) then
            var = me%nc%getVariable('emissions_atmospheric_wetdepo_dissolved')
            allocate(A3(nx,ny,nt)); call var%getData(A3); A3 = zeroFillValue(A3)
            me%emissionsAtmosphericWetDepoDissolvedContaminant = A3
            deallocate(A3)
        end if

        !-----------------------------------------
        ! POINT-SOURCE COORDS (4-D x,y,p,2)
        !-----------------------------------------
        haveCoordVar = .false.
        if (me%nc%hasVariable('emissions_point_water_contaminant_coords')) then
            var = me%nc%getVariable('emissions_point_water_contaminant_coords')
            haveCoordVar = .true.
        else if (me%nc%hasVariable('emissions_point_water_pristine_coords')) then
            var = me%nc%getVariable('emissions_point_water_pristine_coords')
            haveCoordVar = .true.
        else if (me%nc%hasVariable('emissions_point_water_matrixembedded_coords')) then
            var = me%nc%getVariable('emissions_point_water_matrixembedded_coords')
            haveCoordVar = .true.
        end if

        if (allocated(me%emissionsPointWaterCoords)) &
            deallocate(me%emissionsPointWaterCoords)
        allocate(me%emissionsPointWaterCoords(nx, ny, me%maxPointSources, 2))

        if (haveCoordVar) then
            allocate(COORD4(nx, ny, me%maxPointSources, 2))
            call var%getData(COORD4)      ! library returns (x,y,p,d)
            me%emissionsPointWaterCoords = COORD4
            deallocate(COORD4)
        else
            me%emissionsPointWaterCoords = nf90_fill_double
        end if

        !-----------------------------------------
        ! POINT-SOURCE EMISSIONS (4-D p,t,y,x → (x,y,t,p))
        !-----------------------------------------
        if (me%maxPointSources > 0) then

            if (me%nc%hasVariable('emissions_point_water_pristine')) then
                var = me%nc%getVariable('emissions_point_water_pristine')
                allocate(A4(nx,ny,nt,np)); call var%getData(A4); A4 = zeroFillValue(A4)   ! (x,y,t,p)
                do n = 1, nsizes
                    me%emissionsPointWaterContaminant(:,:,:,1:np,n,f_pris,PFAS_AQ) = &
                        A4 * me%defaultDistributionContaminant(n)
                end do
                deallocate(A4)
            end if

            if (me%nc%hasVariable('emissions_point_water_matrixembedded')) then
                var = me%nc%getVariable('emissions_point_water_matrixembedded')
                allocate(A4(nx,ny,nt,np)); call var%getData(A4); A4 = zeroFillValue(A4)   ! (x,y,t,p)
                do n = 1, nsizes
                    me%emissionsPointWaterContaminant(:,:,:,1:np,n,f_mat,PFAS_AQ) = &
                        A4 * me%defaultDistributionContaminant(n)
                end do
                deallocate(A4)
            end if

            ! Optional: only if present in file
            if (me%nc%hasVariable('emissions_point_water_transformed')) then
                var = me%nc%getVariable('emissions_point_water_transformed')
                allocate(A4(nx,ny,nt,np)); call var%getData(A4); A4 = zeroFillValue(A4)   ! (x,y,t,p)
                do n = 1, nsizes
                    me%emissionsPointWaterContaminant(:,:,:,1:np,n,f_tra,PFAS_AQ) = &
                        A4 * me%defaultDistributionContaminant(n)
                end do
                deallocate(A4)
            end if

        end if

        ! -----------------------------------
        ! PFAS SOIL / SURFACE WATER PROPERTIES
        ! -----------------------------------
        ! All optional, (x,y) indexed like deposition_alpha (not transposed). Missing,
        ! fill and negative cell values are stored as -1 so that the parameterisation
        ! falls back cell-by-cell (see PFASParameterisationModule).
        call readPfasMap2D('soil_organic_carbon', me%soilOrganicCarbon)
        call readPfasMap2D('soil_ph', me%soilPh)
        call readPfasMap2D('spm_organic_carbon', me%spmOrganicCarbon)
        call readPfasMap2D('water_ph', me%waterPh)
        call readPfasMap2D('soil_awi_specific_area', me%soilAwiSpecificArea)
        ! Direct Kd maps are given in L kg-1, as in the &pfas group, and stored in m3 kg-1
        call readPfasMap3D('pfas_kd_soil', me%pfasKdSoilMap)
        call readPfasMap3D('pfas_kd_spm', me%pfasKdSpmMap)
        call readPfasMap3D('pfas_kd_sediment', me%pfasKdSedimentMap)
        call readPfasMap2D('sediment_organic_carbon', me%sedimentOrganicCarbon)
        call readPfasMap2D('water_salinity', me%waterSalinity)
        call readPfasMap2D('groundwater_residence_time', me%groundwaterResidenceTime)

        ! -----------------------------------
        ! INITIAL CONCENTRATIONS
        ! -----------------------------------
        if (allocated(me%initialContaminantConcsSoil))      deallocate(me%initialContaminantConcsSoil)
        if (allocated(me%initialContaminantConcsWater))     deallocate(me%initialContaminantConcsWater)
        if (allocated(me%initialContaminantConcsSediment))  deallocate(me%initialContaminantConcsSediment)
        
        allocate(me%initialContaminantConcsSoil(    nx,ny,nsizes,nforms, C%contaminantDim(3)))
        allocate(me%initialContaminantConcsWater(   nx,ny,nsizes,nforms, C%contaminantDim(3)))
        allocate(me%initialContaminantConcsSediment(nx,ny,nsizes,nforms, C%contaminantDim(3)))

        me%initialContaminantConcsSoil     = 0.0_dp
        me%initialContaminantConcsWater    = 0.0_dp
        me%initialContaminantConcsSediment = 0.0_dp

        ! --- 1. SOIL INITIAL CONCENTRATIONS ---
        if (me%nc%hasVariable('initial_contaminant_concs_soil')) then
            ! Try NEW name
            var = me%nc%getVariable('initial_contaminant_concs_soil')
            call var%getData(me%initialContaminantConcsSoil)
            call LOGR%add("DataInput: Read 'initial_contaminant_concs_soil'")
        else if (me%nc%hasVariable('initial_nm_concs_soil')) then
            ! FALLBACK: Try OLD name
            var = me%nc%getVariable('initial_nm_concs_soil')
            call var%getData(me%initialContaminantConcsSoil)
            call LOGR%add("DataInput: Read legacy 'initial_nm_concs_soil'")
        else
            call LOGR%add("DataInput: WARNING - No initial SOIL contaminant data found. Set to 0.0.")
        end if

        ! --- 2. WATER INITIAL CONCENTRATIONS ---
        if (me%nc%hasVariable('initial_contaminant_concs_water')) then
            ! Try NEW name
            var = me%nc%getVariable('initial_contaminant_concs_water')
            call var%getData(me%initialContaminantConcsWater)
            call LOGR%add("DataInput: Read 'initial_contaminant_concs_water'")
        else if (me%nc%hasVariable('initial_nm_concs_water')) then
            ! FALLBACK: Try OLD name
            var = me%nc%getVariable('initial_nm_concs_water')
            call var%getData(me%initialContaminantConcsWater)
            call LOGR%add("DataInput: Read legacy 'initial_nm_concs_water'")
        else
            call LOGR%add("DataInput: WARNING - No initial WATER contaminant data found. Set to 0.0.")
        end if

        ! --- 3. SEDIMENT INITIAL CONCENTRATIONS ---
        if (me%nc%hasVariable('initial_contaminant_concs_sediment')) then
            ! Try NEW name
            var = me%nc%getVariable('initial_contaminant_concs_sediment')
            call var%getData(me%initialContaminantConcsSediment)
            call LOGR%add("DataInput: Read 'initial_contaminant_concs_sediment'")
        else if (me%nc%hasVariable('initial_nm_concs_sediment')) then
            ! FALLBACK: Try OLD name
            var = me%nc%getVariable('initial_nm_concs_sediment')
            call var%getData(me%initialContaminantConcsSediment)
            call LOGR%add("DataInput: Read legacy 'initial_nm_concs_sediment'")
        else
            call LOGR%add("DataInput: WARNING - No initial SEDIMENT contaminant data found. Set to 0.0.")
        end if

        if (me%nc%hasVariable('initial_dissolved_concentrations_soil')) then
            var = me%nc%getVariable('initial_dissolved_concentrations_soil')
            allocate(A2(nx,ny)); call var%getData(A2)
            me%initialDissolvedConcsSoil = A2
            deallocate(A2)
        end if
        if (me%nc%hasVariable('initial_dissolved_concentrations_water')) then
            var = me%nc%getVariable('initial_dissolved_concentrations_water')
            allocate(A2(nx,ny)); call var%getData(A2)
            me%initialDissolvedConcsWater = A2
            deallocate(A2)
        end if
        if (me%nc%hasVariable('initial_dissolved_concentrations_sediment')) then
            var = me%nc%getVariable('initial_dissolved_concentrations_sediment')
            allocate(A2(nx,ny)); call var%getData(A2)
            me%initialDissolvedConcsSediment = A2
            deallocate(A2)
        end if

        ! Count point sources after coords are populated
        call me%calculateNPointSources(me%maxPointSources)

      contains

        !> Read an optional (x,y) PFAS property map. Left unallocated if absent.
        subroutine readPfasMap2D(name, map)
            character(len=*), intent(in)         :: name
            real(dp), allocatable, intent(inout) :: map(:,:)
            type(NcVariable) :: mapVar
            if (allocated(map)) deallocate(map)
            if (.not. me%nc%hasVariable(name)) return
            mapVar = me%nc%getVariable(name)
            call mapVar%getData(map)
            if (size(map, 1) /= nx .or. size(map, 2) /= ny) then
                call LOGR%add("DataInput: WARNING - '" // name // "' is not (x,y) on the model grid; ignored")
                deallocate(map)
                return
            end if
            map = cleanPfasMapValue(map)
            call LOGR%add("DataInput: Read '" // name // "'")
        end subroutine

        !> Read an optional (x,y,species) Kd map given in L kg-1, storing it in m3 kg-1.
        !! A 2D (x,y) variable is accepted when there is a single PFAS species.
        subroutine readPfasMap3D(name, map)
            character(len=*), intent(in)         :: name
            real(dp), allocatable, intent(inout) :: map(:,:,:)
            type(NcVariable) :: mapVar
            real(dp), allocatable :: map2D(:,:)
            if (allocated(map)) deallocate(map)
            if (.not. me%nc%hasVariable(name)) return
            mapVar = me%nc%getVariable(name)
            if (mapVar%getRank() == 2 .and. nsizes == 1) then
                call mapVar%getData(map2D)
                allocate(map(size(map2D, 1), size(map2D, 2), 1))
                map(:,:,1) = map2D
            else if (mapVar%getRank() == 3) then
                call mapVar%getData(map)
            else
                call LOGR%add("DataInput: WARNING - '" // name // "' must be (x,y,species); ignored")
                return
            end if
            if (size(map, 1) /= nx .or. size(map, 2) /= ny .or. size(map, 3) < nsizes) then
                call LOGR%add("DataInput: WARNING - '" // name // "' does not match (x,y,species); ignored")
                deallocate(map)
                return
            end if
            map = cleanPfasMapValue(map)
            where (map >= 0.0_dp) map = map * 1.0e-3_dp
            call LOGR%add("DataInput: Read '" // name // "'")
        end subroutine
    end subroutine readBatchVariablesDatabase

    !> Map fill values, NaNs and negative values to -1, the "no value here" marker
    !! used by the PFAS parameterisation to fall back to the next option.
    elemental function cleanPfasMapValue(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (x /= x .or. x < 0.0_dp .or. x == nf90_fill_double .or. x == real(nf90_fill_real, dp) &
            .or. abs(x) > 1.0e30_dp) then
            y = -1.0_dp
        else
            y = x
        end if
    end function

    !> Hydrology input adapter. The model needs two runoff components per cell and
    !! timestep: total runoff (me%runoff, which generates water inflow to reaches)
    !! and surface/overland runoff (me%quickflow, which drives terrestrial erosion
    !! and material delivery). Rather than requiring the input dataset to provide
    !! exactly these, this routine maps whichever components are available - named
    !! via the &hydrology config group - onto them:
    !!  - total + surface:      both used directly
    !!  - total + subsurface:   surface = max(total - subsurface, 0)
    !!  - surface + subsurface: total = surface + subsurface
    !!  - total only:           surface = surface_runoff_fraction * total
    !!  - surface only:         total = surface (with warning)
    !! All components are in m/timestep. Total runoff is the minimum requirement;
    !! if no component is found, both are set to zero and a warning is logged.
    subroutine parseHydrologyDatabase(me)
        class(Database)     :: me
        type(NcVariable)    :: var
        integer             :: nx, ny, nt
        real, allocatable   :: subsurfaceRunoff(:,:,:)
        logical             :: hasTotal, hasSurface, hasSubsurface

        nx = me%gridShape(1)
        ny = me%gridShape(2)
        nt = C%nTimeSteps

        ! An empty variable name in config means "not provided by this dataset"
        hasTotal = len_trim(C%hydrologyTotalRunoffVariable) > 0
        if (hasTotal) hasTotal = me%nc%hasVariable(trim(C%hydrologyTotalRunoffVariable))
        hasSurface = len_trim(C%hydrologySurfaceRunoffVariable) > 0
        if (hasSurface) hasSurface = me%nc%hasVariable(trim(C%hydrologySurfaceRunoffVariable))
        hasSubsurface = len_trim(C%hydrologySubsurfaceRunoffVariable) > 0
        if (hasSubsurface) hasSubsurface = me%nc%hasVariable(trim(C%hydrologySubsurfaceRunoffVariable))

        if (allocated(me%runoff)) deallocate(me%runoff)
        if (allocated(me%quickflow)) deallocate(me%quickflow)
        allocate(me%runoff(nx,ny,nt))
        allocate(me%quickflow(nx,ny,nt))

        if (hasTotal) then
            var = me%nc%getVariable(trim(C%hydrologyTotalRunoffVariable))    ! (t,y,x) in file
            call var%getData(me%runoff)                                      ! library reverses -> (x,y,t)
        end if
        if (hasSurface) then
            var = me%nc%getVariable(trim(C%hydrologySurfaceRunoffVariable))
            call var%getData(me%quickflow)
        end if
        if (hasSubsurface) then
            allocate(subsurfaceRunoff(nx,ny,nt))
            var = me%nc%getVariable(trim(C%hydrologySubsurfaceRunoffVariable))
            call var%getData(subsurfaceRunoff)
        end if

        if (hasTotal .and. hasSurface) then
            call LOGR%add("DataInput: Hydrology - total runoff '" // trim(C%hydrologyTotalRunoffVariable) // &
                "' and surface runoff '" // trim(C%hydrologySurfaceRunoffVariable) // "' read from input data")
        else if (hasSurface .and. hasSubsurface) then
            me%runoff = me%quickflow + subsurfaceRunoff
            call LOGR%add("DataInput: Hydrology - total runoff derived as surface ('" // &
                trim(C%hydrologySurfaceRunoffVariable) // "') + subsurface ('" // &
                trim(C%hydrologySubsurfaceRunoffVariable) // "') runoff")
        else if (hasTotal .and. hasSubsurface) then
            me%quickflow = max(me%runoff - subsurfaceRunoff, 0.0)
            call LOGR%add("DataInput: Hydrology - surface runoff derived as total ('" // &
                trim(C%hydrologyTotalRunoffVariable) // "') - subsurface ('" // &
                trim(C%hydrologySubsurfaceRunoffVariable) // "') runoff, floored at zero")
        else if (hasTotal) then
            ! Keep the surface component exactly zero (not 0 * runoff) when the
            ! fraction is zero, so fill values in runoff can't leak into quickflow
            if (C%hydrologySurfaceRunoffFraction /= 0.0) then
                me%quickflow = C%hydrologySurfaceRunoffFraction * me%runoff
            else
                me%quickflow = 0.0
            end if
            call LOGR%add("DataInput: Hydrology - only total runoff ('" // &
                trim(C%hydrologyTotalRunoffVariable) // "') found in input data. Surface runoff set using " // &
                "configured surface_runoff_fraction")
        else if (hasSurface) then
            me%runoff = me%quickflow
            call LOGR%add("DataInput: WARNING - Hydrology - only surface runoff ('" // &
                trim(C%hydrologySurfaceRunoffVariable) // "') found in input data. Using it as total runoff too. " // &
                "This is only appropriate for simple tests")
        else
            me%runoff = 0.0
            me%quickflow = 0.0
            call LOGR%add("DataInput: WARNING - Hydrology - no runoff data found in input data. Total and surface " // &
                "runoff set to 0.0. Total runoff is the minimum hydrological requirement for a meaningful run")
        end if
    end subroutine parseHydrologyDatabase

    !> Merge the lake routing network into the main routing arrays. 'lake_outflow'
    !! gives, for each lake/reservoir cell, the cell its segment discharges to -
    !! the next segment of the same lake, or the river/estuary cell at the lake
    !! outlet. These override the general 'outflow' routing for flagged cells, and
    !! the inflow network (inflows, is_headwater, n_waterbodies) is then rebuilt
    !! from the merged outflow field so the whole network stays consistent: a
    !! cell's inflows are exactly the cells whose outflow points at it.
    !! This is only called when 'lake_outflow' is present in the input data, so
    !! datasets without a lake network are completely unaffected.
    subroutine mergeLakeRoutingDatabase(me, lakeOutflow)
        class(Database)     :: me
        integer, intent(in) :: lakeOutflow(:,:,:)   !! Lake outflow targets, (d,x,y) with d = (x,y) of target cell
        integer :: nx, ny, x, y, ox, oy, w
        integer :: nOverridden, nInvalid, maxInflows
        integer, allocatable :: nInflowCells(:,:)
        integer, allocatable :: newInflows(:,:,:,:)
        integer, parameter :: fillInt = -999

        nx = me%gridShape(1)
        ny = me%gridShape(2)

        ! 1. Override the outflow of flagged lake/reservoir cells where a valid
        !    target is given (fill/zero values leave the general routing in place)
        nOverridden = 0
        nInvalid = 0
        do y = 1, ny
            do x = 1, nx
                if ((me%isLake(x, y) .or. me%isReservoir(x, y)) .and. me%nWaterbodies(x, y) > 0) then
                    if (lakeOutflow(1, x, y) >= 1 .and. lakeOutflow(2, x, y) >= 1) then
                        if (lakeOutflow(1, x, y) == x .and. lakeOutflow(2, x, y) == y) then
                            ! A segment can't discharge to itself; keep the general routing
                            nInvalid = nInvalid + 1
                        else
                            me%outflow(1:2, x, y) = lakeOutflow(1:2, x, y)
                            nOverridden = nOverridden + 1
                        end if
                    end if
                end if
            end do
        end do

        ! 2. Count each cell's inflowing cells under the merged outflow field.
        !    Outflows leaving the grid, or pointing at non-waterbody cells, count
        !    as domain outflows and contribute no inflow.
        allocate(nInflowCells(nx, ny))
        nInflowCells = 0
        do y = 1, ny
            do x = 1, nx
                if (me%nWaterbodies(x, y) > 0) then
                    ox = me%outflow(1, x, y)
                    oy = me%outflow(2, x, y)
                    if (ox >= 1 .and. ox <= nx .and. oy >= 1 .and. oy <= ny) then
                        if (me%nWaterbodies(ox, oy) > 0 .and. .not. (ox == x .and. oy == y)) then
                            nInflowCells(ox, oy) = nInflowCells(ox, oy) + 1
                        end if
                    end if
                end if
            end do
        end do

        ! 3. Rebuild the inflows array (a cell's w-th inflow is the w-th cell found
        !    whose outflow points at it, matching the creation-order convention in
        !    Reach%parseInflowsAndOutflow)
        maxInflows = max(1, maxval(nInflowCells))
        allocate(newInflows(2, maxInflows, nx, ny))
        newInflows = fillInt
        nInflowCells = 0
        do y = 1, ny
            do x = 1, nx
                if (me%nWaterbodies(x, y) > 0) then
                    ox = me%outflow(1, x, y)
                    oy = me%outflow(2, x, y)
                    if (ox >= 1 .and. ox <= nx .and. oy >= 1 .and. oy <= ny) then
                        if (me%nWaterbodies(ox, oy) > 0 .and. .not. (ox == x .and. oy == y)) then
                            nInflowCells(ox, oy) = nInflowCells(ox, oy) + 1
                            w = nInflowCells(ox, oy)
                            newInflows(1, w, ox, oy) = x
                            newInflows(2, w, ox, oy) = y
                        end if
                    end if
                end if
            end do
        end do
        call move_alloc(from=newInflows, to=me%inflows)

        ! 4. Update headwater flags and waterbody (branch) counts. Cells with no
        !    waterbodies stay empty; the number of reaches in a cell equals its
        !    number of inflowing cells (min 1, for headwaters).
        do y = 1, ny
            do x = 1, nx
                if (me%nWaterbodies(x, y) > 0) then
                    me%isHeadwater(x, y) = nInflowCells(x, y) == 0
                    me%nWaterbodies(x, y) = max(1, nInflowCells(x, y))
                end if
            end do
        end do
        me%maxNWaterbodies = maxval(me%nWaterbodies)
        me%hasLakeRouting = .true.

        call LOGR%add("DataInput: Merged lake routing network ('lake_outflow'): " &
            // trim(str(nOverridden)) // " cell outflows overridden, inflow network rebuilt")
        if (nInvalid > 0) then
            call LOGR%add("DataInput: WARNING - " // trim(str(nInvalid)) &
                // " 'lake_outflow' entries pointed a cell at itself and were ignored")
        end if
    end subroutine mergeLakeRoutingDatabase

    !> Get the constants from the namelist file
    subroutine parseConstantsDatabase(me, constantsFile)
        class(Database)         :: me
        character(len=*)        :: constantsFile
        integer                 :: nmlIOStat
        character(len=256)      :: nmlIOMsg
        integer :: n_biota, n_contaminant_size_classes, n_default_spm_size_distribution, &
                n_default_matrixembedded_distribution_to_spm, n_vertical_distribution, &
                n_initial_c_org, n_k_growth, n_k_uptake_contaminant, n_k_elim_contaminant, &
                n_name, n_stored_fraction, n_k_uptake_dissolved, n_k_elim_dissolved, &
                n_uptake_from_form, n_harvest_in_month, n_porosity, n_initial_mass, &
                n_spm_density_by_size_class, n_fractional_composition_distribution, &
                n_estuary_mouth_coords, n_compartment, n_default_contaminant_size_distribution, &
                n_default_contaminant_form_distribution
        real :: estuary_mouth_coords(2)
        integer, allocatable :: default_contaminant_size_distribution(:), default_spm_size_distribution(:), &
                            default_matrixembedded_distribution_to_spm(:), vertical_distribution(:), &
                            harvest_in_month(:), prey(:)
        real(dp), allocatable :: k_diet(:), k_death(:)
        real, allocatable :: stored_fraction(:), porosity(:), initial_mass(:), &
                            fractional_composition_distribution(:), spm_density_by_size_class(:), &
                            default_contaminant_form_distribution(:)
        real :: darcy_velocity, default_porosity, particle_density, estuary_tidal_S2, &
                estuary_mean_depth_expA, estuary_mean_depth_expB, estuary_width_expA, &
                estuary_width_expB, estuary_tidal_M2, estuary_meandering_factor, &
                river_meandering_factor, deposition_alpha, deposition_beta, &
                bank_erosion_alpha, bank_erosion_beta
        real(dp) :: hamaker_constant, resuspension_alpha, resuspension_beta, &
                    resuspension_alpha_estuary, resuspension_beta_estuary, k_diss_pristine, &
                    k_diss_transformed, k_transform_pristine, erosivity_a1, erosivity_a2, &
                    erosivity_a3, erosivity_b, contaminant_density, estuary_attachment_efficiency, &
                    soil_attachment_efficiency, river_attachment_efficiency, &
                    sediment_transport_a, sediment_transport_b, sediment_transport_c, &
                    sediment_enrichment_k, sediment_enrichment_a, min_water_temperature, &
                    max_water_temperature, shear_rate
        real(dp), allocatable :: initial_C_org(:), k_growth(:), k_uptake_contaminant(:), &
                                k_elim_contaminant(:), k_uptake_dissolved(:), k_elim_dissolved(:), &
                                contaminant_size_classes(:)
        character(len=100), allocatable :: name(:), compartment(:)
        character(len=17), allocatable :: uptake_from_form(:)
        integer :: min_water_temperature_day_of_year, arable, coniferous, deciduous, grassland, &
                heathland, urban_capped, urban_gardens, urban_parks

        real(dp) :: lake_default_depth, lake_default_area_fraction, &
            lake_stratification_radiation_threshold, lake_stratification_wind_threshold, &
            lake_stratification_temperature_threshold, lake_epilimnion_base_depth, &
            lake_epilimnion_wind_coefficient, lake_epilimnion_min_depth, &
            lake_interlayer_mixing_velocity, lake_wind_resuspension_attenuation_depth, &
            reservoir_default_release_rate

        ! &pfas group (optional). Species arrays are sized from the config's n_pfas_species.
        integer :: nPfas
        character(len=64), allocatable :: species_names(:)
        real(dp), allocatable :: molecular_weight(:), charge(:), pka(:), kd_soil(:), kd_spm(:), &
            koc(:), kd_ph_slope(:), kia(:), k_mt_soil(:), k_mt_spm(:), k_mt_awi(:), &
            reaction_rate(:,:), reaction_yield(:,:), kd_sediment(:), reaction_rate_sediment(:,:), &
            diffusivity(:), kd_salinity_coefficient(:), kia_salinity_coefficient(:), &
            biota_transformation_rate(:,:)
        real(dp) :: kd_ph_reference, soil_organic_carbon, soil_ph, spm_organic_carbon, water_ph, &
            soil_awi_specific_area, awi_geometric_factor, foam_coefficient, volatilisation_rate, &
            sea_spray_aerosol_rate, sediment_organic_carbon, sediment_water_exchange_velocity, &
            estuary_salinity, runoff_mixing_depth, groundwater_residence_time
        namelist /pfas/ species_names, molecular_weight, charge, pka, kd_soil, kd_spm, koc, &
            kd_ph_slope, kd_ph_reference, kia, k_mt_soil, k_mt_spm, k_mt_awi, reaction_rate, &
            reaction_yield, soil_organic_carbon, soil_ph, spm_organic_carbon, water_ph, &
            soil_awi_specific_area, awi_geometric_factor, foam_coefficient, volatilisation_rate, &
            sea_spray_aerosol_rate, kd_sediment, sediment_organic_carbon, reaction_rate_sediment, &
            diffusivity, sediment_water_exchange_velocity, kd_salinity_coefficient, &
            kia_salinity_coefficient, estuary_salinity, runoff_mixing_depth, groundwater_residence_time, &
            biota_transformation_rate
        namelist /allocatable_array_sizes/ n_default_contaminant_size_distribution, &
            n_default_spm_size_distribution, n_default_matrixembedded_distribution_to_spm, &
            n_vertical_distribution, n_initial_c_org, n_k_growth, n_name, &
            n_stored_fraction, n_k_uptake_contaminant, n_k_elim_contaminant, &
            n_compartment, n_k_uptake_dissolved, n_k_elim_dissolved, n_uptake_from_form, &
            n_harvest_in_month, n_porosity, n_spm_density_by_size_class, &
            n_initial_mass, n_fractional_composition_distribution, n_estuary_mouth_coords, &
            n_contaminant_size_classes, n_default_contaminant_form_distribution
        namelist /n_biota_grp/ n_biota
        namelist /contaminant/ contaminant_density, default_contaminant_size_distribution, &
            contaminant_size_classes, k_diss_pristine, k_diss_transformed, k_transform_pristine, &
            default_contaminant_form_distribution
        namelist /biota/ initial_C_org, k_growth, k_elim_contaminant, k_uptake_contaminant, &
            name, stored_fraction, compartment, k_uptake_dissolved, k_elim_dissolved, &
            uptake_from_form, harvest_in_month, prey, k_diet, k_death
        namelist /earthworm_densities/ arable, coniferous, deciduous, grassland, heathland, &
            urban_capped, urban_gardens, urban_parks, vertical_distribution
        namelist /soil/ darcy_velocity, default_porosity, hamaker_constant, particle_density, &
            erosivity_a1, erosivity_a2, erosivity_a3, erosivity_b, soil_attachment_efficiency, &
            sediment_transport_a, sediment_transport_b, sediment_transport_c
        namelist /water/ resuspension_alpha, resuspension_beta, resuspension_alpha_estuary, &
            resuspension_beta_estuary, estuary_tidal_m2, estuary_tidal_s2, estuary_mouth_coords, &
            estuary_mean_depth_expa, estuary_mean_depth_expb, estuary_width_expa, estuary_width_expb, &
            estuary_meandering_factor, river_meandering_factor, river_attachment_efficiency, &
            estuary_attachment_efficiency, deposition_alpha, deposition_beta, bank_erosion_alpha, &
            bank_erosion_beta, shear_rate, min_water_temperature, max_water_temperature, &
            min_water_temperature_day_of_year
        namelist /sediment/ porosity, initial_mass, fractional_composition_distribution, &
            default_spm_size_distribution, default_matrixembedded_distribution_to_spm, &
            sediment_enrichment_a, sediment_enrichment_k, spm_density_by_size_class
        namelist /lake/ lake_default_depth, lake_default_area_fraction, &
            lake_stratification_radiation_threshold, lake_stratification_wind_threshold, &
            lake_stratification_temperature_threshold, lake_epilimnion_base_depth, &
            lake_epilimnion_wind_coefficient, lake_epilimnion_min_depth, &
            lake_interlayer_mixing_velocity, lake_wind_resuspension_attenuation_depth, &
            reservoir_default_release_rate

        ! Initialize variables
        n_biota = 0
        n_default_contaminant_form_distribution = 0
        contaminant_density = default_rho_contaminant
        k_diss_pristine = default_k_diss_pristine
        k_diss_transformed = default_k_diss_transformed
        k_transform_pristine = default_k_transform_pristine
        soil_attachment_efficiency = defaultSoilAttachmentEfficiency
        river_attachment_efficiency = defaultRiverAttachmentEfficiency
        resuspension_alpha_estuary = 0.0_dp
        resuspension_beta_estuary = 0.0_dp
        soil_attachment_efficiency = real(defaultSoilAttachmentEfficiency, dp)
        river_attachment_efficiency = real(defaultRiverAttachmentEfficiency, dp)
        estuary_attachment_efficiency = defaultEstuaryAttachmentEfficiency
        darcy_velocity = defaultSoilDarcyVelocity
        estuary_meandering_factor = 0.0
        river_meandering_factor = 0.0
        shear_rate = defaultShearRate
        min_water_temperature = defaultMinWaterTemperature
        max_water_temperature = defaultMaxWaterTemperature
        min_water_temperature_day_of_year = defaultMinWaterTemperatureDayOfYear
        sediment_transport_a = defaultSedimentTransport_a
        sediment_transport_b = defaultSedimentTransport_b
        sediment_transport_c = defaultSedimentTransport_c
        sediment_enrichment_k = defaultSedimentEnrichment_k
        sediment_enrichment_a = defaultSedimentEnrichment_a
        deposition_alpha = defaultDepositionAlpha
        deposition_beta = defaultDepositionBeta
        bank_erosion_alpha = defaultBankErosionAlpha
        bank_erosion_beta = defaultBankErosionBeta
        lake_default_depth = me%lakeDefaultDepth
        lake_default_area_fraction = me%lakeDefaultAreaFraction
        lake_stratification_radiation_threshold = me%lakeStratificationRadiationThreshold
        lake_stratification_wind_threshold = me%lakeStratificationWindThreshold
        lake_stratification_temperature_threshold = me%lakeStratificationTemperatureThreshold
        lake_epilimnion_base_depth = me%lakeEpilimnionBaseDepth
        lake_epilimnion_wind_coefficient = me%lakeEpilimnionWindCoefficient
        lake_epilimnion_min_depth = me%lakeEpilimnionMinDepth
        lake_interlayer_mixing_velocity = me%lakeInterlayerMixingVelocity
        lake_wind_resuspension_attenuation_depth = me%lakeWindResuspensionAttenuationDepth
        reservoir_default_release_rate = me%reservoirDefaultReleaseRate

        ! &pfas defaults: everything off/zero, so an absent group gives an inert tracer
        nPfas = max(1, C%contaminantDim(1))
        allocate(species_names(nPfas), molecular_weight(nPfas), charge(nPfas), pka(nPfas), &
            kd_soil(nPfas), kd_spm(nPfas), koc(nPfas), kd_ph_slope(nPfas), kia(nPfas), &
            k_mt_soil(nPfas), k_mt_spm(nPfas), k_mt_awi(nPfas), &
            reaction_rate(nPfas, nPfas), reaction_yield(nPfas, nPfas), kd_sediment(nPfas), &
            reaction_rate_sediment(nPfas, nPfas), diffusivity(nPfas), kd_salinity_coefficient(nPfas), &
            kia_salinity_coefficient(nPfas), biota_transformation_rate(nPfas, nPfas))
        species_names = ''
        molecular_weight = 0.0_dp
        charge = 0.0_dp
        pka = -999.0_dp
        kd_soil = 0.0_dp
        kd_spm = 0.0_dp
        koc = 0.0_dp
        kd_ph_slope = 0.0_dp
        kd_ph_reference = me%pfasKdPhReference
        kia = 0.0_dp
        k_mt_soil = 0.0_dp
        k_mt_spm = 0.0_dp
        k_mt_awi = 0.0_dp
        reaction_rate = 0.0_dp
        reaction_yield = 0.0_dp
        soil_organic_carbon = me%soilOrganicCarbonConstant
        soil_ph = me%soilPhConstant
        spm_organic_carbon = me%spmOrganicCarbonConstant
        water_ph = me%waterPhConstant
        soil_awi_specific_area = me%soilAwiSpecificAreaConstant
        awi_geometric_factor = me%soilAwiGeometricFactor
        foam_coefficient = me%pfasFoamCoefficient
        volatilisation_rate = me%pfasVolatilisationRateScalar
        sea_spray_aerosol_rate = me%pfasSeaSprayAerosolRate
        kd_sediment = -1.0_dp
        sediment_organic_carbon = me%sedimentOrganicCarbonConstant
        reaction_rate_sediment = -1.0_dp
        diffusivity = 0.0_dp
        sediment_water_exchange_velocity = me%sedimentWaterExchangeVelocity
        kd_salinity_coefficient = 0.0_dp
        kia_salinity_coefficient = 0.0_dp
        biota_transformation_rate = 0.0_dp
        estuary_salinity = me%estuarySalinityConstant
        runoff_mixing_depth = me%runoffMixingDepth
        groundwater_residence_time = me%groundwaterResidenceTimeConstant

        ! Open and read the NML file
        open(iouConstants, file=constantsFile, status="old", iostat=nmlIOStat)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance( &
                code=200, message="Failed to open constants file: " // trim(constantsFile)))
            return
        end if
        read(iouConstants, nml=allocatable_array_sizes, iostat=nmlIOStat)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance( &
                code=200, message="Failed to read allocatable_array_sizes namelist"))
            close(iouConstants)
            return
        end if
        rewind(iouConstants)

        ! Allocate arrays
        allocate(default_contaminant_size_distribution(n_default_contaminant_size_distribution), &
                default_spm_size_distribution(n_default_spm_size_distribution), &
                default_matrixembedded_distribution_to_spm(n_default_matrixembedded_distribution_to_spm), &
                vertical_distribution(n_vertical_distribution), &
                porosity(n_porosity), &
                initial_mass(n_initial_mass), &
                fractional_composition_distribution(n_fractional_composition_distribution), &
                spm_density_by_size_class(n_spm_density_by_size_class), &
                contaminant_size_classes(n_contaminant_size_classes), &
                default_contaminant_form_distribution(n_default_contaminant_form_distribution))

        ! Allocate class variables
        if (.not. allocated(me%defaultDistributionContaminant)) then
            allocate(me%defaultDistributionContaminant(n_default_contaminant_size_distribution))
        end if
        if (.not. allocated(me%defaultSpmSizeDistribution)) then
            allocate(me%defaultSpmSizeDistribution(n_default_spm_size_distribution))
        end if
        if (.not. allocated(me%defaultMatrixEmbeddedDistributionToSpm)) then
            allocate(me%defaultMatrixEmbeddedDistributionToSpm(n_default_matrixembedded_distribution_to_spm))
        end if
        if (.not. allocated(me%earthwormVerticalDistribution)) then
            allocate(me%earthwormVerticalDistribution(n_vertical_distribution))
        end if
        if (.not. allocated(me%sedimentPorosity)) then
            allocate(me%sedimentPorosity(n_porosity))
        end if
        if (.not. allocated(me%sedimentInitialMass)) then
            allocate(me%sedimentInitialMass(n_initial_mass))
        end if
        if (.not. allocated(me%sedimentFractionalComposition)) then
            allocate(me%sedimentFractionalComposition(n_fractional_composition_distribution))
        end if
        if (.not. allocated(me%contaminantSizeClasses)) then
            allocate(me%contaminantSizeClasses(n_contaminant_size_classes))
        end if
        if (.not. allocated(me%defaultContaminantFormDistribution)) then
            allocate(me%defaultContaminantFormDistribution(n_default_contaminant_form_distribution))
        end if

        ! Read namelists
        read(iouConstants, nml=n_biota_grp, iostat=nmlIOStat)
        rewind(iouConstants)
        me%nBiota = n_biota
        if (nmlIOStat .ge. 0) then
            ! Validate allocation sizes
            if (n_k_uptake_contaminant /= me%nBiota * n_contaminant_size_classes) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in k_uptake_contaminant size: expected " // &
                            trim(str(me%nBiota * n_contaminant_size_classes)) // &
                            ", got " // trim(str(n_k_uptake_contaminant))))
                close(iouConstants)
                return
            end if
            if (n_k_elim_contaminant /= me%nBiota * n_contaminant_size_classes) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in k_elim_contaminant size: expected " // &
                            trim(str(me%nBiota * n_contaminant_size_classes)) // &
                            ", got " // trim(str(n_k_elim_contaminant))))
                close(iouConstants)
                return
            end if
            if (n_initial_c_org /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in initial_C_org size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_initial_c_org))))
                close(iouConstants)
                return
            end if
            if (n_k_growth /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in k_growth size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_k_growth))))
                close(iouConstants)
                return
            end if
            if (n_name /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in name size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_name))))
                close(iouConstants)
                return
            end if
            if (n_stored_fraction /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in stored_fraction size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_stored_fraction))))
                close(iouConstants)
                return
            end if
            if (n_compartment /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in compartment size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_compartment))))
                close(iouConstants)
                return
            end if
            if (n_k_uptake_dissolved /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in k_uptake_dissolved size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_k_uptake_dissolved))))
                close(iouConstants)
                return
            end if
            if (n_k_elim_dissolved /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in k_elim_dissolved size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_k_elim_dissolved))))
                close(iouConstants)
                return
            end if
            if (n_uptake_from_form /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in uptake_from_form size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_uptake_from_form))))
                close(iouConstants)
                return
            end if
            if (n_harvest_in_month /= me%nBiota) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=900, &
                    message="Mismatch in harvest_in_month size: expected " // &
                            trim(str(me%nBiota)) // ", got " // trim(str(n_harvest_in_month))))
                close(iouConstants)
                return
            end if
            rewind(iouConstants)

            ! Allocate biota-related arrays
            allocate(initial_C_org(me%nBiota), k_growth(me%nBiota), &
                    k_uptake_contaminant(me%nBiota * n_contaminant_size_classes), &
                    k_elim_contaminant(me%nBiota * n_contaminant_size_classes), &
                    name(me%nBiota), stored_fraction(me%nBiota), &
                    compartment(me%nBiota), k_uptake_dissolved(me%nBiota), &
                    k_elim_dissolved(me%nBiota), uptake_from_form(me%nBiota), &
                    harvest_in_month(me%nBiota), prey(me%nBiota), k_diet(me%nBiota), k_death(me%nBiota))
            ! Optional food web and loss terms
            prey = 0
            k_diet = 0.0_dp
            k_death = 0.0_dp

            ! Read biota namelist
            read(iouConstants, nml=biota, iostat=nmlIOStat)
            if (nmlIOStat /= 0) then
                call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read biota namelist"))
                close(iouConstants)
                return
            end if
            me%hasBiota = .true.
        end if
        rewind(iouConstants)

        ! Read other namelists
        read(iouConstants, nml=contaminant, iostat=nmlIOStat)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read contaminant namelist, using defaults"))
        end if
        rewind(iouConstants)

        read(iouConstants, nml=earthworm_densities, iostat=nmlIOStat, iomsg=nmlIOMsg)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read earthworm_densities namelist" &
                                                                             // " with message: " // trim(nmlIOMsg)))
            close(iouConstants)
            return
        end if
        rewind(iouConstants)

        read(iouConstants, nml=soil, iostat=nmlIOStat, iomsg=nmlIOMsg)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read soil namelist" &
                                                                             // " with message: " // trim(nmlIOMsg)))
            close(iouConstants)
            return
        end if
        rewind(iouConstants)

        read(iouConstants, nml=water, iostat=nmlIOStat, iomsg=nmlIOMsg)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read water namelist" &
                                                                             // " with message: " // trim(nmlIOMsg)))
            close(iouConstants)
            return
        end if
        rewind(iouConstants)

        read(iouConstants, nml=sediment, iostat=nmlIOStat, iomsg=nmlIOMsg)
        if (nmlIOStat /= 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read sediment namelist" &
                                                                             // " with message: " // trim(nmlIOMsg)))
            close(iouConstants)
            return
        end if
        rewind(iouConstants)

        ! The &lake group is optional; defaults above are kept if it is absent
        read(iouConstants, nml=lake, iostat=nmlIOStat)
        rewind(iouConstants)

        ! The &pfas group is optional (defaults above), but a malformed one is an error
        read(iouConstants, nml=pfas, iostat=nmlIOStat, iomsg=nmlIOMsg)
        if (nmlIOStat > 0) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(code=200, message="Failed to read pfas namelist" &
                // new_line('A') // trim(nmlIOMsg)))
        else if (nmlIOStat < 0) then
            call LOGR%add("DataInput: No &pfas group in the constants file; PFAS partitioning and " // &
                "transformation are disabled")
        end if
        close(iouConstants)

        ! Save to class variables
        me%contaminantDensity = contaminant_density
        me%contaminantSizeClasses = contaminant_size_classes
        if (size(default_contaminant_size_distribution) /= n_contaminant_size_classes) then
            call ERROR_HANDLER%trigger(error=ErrorInstance( &
                code=900, message="Mismatch in default_contaminant_size_distribution size: expected " // &
                trim(str(n_contaminant_size_classes)) // ", got " // &
                trim(str(size(default_contaminant_size_distribution)))))
            return
        end if
        me%defaultDistributionContaminant = default_contaminant_size_distribution / 100.0
        if (size(default_contaminant_form_distribution) /= C%contaminantDim(2)+1) then
            call ERROR_HANDLER%trigger(error=ErrorInstance( &
                code=900, message="Mismatch in default_contaminant_form_distribution size: expected " // &
                trim(str(C%contaminantDim(2)+1)) // ", got " // &
                trim(str(size(default_contaminant_form_distribution)))))
            return
        end if
        me%defaultContaminantFormDistribution = default_contaminant_form_distribution / 100.0
        me%nContaminantSizeClasses = n_contaminant_size_classes
        me%defaultSpmSizeDistribution = default_spm_size_distribution / 100.0
        me%defaultMatrixEmbeddedDistributionToSpm = default_matrixembedded_distribution_to_spm / 100.0
        me%soilDarcyVelocity = darcy_velocity
        me%soilDefaultPorosity = default_porosity
        me%soilHamakerConstant = hamaker_constant
        me%soilParticleDensity = particle_density
        me%soilAttachmentEfficiencyConstant = soil_attachment_efficiency
        me%soilErosivity_a1 = erosivity_a1
        me%soilErosivity_a2 = erosivity_a2
        me%soilErosivity_a3 = erosivity_a3
        me%soilErosivity_b = erosivity_b
        me%sedimentTransport_aConstant = sediment_transport_a
        me%sedimentTransport_bConstant = sediment_transport_b
        me%sedimentTransport_cConstant = sediment_transport_c
        me%earthwormDensityArable = arable
        me%earthwormDensityConiferous = coniferous
        me%earthwormDensityDeciduous = deciduous
        me%earthwormDensityGrassland = grassland
        me%earthwormDensityHeathland = heathland
        me%earthwormDensityUrbanCapped = urban_capped
        me%earthwormDensityUrbanGardens = urban_gardens
        me%earthwormDensityUrbanParks = urban_parks
        me%earthwormVerticalDistribution = vertical_distribution / 100.0
        if (me%hasBiota) then
            me%biotaName = name
            me%biotaInitial_C_org = initial_C_org
            me%biota_k_growth = k_growth
            me%biotaStoredFraction = stored_fraction
            me%biotaCompartment = compartment
            me%biotaUptakeFromForm = uptake_from_form
            me%biotaHarvestInMonth = harvest_in_month
            me%biotaPrey = prey
            where (me%biotaPrey < 0 .or. me%biotaPrey > me%nBiota) me%biotaPrey = 0
            me%biota_k_diet = max(0.0_dp, k_diet)
            me%biota_k_death = max(0.0_dp, k_death)
            if (allocated(me%biota_k_uptake_contaminant)) deallocate(me%biota_k_uptake_contaminant)
            if (allocated(me%biota_k_elim_contaminant)) deallocate(me%biota_k_elim_contaminant)
            allocate(me%biota_k_uptake_contaminant(me%nBiota, C%contaminantDim(2)))
            allocate(me%biota_k_elim_contaminant(me%nBiota, C%contaminantDim(2)))
            me%biota_k_uptake_contaminant(:,1) = k_uptake_contaminant(1:me%nBiota)
            me%biota_k_elim_contaminant(:,1) = k_elim_contaminant(1:me%nBiota)
            if (C%contaminantDim(2) > 1) then
                me%biota_k_uptake_contaminant(:,2) = k_uptake_contaminant(me%nBiota+1:2*me%nBiota)
                me%biota_k_elim_contaminant(:,2) = k_elim_contaminant(me%nBiota+1:2*me%nBiota)
            end if
            if (.not. allocated(me%biota_k_uptake_dissolved)) then
                allocate(me%biota_k_uptake_dissolved(me%nBiota))
                me%biota_k_uptake_dissolved = k_uptake_dissolved
            end if
            if (.not. allocated(me%biota_k_elim_dissolved)) then
                allocate(me%biota_k_elim_dissolved(me%nBiota))
                me%biota_k_elim_dissolved = k_elim_dissolved
            end if
        end if
        me%riverMeanderingFactor = river_meandering_factor
        me%waterResuspensionAlpha = resuspension_alpha
        me%waterResuspensionBeta = resuspension_beta
        me%depositionAlphaConstant = deposition_alpha
        me%depositionBetaConstant = deposition_beta
        me%bankErosionAlphaConstant = bank_erosion_alpha
        me%bankErosionBetaConstant = bank_erosion_beta
        me%contaminant_k_diss_pristine = k_diss_pristine
        me%contaminant_k_diss_transformed = k_diss_transformed
        me%contaminant_k_transform_pristine = k_transform_pristine
        me%riverAttachmentEfficiency = river_attachment_efficiency
        me%estuaryAttachmentEfficiency = estuary_attachment_efficiency
        me%waterResuspensionAlphaEstuary = merge(resuspension_alpha_estuary, me%waterResuspensionAlpha, &
                                                resuspension_alpha_estuary /= 0.0_dp)
        me%waterResuspensionBetaEstuary = merge(resuspension_beta_estuary, me%waterResuspensionBeta, &
                                                resuspension_beta_estuary /= 0.0_dp)
        me%shearRate = shear_rate
        me%waterTemperature = me%calculateWaterTemperatureTimeSeries(min_water_temperature, &
                                                                    max_water_temperature, &
                                                                    min_water_temperature_day_of_year)
        me%estuaryTidalM2 = estuary_tidal_M2
        me%estuaryTidalS2 = estuary_tidal_S2
        me%estuaryMeanDepthExpA = estuary_mean_depth_expA
        me%estuaryMeanDepthExpB = estuary_mean_depth_expB
        me%estuaryWidthExpA = estuary_width_expA
        me%estuaryWidthExpB = estuary_width_expb
        me%estuaryMeanderingFactor = estuary_meandering_factor
        me%estuaryMouthCoords = estuary_mouth_coords
        me%sedimentInitialMass = initial_mass
        me%sedimentPorosity = porosity
        me%sedimentFractionalComposition = fractional_composition_distribution
        me%sedimentEnrichment_k = sediment_enrichment_k
        me%sedimentEnrichment_a = sediment_enrichment_a
        me%spmDensityBySizeClass = spm_density_by_size_class
        me%lakeDefaultDepth = lake_default_depth
        me%lakeDefaultAreaFraction = lake_default_area_fraction
        me%lakeStratificationRadiationThreshold = lake_stratification_radiation_threshold
        me%lakeStratificationWindThreshold = lake_stratification_wind_threshold
        me%lakeStratificationTemperatureThreshold = lake_stratification_temperature_threshold
        me%lakeEpilimnionBaseDepth = lake_epilimnion_base_depth
        me%lakeEpilimnionWindCoefficient = lake_epilimnion_wind_coefficient
        me%lakeEpilimnionMinDepth = lake_epilimnion_min_depth
        me%lakeInterlayerMixingVelocity = lake_interlayer_mixing_velocity
        me%lakeWindResuspensionAttenuationDepth = lake_wind_resuspension_attenuation_depth
        me%reservoirDefaultReleaseRate = reservoir_default_release_rate

        ! PFAS. Partition coefficients are converted from L kg-1 to m3 kg-1
        me%pfasSpeciesNames = species_names
        me%pfasMolecularWeight = molecular_weight
        me%pfasCharge = charge
        me%pfasPka = pka
        me%pfasKdSoil = max(0.0_dp, kd_soil) * 1.0e-3_dp
        me%pfasKdSpm = max(0.0_dp, kd_spm) * 1.0e-3_dp
        me%pfasKoc = max(0.0_dp, koc) * 1.0e-3_dp
        me%pfasKdPhSlope = kd_ph_slope
        me%pfasKdPhReference = kd_ph_reference
        me%pfasKia = max(0.0_dp, kia)
        me%pfasKMtSoil = max(0.0_dp, k_mt_soil)
        me%pfasKMtSpm = max(0.0_dp, k_mt_spm)
        me%pfasKMtAwi = max(0.0_dp, k_mt_awi)
        me%pfasReactionRate = max(0.0_dp, reaction_rate)
        me%pfasReactionYield = max(0.0_dp, reaction_yield)
        me%soilOrganicCarbonConstant = soil_organic_carbon
        me%soilPhConstant = soil_ph
        me%spmOrganicCarbonConstant = spm_organic_carbon
        me%waterPhConstant = water_ph
        me%soilAwiSpecificAreaConstant = soil_awi_specific_area
        me%soilAwiGeometricFactor = awi_geometric_factor
        me%pfasFoamCoefficient = foam_coefficient
        me%pfasVolatilisationRateScalar = volatilisation_rate
        me%pfasSeaSprayAerosolRate = sea_spray_aerosol_rate
        ! Negative kd_sediment = not given, so the bed uses the SPM Kd
        me%pfasKdSediment = merge(kd_sediment * 1.0e-3_dp, -1.0_dp, kd_sediment >= 0.0_dp)
        me%sedimentOrganicCarbonConstant = sediment_organic_carbon
        ! If any reaction_rate_sediment element was given, the matrix is the bed's own
        ! (unset elements = 0); otherwise the bed uses reaction_rate
        if (any(reaction_rate_sediment >= 0.0_dp)) then
            me%pfasReactionRateSediment = max(0.0_dp, reaction_rate_sediment)
        else
            me%pfasReactionRateSediment = me%pfasReactionRate
        end if
        me%pfasDiffusivity = max(0.0_dp, diffusivity)
        me%sedimentWaterExchangeVelocity = sediment_water_exchange_velocity
        me%pfasKdSalinityCoefficient = kd_salinity_coefficient
        me%pfasKiaSalinityCoefficient = kia_salinity_coefficient
        me%pfasBiotaTransformationRate = max(0.0_dp, biota_transformation_rate)
        me%estuarySalinityConstant = max(0.0_dp, estuary_salinity)
        me%runoffMixingDepth = runoff_mixing_depth
        me%groundwaterResidenceTimeConstant = groundwater_residence_time

    end subroutine

    !> Replace NetCDF fill values with zero. Emissions variables use the fill value to
    !! mean "nothing emitted here" (e.g. cells with no point source, or a point source
    !! slot that isn't occupied in this cell), so it must be zeroed before the value is
    !! used in arithmetic. 
    elemental function zeroFillValue(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (x == nf90_fill_double .or. x == real(nf90_fill_real, dp)) then
            y = 0.0_dp
        else
            y = x
        end if
    end function

    !> Elemental function for getting a mask from an int2 array, where the NetCDF
    !! fill value nf90_fill_int2 is used to mask values
    elemental function maskDatabase(me, int) result(mask)
        class(Database), intent(in) :: me
        integer, intent(in) :: int
        logical :: mask
        if (int == nf90_fill_int2) then
            mask = .true.
        else
            mask = .false.
        end if
    end function

    ! Compute number of point sources per (x,y) cell by inspecting coordinates.
    subroutine calculateNPointSourcesDatabase(me, maxPointSources)
        class(Database) :: me
        integer, intent(in) :: maxPointSources
        integer :: i, j, p
        real(dp) :: px, py

        if (allocated(me%nPointSources)) deallocate(me%nPointSources)
        allocate(me%nPointSources(me%gridShape(1), me%gridShape(2)))
        me%nPointSources = 0

        if (.not. allocated(me%emissionsPointWaterCoords)) return
        if (maxPointSources <= 0) return

        do j = 1, me%gridShape(2)
            do i = 1, me%gridShape(1)
                do p = 1, maxPointSources
                    px = me%emissionsPointWaterCoords(i, j, p, 1)
                    py = me%emissionsPointWaterCoords(i, j, p, 2)
                    if (px /= nf90_fill_double .and. py /= nf90_fill_double) then
                        if ((abs(px) > C%epsilon .or. abs(py) > C%epsilon) .and. &
                            .not. (px < -9.9e8_dp .and. py < -9.9e8_dp)) then
                            me%nPointSources(i, j) = me%nPointSources(i, j) + 1
                        end if
                    end if
                end do
            end do
        end do
    end subroutine calculateNPointSourcesDatabase


    !> Check whether a set of coordinates (x,y) is in the model domain
    function inModelDomainDatabase(me, x, y) result(inModelDomain)
        class(Database), intent(in) :: me
        integer, intent(in) :: x, y
        logical :: inModelDomain
        logical :: xInDomain
        logical :: yInDomain

        xInDomain = x .ge. 1 .and. x .le. me%gridShape(1)
        yInDomain = y .ge. 1 .and. y .le. me%gridShape(2)

        if (xInDomain .and. yInDomain) then
            if (.not. me%gridMask(x, y)) then
                inModelDomain = .true.
            else
                inModelDomain = .false.
            end if
        else
            inModelDomain = .false.
        end if
    end function

    !> Convert easting and northing coordinates to cell indicies
    function coordsToCellIndexDatabase(me, easts, norths) result(indicies)
        class(Database) :: me
        real :: easts
        real :: norths
        integer :: indicies(2)
        integer :: x, y
        x = (int(easts) - mod(int(easts), int(me%gridRes(1))) - int(me%gridBounds(1)))/int(me%gridRes(1)) + 1
        y = (int(me%gridBounds(4)) - (int(norths) - mod(int(norths), int(me%gridRes(2)))))/int(me%gridRes(2))
        indicies = [x, y]
    end function

    !> Convert easting and northing coordinates to fractional cell indicies,
    !! i.e. as a fraction of the position of the coords within the cell
    function coordsToFractionalCellIndexDatabase(me, easts, norths) result(fracIndicies)
        class(Database) :: me
        real :: easts
        real :: norths
        real :: fracIndicies(2)
        integer :: indicies(2)
        indicies = me%coordsToCellIndex(easts, norths)
        fracIndicies(1) = indicies(1) + mod(easts, me%gridRes(1)) / me%gridRes(1) 
        fracIndicies(2) = indicies(2) + 1 - mod(norths, me%gridRes(2)) / me%gridRes(2)
    end function

    !> Calculate the meandering factor from cell size, based on Fekete et al. 2001 (https://doi.org/10.1029/2001WR900024).
    !! This method is empirical and was used in GWAVA. It assumes square cells, so how we'll average grid size to
    !! account for that, in case the cells aren't square. This value is used if meandering factors aren't specified in data.
    function calculateMeanderingFactorFromCellSizeDatabase(me) result(f_m)
        class(Database) :: me
        real            :: f_m
        f_m = 1.024 - 0.077 * log(100 / ((me%gridRes(1) + me%gridRes(2)) / 2))
    end function

    !> Calculate the water temperature timeseries using min and max water temperatures, and the Julian day
    !! (day of year) on which the min water temperature occurs. This function using those parameters to construct
    !! a cosine function representing the time series:
    !! $$
    !!  T_t = 0.5 (T_\text{max} - T_\text{min}) \cos\( \frac{2 D_t \pi}{366} - D_\text{min} \) + 0.5 (T_\text{max} + T_\text{min}) 
    !! $$
    function calculateWaterTemperatureTimeSeriesWaterBody(me, minTemp, maxTemp, minTempDay) result(waterTemperature)
        class(Database) :: me
        real(dp)        :: minTemp, maxTemp
        integer         :: minTempDay
        real(dp)        :: waterTemperature(366)
        integer         :: i
        real(dp)        :: angle(366)

        ! angle = 2*pi*(day - day_min)/366
        do i = 1, 366
            angle(i) = 2.0_dp*C%pi * real(i - minTempDay, dp) / 366.0_dp
        end do

        waterTemperature = 0.5_dp*(maxTemp - minTemp) * cos(angle) + 0.5_dp*(maxTemp + minTemp)
    end function calculateWaterTemperatureTimeSeriesWaterBody


    !> Audit the database
    function auditDatabase(me) result(rslt)
        class(Database) :: me
        type(Result)    :: rslt
        integer         :: x, y, i
        integer         :: xi, yi
        integer         :: nx, ny
        logical         :: simulationMaskError

        simulationMaskError = .false.

        ! Is the simulation mask self-contained (no inflows to area to simulate)?
        if (C%hasSimulationMask) then
            nx = me%gridShape(1)
            ny = me%gridShape(2)
            do y = 1, ny
                do x = 1, nx
                    if (me%simulationMask(x, y)) then
                        ! me%inflows is (d,w,x,y) where d=2 holds (x,y) origin indices
                        do i = 1, size(me%inflows, dim=2)
                            xi = me%inflows(1, i, x, y)
                            yi = me%inflows(2, i, x, y)
                            ! Skip invalid/out-of-domain inflow indices
                            if (xi >= 1 .and. xi <= nx .and. yi >= 1 .and. yi <= ny) then
                                if (.not. me%simulationMask(xi, yi)) simulationMaskError = .true.
                            end if
                        end do
                    end if
                end do
            end do
        end if

        if (simulationMaskError) then
            call rslt%addError(ErrorInstance( &
                message="Simulation mask provided has inflows from outside the area to " // &
                        "simulate. Please provide a simulation mask that is self-contained." ))
        end if

        ! Bounds checks for sediment calibration parameters (arrays may be unallocated)
        if (allocated(me%depositionAlpha)) then
            if (any(me%depositionAlpha < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="Value provided for deposition_alpha must be >= 0. At least one < 0." ))
            end if
        end if
        if (allocated(me%resuspensionAlpha)) then
            if (any(me%resuspensionAlpha < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="Value provided for resuspension_alpha must be >= 0. At least one < 0." ))
            end if
        end if
        if (allocated(me%resuspensionBeta)) then
            if (any(me%resuspensionBeta < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="Value provided for resuspension_beta must be >= 0. At least one < 0." ))
            end if
        end if

        ! Does sediment fractional composition sum to unity?
        if (.not. isZero(1.0_dp - sum(me%sedimentFractionalComposition))) then
            call rslt%addError(ErrorInstance( &
                message="sedimentFractionalComposition must sum to 1. Found: " // &
                        str(sum(me%sedimentFractionalComposition)) ))
        end if

        if (any(me%defaultDistributionContaminant < 0.0)) then
            call rslt%addError(ErrorInstance( &
                message="defaultDistributionContaminant must be non-negative." ))
        end if
        if (.not. isZero(1.0_dp - sum(me%defaultDistributionContaminant))) then
            call rslt%addError(ErrorInstance( &
                message="defaultDistributionContaminant must sum to 1. Found: " // &
                        str(sum(me%defaultDistributionContaminant)) ))
        end if

        if (me%contaminantDensity <= 0.0) then
            call rslt%addError(ErrorInstance( &
                message="contaminantDensity must be positive." ))
        end if
        if (any(me%contaminantSizeClasses <= 0.0)) then
            call rslt%addError(ErrorInstance( &
                message="contaminantSizeClasses must be positive." ))
        end if

        ! Bounds checks for initial concentrations
        if (allocated(me%initialContaminantConcsSoil)) then
            if (any(me%initialContaminantConcsSoil < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="initialContaminantConcsSoil must be non-negative. At least one < 0." ))
            end if
        end if
        if (allocated(me%initialContaminantConcsWater)) then
            if (any(me%initialContaminantConcsWater < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="initialContaminantConcsWater must be non-negative. At least one < 0." ))
            end if
        end if
        if (allocated(me%initialContaminantConcsSediment)) then
            if (any(me%initialContaminantConcsSediment < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="initialContaminantConcsSediment must be non-negative. At least one < 0." ))
            end if
        end if
        if (allocated(me%initialDissolvedConcsSoil)) then
            if (any(me%initialDissolvedConcsSoil < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="initialDissolvedConcsSoil must be non-negative. At least one < 0." ))
            end if
        end if
        if (allocated(me%initialDissolvedConcsWater)) then
            if (any(me%initialDissolvedConcsWater < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="initialDissolvedConcsWater must be non-negative. At least one < 0." ))
            end if
        end if
        if (allocated(me%initialDissolvedConcsSediment)) then
            if (any(me%initialDissolvedConcsSediment < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="initialDissolvedConcsSediment must be non-negative. At least one < 0." ))
            end if
        end if

        ! Bounds checks for dissolved emissions
        if (allocated(me%emissionsArealSoilDissolvedContaminant)) then
            if (any(me%emissionsArealSoilDissolvedContaminant < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="emissionsArealSoilDissolvedContaminant must be >= 0. At least one < 0." ))
            end if
        end if
        if (allocated(me%emissionsArealWaterDissolvedContaminant)) then
            if (any(me%emissionsArealWaterDissolvedContaminant < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="emissionsArealWaterDissolvedContaminant must be >= 0. At least one < 0." ))
            end if
        end if
        if (allocated(me%emissionsAtmosphericDryDepoDissolvedContaminant)) then
            if (any(me%emissionsAtmosphericDryDepoDissolvedContaminant < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="emissionsAtmosphericDryDepoDissolvedContaminant must be >= 0. " // &
                            "At least one < 0." ))
            end if
        end if
        if (allocated(me%emissionsAtmosphericWetDepoDissolvedContaminant)) then
            if (any(me%emissionsAtmosphericWetDepoDissolvedContaminant < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="emissionsAtmosphericWetDepoDissolvedContaminant must be >= 0. " // &
                            "At least one < 0." ))
            end if
        end if
        if (allocated(me%emissionsPointWaterDissolvedContaminant)) then
            if (any(me%emissionsPointWaterDissolvedContaminant < 0.0_dp)) then
                call rslt%addError(ErrorInstance( &
                    message="emissionsPointWaterDissolvedContaminant must be >= 0. At least one < 0." ))
            end if
        end if
    end function

end module


! AUDITING TO DO:
!   array size checks, particularly variables for size classes of NM and SPM
!   sedimentPorosity: 0 <= x <= 1
!   water content at saturation is greater than water content at field capacity
!   are batch runs contiguous (dates follow on from each other for different chunks)