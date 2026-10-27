!> Spatially explicit (per grid cell) PFAS partitioning parameters, derived from
!! the soil and surface water properties in DATASET.
!!
!! Every parameter is resolved per cell with the same order of precedence:
!!   1. a value from the NetCDF input map, where the map exists and the cell has a value;
!!   2. a value computed by the model from soil/water properties, where the parameters
!!      that computation needs have been supplied;
!!   3. the constant from the &pfas group of the constants file.
!! Soil/water properties used by step 2 (organic carbon, pH) follow the same map-then-
!! constant fallback. Steps 1 and 3 are static; the air-water interfacial area also
!! responds to the soil water content on every time step (pfas_soil_awi_area).
module PFASParameterisationModule
    use GlobalsModule, only: dp, C
    use DataInputModule, only: DATASET
    implicit none
    private

    public :: pfas_kd_soil
    public :: pfas_kd_spm
    public :: pfas_kd_sediment
    public :: pfas_water_ph
    public :: pfas_water_salinity
    public :: pfas_kia_water
    public :: pfas_groundwater_residence_time
    public :: pfas_soil_awi_specific_area
    public :: pfas_soil_awi_area

contains

    !> Soil (and bed) Kd for each species in cell (x,y) [m3 kg-1].
    !! Computed route: Kd = Koc * foc * 10**(slope * (pH - pH_ref)), with foc and pH
    !! from the 'soil_organic_carbon'/'soil_ph' maps, else the &pfas constants.
    function pfas_kd_soil(x, y) result(kd)
        integer, intent(in) :: x, y
        real(dp) :: kd(C%contaminantDim(1))
        kd = resolve_kd(x, y, DATASET%pfasKdSoilMap, species_array(DATASET%pfasKdSoil), &
                        cell_value(DATASET%soilOrganicCarbon, x, y, DATASET%soilOrganicCarbonConstant), &
                        cell_value(DATASET%soilPh, x, y, DATASET%soilPhConstant))
    end function

    !> Suspended particulate matter Kd for each species in cell (x,y) [m3 kg-1].
    !! Computed route as for soil, with foc and pH from 'spm_organic_carbon'/'water_ph'.
    !! Kd increases with salinity (salting out): the result is scaled by
    !! 10**(kd_salinity_coefficient * salinity), salinity from pfas_water_salinity.
    function pfas_kd_spm(x, y) result(kd)
        integer, intent(in) :: x, y
        real(dp) :: kd(C%contaminantDim(1))
        kd = kd_spm_fresh(x, y) * salinity_factor(DATASET%pfasKdSalinityCoefficient, pfas_water_salinity(x, y))
    end function

    !> SPM Kd before the salinity adjustment
    function kd_spm_fresh(x, y) result(kd)
        integer, intent(in) :: x, y
        real(dp) :: kd(C%contaminantDim(1))
        kd = resolve_kd(x, y, DATASET%pfasKdSpmMap, species_array(DATASET%pfasKdSpm), &
                        cell_value(DATASET%spmOrganicCarbon, x, y, DATASET%spmOrganicCarbonConstant), &
                        cell_value(DATASET%waterPh, x, y, DATASET%waterPhConstant))
    end function

    !> Bed sediment Kd for each species in cell (x,y) [m3 kg-1]. Map 'pfas_kd_sediment',
    !! else Koc * foc from 'sediment_organic_carbon' (map, else constant) and the water
    !! pH, else the constant kd_sediment, else - as bed sediment is deposited SPM - the
    !! SPM Kd of the cell; then scaled for salinity as for SPM.
    function pfas_kd_sediment(x, y) result(kd)
        integer, intent(in) :: x, y
        real(dp) :: kd(C%contaminantDim(1))
        real(dp) :: fallback(C%contaminantDim(1))
        fallback = kd_spm_fresh(x, y)
        if (allocated(DATASET%pfasKdSediment)) then
            where (DATASET%pfasKdSediment(:size(fallback)) >= 0.0_dp) &
                fallback = DATASET%pfasKdSediment(:size(fallback))
        end if
        kd = resolve_kd(x, y, DATASET%pfasKdSedimentMap, fallback, &
                        cell_value(DATASET%sedimentOrganicCarbon, x, y, DATASET%sedimentOrganicCarbonConstant), &
                        pfas_water_ph(x, y))
        kd = kd * salinity_factor(DATASET%pfasKdSalinityCoefficient, pfas_water_salinity(x, y))
    end function

    !> Surface water salinity in cell (x,y) [psu]: the 'water_salinity' map, else the
    !! pfas estuary_salinity in estuary cells, else 0 (fresh water). 
    function pfas_water_salinity(x, y) result(salinity)
        integer, intent(in) :: x, y
        real(dp) :: salinity
        logical :: estuary
        estuary = .false.
        if (allocated(DATASET%isEstuary)) then
            if (in_grid(x, y, size(DATASET%isEstuary, 1), size(DATASET%isEstuary, 2))) estuary = DATASET%isEstuary(x, y)
        end if
        salinity = cell_value(DATASET%waterSalinity, x, y, merge(DATASET%estuarySalinityConstant, 0.0_dp, estuary))
    end function

    !> Air-water interfacial adsorption coefficient in the surface water of cell (x,y) [m].
    !! Interfacial adsorption increases with ionic strength, so Kia is scaled by
    !! 10**(kia_salinity_coefficient * salinity).
    function pfas_kia_water(x, y) result(kia)
        integer, intent(in) :: x, y
        real(dp) :: kia(C%contaminantDim(1))
        kia = species_array(DATASET%pfasKia) &
              * salinity_factor(DATASET%pfasKiaSalinityCoefficient, pfas_water_salinity(x, y))
    end function

    !> Residence time of the groundwater reservoir below the soil in cell (x,y) [s]:
    !! the 'groundwater_residence_time' map [days], else the constant [days]; <= 0 = none.
    function pfas_groundwater_residence_time(x, y) result(tau)
        integer, intent(in) :: x, y
        real(dp) :: tau
        tau = cell_value(DATASET%groundwaterResidenceTime, x, y, DATASET%groundwaterResidenceTimeConstant)
        tau = tau * 86400.0_dp
    end function

    !> Salting-out factor 10**(coefficient * salinity) per species (1 when not configured)
    function salinity_factor(coefficient, salinity) result(f)
        real(dp), allocatable, intent(in) :: coefficient(:)
        real(dp), intent(in) :: salinity
        real(dp) :: f(C%contaminantDim(1))
        f = 10.0_dp ** min(30.0_dp, species_array(coefficient) * max(0.0_dp, salinity))
    end function

    !> Surface water pH in cell (x,y), from the 'water_ph' map, else the constant.
    function pfas_water_ph(x, y) result(ph)
        integer, intent(in) :: x, y
        real(dp) :: ph
        ph = cell_value(DATASET%waterPh, x, y, DATASET%waterPhConstant)
    end function

    !> Air-water interfacial area per unit bulk soil volume when the soil is dry
    !! (zero water saturation) [m2 m-3].
    function pfas_soil_awi_specific_area(x, y, porosity, d_grain) result(a_aw)
        integer, intent(in)  :: x, y
        real(dp), intent(in) :: porosity        !! Soil porosity [-]
        real(dp), intent(in) :: d_grain         !! Mean grain diameter [m]
        real(dp) :: a_aw
        a_aw = cell_value(DATASET%soilAwiSpecificArea, x, y, -1.0_dp)
        if (a_aw >= 0.0_dp) return
        if (DATASET%soilAwiGeometricFactor > 0.0_dp .and. d_grain > C%epsilon) then
            a_aw = DATASET%soilAwiGeometricFactor * 6.0_dp * max(0.0_dp, 1.0_dp - porosity) / d_grain
        else
            a_aw = max(0.0_dp, DATASET%soilAwiSpecificAreaConstant)
        end if
    end function

    !> Air-water interfacial area in a soil layer on this time step [m2]. The area
    !! available for PFAS adsorption shrinks as the pores fill with water, so it is
    !! taken as linear in the air saturation: A_aw = a_aw * (1 - S_w) * V_bulk.
    pure function pfas_soil_awi_area(a_aw, V_w, V_sat, bulk_volume) result(area)
        real(dp), intent(in) :: a_aw            !! Specific area at zero saturation [m2 m-3]
        real(dp), intent(in) :: V_w, V_sat      !! Current and saturated water content [m3 m-2]
        real(dp), intent(in) :: bulk_volume     !! Bulk volume of the layer [m3]
        real(dp) :: area, S_w
        if (V_sat > 0.0_dp) then
            S_w = max(0.0_dp, min(1.0_dp, V_w / V_sat))
        else
            S_w = 1.0_dp
        end if
        area = max(0.0_dp, a_aw) * (1.0_dp - S_w) * max(0.0_dp, bulk_volume)
    end function

    !> Kd per species: direct map > Koc*foc (pH-adjusted) > fallback.
    function resolve_kd(x, y, kd_map, kd_fallback, foc, ph) result(kd)
        integer, intent(in)  :: x, y
        real(dp), allocatable, intent(in) :: kd_map(:,:,:)
        real(dp), intent(in) :: kd_fallback(:)
        real(dp), intent(in) :: foc, ph         !! Negative = not available
        real(dp) :: kd(C%contaminantDim(1))
        integer :: i
        kd = 0.0_dp
        do i = 1, size(kd)
            if (allocated(kd_map)) then
                if (in_grid(x, y, size(kd_map, 1), size(kd_map, 2)) .and. i <= size(kd_map, 3)) then
                    if (kd_map(x, y, i) >= 0.0_dp) then
                        kd(i) = kd_map(x, y, i)
                        cycle
                    end if
                end if
            end if
            if (allocated(DATASET%pfasKoc) .and. foc >= 0.0_dp) then
                if (i <= size(DATASET%pfasKoc)) then
                    if (DATASET%pfasKoc(i) > 0.0_dp) then
                        kd(i) = DATASET%pfasKoc(i) * min(1.0_dp, foc)
                        if (ph >= 0.0_dp) kd(i) = kd(i) * 10.0_dp ** &
                            (DATASET%pfasKdPhSlope(i) * (ph - DATASET%pfasKdPhReference))
                        cycle
                    end if
                end if
            end if
            if (i <= size(kd_fallback)) kd(i) = kd_fallback(i)
        end do
    end function

    !> A species-sized copy of an optional constants array (zeros where absent).
    function species_array(a) result(v)
        real(dp), allocatable, intent(in) :: a(:)
        real(dp) :: v(C%contaminantDim(1))
        integer :: n
        v = 0.0_dp
        if (.not. allocated(a)) return
        n = min(size(a), size(v))
        v(:n) = a(:n)
    end function

    !> Value of a (x,y) property map in a cell, or the fallback where the map is
    !! absent or has no value (negative) in that cell.
    function cell_value(map, x, y, fallback) result(v)
        real(dp), allocatable, intent(in) :: map(:,:)
        integer, intent(in)  :: x, y
        real(dp), intent(in) :: fallback
        real(dp) :: v
        v = fallback
        if (.not. allocated(map)) return
        if (.not. in_grid(x, y, size(map, 1), size(map, 2))) return
        if (map(x, y) >= 0.0_dp) v = map(x, y)
    end function

    pure logical function in_grid(x, y, nx, ny)
        integer, intent(in) :: x, y, nx, ny
        in_grid = x >= 1 .and. y >= 1 .and. x <= nx .and. y <= ny
    end function

end module PFASParameterisationModule
