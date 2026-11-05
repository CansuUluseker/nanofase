# Input data

The model requires a NetCDF file of spatio(temporal) input data, and a Fortran namelist file for constants. The config file is responsible for telling the model where these input data are located (via the `&data` group). There are two options for generating these input data files:
* *Recommended*: Use the [NanoFASE data module](https://github.com/nerc-ceh/nanofase-data) - see [](nanofase-data). This is a Python library that compiles multiple spatio(temporal) input files into the main NetCDF file, at the same time as deriving secondary variables from these data, which are required by the model.
* Manually creating (or editing a copy of) the NetCDF and namelist files. This is not recommended for compiling new data (though can be useful for editing existing data) due to the requirement for a range of [secondary derived variables](netcdf-namelist-input:secondary-derived-variables) that would have to be calculated manually. If you wish to go down this route, see [](netcdf-namelist-input) for details of the variables required in the NetCDF and constants namelist file.

The complete list of parameters required by the model is givin in [](parameter-reference). Example NetCDF and constants files are given in the [data.example](https://github.com/NERC-CEH/nanofase/tree/develop/data.example) directory.

% TODO draw diagram of data input components

(input-data:hydrology)=
## Hydrological input data

The model routes runoff provided in the input data, rather than generating it. How these inputs are used is described in [](../theory/hydrology.md). This section summarises the hydrological inputs and the options that control how they are read. All water fluxes are depths per timestep (m/timestep), so they must be recalculated if the timestep length changes.

### Runoff

The model needs the total runoff and the surface runoff (quickflow) for each grid cell and timestep, as (t, y, x) variables in the NetCDF file. Datasets provide different combinations of runoff components, so the model can read any of total, surface and sub-surface runoff, and derive the missing components (see [](hydrology:runoff-inputs) for the combinations supported). Which NetCDF variables hold each component is set in the optional `&hydrology` group of the model config file:

`total_runoff_variable`
: *Optional, defaults to `'runoff'`.* \
Name of the NetCDF variable containing total (surface plus sub-surface) runoff [m/timestep].

`surface_runoff_variable`
: *Optional, defaults to `'quickflow'`.* \
Name of the NetCDF variable containing surface (overland) runoff [m/timestep].

`subsurface_runoff_variable`
: *Optional, defaults to `''` (not provided).* \
Name of the NetCDF variable containing sub-surface runoff [m/timestep].

`surface_runoff_fraction`
: *Optional, defaults to 0. Must be between 0 and 1.* \
Fraction of total runoff treated as surface runoff, used only when total runoff is the only component found in the input data.

Set a variable name to `''` if the dataset doesn't provide that component. At least one of `total_runoff_variable` and `surface_runoff_variable` must be named. A component is also treated as not provided if the named variable isn't in the NetCDF file. If the group is left out, the defaults are used, so the model reads total runoff from `runoff` and surface runoff from `quickflow`. For example, for a dataset that provides total and sub-surface runoff:

```fortran
&hydrology
    total_runoff_variable = 'runoff'
    surface_runoff_variable = ''
    subsurface_runoff_variable = 'subsurface_runoff'
/
```

### Precipitation and evapotranspiration

`precip(t, y, x)` and `evap(t, y, x)` [m/timestep] drive percolation through the soil profile, and precipitation also drives soil erosion. Both are optional and default to zero. See [](netcdf-namelist-input) for details.

### River network and topography

The surface water network is given by the `outflow`, `inflows`, `n_waterbodies` and `is_headwater` variables (see [](netcdf-namelist-input:secondary-derived-variables)). The optional `dem` variable is used to calculate the slope of each reach, which is floored at the `min_stream_slope` option in the `&water` group of the config file (default 0.0001 m/m). If `dem` isn't provided, a slope of 0.0005 m/m is used for every reach.

(input-data:groundwater)=
### Groundwater store

The optional groundwater store delays the sub-surface part of the runoff before it reaches the water bodies (see [](hydrology:groundwater-store)). It is configured in the optional `&groundwater` group of the constants namelist file, and is switched off if the group is left out:

```fortran
&groundwater
    residence_time = 180.0              ! [days], <= 0 switches the store off (default -1)
    subsurface_recharge_fraction = 1.0  ! [-], default 1.0
    surface_recharge_fraction = 0.0     ! [-], default 0.0
    loss_fraction = 0.0                 ! [-], default 0.0
    initial_storage = -1.0              ! [m], < 0 starts at steady state (default -1)
/
```

`residence_time`
: *Optional, defaults to -1 (store off). Units: days.* \
Residence time of water in the store. The store is switched off where this is zero or negative, and the runoff passes to the water bodies unchanged.

`subsurface_recharge_fraction`
: *Optional, defaults to 1. Unitless, between 0 and 1.* \
Fraction of the sub-surface runoff (total minus surface runoff) that recharges the store. The rest goes directly to the water bodies.

`surface_recharge_fraction`
: *Optional, defaults to 0. Unitless, between 0 and 1.* \
Fraction of the surface runoff that recharges the store. The rest goes directly to the water bodies. This does not change the surface runoff used for soil erosion delivery and contaminant wash-off.

`loss_fraction`
: *Optional, defaults to 0. Unitless, between 0 and 1.* \
Fraction of the recharge that leaves the model without reaching the water bodies, for example to a deep aquifer or by abstraction.

`initial_storage`
: *Optional, defaults to -1 (steady state). Units: m.* \
Water in the store at the start of the run, as a depth over the grid cell. If negative, the store starts at steady state with the mean recharge, so that baseflow starts at the mean recharge rate.

The fractions are limited to between 0 and 1. All four parameters can also vary in space, using the following optional (y, x) variables in the NetCDF file:

`float gw_residence_time(y, x)`
: *Optional, defaults to `residence_time` in the constants file. Units: days.* \
Residence time of water in each grid cell's store.

`float gw_loss_fraction(y, x)`
: *Optional, defaults to `loss_fraction` in the constants file. Unitless.* \
Fraction of each grid cell's recharge that leaves the model.

`float gw_subsurface_recharge_fraction(y, x)`
: *Optional, defaults to `subsurface_recharge_fraction` in the constants file. Unitless.* \
Fraction of each grid cell's sub-surface runoff that recharges the store.

`float gw_surface_recharge_fraction(y, x)`
: *Optional, defaults to `surface_recharge_fraction` in the constants file. Unitless.* \
Fraction of each grid cell's surface runoff that recharges the store.

Where these variables are missing or a fill value, the constant from the `&groundwater` group is used. The same applies to residence times that are zero or negative, and to negative fractions.

```{note}
The `groundwater_residence_time` parameter in the `&pfas` constants group (and the `groundwater_residence_time` NetCDF variable) is not part of this water store. It sets the residence time of a separate reservoir that holds only the PFAS leached from the soil profile.
```

(input-data:routing)=
### Channel routing

Muskingum-Cunge channel routing (see [](../theory/hydrology.md)) is off by default. It is switched on, and its sub-timestep set, in the `&water` group of the model config file:

```fortran
&water
    include_routing = .true.        ! Route water through river reaches at a sub-timestep? Default .false.
    routing_timestep = 3600         ! Routing sub-step [s]; the timestep is split into a whole number of these. Default 3600
/
```

The routing parameters are physical, and set in the `&water` group of the constants file:

`river_manning_n`
: *Optional, defaults to 0.035. Units: s/m<sup>1/3</sup>.* \
Manning's roughness coefficient of river channels. It sets the river depth (with or without routing) and, through the flow velocity, the speed at which flood waves travel. It can also vary in space, using the optional `float river_manning_n(y, x)` NetCDF variable; where that is missing, zero or negative, the constant is used.

`routing_celerity_factor`
: *Optional, defaults to 1.0. Unitless.* \
Multiplier on the kinematic wave celerity used for routing. Values below 1 slow flood waves down, values above 1 speed them up. It is a calibration parameter, to be left at 1 unless routed flows can be compared with gauged flows.

The routing state (the flows at the end of the last sub-step and the water held in each reach) is saved in checkpoints (see [](checkpointing.md)).
