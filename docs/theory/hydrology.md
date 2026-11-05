# Hydrology

The model does not generate runoff itself. Instead, it takes runoff from the input data, which might come from a separate hydrological model or be derived from observations, and routes this runoff through the surface water network to give flows in each water body. Optionally, part of the runoff can pass through a simple groundwater store before it reaches the water bodies, which delays and smooths the sub-surface contribution to river flow (see [](hydrology:groundwater-store)).

The hydrological inputs play three separate roles:

* **Total runoff** sets the volume of water entering the water bodies in each grid cell, and therefore river flow, depth, width and velocity.
* **Surface runoff** (also called *quickflow* or overland flow) drives the delivery of eroded soil to water bodies (via the sediment transport capacity) and the wash-off of dissolved contaminant from the soil surface.
* **Precipitation** and **evapotranspiration** drive the percolation of water down through the soil profile, which controls the leaching of contaminant. Precipitation also sets the rainfall erosivity used to calculate soil erosion.

These roles are not coupled to each other by a water balance within the model. For example, the water that percolates out of the bottom of the soil profile is not added to the runoff, and the runoff is not reduced by the water held in the soil. This means the runoff, surface runoff, precipitation and evapotranspiration inputs should be consistent with each other when they are prepared.

The figure below shows how runoff moves from the input data to the water bodies.

```{mermaid}
flowchart TD
    IN["Runoff input data<br/>(total, surface and/or sub-surface)"] --> AD["Hydrology input adapter"]
    AD --> QS["Surface runoff (quickflow)"]
    AD --> QB["Sub-surface runoff"]
    QS -- "surface_recharge_fraction" --> REC["Recharge"]
    QB -- "subsurface_recharge_fraction" --> REC
    QS -- "remainder" --> DIR["Direct runoff"]
    QB -- "remainder" --> DIR
    REC -- "loss_fraction" --> LOSS["Lost from the model"]
    REC -- "remainder" --> GW[("Groundwater store")]
    GW -- "baseflow" --> CR["Grid cell runoff"]
    DIR --> CR
    CR -- "shared by reach length" --> R["Water bodies in the grid cell"]
    UP["Upstream water bodies"] -- "inflow" --> R
    R -- "outflow" --> DOWN["Downstream water body"]
    QS -. "soil erosion delivery and<br/>dissolved contaminant wash-off" .-> R
```

When the groundwater store is switched off (the default), all of the total runoff is direct runoff and reaches the water bodies on the same timestep. All runoff inputs are depths of water per unit grid cell area per timestep (m/timestep). Missing (fill), NaN or negative runoff is treated as zero.

(hydrology:runoff-inputs)=
## Runoff inputs

The model needs two runoff components for every grid cell and timestep: the total runoff $q_T$ and the surface runoff $q_S$. The sub-surface runoff is the difference between them, $q_B = q_T - q_S$. Hydrological datasets provide different combinations of these components, so the model reads whichever are available and derives the rest. The names of the variables to read from the input NetCDF file are set in the `&hydrology` group of the config file (see [](input-data:hydrology)). The following combinations are supported, in order of preference:

```{list-table}
:header-rows: 1

* - Components provided
  - Total runoff $q_T$
  - Surface runoff $q_S$
* - Total and surface
  - As provided
  - As provided
* - Surface and sub-surface
  - $q_S + q_B$
  - As provided
* - Total and sub-surface
  - As provided
  - $\max(q_T - q_B, 0)$
* - Total only
  - As provided
  - $f_S \, q_T$, where $f_S$ is the `surface_runoff_fraction` config option (default 0)
* - Surface only
  - $q_S$ (a warning is given)
  - As provided
```

If all three components are provided, the total and surface runoff are used and the sub-surface runoff is ignored. Using surface runoff alone as total runoff is only suitable for simple tests, because it leaves out the sub-surface contribution to river flow. If no runoff component can be found, both components are set to zero and a warning is logged. The model will run, but there will be no flow in the water bodies.

Total runoff is the minimum requirement for a meaningful simulation. If only total runoff is available, the `surface_runoff_fraction` option gives a simple way to estimate the surface runoff that drives soil erosion delivery. With the default fraction of zero, no eroded soil reaches the water bodies.

(hydrology:soil-water)=
## Precipitation and evapotranspiration

On each timestep, the water entering the top of each (non-urban) soil profile is the precipitation $P$ minus the evapotranspiration $E$, floored at zero:

$$
q_{in} = \max(P - E, 0)
$$

This water percolates down through the soil layers, and water that exceeds a layer's capacity is pushed back up to the layer above. Water that percolates out of the bottom soil layer carries leached contaminant with it, but the water itself leaves the model: it is not added to the runoff. Precipitation also sets the rainfall kinetic energy used to calculate soil erosion. Both inputs are optional and default to zero, but without precipitation there is no percolation and no soil erosion.

(hydrology:groundwater-store)=
## Groundwater store

The groundwater store represents the slow sub-surface pathway from the land to the river: water that recharges an aquifer and discharges to the river as baseflow over days to years. Runoff inputs often route sub-surface water to the river on the same timestep that it is generated, which gives river flows that recede too quickly after rain and are too low in dry periods. The store adds this delay. It is optional and switched off by default.

Each grid cell has its own *linear reservoir*: a store whose outflow is proportional to the amount of water it holds. The store has no feedback on the rest of the model, so it is applied once to each grid cell's runoff time series when the input data are read, rather than on every timestep of the simulation. In batch runs (see [](../users/batch.md)), the state of the store is carried over from one batch to the next.

### Recharge

On each timestep $t$, the total runoff is first split into its surface and sub-surface parts:

$$
q_{T,t} = \max(\text{runoff}_t, 0), \qquad
q_{S,t} = \min\left(\max(\text{quickflow}_t, 0), \, q_{T,t}\right), \qquad
q_{B,t} = q_{T,t} - q_{S,t}
$$

Missing values in the runoff are treated as zero runoff. A fraction of each part recharges the store, and the rest goes directly to the water bodies:

$$
R^*_t = f_{sub} \, q_{B,t} + f_{surf} \, q_{S,t}
$$

$$
q_{D,t} = q_{T,t} - R^*_t
$$

where $R^*_t$ is the gross recharge [m/timestep], $q_{D,t}$ is the direct runoff [m/timestep], and $f_{sub}$ and $f_{surf}$ are the `subsurface_recharge_fraction` (default 1) and `surface_recharge_fraction` (default 0) parameters. By default, all sub-surface runoff and no surface runoff passes through the store.

A fraction $f_{loss}$ (the `loss_fraction` parameter, default 0) of the gross recharge leaves the model and never reaches the water bodies. This can represent, for example, water that flows to a deep aquifer that discharges outside the model domain, or groundwater abstraction. The net recharge to the store is:

$$
R_t = (1 - f_{loss}) \, R^*_t
$$

### Baseflow

The store gains the net recharge, and then drains a fraction $k$ of its contents as baseflow $B_t$:

$$
k = 1 - \exp\left(-\frac{\Delta t}{\tau}\right)
$$

$$
B_t = k \left(S_{t-1} + R_t\right)
$$

$$
S_t = (1 - k)\left(S_{t-1} + R_t\right)
$$

Here $S$ is the water in the store, as a depth over the grid cell [m], $\Delta t$ is the model timestep, and $\tau$ is the residence time of water in the store (the `residence_time` parameter or the `gw_residence_time` input variable, given in days). This is the exact solution for a linear reservoir draining over one timestep, so the results do not depend on how $\tau$ compares with the timestep. The baseflow is added to the direct runoff to give the grid cell's runoff to its water bodies:

$$
q_{runoff,t} = q_{D,t} + B_t
$$

Over a long period, the baseflow is equal to the net recharge, so the store changes when water reaches the river but not how much (apart from the loss fraction). A longer residence time gives a slower, smoother baseflow response.

### Initial storage

The amount of water in the store at the start of the simulation can either be given (the `initial_storage` parameter, in metres), or calculated as the steady state for the mean net recharge $\bar{R}$ over the input time series (the default):

$$
S_0 = \bar{R} \, \frac{1 - k}{k}
$$

With this initial storage, baseflow on the first timestep is equal to the mean recharge rate, which avoids a spin-up period with no baseflow. In batch runs, $\bar{R}$ is the mean over the first batch. After a warm up period, the store continues from its state at the end of the warm up, like the rest of the model. The store's contents are saved in [checkpoints](../users/checkpointing.md), and a run started from a checkpoint continues the store from the saved storage.

### Parameters

The store is configured in the optional `&groundwater` group of the constants file, and all four parameters (residence time, the two recharge fractions and the loss fraction) can also be given per grid cell in the input NetCDF file (see [](input-data:groundwater)). A residence time map value is used where it is positive, and a fraction map value where it is zero or more; elsewhere the constant applies. The store is switched off in any grid cell where the resulting residence time is zero or negative, and in that cell the runoff passes to the water bodies unchanged.

### What the store does not change

The store only changes the grid cell runoff that generates flow in the water bodies. The surface runoff (quickflow) that drives the delivery of eroded soil to the water bodies and the wash-off of dissolved contaminant from the soil is *not* changed by the store, even if part of it recharges the store (`surface_recharge_fraction` > 0).

```{note}
The water groundwater store described here is separate from the contaminant groundwater reservoir set by `groundwater_residence_time` in the `&pfas` constants group. That reservoir holds only the PFAS leached from the bottom of the soil profile, and releases it to the water bodies as dissolved PFAS with its own residence time. It does not hold or release any water, and the two can be used independently.
```

(hydrology:routing)=
## Routing runoff to water bodies

### Surface water network

The surface water network is defined by the input variables `outflow`, `inflows`, `n_waterbodies` and `is_headwater`, which are usually derived from flow directions calculated from a digital elevation model (see [](../users/netcdf-namelist-input.md)). Each grid cell drains to one outflow cell and can have up to seven inflow cells. The number of water bodies in each grid cell is given by `n_waterbodies`. Water bodies are updated in routing order, from the headwaters downstream, so that each river reach receives the outflow from its upstream reaches on the same timestep. Without channel routing, this means that water (and the sediment and contaminant it carries) can travel through the whole river network within a single timestep. With channel routing switched on (see [](theory:channel-routing)), each reach receives its upstream reaches' outflow as a sub-timestep hydrograph and delays and attenuates it, so travel times build up with distance down the network.

### Runoff to each water body

A grid cell's runoff is shared between its water bodies in proportion to their length. For a reach $r$ of length $L_r$ in a grid cell of area $A$ with total reach length $\sum L$, the runoff volume to the reach is:

$$
Q_{runoff,r} = q_{runoff} \, A \, \frac{L_r}{\sum L}
$$

and the total water inflow to the reach on the timestep is the sum of the runoff and the outflows from its upstream reaches:

$$
Q_{in} = Q_{inflow} + Q_{runoff}
$$

Eroded soil and dissolved contaminant from the land are shared between the reaches in the same way.

### River reach hydraulics

Without channel routing (the default), river reaches do not store water from one timestep to the next. The flow through the reach is $Q = Q_{in} / \Delta t$ [m<sup>3</sup>/s], and the whole of $Q_{in}$ leaves the reach as outflow on the same timestep. The reach dimensions are calculated from the flow on each timestep. The width $W$ [m] is:

$$
W = 1.22 \, Q^{0.557}
$$

The depth $D$ [m] is found from Manning's equation for a rectangular channel, solved iteratively with Newton's method:

$$
Q = \frac{\sqrt{S}}{n} \, W D \left( \frac{W D}{W + 2D} \right)^{2/3}
$$

where $S$ is the reach slope and $n$ is Manning's roughness coefficient, given by the `river_manning_n` constant (default 0.035) or the optional `river_manning_n` input map. The velocity is $v = Q / (W D)$, and the volume of water in the reach is $V = D W L f_m$, where $L$ is the reach length and $f_m$ is the meandering factor (see the `river_meandering_factor` constant).

The reach length is the straight-line distance through the grid cell from where the reach enters (or the cell centre, for headwaters) to where it leaves. The slope is the drop in elevation along this distance, calculated from the `dem` input variable, with a minimum given by the `min_stream_slope` config option (default 0.0001 m/m). If no `dem` is provided, every reach has a slope of 0.0005 m/m.

To keep sediment and contaminant transport stable when a lot of water passes through a small reach, each timestep is split into $n_{disp}$ equal *displacements*:

$$
n_{disp} = \left\lceil \frac{Q_{in}}{V} \right\rceil
$$

so that no more than one reach volume of water passes through the reach in each displacement. Each displacement carries $Q_{in} / n_{disp}$ of water, and sediment and contaminant are moved out of the reach in proportion to the fraction of the reach volume that flows out.

The outflow from a grid cell is the sum of the outflows from those of its reaches that flow out of the cell.

(theory:channel-routing)=
### Channel routing

Channel routing is switched on with the `include_routing` config option (in the `&water` group). River reaches then hold water between timesteps, and the time it takes a flood wave to travel through each reach is represented using the variable-parameter Muskingum-Cunge method {cite:p}`cunge_1969`, at a sub-timestep given by the `routing_timestep` config option (default 3600 s). Each timestep is split into $N = \Delta t / \Delta t_r$ routing sub-steps, rounded to a whole number.

On each timestep, each river reach, in routing order:

1. Builds its inflow hydrograph for each sub-step $k$: the sub-step outflows of its upstream river reaches, plus its runoff spread evenly over the timestep. Upstream water bodies that aren't routed river reaches (lakes, estuaries) contribute their outflow for the timestep, spread evenly. Sediment and contaminant inflows are built in the same way, with eroded soil, bank erosion and contaminant sources spread evenly over the sub-steps.
2. Routes the water through the reach (below), giving the outflow at the end of each sub-step.
3. Sets its dimensions (width, depth, volume) and the resuspension and settling rates from the routed flow, using the mean of the inflow and outflow over the timestep.
4. Moves sediment and contaminant through the reach sub-step by sub-step, with the same displacement physics as without routing, but driven by the routed outflow. Each sub-step is split into as many displacements as needed for no more than one reach volume to flow in or out per displacement.
5. Stores its outflows of water, sediment and contaminant for each sub-step, to be the inflow hydrograph of the water body downstream.

Sediment and contaminant therefore travel down the network with the water. Within a reach, they mix through the reach volume and leave in proportion to the routed outflow, so they move at about the flow velocity, while the flood wave moves at the (faster) wave celerity.

```{mermaid}
flowchart LR
    UP["Upstream reaches:<br/>outflow per sub-step"] --> IN["Inflow hydrograph<br/>(sub-steps k = 1..N)"]
    RO["Runoff, erosion, sources<br/>(spread evenly)"] --> IN
    IN --> MC["Muskingum-Cunge<br/>water routing"]
    MC --> DIM["Reach dimensions and rates<br/>from routed flow"]
    DIM --> TR["Sediment and contaminant<br/>transport per sub-step"]
    MC --> TR
    TR --> OUT["Outflow per sub-step<br/>to the reach downstream"]
```

#### Muskingum-Cunge routing

The outflow at the end of a routing step is a weighted combination of the inflow at the start and end of the step and the outflow at the start of the step:

$$
O_{k} = C_0 I_{k} + C_1 I_{k-1} + C_2 O_{k-1}
$$

with

$$
C_0 = \frac{\tfrac{1}{2}\Delta t - K X}{K(1-X) + \tfrac{1}{2}\Delta t}, \qquad
C_1 = \frac{\tfrac{1}{2}\Delta t + K X}{K(1-X) + \tfrac{1}{2}\Delta t}, \qquad
C_2 = \frac{K(1-X) - \tfrac{1}{2}\Delta t}{K(1-X) + \tfrac{1}{2}\Delta t}
$$

where $K$ [s] is the travel time of the flood wave through the reach and $X$ [-] weights the inflow against the outflow. In the Cunge form, $K$ and $X$ are not calibrated per reach but calculated from the reach properties, at a reference flow $Q_0$ (the mean of the current inflow and the previous outflow):

$$
K = \frac{L f_m}{c}, \qquad
X = \frac{1}{2}\left(1 - \frac{Q_0}{W S c L f_m}\right)
$$

with $X$ limited to between 0 and 0.5. The wave celerity $c$ [m/s] is that of a kinematic wave in a wide channel, $c = \tfrac{5}{3} v$, multiplied by the `routing_celerity_factor` constant (default 1). Here $v$ is the flow velocity from Manning's equation in its wide-channel form (hydraulic radius equal to depth), which is explicit and cheap enough to evaluate on every routing step:

$$
D = \left( \frac{n \, Q_0}{W \sqrt{S}} \right)^{3/5}, \qquad v = \frac{Q_0}{W D}
$$

and the width $W$ is that given by the width-flow relation above. The diffusion of the flood wave that the method represents comes from the $Q_0 / (W S)$ term in $X$, so flood peaks are attenuated more in wide, gently sloping reaches.

To keep the routing stable and free of negative flows, each sub-step is split into internal steps no longer than $K$ (so that $C_2 \geq 0$), with the inflow interpolated linearly between the sub-steps, and $X$ is limited to at most $\Delta t / (2K)$ (so that $C_0 \geq 0$). Water volumes over each step are calculated from the trapezoidal mean of the flows at its start and end, and the water stored in the reach (channel storage) is updated from the inflow and outflow volumes, so water is conserved exactly. The outflow on a step is limited to the water stored in the reach plus its inflow.

On the first timestep, each reach starts at steady state with the inflow on its first sub-step (outflow equal to inflow). The routing state (the inflow and outflow at the end of the last sub-step, and the channel storage) carries over between timesteps and batches, and is saved in checkpoints.

#### What routing does and doesn't represent

Routing represents the time water takes to travel along the river network, and the attenuation of flood peaks on the way. It doesn't represent delays before water reaches the river, such as the response of hillslopes and soils, which are part of the runoff inputs (and, for the slow component, of the groundwater store). Because the delay through each reach is about $L/c$, routing matters most in large river networks; in small catchments the travel time through the network is often much shorter than a day.

The routing timestep should be no longer than the travel time through a typical reach (for example, about an hour for 5 km reaches at 1.5 m/s), so that a flood wave doesn't cross several reaches within one sub-step. Routing makes each timestep more costly, as the transport of sediment and contaminant runs for each sub-step.

### Estuaries and lakes

Estuary reaches and lake and reservoir segments receive runoff and upstream inflows in the same way as river reaches, but use their own hydraulics. Unlike river reaches, lake and reservoir segments hold a volume of water and release it according to their own outflow rules.
