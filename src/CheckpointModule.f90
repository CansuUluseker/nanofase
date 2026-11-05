module CheckpointModule
    use AbstractEnvironmentModule, only: EnvironmentPointer
    use EnvironmentModule
    use DefaultsModule, only: iouCheckpoint
    use GlobalsModule, only: dp, C, ERROR_HANDLER
    use DataInputModule, only: DATASET
    use LoggerModule, only: LOGR
    use FlowModule
    use ErrorInstanceModule
    use ResultModule
    use UtilModule
    use ContaminantModule
    use WaterBodyModule
    use SoilProfileModule
    use BedSedimentModule
    use GridCellModule, only: GridCell
    implicit none
    private

    type, public :: Checkpoint
        type(EnvironmentPointer) :: env
        character(len=256) :: checkpointFile
    contains
        procedure, public :: init => initCheckpoint
        procedure, public :: save => saveCheckpoint
        procedure, public :: reinstate => reinstateCheckpoint
    end type

contains
    subroutine initCheckpoint(me, env, checkpointFile)
        class(Checkpoint) :: me
        type(Environment), target :: env
        character(len=*) :: checkpointFile
        me%env%item => env
        me%checkpointFile = checkpointFile
    end subroutine

    subroutine saveCheckpoint(me, t)
        class(Checkpoint) :: me
        integer, intent(in) :: t
        integer :: i, j, k, l, m, alloc_stat
        real(dp), allocatable :: soilProfile_contaminant(:,:,:,:,:,:)
        real(dp), allocatable :: soilProfile_m_dissolved(:,:,:)
        real(dp), allocatable :: soilLayer_contaminant(:,:,:,:,:,:,:)
        real(dp), allocatable :: soilLayer_m_dissolved(:,:,:,:)
        real(dp), allocatable :: soilLayer_V_w(:,:,:,:)
        real(dp), allocatable :: waterBody_contaminant(:,:,:,:,:,:)
        real(dp), allocatable :: waterBody_m_dissolved(:,:,:)
        real(dp), allocatable :: waterBody_volume(:,:,:)
        real(dp), allocatable :: waterBody_bedArea(:,:,:)
        real(dp), allocatable :: waterBody_Q(:,:,:,:)
        real(dp), allocatable :: waterBody_Q_final(:,:,:,:)
        real(dp), allocatable :: waterBody_j_spm(:,:,:,:,:)
        real(dp), allocatable :: waterBody_j_spm_final(:,:,:,:,:)
        real(dp), allocatable :: bedSediment_contaminant(:,:,:,:,:,:,:)
        real(dp), allocatable :: bedSediment_m_dissolved(:,:,:,:)
        real(dp), allocatable :: sedimentLayer_M_f(:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_M_f_backup(:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_V_w(:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_f_comp(:,:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_pd_comp(:,:,:,:,:,:)
        real(dp), allocatable :: gridCell_gwStorage(:,:)
        real(dp), allocatable :: waterBody_routing(:,:,:,:)
        real(dp), allocatable :: waterBody_m_spm(:,:,:,:)
        character(len=256) :: tr
        tr = "Checkpoint%saveCheckpoint"

        allocate(soilProfile_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), 1, &
                                        C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 soilProfile_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), 1), &
                 soilLayer_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), 1, C%nSoilLayers, &
                                       C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 soilLayer_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), 1, C%nSoilLayers), &
                 soilLayer_V_w(DATASET%gridShape(1), DATASET%gridShape(2), 1, C%nSoilLayers), &
                 waterBody_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                       C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 waterBody_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies)), &
                 waterBody_volume(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies)), &
                 waterBody_bedArea(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies)), &
                 waterBody_Q(5, maxval(DATASET%nWaterbodies), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_Q_final(5, maxval(DATASET%nWaterbodies), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_j_spm(8, C%nSizeClassesSpm, maxval(DATASET%nWaterbodies), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_j_spm_final(8, C%nSizeClassesSpm, maxval(DATASET%nWaterbodies), &
                 DATASET%gridShape(1), DATASET%gridShape(2)), &
                 bedSediment_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                         C%nSedimentLayers+3, C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 bedSediment_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                         C%nSedimentLayers+3), &
                 sedimentLayer_M_f(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                   C%nSedimentLayers, C%nSizeClassesSpm), &
                 sedimentLayer_M_f_backup(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                          C%nSedimentLayers, C%nSizeClassesSpm), &
                 sedimentLayer_V_w(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                   C%nSedimentLayers, C%nSizeClassesSpm), &
                 sedimentLayer_f_comp(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                      C%nSedimentLayers, C%nSizeClassesSpm, C%nFracCompsSpm), &
                 sedimentLayer_pd_comp(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                       C%nSedimentLayers, C%nSizeClassesSpm, C%nFracCompsSpm), &
                 stat=alloc_stat)
        if (alloc_stat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Checkpoint allocation failed", trace=[tr]))

        soilProfile_contaminant = 0.0_dp
        soilProfile_m_dissolved = 0.0_dp
        soilLayer_contaminant = 0.0_dp
        soilLayer_m_dissolved = 0.0_dp
        soilLayer_V_w = 0.0_dp
        waterBody_contaminant = 0.0_dp
        waterBody_m_dissolved = 0.0_dp
        waterBody_volume = 0.0_dp
        waterBody_bedArea = 0.0_dp
        waterBody_Q = 0.0_dp
        waterBody_Q_final = 0.0_dp
        waterBody_j_spm = 0.0_dp
        waterBody_j_spm_final = 0.0_dp
        bedSediment_contaminant = 0.0_dp
        bedSediment_m_dissolved = 0.0_dp
        sedimentLayer_M_f = 0.0_dp
        sedimentLayer_M_f_backup = 0.0_dp
        sedimentLayer_V_w = 0.0_dp
        sedimentLayer_f_comp = 0.0_dp
        sedimentLayer_pd_comp = 0.0_dp

        do j = 1, DATASET%gridShape(2)
            do i = 1, DATASET%gridShape(1)
                associate (cell => me%env%item%colGridCells(i,j)%item)
                    do k = 1, cell%nSoilProfiles
                        associate (profile => cell%colSoilProfiles(k)%item)
                            soilProfile_contaminant(i,j,k,:,:,:) = profile%m_contaminant%c
                            soilProfile_m_dissolved(i,j,k) = profile%m_contaminant%m_dissolved
                            do l = 1, C%nSoilLayers
                                associate (layer => profile%colSoilLayers(l)%item)
                                    soilLayer_contaminant(i,j,k,l,:,:,:) = layer%m_contaminant%c
                                    soilLayer_m_dissolved(i,j,k,l) = layer%m_contaminant%m_dissolved
                                    soilLayer_V_w(i,j,k,l) = layer%V_w
                                end associate
                            end do
                        end associate
                    end do
                    do k = 1, cell%nReaches
                        associate (water => cell%colRiverReaches(k)%item)
                            waterBody_contaminant(i,j,k,:,:,:) = water%reactor%contaminant%c
                            waterBody_m_dissolved   (i,j,k)     = water%reactor%contaminant%m_dissolved
                            waterBody_volume(i,j,k) = water%volume
                            waterBody_bedArea(i,j,k) = water%bedArea
                            waterBody_Q(:,k,i,j) = water%Q%asArray()
                            waterBody_Q_final(:,k,i,j) = water%Q_final%asArray()
                            waterBody_j_spm(:,:,k,i,j) = water%j_spm%asArray()
                            waterBody_j_spm_final(:,:,k,i,j) = water%j_spm_final%asArray()
                            associate (sediment => water%bedSediment)
                                do l = 1, C%nSedimentLayers+3
                                    bedSediment_contaminant(i,j,k,l,:,:,:) = sediment%m_contaminant(l)%c
                                    bedSediment_m_dissolved(i,j,k,l) = sediment%m_contaminant(l)%m_dissolved
                                end do
                                do l = 1, C%nSedimentLayers
                                    associate (layer => sediment%colBedSedimentLayers(l)%item)
                                        do m = 1, C%nSizeClassesSpm
                                            sedimentLayer_M_f(i,j,k,l,m) = layer%colFineSediment(m)%M_f()
                                            sedimentLayer_M_f_backup(i,j,k,l,m) = layer%colFineSediment(m)%M_f_backup()
                                            sedimentLayer_V_w(i,j,k,l,m) = layer%colFineSediment(m)%V_w()
                                            sedimentLayer_f_comp(i,j,k,l,m,:) = layer%colFineSediment(m)%f_comp
                                            sedimentLayer_pd_comp(i,j,k,l,m,:) = layer%colFineSediment(m)%pd_comp
                                        end do
                                    end associate
                                end do
                            end associate
                        end associate
                    end do
                end associate
            end do
        end do

        ! State added with the groundwater store and channel routing (an extra record at
        ! the end, so that older checkpoint files can still be read): groundwater storage
        ! per cell, routing state (inflow and outflow rates, channel storage) and SPM
        ! mass per water body
        call collectExtraState(me, t, gridCell_gwStorage, waterBody_routing, waterBody_m_spm)

        open(iouCheckpoint, file=trim(me%checkpointFile), form='unformatted', status='replace')
        write(iouCheckpoint) DATASET%gridBounds, DATASET%gridRes
        write(iouCheckpoint) t
        write(iouCheckpoint) soilProfile_contaminant, soilProfile_m_dissolved
        write(iouCheckpoint) soilLayer_contaminant, soilLayer_m_dissolved, soilLayer_V_w
        write(iouCheckpoint) waterBody_contaminant, waterBody_m_dissolved, waterBody_volume, waterBody_bedArea, &
                             waterBody_Q, waterBody_Q_final, waterBody_j_spm, waterBody_j_spm_final
        write(iouCheckpoint) bedSediment_contaminant, bedSediment_m_dissolved, sedimentLayer_M_f, &
                             sedimentLayer_M_f_backup, sedimentLayer_V_w, sedimentLayer_f_comp, sedimentLayer_pd_comp
        write(iouCheckpoint) gridCell_gwStorage, waterBody_routing, waterBody_m_spm
        close(iouCheckpoint)

        call LOGR%toConsole('Saving checkpoint to '//trim(me%checkpointFile)//': '//COLOR_GREEN//'success'//COLOR_RESET)
        call LOGR%toFile('Saving checkpoint to '//trim(me%checkpointFile)//': success')
    end subroutine

    subroutine reinstateCheckpoint(me, preserve_timestep)
        class(Checkpoint) :: me
        logical, optional :: preserve_timestep
        integer :: i, j, k, l, m, ioStat, alloc_stat
        real :: gridRes(2), gridBounds(4)
        real(dp), allocatable :: soilProfile_contaminant(:,:,:,:,:,:)
        real(dp), allocatable :: soilProfile_m_dissolved(:,:,:)
        real(dp), allocatable :: soilLayer_contaminant(:,:,:,:,:,:,:)
        real(dp), allocatable :: soilLayer_m_dissolved(:,:,:,:)
        real(dp), allocatable :: soilLayer_V_w(:,:,:,:)
        real(dp), allocatable :: waterBody_contaminant(:,:,:,:,:,:)
        real(dp), allocatable :: waterBody_m_dissolved(:,:,:)
        real(dp), allocatable :: waterBody_volume(:,:,:)
        real(dp), allocatable :: waterBody_bedArea(:,:,:)
        real(dp), allocatable :: waterBody_Q(:,:,:,:)
        real(dp), allocatable :: waterBody_Q_final(:,:,:,:)
        real(dp), allocatable :: waterBody_j_spm(:,:,:,:,:)
        real(dp), allocatable :: waterBody_j_spm_final(:,:,:,:,:)
        real(dp), allocatable :: bedSediment_contaminant(:,:,:,:,:,:,:)
        real(dp), allocatable :: bedSediment_m_dissolved(:,:,:,:)
        real(dp), allocatable :: sedimentLayer_M_f(:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_M_f_backup(:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_V_w(:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_f_comp(:,:,:,:,:,:)
        real(dp), allocatable :: sedimentLayer_pd_comp(:,:,:,:,:,:)
        integer :: t
        real(dp), allocatable :: gridCell_gwStorage(:,:)
        real(dp), allocatable :: waterBody_routing(:,:,:,:)
        real(dp), allocatable :: waterBody_m_spm(:,:,:,:)
        logical :: hasExtraState
        character(len=256) :: tr
        tr = "Checkpoint%reinstateCheckpoint"

        allocate(soilProfile_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), 1, &
                                        C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 soilProfile_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), 1), &
                 soilLayer_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), 1, C%nSoilLayers, &
                                       C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 soilLayer_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), 1, C%nSoilLayers), &
                 soilLayer_V_w(DATASET%gridShape(1), DATASET%gridShape(2), 1, C%nSoilLayers), &
                 waterBody_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                       C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 waterBody_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies)), &
                 waterBody_volume(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies)), &
                 waterBody_bedArea(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies)), &
                 waterBody_Q(5, maxval(DATASET%nWaterbodies), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_Q_final(5, maxval(DATASET%nWaterbodies), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_j_spm(8, C%nSizeClassesSpm, maxval(DATASET%nWaterbodies), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_j_spm_final(8, C%nSizeClassesSpm, maxval(DATASET%nWaterbodies), &
                 DATASET%gridShape(1), DATASET%gridShape(2)), &
                 bedSediment_contaminant(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                         C%nSedimentLayers+3, C%contaminantDim(1), C%contaminantDim(2), C%contaminantDim(3)), &
                 bedSediment_m_dissolved(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                         C%nSedimentLayers+3), &
                 sedimentLayer_M_f(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                   C%nSedimentLayers, C%nSizeClassesSpm), &
                 sedimentLayer_M_f_backup(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                          C%nSedimentLayers, C%nSizeClassesSpm), &
                 sedimentLayer_V_w(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                   C%nSedimentLayers, C%nSizeClassesSpm), &
                 sedimentLayer_f_comp(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                      C%nSedimentLayers, C%nSizeClassesSpm, C%nFracCompsSpm), &
                 sedimentLayer_pd_comp(DATASET%gridShape(1), DATASET%gridShape(2), maxval(DATASET%nWaterbodies), &
                                       C%nSedimentLayers, C%nSizeClassesSpm, C%nFracCompsSpm), &
                 stat=alloc_stat)
        if (alloc_stat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Checkpoint allocation failed", trace=[tr]))

        open(iouCheckpoint, file=trim(me%checkpointFile), form='unformatted', status='old')
        read(iouCheckpoint, iostat=ioStat) gridBounds, gridRes
        if (ioStat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Error reading grid properties", trace=[tr]))
        if (any(abs(gridBounds - DATASET%gridBounds) > C%epsilon)) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(message="Grid bounds mismatch. Checkpoint: " // &
                trim(adjustl(str(gridBounds))) // ". Simulation: " // trim(adjustl(str(DATASET%gridBounds))), trace=[tr]))
        end if
        if (any(abs(gridRes - DATASET%gridRes) > C%epsilon)) then
            call ERROR_HANDLER%trigger(error=ErrorInstance(message="Grid resolution mismatch. Checkpoint: " // &
                trim(adjustl(str(gridRes))) // ". Simulation: " // trim(adjustl(str(DATASET%gridRes))), trace=[tr]))
        end if
        read(iouCheckpoint, iostat=ioStat) t
        if (ioStat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Error reading timestep", trace=[tr]))
        read(iouCheckpoint, iostat=ioStat) soilProfile_contaminant, soilProfile_m_dissolved
        if (ioStat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Error reading soil profile data", trace=[tr]))
        read(iouCheckpoint, iostat=ioStat) soilLayer_contaminant, soilLayer_m_dissolved, soilLayer_V_w
        if (ioStat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Error reading soil layer data", trace=[tr]))
        read(iouCheckpoint, iostat=ioStat) waterBody_contaminant, waterBody_m_dissolved, waterBody_volume, waterBody_bedArea, &
                                           waterBody_Q, waterBody_Q_final, waterBody_j_spm, waterBody_j_spm_final
        if (ioStat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Error reading water body data", trace=[tr]))
        read(iouCheckpoint, iostat=ioStat) bedSediment_contaminant, bedSediment_m_dissolved, sedimentLayer_M_f, &
                                           sedimentLayer_M_f_backup, sedimentLayer_V_w, sedimentLayer_f_comp, sedimentLayer_pd_comp
        if (ioStat /= 0) call ERROR_HANDLER%trigger(error=ErrorInstance(message="Error reading bed sediment data", trace=[tr]))
        ! Groundwater store, routing and SPM state (absent from older checkpoint files)
        allocate(gridCell_gwStorage(DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_routing(3, max(1, maxval(DATASET%nWaterbodies)), DATASET%gridShape(1), DATASET%gridShape(2)), &
                 waterBody_m_spm(C%nSizeClassesSpm, max(1, maxval(DATASET%nWaterbodies)), DATASET%gridShape(1), &
                                 DATASET%gridShape(2)))
        read(iouCheckpoint, iostat=ioStat) gridCell_gwStorage, waterBody_routing, waterBody_m_spm
        hasExtraState = (ioStat == 0)
        if (.not. hasExtraState) then
            call LOGR%add("Checkpoint: WARNING - no groundwater store, routing or SPM state in " // &
                trim(me%checkpointFile) // " (older format); these start from their initial state")
        end if
        close(iouCheckpoint)

        do j = 1, DATASET%gridShape(2)
            do i = 1, DATASET%gridShape(1)
                associate (cell => me%env%item%colGridCells(i,j)%item)
                    do k = 1, cell%nSoilProfiles
                        associate (profile => cell%colSoilProfiles(k)%item)
                            profile%m_contaminant%c = soilProfile_contaminant(i,j,k,:,:,:)
                            profile%m_contaminant%m_dissolved = soilProfile_m_dissolved(i,j,k)
                            do l = 1, C%nSoilLayers
                                associate (layer => profile%colSoilLayers(l)%item)
                                    layer%m_contaminant%c = soilLayer_contaminant(i,j,k,l,:,:,:)
                                    layer%m_contaminant%m_dissolved = soilLayer_m_dissolved(i,j,k,l)
                                    layer%V_w = soilLayer_V_w(i,j,k,l)
                                end associate
                            end do
                        end associate
                    end do
                    do k = 1, cell%nReaches
                        associate (water => cell%colRiverReaches(k)%item)
                            water%reactor%contaminant%c        = waterBody_contaminant(i,j,k,:,:,:)
                            water%reactor%contaminant%m_dissolved = waterBody_m_dissolved(i,j,k)
                            water%volume = waterBody_volume(i,j,k)
                            water%bedArea = waterBody_bedArea(i,j,k)
                            water%Q = waterBody_Q(:,k,i,j)
                            water%Q_final = waterBody_Q_final(:,k,i,j)
                            water%j_spm = waterBody_j_spm(:,:,k,i,j)
                            water%j_spm_final = waterBody_j_spm_final(:,:,k,i,j)
                            associate (sediment => water%bedSediment)
                                do l = 1, C%nSedimentLayers+3
                                    sediment%m_contaminant(l)%c = bedSediment_contaminant(i,j,k,l,:,:,:)
                                    sediment%m_contaminant(l)%m_dissolved = bedSediment_m_dissolved(i,j,k,l)
                                end do
                                do l = 1, C%nSedimentLayers
                                    associate (layer => sediment%colBedSedimentLayers(l)%item)
                                        do m = 1, C%nSizeClassesSpm
                                            ! Composition first: set() derives the sediment volume from it
                                            layer%colFineSediment(m)%f_comp = sedimentLayer_f_comp(i,j,k,l,m,:)
                                            layer%colFineSediment(m)%pd_comp = sedimentLayer_pd_comp(i,j,k,l,m,:)
                                            ! The saved backup mass (there's no direct setter: set the mass
                                            ! to it, take the backup, then set the actual mass)
                                            call layer%colFineSediment(m)%set(Mf_in=sedimentLayer_M_f_backup(i,j,k,l,m), &
                                                                              Vw_in=sedimentLayer_V_w(i,j,k,l,m))
                                            call layer%colFineSediment(m)%backup_M_f()
                                            call layer%colFineSediment(m)%set(Mf_in=sedimentLayer_M_f(i,j,k,l,m), &
                                                                              Vw_in=sedimentLayer_V_w(i,j,k,l,m))
                                        end do
                                    end associate
                                end do
                            end associate
                        end associate
                    end do
                end associate
            end do
        end do

        if (hasExtraState) call reinstateExtraState(me, gridCell_gwStorage, waterBody_routing, waterBody_m_spm)
        if (present(preserve_timestep) .and. preserve_timestep) then
            C%t0 = t
        end if

        call LOGR%toConsole('Reinstating checkpoint from '//trim(me%checkpointFile)//': '//COLOR_GREEN//'success'//COLOR_RESET)
        call LOGR%toFile('Reinstating checkpoint from '//trim(me%checkpointFile)//': success')
    end subroutine
    !> Gather the groundwater store, routing and water-column SPM state. At t = 0 (a
    !! checkpoint saved after the warm up period) the groundwater storage is that at the end
    !! of the warm up; otherwise it is the current (end of run) storage.
    subroutine collectExtraState(me, t, gw, routing, m_spm)
        class(Checkpoint) :: me
        integer, intent(in) :: t
        real(dp), allocatable, intent(out) :: gw(:,:), routing(:,:,:,:), m_spm(:,:,:,:)
        integer :: i, j, k, nwb
        nwb = max(1, maxval(DATASET%nWaterbodies))
        allocate(gw(DATASET%gridShape(1), DATASET%gridShape(2)), routing(3, nwb, DATASET%gridShape(1), &
                 DATASET%gridShape(2)), m_spm(C%nSizeClassesSpm, nwb, DATASET%gridShape(1), DATASET%gridShape(2)))
        gw = -1.0_dp
        routing = -1.0_dp
        m_spm = 0.0_dp
        do j = 1, DATASET%gridShape(2)
            do i = 1, DATASET%gridShape(1)
                associate (cell => me%env%item%colGridCells(i,j)%item)
                    gw(i,j) = cell%gwStorage
                    if (t == 0 .and. C%warmUpPeriod > 0 .and. allocated(cell%gwStorage_t)) then
                        if (size(cell%gwStorage_t) >= C%warmUpPeriod) gw(i,j) = cell%gwStorage_t(C%warmUpPeriod)
                    end if
                    do k = 1, cell%nReaches
                        associate (water => cell%colRiverReaches(k)%item)
                            routing(:,k,i,j) = [water%mc_I, water%mc_O, water%mc_storage]
                            if (allocated(water%m_spm)) m_spm(:,k,i,j) = water%m_spm
                        end associate
                    end do
                end associate
            end do
        end do
    end subroutine

    !> Restore the groundwater store (recomputing each cell's runoff series from the saved
    !! storage), the routing state and the water-column SPM
    subroutine reinstateExtraState(me, gw, routing, m_spm)
        class(Checkpoint) :: me
        real(dp), intent(in) :: gw(:,:), routing(:,:,:,:), m_spm(:,:,:,:)
        integer :: i, j, k
        do j = 1, DATASET%gridShape(2)
            do i = 1, DATASET%gridShape(1)
                associate (cell => me%env%item%colGridCells(i,j)%item)
                    select type (cell)
                        type is (GridCell)
                            if (gw(i,j) >= 0.0_dp) call cell%reapplyGroundwaterStore(gw(i,j))
                    end select
                    do k = 1, cell%nReaches
                        associate (water => cell%colRiverReaches(k)%item)
                            water%mc_I = routing(1,k,i,j)
                            water%mc_O = routing(2,k,i,j)
                            water%mc_storage = max(routing(3,k,i,j), 0.0_dp)
                            if (allocated(water%m_spm)) water%m_spm = m_spm(:,k,i,j)
                        end associate
                    end do
                end associate
            end do
        end do
    end subroutine

end module