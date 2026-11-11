module PFASProcessModule
    !! Generic PFAS mass-transfer operators shared across compartments: fractional
    !! removal, three-phase (aqueous/sorbed/air-water interface) partitioning with
    !! optional rate limitation, and a first-order multi-species transformation network.
    use GlobalsModule, only: dp, C
    use PFASConstantsModule
    implicit none

contains

    pure function pfas_clamp01(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = max(0.0_dp, min(1.0_dp, x))
    end function pfas_clamp01


    subroutine pfas_remove_by_fraction(pool, fraction, removed)
        real(dp), intent(inout) :: pool(:)
        real(dp), intent(in)    :: fraction
        real(dp), intent(out)   :: removed(:)
        real(dp) :: f

        f = pfas_clamp01(fraction)
        removed = min(max(0.0_dp, pool), max(0.0_dp, pool) * f)
        pool = max(0.0_dp, pool - removed)
    end subroutine pfas_remove_by_fraction


    subroutine pfas_partition_phases(pfas, dt, volume_water, sorbed_phase, mass_sorbent, kd, k_mt_sorbed, &
                                     awi_area, kia, k_mt_awi)
        !! Partitioning of each species between the aqueous phase, one sorbed phase
        !! (soil/bed solids or SPM) and the air-water interface (AWI), solved together.
        !! At equilibrium the three pools hold the total in proportion to their capacities
        !!   V_w : Kd*M_s : Kia*A_aw   [all m3],
        !! i.e. linear sorption (Kd [m3 kg-1]) and linear interfacial adsorption
        !! (Gamma = Kia*C_w, Kia [m], A_aw [m2]). With a mass-transfer rate k > 0 [s-1]
        !! a pool relaxes towards its equilibrium mass as exp(-k*dt) (rate-limited,
        !! first-order sorption/desorption); k = 0 means instantaneous equilibrium.
        !! Mass is conserved per species and form; other phases are not touched.
        !! With no water (a dry soil or bed layer) nothing can be dissolved: the
        !! equilibrium is then all on the sorbent and the interface.
        real(dp), intent(inout) :: pfas(:,:,:)
        real(dp), intent(in)    :: dt
        real(dp), intent(in)    :: volume_water     !! Water volume [m3]
        integer, intent(in)     :: sorbed_phase     !! PFAS_SOL or PFAS_SPM
        real(dp), intent(in)    :: mass_sorbent     !! Sorbent (solid or SPM) mass [kg]
        real(dp), intent(in)    :: kd(:), k_mt_sorbed(:)
        real(dp), intent(in)    :: awi_area         !! Air-water interfacial area [m2]
        real(dp), intent(in)    :: kia(:), k_mt_awi(:)

        integer :: i, f, ns
        real(dp) :: cap_s, cap_aw, denom, total, eq_s, eq_aw, m_s, m_aw

        if (sorbed_phase < 1 .or. sorbed_phase > size(pfas,3) .or. sorbed_phase == PFAS_AQ &
            .or. sorbed_phase == PFAS_AWI) return

        ns = min(size(pfas,1), size(kd), size(k_mt_sorbed), size(kia), size(k_mt_awi))
        do i = 1, ns
            cap_s = max(0.0_dp, kd(i)) * max(0.0_dp, mass_sorbent)
            cap_aw = max(0.0_dp, kia(i)) * max(0.0_dp, awi_area)
            denom = max(0.0_dp, volume_water) + cap_s + cap_aw
            if (denom <= C%epsilon) cycle     ! No water, sorbent or interface: nothing to partition into
            do f = 1, size(pfas,2)
                total = max(0.0_dp, pfas(i,f,PFAS_AQ)) + max(0.0_dp, pfas(i,f,sorbed_phase)) &
                        + max(0.0_dp, pfas(i,f,PFAS_AWI))
                if (total <= 0.0_dp) cycle
                eq_s = total * cap_s / denom
                eq_aw = total * cap_aw / denom
                m_s = relax(max(0.0_dp, pfas(i,f,sorbed_phase)), eq_s, k_mt_sorbed(i))
                m_aw = relax(max(0.0_dp, pfas(i,f,PFAS_AWI)), eq_aw, k_mt_awi(i))
                ! Pools relaxing independently can together ask for more than the total
                if (m_s + m_aw > total) then
                    m_s = m_s * total / (m_s + m_aw)
                    m_aw = total - m_s
                end if
                pfas(i,f,sorbed_phase) = m_s
                pfas(i,f,PFAS_AWI) = m_aw
                pfas(i,f,PFAS_AQ) = max(0.0_dp, total - m_s - m_aw)
            end do
        end do

    contains

        pure function relax(m, m_eq, k) result(m_new)
            real(dp), intent(in) :: m, m_eq, k
            real(dp) :: m_new
            if (k > 0.0_dp) then
                m_new = m_eq + (m - m_eq) * exp(-k * dt)
            else
                m_new = m_eq
            end if
        end function
    end subroutine pfas_partition_phases


    subroutine pfas_transform_network(pfas, dt, k_react, reaction_yield)
        !! First-order precursor/intermediate/terminal network.
        !! k_react(i,j) = rate from species i to species j [s-1].
        !! Competing outgoing reactions are handled as one total loss with branching.
        real(dp), intent(inout) :: pfas(:,:,:)
        real(dp), intent(in) :: dt, k_react(:,:), reaction_yield(:,:)

        integer :: i, j, f, p, ns
        real(dp), allocatable :: delta(:,:,:)
        real(dp) :: k_total, frac_loss, branch, loss

        if (dt <= 0.0_dp) return

        ns = min( &
            size(pfas,1), size(k_react,1), size(k_react,2), &
            size(reaction_yield,1), size(reaction_yield,2) &
        )

        allocate(delta(size(pfas,1), size(pfas,2), size(pfas,3)))
        delta = 0.0_dp

        do i = 1, ns
            k_total = 0.0_dp

            do j = 1, ns
                if (j /= i) k_total = k_total + max(0.0_dp, k_react(i,j))
            end do

            if (k_total <= 0.0_dp) cycle

            frac_loss = pfas_clamp01(1.0_dp - exp(-k_total * dt))

            do f = 1, size(pfas,2)
                do p = 1, size(pfas,3)
                    loss = min(max(0.0_dp, pfas(i,f,p)), max(0.0_dp, pfas(i,f,p)) * frac_loss)
                    delta(i,f,p) = delta(i,f,p) - loss

                    do j = 1, ns
                        if (j == i) cycle
                        if (k_react(i,j) <= 0.0_dp) cycle

                        branch = max(0.0_dp, k_react(i,j)) / k_total
                        delta(j,f,p) = delta(j,f,p) + loss * branch * max(0.0_dp, reaction_yield(i,j))
                    end do
                end do
            end do
        end do

        pfas = max(0.0_dp, pfas + delta)
        deallocate(delta)
    end subroutine pfas_transform_network

end module PFASProcessModule
