  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2016-2018 Samuel Ponce'
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE ibte( totq, etf_all, vkk_all, wkf_all, trans_prob, ef0, trans_probcb, efcb ) 
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine computes the scattering rate with the iterative BTE
  !!  (inv_tau).
  !!  The fine k-point and q-point grid have to be commensurate. 
  !!  The k-point grid uses crystal symmetry to decrease computational cost.
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : alat, at, omega, bg
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : fsthick, & 
                            eps_acustic, degaussw, nstemp, & 
                            system_2d, int_mob, ncarrier, restart, restart_freq,&
                            mp_mesh_k, nkf1, nkf2, nkf3, vme, broyden_beta
  USE pwcom,         ONLY : ef 
  USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, wkf, dmef, vmef, & 
                            wf, wqf, xkf, epf17, nqtotf, nkqtotf, inv_tau_all, xqf, & 
                            F_current, Fi_all, F_currentcb, Fi_allcb, &
                            inv_tau_allcb, map_rebal
  USE transportcom,  ONLY : transp_temp, mobilityh_save, mobilityel_save, lower_bnd, &
                            ixkqf_tr, s_BZtoIBZ_full
  USE constants_epw, ONLY : zero, one, two, pi, kelvin2eV, ryd2ev, & 
                            electron_SI, bohr2ang, ang2cm, hbarJ, eps6, eps8, eps10
  USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_global,     ONLY : inter_pool_comm, world_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id
  USE symm_base,     ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
  USE superconductivity, ONLY : kpmq_map
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: totq
  !! Total number of q-points in the fsthick
  REAL(KIND=DP), INTENT(IN) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
  !! Eigenenergies
  REAL(KIND=DP), INTENT(IN) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
  !! Velocity of k
  REAL(KIND=DP), INTENT(IN) :: wkf_all(nkqtotf/2)
  !! Weight of k
  REAL(KIND=DP), INTENT(IN) :: trans_prob(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nstemp, nkqtotf/2, totq)
  !! Transition probability
  REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
  !! The Fermi level 
  REAL(KIND=DP), INTENT(IN) :: trans_probcb(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nstemp, nkqtotf/2, totq)
  !! Transition probability for the electron
  REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
  !! The Fermi level for the electron
  !
  ! Local variables
  INTEGER :: i, iiq, iq
  !! Cartesian direction index 
  INTEGER :: j
  !! Cartesian direction index 
  INTEGER :: ij
  !! Cartesian direction index 
  INTEGER :: ik
  !! K-point index
  INTEGER :: ikk
  !! Odd index to read etf
  INTEGER :: ikq
  !! Even k+q index to read etf
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: jbnd
  !! Local band index
  INTEGER :: imode
  !! Local mode index
  INTEGER :: itemp
  !! Temperature index
  INTEGER :: ipool
  !! Index of the pool
  INTEGER :: nkq
  !! Index of the pool the the k+q point is
  INTEGER :: nkq_abs
  !! Index of the k+q point from the full grid. 
  INTEGER :: nkqtotf_tmp
  !! Temporary k-q points.
  INTEGER :: ikbz
  !! k-point index that run on the full BZ
  INTEGER :: nb
  !! Number of points in the BZ corresponding to a point in IBZ 
  INTEGER :: BZtoIBZ_tmp(nkf1*nkf2*nkf3)
  !! Temporary mapping
  INTEGER :: BZtoIBZ(nkf1*nkf2*nkf3)
  !! BZ to IBZ mapping
  INTEGER :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
  !! symmetry 
  ! 
  REAL(KIND=DP) :: tau
  !! Relaxation time
  REAL(KIND=DP) :: ekk
  !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
  REAL(KIND=DP) :: ekq
  !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
  REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
  !! Electronic velocity $$v_{n\mathbf{k}}$$
  REAL(kind=DP) :: xkf_all(3,nkqtotf)
  !! Collect k-point coordinate (and k+q) from all pools in parallel case
  REAL(kind=DP) :: xkf_red(3,nkqtotf/2)
  !! Collect k-point coordinate from all pools in parallel case
  REAL(kind=DP) :: F_SERTA(3, ibndmax-ibndmin+1, nstemp, nkqtotf/2)
  !! SERTA solution
  REAL(kind=DP) :: F_SERTAcb(3, ibndmax-ibndmin+1, nstemp, nkqtotf/2)
  !! SERTA solution
  REAL(kind=DP) :: tmp
  !
  REAL(kind=DP) :: xkf_tmp (3, nkqtotf)
  !! Temporary k-point coordinate (dummy variable)
  REAL(kind=DP) :: wkf_tmp(nkqtotf)
  !! Temporary k-weights (dummy variable)
  ! 
  ! Gather all the k-point coordinate from all the pools
  xkf_all(:,:) = zero
  xkf_red(:,:) = zero
#ifdef __MPI
  CALL poolgather2 ( 3, nkqtotf, nkqf, xkf, xkf_all)
#else
  xkf_all = xkf
#endif
  !
  ! Deal with symmetries
  IF (mp_mesh_k) THEN
    ALLOCATE(ixkqf_tr(nkf,nqtotf))
    ALLOCATE(s_BZtoIBZ_full(3,3,nkf,nqtotf))
    BZtoIBZ(:) = 0
    s_BZtoIBZ(:,:,:) = 0
    ixkqf_tr(:,:) = 0
    s_BZtoIBZ_full(:,:,:,:) = 0
    ! 
    IF ( mpime .eq. ionode_id ) THEN
      ! 
      CALL set_sym_bl( )
      !
      ! What we get from this call is BZtoIBZ
      CALL kpoint_grid_epw ( nrot, time_reversal, .false., s, t_rev, bg, nkf1*nkf2*nkf3, &
                 nkf1,nkf2,nkf3, nkqtotf_tmp, xkf_tmp, wkf_tmp,BZtoIBZ,s_BZtoIBZ)
      ! 
      DO ik = 1, nkqtotf/2
        ikk = 2 * ik - 1
        xkf_red(:,ik) = xkf_all(:,ikk)
      ENDDO
      ! 
      BZtoIBZ_tmp(:) = 0
      DO ikbz=1, nkf1*nkf2*nkf3
        BZtoIBZ_tmp(ikbz) = map_rebal( BZtoIBZ( ikbz ) )
      ENDDO
      BZtoIBZ(:) = BZtoIBZ_tmp(:)
      ! 
    ENDIF ! mpime
    CALL mp_bcast( xkf_red, ionode_id, inter_pool_comm )
    CALL mp_bcast( s_BZtoIBZ, ionode_id, inter_pool_comm )
    CALL mp_bcast( BZtoIBZ, ionode_id, inter_pool_comm )
    ! 
    DO ik = 1, nkf
      !
      DO iiq=1, nqtotf
        ! 
        CALL kpmq_map( xkf_red(:,ik+lower_bnd-1), xqf (:, iiq), +1, nkq_abs )
        ! 
        ! We want to map k+q onto the full fine k and keep the symm that bring
        ! that point onto the IBZ one.
        s_BZtoIBZ_full(:,:,ik,iiq) = s_BZtoIBZ(:,:,nkq_abs)
        !
        ixkqf_tr(ik,iiq) = BZtoIBZ(nkq_abs)
        ! 
      ENDDO ! q-loop
    ENDDO ! k-loop
  ENDIF
  !
  ! First computes the SERTA solution as the first step of the IBTE
  F_SERTA(:,:,:,:) = 0.0d0
  F_SERTAcb(:,:,:,:) = 0.0d0
  ! 
  DO ik = 1, nkf
    DO itemp = 1, nstemp
      DO ibnd = 1, ibndmax-ibndmin+1
        ! 
        ! Sum over jbnd and iq
        tmp = SUM(trans_prob(:, ibnd, itemp, ik+lower_bnd-1, :))
        IF ( ABS(tmp) > eps10 ) THEN
          F_SERTA(:, ibnd, itemp, ik+lower_bnd-1) = F_SERTA(:, ibnd, itemp,ik+lower_bnd-1) + &
                     vkk_all(:,ibnd,ik+lower_bnd-1) / ( 2.0d0 * tmp )
        ENDIF
        ! 
        tmp = SUM(trans_probcb(:, ibnd, itemp, ik+lower_bnd-1, :))
        IF ( ABS(tmp) > eps10 ) THEN
          F_SERTAcb(:, ibnd, itemp, ik+lower_bnd-1) = F_SERTAcb(:, ibnd, itemp,ik+lower_bnd-1) + &
                     vkk_all(:,ibnd,ik+lower_bnd-1) / ( 2.0d0 * tmp )
        ENDIF
        !
 
      ENDDO
    ENDDO
  ENDDO
  CALL mp_sum(F_SERTA, world_comm)
  CALL mp_sum(F_SERTAcb, world_comm)
  ! Now compute and print the hole mobility
  CALL print_mob(F_SERTA, F_SERTAcb, BZtoIBZ, s_BZtoIBZ, vkk_all, etf_all, wkf_all, ef0, efcb)
  ! 
  RETURN
  !
  ! ---------------------------------------------------------------------------
  END SUBROUTINE ibte
  !----------------------------------------------------------------------------
