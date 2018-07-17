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
  SUBROUTINE print_ibte( iq, ef0, efcb, first_cycle, totq, lrepmatw2 ) 
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine computes the scattering rate (inv_tau)
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, eps_acustic, degaussw, & 
                            nstemp, scattering_serta, scattering_0rta, shortrange,&
                            restart, restart_freq, restart_filq, vme
  USE pwcom,         ONLY : ef
  USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, dmef, vmef, wf, wqf, & 
                            epf17, nqtotf, nkqtotf, inv_tau_all, inv_tau_allcb, &
                            xqf, zi_allvb, zi_allcb
  USE transportcom,  ONLY : transp_temp, lower_bnd
  USE constants_epw, ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev, & 
                            eps6
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : world_comm, my_pool_id, npool
  USE io_epw,        ONLY : iunepmat, iunepmatcb
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, &
                               MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (INOUT) :: first_cycle
  !! Use to determine weather this is the first cycle after restart 
  INTEGER, INTENT(IN) :: iq
  !! Q-point inde
  REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp. Could be unused (0).
  INTEGER , INTENT(INOUT) :: totq
#if defined(__MPI)  
  INTEGER (kind=MPI_OFFSET_KIND), INTENT(inout) :: lrepmatw2
#else
  INTEGER (kind=8), INTENT (inout) :: lrepmatw2
#endif   
  !
  ! Local variables
  INTEGER :: n
  !! Integer for the degenerate average over eigenstates  
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
  !! Index over temperature range
  INTEGER :: nqtotf_new
  !! Number of q-point in the new dataset
  INTEGER :: k_all(npool)
#if defined(__MPI)
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw
  !! 
  INTEGER (kind=MPI_OFFSET_KIND) :: lsize
  !! Offset to tell where to start reading the file
#else
  INTEGER(kind=8) :: lrepmatw
  INTEGER(kind=8) :: lsize
  !! Offset to tell where to start reading the file
#endif
  INTEGER :: ierr

  !
  REAL(kind=DP) :: tmp
  !! Temporary variable to store real part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp2
  !! Temporary variable for zi_all
  REAL(kind=DP) :: ekk2
  !! Temporary variable to the eigenenergies for the degenerate average  
  REAL(KIND=DP) :: ekk
  !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
  REAL(KIND=DP) :: ekq
  !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
  REAL(KIND=DP) :: g2
  !! Electron-phonon matrix elements squared (g2 is Ry^2) 
  REAL(KIND=DP) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND=DP) :: w0g1
  !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q}] $$ 
  REAL(KIND=DP) :: w0g2 
  !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q}] $$
  REAL(KIND=DP) :: inv_wq 
  !! Inverse phonon frequency. Defined for efficiency reasons.
  REAL(KIND=DP) :: inv_etemp
  !! Invese temperature inv_etemp = 1/etemp. Defined for efficiency reasons.
  REAL(KIND=DP) :: g2_tmp 
  !! Used to set component to 0 if the phonon freq. is too low. This is defined
  !! for efficiency reasons as if statement should be avoided in inner-most loops.
  REAL(KIND=DP) :: inv_degaussw
  !! 1.0/degaussw. Defined for efficiency reasons. 
  REAL(KIND=DP) :: wq
  !! Phonon frequency $$\omega_{q\nu}$$ on the fine grid.  
  REAL(KIND=DP) :: wgq
  !! Bose-Einstein occupation function $$n_{q\nu}$$
  REAL(kind=DP) :: weight
  !! Self-energy factor 
  REAL(KIND=DP) :: fmkq
  !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
  REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
  !! Electronic velocity $$v_{n\mathbf{k}}$$
  REAL(kind=DP) :: trans_prob(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nkf, nstemp)
  !! Temporary array to store the scattering rates
  REAL(kind=DP) :: trans_probcb(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nkf, nstemp)
  !! Temporary array to store the scattering rates
  REAL(kind=DP) :: zi_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the zi
  REAL(KIND=DP), ALLOCATABLE :: inv_tau_all_new (:,:,:)
  !! New scattering rates to be merged
  !
  REAL(KIND=DP), ALLOCATABLE :: etf_all(:,:)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  REAL(KIND=DP), EXTERNAL :: DDOT
  !! Dot product function
  REAL(KIND=DP), EXTERNAL :: efermig
  !! Function that returns the Fermi energy
  REAL(KIND=DP), EXTERNAL :: wgauss
  !! Compute the approximate theta function. Here computes Fermi-Dirac 
  REAL(KIND=DP), EXTERNAL :: w0gauss
  !! The derivative of wgauss:  an approximation to the delta function  
  REAL(kind=DP), PARAMETER :: eps = 1.d-4
  !! Tolerence parameter for the velocity
  REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
  !! Tolerence
  ! 

  ! 
  IF ( iq .eq. 1 ) THEN
    !
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Scattering rate")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick .lt. 1.d3 ) &
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(5x,a,f10.6,a)' ) 'This is computed with respect to the fine Fermi level ',ef * ryd2ev, ' eV'
      WRITE(stdout, '(5x,a,f10.6,a,f10.6,a)' ) 'Only states between ',(ef-fsthick) * ryd2ev, ' eV and ',&
              (ef+fsthick) * ryd2ev, ' eV will be included'
      WRITE(stdout,'(5x,a/)')
    !
  ENDIF
  ! 
  ! In the case of a restart do not add the first step
  IF (first_cycle) THEN
    first_cycle = .FALSE.
  ELSE

    k_all(:) = 0 
    ! loop over temperatures
    DO itemp = 1, nstemp
      !
      etemp = transp_temp(itemp)
      !
      ! SP: Define the inverse so that we can efficiently multiply instead of
      ! dividing
      !
      inv_etemp = 1.0/etemp
      inv_degaussw = 1.0/degaussw
      !  
      !
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! We are not consistent with ef from ephwann_shuffle but it should not 
        ! matter if fstick is large enough.
        IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
          
          IF ( itemp == 1 ) THEN   
            ! Number of k-points for that cpu for that q-points
            ! This is temperature indepdendent
            k_all(my_pool_id+1) = k_all(my_pool_id+1) + 1 
          ENDIF
          !
          DO imode = 1, nmodes
            !
            ! the phonon frequency and bose occupation
            wq = wf (imode, iq)
            !
            ! SP : Avoid if statement in inner loops
            ! the coupling from Gamma acoustic phonons is negligible
            IF ( wq .gt. eps_acustic ) THEN
              g2_tmp = 1.0
              wgq = wgauss( -wq*inv_etemp, -99)
              wgq = wgq / ( one - two * wgq )
              ! SP : Define the inverse for efficiency
              inv_wq =  1.0/(two * wq) 
            ELSE
              g2_tmp = 0.0
              wgq = 0.0
              inv_wq = 0.0
            ENDIF
            !
            DO ibnd = 1, ibndmax-ibndmin+1
              !
              !  energy at k (relative to Ef)
              ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
              !
              DO jbnd = 1, ibndmax-ibndmin+1
                !
                !  energy and fermi occupation at k+q
                ekq = etf (ibndmin-1+jbnd, ikq) - ef0(itemp)
                fmkq = wgauss( -ekq*inv_etemp, -99)
                !
                ! here we take into account the zero-point sqrt(hbar/2M\omega)
                ! with hbar = 1 and M already contained in the eigenmodes
                ! g2 is Ry^2, wkf must already account for the spin factor
                !
                ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                ! need to treat it like the normal g with abs(g).
                IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                   .OR. abs(xqf (3, iq))> eps2 )) THEN
                  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                  !     number, in which case its square will be a negative number. 
                  g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp, KIND=DP )
                ELSE
                  g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                ENDIF
                !
                ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                !
                ! transition probability 
                ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                !
                ! This is summed over modes
                trans_prob(jbnd,ibnd,k_all(my_pool_id+1),itemp) = trans_prob(jbnd,ibnd,k_all(my_pool_id+1),itemp) &
                                    + pi * wqf(iq) * g2 * ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                !
                !  
              ENDDO !jbnd
              !
            ENDDO !ibnd
            !
            ! In this case we are also computing the scattering rate for another Fermi level position
            ! This is used to compute both the electron and hole mobility at the same time.  
            IF ( ABS(efcb(itemp)) > eps ) THEN
              ! 
              DO ibnd = 1, ibndmax-ibndmin+1
                !
                !  energy at k (relative to Ef)
                ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                !
                DO jbnd = 1, ibndmax-ibndmin+1
                  !
                  !  energy and fermi occupation at k+q
                  ekq = etf (ibndmin-1+jbnd, ikq) - efcb(itemp)
                  fmkq = wgauss( -ekq*inv_etemp, -99)
                  !
                  ! here we take into account the zero-point sqrt(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                  ! need to treat it like the normal g with abs(g).
                  IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                     .OR. abs(xqf (3, iq))> eps2 )) THEN
                    ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                    !     number, in which case its square will be a negative number. 
                    g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp, KIND=DP)
                  ELSE
                    g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                  ENDIF
                  !
                  ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                  w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                  w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                  !
                  ! transition probability 
                  ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                  ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                  !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                  !
                  trans_probcb(jbnd,ibnd,k_all(my_pool_id+1),itemp) = trans_probcb(jbnd,ibnd,k_all(my_pool_id+1),itemp) &
                                    + pi * wqf(iq)* g2 * ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )

                  ! 
                ENDDO !jbnd
                !
              ENDDO !ibnd
              ! 
            ENDIF ! ABS(efcb) < eps
            !
          ENDDO !imode
          !
        ENDIF ! endif  fsthick
        !
      ENDDO ! end loop on k
    ENDDO ! itemp
    !
    ! If the q-point is taken, write on file
    CALL mp_sum( k_all,    world_comm )
    IF ( sum(k_all) > 0 ) THEN
      !
      totq = totq + 1
      ! 
      ! Offset where we need to start writing (we increment for each q-points)
      lrepmatw = lrepmatw2 + 2_MPI_OFFSET_KIND * INT( ibndmax-ibndmin+1  , kind = MPI_OFFSET_KIND ) * &
                    INT( ibndmax-ibndmin+1  , kind = MPI_OFFSET_KIND ) * &
                    INT( nstemp , kind = MPI_OFFSET_KIND ) * &
                    INT( SUM( k_all(1:my_pool_id+1) ), kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND 
      ! 
      lrepmatw2 = lrepmatw2 + 2_MPI_OFFSET_KIND * INT( ibndmax-ibndmin+1  , kind= MPI_OFFSET_KIND ) * &
                    INT( ibndmax-ibndmin+1  , kind = MPI_OFFSET_KIND ) * &
                    INT( nstemp , kind = MPI_OFFSET_KIND ) * &
                    INT( SUM( k_all(:) ), kind = MPI_OFFSET_KIND )* 8_MPI_OFFSET_KIND 
      ! 
      ! Size of what we write
      lsize = 2_MPI_OFFSET_KIND * INT( ibndmax-ibndmin+1  , kind =MPI_OFFSET_KIND ) * &
                    INT( ibndmax-ibndmin+1  , kind = MPI_OFFSET_KIND ) * &
                    INT( nstemp , kind = MPI_OFFSET_KIND ) * &
                    INT( k_all(my_pool_id+1) , kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      
      !print*,'my_pool_id ',my_pool_id
      !print*,'k_all ',k_all 
      !print*,'lrepmatw ',lrepmatw
      !print*,'lsize ',lsize
      !print*,'iq totq size k_all ',iq, totq, lrepmatw, k_all
      ! 
      CALL MPI_FILE_SEEK(iunepmat,lrepmatw,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      ! 
      CALL MPI_FILE_WRITE(iunepmat, trans_prob, lsize, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      ! 
      CALL MPI_FILE_SEEK(iunepmatcb,lrepmatw,MPI_SEEK_SET,ierr) 
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_SEEK',1 )
      ! 
      CALL MPI_FILE_WRITE(iunepmatcb, trans_probcb, lsize, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'print_ibte', 'error in MPI_FILE_WRITE',1 )      
      !  
    ENDIF 
    ! 
  ENDIF ! first_cycle
  ! 
  !
  ! The k points are distributed among pools: here we collect them
  !
  !IF ( iq .eq. nqtotf ) THEN
    !
    ! 
  !ENDIF ! iq 
  !
  RETURN
  !
  END SUBROUTINE print_ibte
  !-----------------------------------------------------------------------
