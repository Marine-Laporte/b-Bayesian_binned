subroutine Likelihood_FMD_gridsearch(mag_bins, delta_bin, min_bin, Beta_prior, Mu_prior, Sigma_prior, &
					Likelihood_FMD_3D, Maxloglike) 
! C ----------------------------------------------------------------------
! C   Compute the 3D posterior distribution of the Frequency-Magnitude Distribution 
! C   considering a detection function (Error function) and a Gutenberg Richter law
! C   using a GRIDSEARCH exploration over the model prior 
! C   Inputs : 
! C    . mag_bins    : number of magnitudes/bin (/subset for subfunctions) {array or list}
! C    . delta_bin   : size of the magnitude binning (usually 0.1 for traditional catalogues) {float}
! C    . min_bin     : minimum magnitude bin (eq. minimum magnitude) for integration bounds {float}
! C    . n_draws     : number of random draw for Beta/Mu/Sigma exploration of the prior {integer}
! C    . Beta_prior    : exploration space of Beta {array or list} 
! C    . Mu_prior      : exploration space of Mu {array or list}
! C    . Sigma_prior   : exploration space of Sigma {array or list}
! C   Outputs :
! C    . Likelihood_FMD_3D  : 3D posterior density distribution for Beta/Mu/Sigma
! C    . Maxloglike    : Maximum Loglikelihoodof the 3D posterior distribution 
! C ----------------------------------------------------------------------
  use ieee_arithmetic

  integer, intent(in) :: mag_bins(:)
  real*8, intent(in) :: delta_bin
  real*8, intent(in) :: min_bin
  real*8, intent(in) :: Mu_prior(:)
  real*8, intent(in) :: Beta_prior(:)  
  real*8, intent(in) :: Sigma_prior(:)
  real*8, intent(out) :: Likelihood_FMD_3D(size(Mu_prior), size(Beta_prior),size(Sigma_prior))
  real*8, intent(out) :: Maxloglike
  
  ! Priors
  integer :: N_bins, N_b, N_m, N_s
  real*8 :: Mmin
  real*8 :: sigma_i, beta_i, mu_i
  
  ! Detection function 
  real*8 :: Q_Mmin, Q_Mc
  real*8 :: Phi_cste, Phi_bin
  real*8 :: terme_m1, terme_m2, terme_3
  real*8 :: m1, m2, q_m1, q_m2, m1_t3, m2_t3, q_m1_t3, q_m2_t3 
  
  ! Likelihood/Loglikelihood
  real*8 :: Likelihood, Loglikelihood

  integer :: i,  k, ib, im, is ! iteratives

  ! ----Init params for gridsearch exploration----- !
  !----------------------------------------------
  
  N_b = size(Beta_prior)
  N_m = size(Mu_prior)
  N_s = size(Sigma_prior)
  N_bins = size(mag_bins)
  Mmin = min_bin
  
  ! Init -Posterior (outputs)
  Likelihood_FMD_3D = 0.0
  
  ! Terms of the detection function
  Q_Mmin = 0
  Q_Mc = 0
  Phi_bin = 0.0
  Phi_cste = 0.0

!C    # --------------------------
!C    # 3D likelihood function
!C    # --------------------------   

!C  Step 1 - Sample over the 3 directions
    do im = 1, N_m
        do ib = 1, N_b
        	do is = 1, N_s
        		!%----Values of beta mu sigma for each sample of the priors in the gridsearch
		       beta_i  = Beta_prior(ib)
		       mu_i    = Mu_prior(im)
		       sigma_i = Sigma_prior(is)

	    		!%----Integral Bounds for q error detection function
	    		Q_Mmin = 0.5 + 0.5 * erf( (Mmin-mu_erf) / (sqrt(2.0)*sigma_erf) )
	    		Q_Mc = 0.5 + 0.5 * erf( ((Mmin+sigma_erf**2 * beta_erf)- mu_erf) / (sqrt(2.0)*sigma_erf) )
	    		
	    		Phi_cste = Q_Mmin + exp( ((sigma_erf**2) * (beta_erf**2) /2) - (beta_erf*(mu_erf-Mmin)))*(1-Q_Mc)
			
			Likelihood = 0.0
			LogLikelihood = 0.0
			
			! Step 2 - Loop over the number of bins
	    		do k = 1, N_bins
				!%---- Likelihood of the frequency-magnitude distribution
				!%---- GR(m) x Q(m)
	        		
	        		if (mag_bins(k) .ne. 0) then
	        			!%- terme M1
					m1 = Mmin + delta_bin * (k-1)
					q_m1 = 0.5 + 0.5 * erf( (m1-mu_erf)/(sqrt(2.0)*sigma_erf))
					terme_m1 =  exp(-beta_erf*(m1-Mmin))*q_m1
					
					!%- terme M2
					m2 = Mmin + delta_bin * k
					q_m2 = 0.5 + 0.5 * erf( (m2-mu_erf)/(sqrt(2.0)*sigma_erf))
					terme_m2 =  exp(-beta_erf*(m2-Mmin))*q_m2
					
					!%- terme 3		
					m1_t3 = m1 + (sigma_erf**2 * beta_erf)
					q_m1_t3 = 0.5 + 0.5 * erf( (m1_t3-mu_erf)/(sqrt(2.0)*sigma_erf))
					m2_t3 = m2 + (sigma_erf**2 * beta_erf)
					q_m2_t3 = 0.5 + 0.5 * erf( (m2_t3-mu_erf)/(sqrt(2.0)*sigma_erf))
					
					terme_3 = exp(((sigma_erf**2)*(beta_erf**2))/2 - (beta_erf*(mu_erf-Mmin))) *( q_m2_t3 - q_m1_t3)
	    	        		
	    	        		!%---- Likelihood formula (from Marsan private comm see Laporte et al.,2025)
	    	        		Phi_bin =  - terme_m2 + terme_m1 + terme_3
			    		Likelihood = Phi_bin/(Phi_cste*delta_bin)
			    		!%--- Loglikelihood = sum of the likelihood in each bin multiplied by the number of magnitudes in the bin
					!%--- Consider all cases where Likelihood = 0 
			    		if (Phi_bin .lt. 0) then
			    			LogLikelihood = -Infinity
			    		else 
			    			LogLikelihood = LogLikelihood + mag_bins(k)* dlog(Likelihood) 
			    			
					end if			    		
				else 
					LogLikelihood = LogLikelihood 
				end if

	    		end do
	    		
	    		Likelihood_FMD_3D(im, ib, is) = LogLikelihood

		enddo
     	enddo
    enddo
    Maxloglike = maxval(Likelihood_FMD_3D) 
    write(*,*) "Maximum loglikelihood (f)", Maxloglike
    Likelihood_FMD_3D = exp(Likelihood_FMD_3D-Maxloglike)
    write(*,*) "Likelihood (sum of 3D mat) (f)", sum(Likelihood_FMD_3D)
    return 
end subroutine   Likelihood_FMD_gridsearch







