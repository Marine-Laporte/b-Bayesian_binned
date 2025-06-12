subroutine Likelihood_FMD_random(mag_bins, delta_bin, min_bin, n_draws, Beta_prior , Mu_prior, Sigma_prior, &
					Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
! C ----------------------------------------------------------------------
! C   Compute the 3D posterior distribution of the Frequency-Magnitude Distribution 
! C   considering a detection function (Error function) and a Gutenberg Richter law
! C   using a random exploration over the model prior (Monte-Carlo)
! C   Inputs : 
! C    . mag_bins    : number of magnitudes/bin (/subset for subfunctions) {array or list}
! C    . delta_bin   : size of the magnitude binning (usually 0.1 for traditional catalogues) {float}
! C    . min_bin     : minimum magnitude bin (eq. minimum magnitude) for integration bounds {float}
! C    . n_draws     : number of random draw for Beta/Mu/Sigma exploration of the prior {integer}
! C    . Beta_prior    : exploration space of Beta {array or list} (careful : here we only retrieve the bounds of this space)
! C    . Mu_prior      : exploration space of Mu {array or list}
! C    . Sigma_prior   : exploration space of Sigma {array or list}
! C   Outputs :
! C    . Likelihood_FMD_MC_3D  : 3D posterior density distribution for Beta/Mu/Sigma
! C    . Maxloglike    : Maximum Loglikelihoodof the 3D posterior distribution 
! C    . Beta_posterior: Marginal density distribution for Beta
! C    . Mu_posterior  : Marginal density distribution for Mu
! C    . Sigma_posterior  : Marginal density distribution for Sigma
! C ----------------------------------------------------------------------
  use ieee_arithmetic

  integer, intent(in) :: mag_bins(:)
  real*8, intent(in) :: delta_bin
  real*8, intent(in) :: min_bin
  integer, intent(in):: n_draws
  real*8, intent(in) :: Beta_prior(:), Mu_prior(:), Sigma_prior(:)
  real*8, intent(out) :: Likelihood_FMD_MC_3D(n_draws)
  real*8, intent(out) :: Maxloglike
  real*8, intent(out) :: Beta_posterior(size(Beta_prior)), Mu_posterior(size(Mu_prior)), Sigma_posterior(size(Sigma_prior))
  
  
  ! Priors
  integer :: N_bins, N_b, N_m, N_s
  real*8 :: Mmin
  real*8 :: Beta_min, Beta_max, Mu_min, Mu_max, Sigma_min, Sigma_max
  real*8 :: Beta_int, Mu_int, Sigma_int
  real*8 :: random_samples(n_draws, 3)
  real*8 :: sigma_i, beta_i, mu_i
  
  ! Detection function 
  real*8 :: Q_Mmin, Q_Mc
  real*8 :: Phi_cste, Phi_bin
  real*8 :: terme_m1, terme_m2, terme_3
  real*8 :: m1, m2, q_m1, q_m2, m1_t3, m2_t3, q_m1_t3, q_m2_t3 
  
  ! Likelihood/Loglikelihood
  real*8 :: Likelihood, Loglikelihood
  
  ! new variables for binning
  integer :: index_beta(n_draws), index_mu(n_draws), index_sigma(n_draws)
  integer :: counts_beta(size(Beta_prior)), counts_mu(size(Mu_prior)), counts_sigma(size(Sigma_prior))
  integer :: i,  k, ib, im, is ! iteratives

  ! ----Init params for random exploration----- !
  !----------------------------------------------
  
  N_b = size(Beta_prior)
  N_m = size(Mu_prior)
  N_s = size(Sigma_prior)
  N_bins = size(mag_bins)
  Mmin = 0
  
  Beta_min = minval(Beta_prior)
  Beta_max = maxval(Beta_prior)
  Mu_min = minval(Mu_prior)
  Mu_max = maxval(Mu_prior)
  Sigma_min = minval(Sigma_prior)
  Sigma_max = maxval(Sigma_prior)
  
  ! Compute prior integrals
  Beta_int = (Beta_max-Beta_min)/n_draws
  Mu_int = (Mu_max-Mu_min)/n_draws
  Sigma_int = (Sigma_max-Sigma_min)/n_draws
  
  ! Init -Posterior (outputs)
  Likelihood_FMD_MC_3D = 0.0
  Beta_posterior = 0.0
  Mu_posterior = 0.0
  Sigma_posterior = 0.0
  
  ! terms of th detection function
  Q_Mmin = 0
  Q_Mc = 0
  Phi_bin = 0.0
  
  ! for reconstructing marginal posteriors
  counts_beta = 0
  index_beta = 0 
  counts_mu = 0 
  index_mu = 0
  counts_sigma = 0
  index_sigma = 0

  !      Computing likelihood over the prior distributions = sampling the posterior !
  !----------------------------------------------------------------------------------
  ! Step 1 - Sample over the 3 directions
  call random_number(random_samples)
  ! Step 2- Loop over the number of random draws
  do i = 1, n_draws
  
  	!%----Values of beta mu sigma for each randow draw over the prior distribution
	beta_i  = random_samples(i,1)* (Beta_max-Beta_min) + Beta_min
	mu_i    = random_samples(i,2)* (Mu_max-Mu_min) + Mu_min
	sigma_i = random_samples(i,3)* (Sigma_max-Sigma_min) + Sigma_min
	
	!%----Integral Bounds for q error detection function
	Q_Mmin = 0.5 + 0.5 * erf( (Mmin-mu_i) / (sqrt(2.0)*sigma_i) )
	Q_Mc = 0.5 + 0.5 * erf( ((Mmin+sigma_i**2 * beta_i)- mu_i) / (sqrt(2.0)*sigma_i) )
	
	Phi_cste = Q_Mmin + exp( ((sigma_i**2) * (beta_i**2) /2) - (beta_i*(mu_i-Mmin)))*(1-Q_Mc)
	
	Likelihood = 0.0
	LogLikelihood = 0.0

	! Step 3 - Loop over the number of bins
	do k = 1, N_bins
		!%---- Likelihood of the frequency-magnitude distribution
		!%---- GR(m) x Q(m)
		if (mag_bins(k) .ne. 0) then
		
			!%- terme M1
			m1 = Mmin + delta_bin * (k-1)
			q_m1 = 0.5 + 0.5 * erf( (m1-mu_i)/(sqrt(2.0)*sigma_i))
			terme_m1 =  exp(-beta_i*(m1-Mmin))*q_m1
			
			!%- terme M2
			m2 = Mmin + delta_bin * k	
			q_m2 = 0.5 + 0.5 * erf( (m2-mu_i)/(sqrt(2.0)*sigma_i))
			terme_m2 =  exp(-beta_i*(m2-Mmin))*q_m2
			
			!%- terme 3	
			m1_t3 = m1 + (sigma_i**2 * beta_i)
			q_m1_t3 = 0.5 + 0.5 * erf( (m1_t3-mu_i)/(sqrt(2.0)*sigma_i))
			m2_t3 = m2 + (sigma_i**2 * beta_i)
			q_m2_t3 = 0.5 + 0.5 * erf( (m2_t3-mu_i)/(sqrt(2.0)*sigma_i))
			terme_3 = exp(((sigma_i**2)*(beta_i**2))/2 - (beta_i*(mu_i-Mmin))) *( q_m2_t3 - q_m1_t3)
			
			!%---- Likelihood formula (from Marsan private comm see Laporte et al.,2025)
        		Phi_bin =  - terme_m2 + terme_m1 + terme_3  ! Likelihood of the bin [M1-M2]
	    		Likelihood = Phi_bin/(Phi_cste*delta_bin)
	    		
	    		!%--- Loglikelihood = sum of the likelihood in each bin multiplied by the number of magnitudes in the bin
			!%--- Consider all cases where Likelihood = 0 
			if (Phi_bin .lt. 0) then
			    	LogLikelihood = -99999999
			else 
			    	LogLikelihood = LogLikelihood + mag_bins(k)* dlog(Likelihood) 
	    	        end if 
		else 
			LogLikelihood = LogLikelihood
		end if
		if ((k.eq.(N_bins-1)).and.(LogLikelihood.eq.0)) then
			LogLikelihood = -99999999 ! if not working put back -999999
		end if
  	end do
  	
  	Likelihood_FMD_MC_3D(i) = LogLikelihood
	
	! For reconstructing the marginals afterward
	! Keeping in mind the values explored
	index_beta(i) = ceiling(random_samples(i,1)*N_b)
	index_mu(i) = ceiling(random_samples(i,2)*N_m)
	index_sigma(i) = ceiling(random_samples(i,3)*N_s)
    end do
    

    Maxloglike = maxval(Likelihood_FMD_MC_3D)
    Likelihood_FMD_MC_3D = exp(Likelihood_FMD_MC_3D - Maxloglike)
    !      Interpolation of B-Mu-Sigma marginals - projection into the respective arrays !
    !----------------------------------------------------------------------------------
    
    do i = 1, n_draws
    	    ib = index_beta(i) 
    	    im = index_mu(i) 
    	    is = index_sigma(i) 
	    Beta_posterior(ib) =  Beta_posterior(ib) + Likelihood_FMD_MC_3D(i)
	    Mu_posterior(im) = Mu_posterior(im) + Likelihood_FMD_MC_3D(i)
	    Sigma_posterior(is) = Sigma_posterior(is) + Likelihood_FMD_MC_3D(i)
    end do

    ! Double step of normalisation of marginals
    ! Divides by the number of occurrence in each bin
    do ib =  1, size(Beta_prior) 
    	counts_beta(ib) = count(index_beta == ib)
    	Beta_posterior(ib) =  Beta_posterior(ib)/counts_beta(ib)
    end do 
    
    do im = 1, size(Mu_prior)
    	counts_mu(im) = count(index_mu == im)
    	Mu_posterior(im) = Mu_posterior(im)/counts_mu(im)
    end do
    
    do is = 1, size(Sigma_prior)
    	counts_sigma(is) = count(index_sigma == is)
    	Sigma_posterior(is) = Sigma_posterior(is)/counts_sigma(is)
    end do
    
    ! Divides by the integral of the prior
    Beta_posterior = Beta_posterior / (Beta_int * sum(Beta_posterior))
    Mu_posterior = Mu_posterior/(Mu_int * sum(Mu_posterior))
    Sigma_posterior = Sigma_posterior/(Sigma_int * sum(Sigma_posterior))
    return 
end subroutine Likelihood_FMD_random





























subroutine Likelihood_FMD_random_DEATH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
					Beta_prior, Mu_prior, Sigma_prior, cpert, New_Nc, T_bins, &
                                        New_model_idata, New_model_iprop, & 
                                        Acc_Likelihood, Acc_Maxloglike, &
                                        Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        New_Likelihood, New_Maxloglike, &
                                        New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)

! C ----------------------------------------------------------------------
! C   Compute the 3D posterior distribution of the Frequency-Magnitude Distribution 
! C   considering a detection function (Error function) and a Gutenberg Richter law
! C   ONLY in the segments that have changed during a "death" of a discontinuity
! C
! C   Inputs : 
! C    . mag_obs    : list of magnitudes in the entire input catalog {array or list}
! C    . mag_bins    : number of magnitudes/bin (/subset for subfunctions) {array or list}
! C    . delta_bin   : size of the magnitude binning (usually 0.1 for traditional catalogues) {float}
! C    . min_bin     : minimum magnitude bin (eq. minimum magnitude) for integration bounds {float}
! C    . n_draws     : number of random draw for Beta/Mu/Sigma exploration of the prior {integer}
! C    . Beta_prior    : exploration space of Beta {array or list} (careful : here we only retrieve the bounds of this space)
! C    . Mu_prior      : exploration space of Mu {array or list}
! C    . Sigma_prior   : exploration space of Sigma {array or list}
! C
! C    . cpert   : [Proposed model] layer perturbed (removed/added or moved) {integer}
! C    . New_Nc  : [Proposed model] Number of layers for the proposed model of discontinuities {integer}
! C    . T_bins  : Number of bins in the discretized time vector
! C
! C    . New_model_idata  : [Proposed model]
! C    . New_model_iprop  : [Proposed model]
! C
! C    . Acc_Likelihood  : [Current model] Sum of the 3D likelihood for Beta/Mu/Sigma {float}
! C    . Acc_maxloglike   : [Current model] Maximum Loglikelihoodof the 3D posterior distribution {float}
! C
! C    . Acc_Beta_T_posterior  : [Current model] Beta posterior for each time of the discrete time vector {N_b x T_bins array}
! C    . Acc_Mu_T_posterior    : [Current model] Mu posterior for each time of the discrete time vector   {N_m x T_bins array}
! C    . Acc_Sigma_T_posterior : [Current model] Sigma posterior for each time of the discrete time vector{N_s x T_bins array}
! C
! C   Outputs :
! C
! C    . New_Likelihood : [Proposed model] Sum of the 3D likelihood for Beta/Mu/Sigma {float}
! C    . New_Maxloglike : [Proposed model] Maximum Loglikelihoodof the 3D posterior distribution {float} 
! C
! C    . New_Beta_T_posterior  : [Proposed model] Beta posterior for each time of the discrete time vector {N_b x T_bins array}
! C    . New_Mu_T_posterior    : [Proposed model] Mu posterior for each time of the discrete time vector   {N_m x T_bins array}
! C    . New_Sigma_T_posterior : [Proposed model] Sigma posterior for each time of the discrete time vector{N_s x T_bins array}
! C ----------------------------------------------------------------------
    interface
        subroutine Likelihood_FMD_random(mag_bins, delta_bin, min_bin, n_draws, Beta_prior , Mu_prior, Sigma_prior, &
					Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
		  integer, intent(in) :: mag_bins(:)
		  real*8, intent(in)  :: delta_bin
		  real*8, intent(in)  :: min_bin
		  integer, intent(in) :: n_draws
		  real*8, intent(in)  :: Beta_prior(:), Mu_prior(:), Sigma_prior(:)
		  real*8, intent(out) :: Likelihood_FMD_MC_3D(n_draws)
		  real*8, intent(out) :: Beta_posterior(size(Beta_prior)), Mu_posterior(size(Mu_prior)), Sigma_posterior(size(Sigma_prior))
		  real*8, intent(out) :: Maxloglike
        end subroutine Likelihood_FMD_random
    end interface

    real*8, intent(in) :: mag_obs(:)
    integer,intent(in) :: mag_bins(:)
    real*8, intent(in) :: delta_bin
    real*8, intent(in) :: min_bin
    integer,intent(in):: n_draws
    real*8, intent(in):: Beta_prior(:), Mu_prior(:), Sigma_prior(:)
    integer,intent(in):: cpert, T_bins
    integer,intent(in):: New_Nc
    integer,intent(in):: New_model_idata(:), New_model_iprop(:)
    real*8, intent(in):: Acc_Likelihood(:), Acc_Maxloglike(:)
    real*8, intent(in):: Acc_Beta_T_posterior(:,:), Acc_Mu_T_posterior(:,:), Acc_Sigma_T_posterior(:,:)
    
    
    real*8,intent(out):: New_Likelihood(:), New_Maxloglike(:)
    real*8,intent(out):: New_Beta_T_posterior(:,:), New_Mu_T_posterior(:,:), New_Sigma_T_posterior(:,:)

    integer :: j, k, ma
    integer :: N_m, N_b, N_s
    integer :: N_data_subset, N_bins
    integer :: flag_ndata(New_Nc + 1)
    integer :: a_death, b_death, a_prop, b_prop
    
    integer :: bin_index
    integer :: mag_bins_subset(size(mag_bins))
    
    real*8 :: Likelihood_FMD_MC_3D(n_draws)
    real*8 :: Maxloglike
    real*8 :: Beta_posterior(size(Beta_prior)), Mu_posterior(size(Mu_prior)), Sigma_posterior(size(Sigma_prior))

    ! Initialize key variables
    N_bins = size(mag_bins)
    N_m = size(Mu_prior)
    N_s = size(Sigma_prior)
    N_b = size(Beta_prior)

    flag_ndata = 0
    
    New_Likelihood = 0.0
    New_Maxloglike = 0.0
    write(*,*) ">"
    write(*,*) ". Proposed model idata:", New_model_idata(1:New_Nc)
    write(*,*) ". Poposed model prop :", New_model_iprop(1:New_Nc)
    write(*,*) ">"
    write(*,*) ".. Death"
    write(*,*) ".. Id of the concerned discontinuity:", cpert
    
    ! Loop over the proposed number of temporal subsets (number of proposed discontinuities +1)
    do j=1, New_Nc+1
		! First get the indexes of each discontinuity in the data and in the discretized time vector									
	      	if (j==1) then								
	    		flag_ndata(1) = New_model_idata(1) 
	    		
	    		a_death = 1 
	      		b_death = a_death + flag_ndata(1) - 1
	      		
	      		a_prop = 1
	      		b_prop = New_model_iprop(j)

    		else if (j == New_Nc+1) then
    			flag_ndata(j) = size(mag_obs) - New_model_idata(j-1)  
    		
    			a_death = a_death + flag_ndata(j-1)
	    		b_death = a_death + flag_ndata(j) - 1
    		
    			a_prop = New_model_iprop(j-1)
      			b_prop = T_bins

    		else
    			flag_ndata(j) = New_model_idata(j) - New_model_idata(j-1) 
    	
    			a_death = a_death + flag_ndata(j-1)
	    		b_death = a_death + flag_ndata(j) - 1

	      		a_prop = New_model_iprop(j-1)
	      		b_prop = New_model_iprop(j)
    		end if
    		
		N_data_subset = flag_ndata(j)
		! 1 - Before the dead discontinuity : CC<CPERT
		! subset unchanged : keeping previous likelihood
		if (j .lt. cpert) then  				
			New_Likelihood(j) = Acc_Likelihood(j)		
			New_Maxloglike(j) = Acc_Maxloglike(j)
			
			do k = a_prop, b_prop
					New_Beta_T_posterior(k, :) = Acc_Beta_T_posterior(k,  :)
					New_Mu_T_posterior(k,  :) = Acc_Mu_T_posterior(k,  :)
					New_Sigma_T_posterior(k , :) = Acc_Sigma_T_posterior(k,  :)
			end do
		! 2 - Dead discontinuity CC = CPERT
		else if (j==cpert) then 
			! 2.a - If no events in the new/largest subset							
			if (flag_ndata(j)==0) then 						
    	    			!write(*,*) '.. No events in the new partition'
    	    			New_Likelihood(j) = 1.0						
    				New_Maxloglike(j) = 0.0
    				do k = a_prop, b_prop
						New_Beta_T_posterior(k,  :) = 1.0/N_b
						New_Mu_T_posterior(k, :) = 1.0/N_m
						New_Sigma_T_posterior(k, :) = 1.0/N_s
				end do
			! 2.b - If events in the new/largest subset 
			else 									! 2.b - If events in the new/largest subset
				!write(*,*) "... We recompute the likelihood on the segment (Timegrid)",  a_death, b_death	

				! 2.c Count events in magnitude bins
				mag_bins_subset = 0 
				do ma = 1, N_data_subset
					bin_index = int((mag_obs(a_death+ma-1)-min_bin) / delta_bin) + 1
					if (bin_index >= 1 .and. bin_index <= N_bins) then
						mag_bins_subset(bin_index) = mag_bins_subset(bin_index) + 1
					end if
				end do
				! 2.d Recompute Likelihood in the new subset
				call Likelihood_FMD_random(mag_bins_subset, delta_bin, min_bin, n_draws, Beta_prior, Mu_prior, Sigma_prior, &
						Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
				
				New_Likelihood(j) = sum(Likelihood_FMD_MC_3D)
				New_Maxloglike(j) = Maxloglike

				do k = a_prop, b_prop
					New_Beta_T_posterior(k, :) = Beta_posterior
					New_Mu_T_posterior(k,  :) =  Mu_posterior
					New_Sigma_T_posterior(k, :) = Sigma_posterior
				end do
			end if
		! 3 - After the dead discontinuity
		! subset cc = old subset cc+1
		else if (j .gt.cpert) then							
			New_Likelihood(j) = Acc_Likelihood(j+1)						
			New_Maxloglike(j) = Acc_Maxloglike(j+1)
			do k = a_prop, b_prop
					New_Beta_T_posterior(k, :) = Acc_Beta_T_posterior(k,  :)
					New_Mu_T_posterior(k,  :) = Acc_Mu_T_posterior(k,  :)
					New_Sigma_T_posterior(k , :) = Acc_Sigma_T_posterior(k,  :)
			end do
		end if
		write(*,*) "... Likelihood of the partition", j, a_prop, b_prop, N_data_subset, &
								New_Beta_T_posterior(a_prop, 2), New_Maxloglike(j)
	end do
	write(*,*)".. Likelihood :       id        begin       end       n_data	    Like_b(sample)     MaxLogLike      "
	New_Likelihood(New_Nc+2) = Acc_Likelihood(New_Nc+3) ! We put -99999 on the last segment
	return

end subroutine Likelihood_FMD_random_DEATH








subroutine Likelihood_FMD_random_BIRTH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
					Beta_prior, Mu_prior, Sigma_prior, cpert, New_Nc, T_bins, &
                                        New_model_idata, New_model_iprop, & 
                                        Acc_Likelihood, Acc_Maxloglike, &
                                        Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        New_Likelihood, New_Maxloglike, &
                                        New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)

! C ----------------------------------------------------------------------
! C   Compute the 3D posterior distribution of the Frequency-Magnitude Distribution 
! C   considering a detection function (Error function) and a Gutenberg Richter law
! C   ONLY in the segments that have changed during a "birth" of a discontinuity
! C
! C   Inputs : 
! C    . mag_obs    : list of magnitudes in the entire input catalog {array or list}
! C    . mag_bins    : number of magnitudes/bin (/subset for subfunctions) {array or list}
! C    . delta_bin   : size of the magnitude binning (usually 0.1 for traditional catalogues) {float}
! C    . min_bin     : minimum magnitude bin (eq. minimum magnitude) for integration bounds {float}
! C    . n_draws     : number of random draw for Beta/Mu/Sigma exploration of the prior {integer}
! C    . Beta_prior    : exploration space of Beta {array or list} (careful : here we only retrieve the bounds of this space)
! C    . Mu_prior      : exploration space of Mu {array or list}
! C    . Sigma_prior   : exploration space of Sigma {array or list}
! C
! C    . cpert   : [Proposed model] layer perturbed (removed/added or moved) {integer}
! C    . New_Nc  : [Proposed model] Number of layers for the proposed model of discontinuities {integer}
! C    . T_bins  : Number of bins in the discretized time vector
! C
! C    . New_model_idata  : [Proposed model]
! C    . New_model_iprop  : [Proposed model]
! C
! C    . Acc_Likelihood  : [Current model] Sum of the 3D likelihood for Beta/Mu/Sigma {float}
! C    . Acc_maxloglike   : [Current model] Maximum Loglikelihoodof the 3D posterior distribution {float}
! C
! C    . Acc_Beta_T_posterior  : [Current model] Beta posterior for each time of the discrete time vector {N_b x T_bins array}
! C    . Acc_Mu_T_posterior    : [Current model] Mu posterior for each time of the discrete time vector   {N_m x T_bins array}
! C    . Acc_Sigma_T_posterior : [Current model] Sigma posterior for each time of the discrete time vector{N_s x T_bins array}
! C
! C   Outputs :
! C
! C    . New_Likelihood : [Proposed model] Sum of the 3D likelihood for Beta/Mu/Sigma {float}
! C    . New_Maxloglike : [Proposed model] Maximum Loglikelihoodof the 3D posterior distribution {float} 
! C
! C    . New_Beta_T_posterior  : [Proposed model] Beta posterior for each time of the discrete time vector {N_b x T_bins array}
! C    . New_Mu_T_posterior    : [Proposed model] Mu posterior for each time of the discrete time vector   {N_m x T_bins array}
! C    . New_Sigma_T_posterior : [Proposed model] Sigma posterior for each time of the discrete time vector{N_s x T_bins array}
! C ----------------------------------------------------------------------
    interface
        subroutine Likelihood_FMD_random(mag_bins, delta_bin, min_bin, n_draws, Beta_prior , Mu_prior, Sigma_prior, &
					Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
		  integer, intent(in) :: mag_bins(:)
		  real*8, intent(in)  :: delta_bin
		  real*8, intent(in)  :: min_bin
		  integer, intent(in) :: n_draws
		  real*8, intent(in)  :: Beta_prior(:), Mu_prior(:), Sigma_prior(:)
		  real*8, intent(out) :: Likelihood_FMD_MC_3D(n_draws)
		  real*8, intent(out) :: Beta_posterior(size(Beta_prior)), Mu_posterior(size(Mu_prior)), Sigma_posterior(size(Sigma_prior))
		  real*8, intent(out) :: Maxloglike
        end subroutine Likelihood_FMD_random
    end interface

    real*8, intent(in) :: mag_obs(:)
    integer,intent(in) :: mag_bins(:)
    real*8, intent(in) :: delta_bin
    real*8, intent(in) :: min_bin
    integer,intent(in):: n_draws
    real*8, intent(in):: Beta_prior(:), Mu_prior(:), Sigma_prior(:)
    integer,intent(in):: cpert, T_bins
    integer,intent(in):: New_Nc
    integer,intent(in):: New_model_idata(:), New_model_iprop(:)
    real*8, intent(in):: Acc_Likelihood(:), Acc_Maxloglike(:)
    real*8, intent(in):: Acc_Beta_T_posterior(:,:), Acc_Mu_T_posterior(:,:), Acc_Sigma_T_posterior(:,:)
    
    
    real*8,intent(out):: New_Likelihood(:), New_Maxloglike(:)
    real*8,intent(out):: New_Beta_T_posterior(:,:), New_Mu_T_posterior(:,:), New_Sigma_T_posterior(:,:)

    integer :: j,k , ma
    integer :: N_m, N_b, N_s
    integer :: N_data_subset, N_bins
    integer :: flag_ndata(New_Nc + 1)
    integer :: a_birth, b_birth, a_prop, b_prop
    
    integer :: bin_index
    integer :: mag_bins_subset(size(mag_bins))
    
    real*8 :: Likelihood_FMD_MC_3D(n_draws)
    real*8 :: Maxloglike
    real*8 :: Beta_posterior(size(Beta_prior)), Mu_posterior(size(Mu_prior)), Sigma_posterior(size(Sigma_prior))

    ! Initialize key variables
    N_bins = size(mag_bins)
    N_m = size(Mu_prior)
    N_s = size(Sigma_prior)
    N_b = size(Beta_prior)

    flag_ndata = 0
    
    New_Likelihood = 0.0
    New_Maxloglike = 0.0
    write(*,*) ">"
    write(*,*) ". Proposed model idata:", New_model_idata(1:New_Nc)
    write(*,*) ". Proposed model prop :", New_model_iprop(1:New_Nc)
    write(*,*) ">"
    write(*,*) ".. Birth"
    write(*,*) ".. Id of the concerned discontinuity:", cpert
    
    ! Loop over the proposed number of temporal subsets (number of proposed discontinuities +1)
    do j = 1, New_Nc+1
		! First get the indexes of each discontinuity in the data and in the discretized time vector									
	      	if (j==1) then								
	    		flag_ndata(1) = New_model_idata(1) 
	    		
	    		a_birth = 1 
	      		b_birth = a_birth + flag_ndata(1) - 1
	      		
	      		a_prop = 1
	      		b_prop = New_model_iprop(j)

    		else if (j == New_Nc+1) then
    			flag_ndata(j) = size(mag_obs) - New_model_idata(j-1)  
    		
    			a_birth = a_birth + flag_ndata(j-1)
	    		b_birth = a_birth + flag_ndata(j) - 1
    		
    			a_prop = New_model_iprop(j-1)
      			b_prop = T_bins

    		else
    			flag_ndata(j) = New_model_idata(j) - New_model_idata(j-1) 
    	
    			a_birth = a_birth + flag_ndata(j-1)
	    		b_birth = a_birth + flag_ndata(j) - 1

	      		a_prop = New_model_iprop(j-1)
	      		b_prop = New_model_iprop(j)
    		end if
    		
		N_data_subset = flag_ndata(j)
		! 1 - Before the new discontinuity : CC<CPERT
		! subset unchanged : keeping previous likelihood
		if (j .lt. cpert) then  				
			New_Likelihood(j) = Acc_Likelihood(j)		
			New_Maxloglike(j) = Acc_Maxloglike(j)
			
			do k = a_prop, b_prop
					New_Beta_T_posterior(k , :) = Acc_Beta_T_posterior(k,  :)
					New_Mu_T_posterior(k,  :) = Acc_Mu_T_posterior(k,  :)
					New_Sigma_T_posterior( k , :) = Acc_Sigma_T_posterior(k,  :)
			end do
		! 2- Subset concerned by a change CC=CPERT
		else if (j==cpert) then 
			! 2.a - If no events in the new/smallest subset							
			if (flag_ndata(j)==0) then 						
    	    			!write(*,*) '.. No events in the new partition'
    	    			New_Likelihood(j) = 1.0						
    				New_Maxloglike(j) = 0.0
    				N_data_subset = 0
    				do k = a_prop, b_prop
						New_Beta_T_posterior(k,  :) = 1./N_b
						New_Mu_T_posterior(k, :) = 1./N_m
						New_Sigma_T_posterior(k, :) = 1./N_s
				end do
			! 2.b - If events in the new/smallest subset
			else 									
				!write(*,*) "... We recompute the likelihood on the segment (Timegrid)",  a_birth, b_birth	

				! 2.c Count events in magnitude bins
				mag_bins_subset = 0 
				do ma = 1, N_data_subset
					bin_index = int((mag_obs(a_birth+ma-1)-min_bin) / delta_bin) + 1
					if (bin_index >= 1 .and. bin_index <= N_bins) then
						mag_bins_subset(bin_index) = mag_bins_subset(bin_index) + 1
					end if
				end do
				! 2.d Recompute Likelihood in the new subset
				call Likelihood_FMD_random(mag_bins_subset, delta_bin, min_bin, n_draws, Beta_prior, Mu_prior, Sigma_prior, &
						Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
				
				New_Likelihood(j) = sum(Likelihood_FMD_MC_3D)
				New_Maxloglike(j) = Maxloglike

				do k = a_prop, b_prop
					New_Beta_T_posterior(k , :) = Beta_posterior
					New_Mu_T_posterior(k,  :) =  Mu_posterior
					New_Sigma_T_posterior(k , :) = Sigma_posterior
				end do
			end if
			
		! 3- Subset concerned by a change CC=CPERT+1
		else if (j==cpert+1) then 
			! 3.a - No events in the new subset									
			if (flag_ndata(j)==0) then 								
    				!write(*,*) '.. No events in the new partition'
    	    			New_Likelihood(j) = 1.0						
    				New_Maxloglike(j) = 0.0
    				N_data_subset = 0
    				do k = a_prop, b_prop
						New_Beta_T_posterior(k,  :) = 1./N_b
						New_Mu_T_posterior(k, :) = 1./N_m
						New_Sigma_T_posterior(k, :) = 1./N_s
				end do
			! 3.b - If events in the new/smallest subset
    			else 
    
				!write(*,*) "... We recompute the likelihood on the segment (Timegrid) ", a_birth, b_birth ! 3.b - At least one event in the new subset

				
														! 3.c Count events in magnitude bins
				! 3.c Count events in magnitude bins
				mag_bins_subset = 0 
				do ma = 1, N_data_subset
					bin_index = int((mag_obs(a_birth + ma-1)-min_bin) / delta_bin)+1
					if (bin_index >= 1 .and. bin_index <= n_bins) then
						mag_bins_subset(bin_index) = mag_bins_subset(bin_index) + 1
					end if
				end do	
				! 3.d Recompute Likelihood in the new subset
				call Likelihood_FMD_random(mag_bins_subset, delta_bin, min_bin, n_draws, Beta_prior, Mu_prior, Sigma_prior, &
							Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
				
				New_Likelihood(j) = sum(Likelihood_FMD_MC_3D)
				New_Maxloglike(j) = Maxloglike
				
				do k = a_prop, b_prop
					New_Beta_T_posterior(k, :) = Beta_posterior
					New_Mu_T_posterior(k, :) =  Mu_posterior
					New_Sigma_T_posterior(k, :) = Sigma_posterior
				end do
			end if
		! 3 - After the dead discontinuity
		! subset cc = old subset cc-1
		else if (j .gt.cpert+1) then							
			New_Likelihood(j) = Acc_Likelihood(j-1)						
			New_Maxloglike(j) = Acc_Maxloglike(j-1)
			do k = a_prop, b_prop
					New_Beta_T_posterior(k , :) = Acc_Beta_T_posterior(k,  :)
					New_Mu_T_posterior(k,  :) = Acc_Mu_T_posterior(k,  :)
					New_Sigma_T_posterior(k, :) = Acc_Sigma_T_posterior(k,  :)
			end do
		end if
		write(*,*) "... Likelihood of the partition", j, a_prop, b_prop, N_data_subset, &
								New_Beta_T_posterior(a_prop, 2), New_Maxloglike(j)
	end do
	write(*,*)".. Likelihood :       id        begin       end       n_data	    Like_b(sample)     MaxLogLike      "
	return

end subroutine Likelihood_FMD_random_BIRTH















function gammln(xx)
  !use par_mod, only: dp
  implicit none
  integer :: j
  real*8 :: x,tmp,ser,xx,gammln
  real*8 :: cof(6) = (/ &
       76.18009173, -86.50532033, 24.01409822, &
       -1.231739516, .120858003e-2, -.536382e-5    /)
  real*8 :: stp = 2.50662827465
  real*8 :: half = 0.5, one = 1.0, fpf = 5.5
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one
  do j=1,6
    x=x+one
    ser=ser+cof(j)/x
  end do
  gammln=tmp+log(stp*ser)
end function gammln



function gammp(a,x)
  implicit none
  real*8 :: a, x, gln, gamser, gammp, gammcf
  if(x .lt. 0. .or. a .le. 0.) then
     print*, 'gammp'
     stop
  end if
  if(x.lt.a+1.)then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=1.-gammcf
  endif
end function gammp



function gammq(a,x)
  implicit none
  real*8 :: a, x, gln, gamser, gammq, gammcf
  if(x.lt.0..or.a.le.0.) then
     print*, 'gammq'
     stop
  end if
  if(x.lt.a+1.)then
    call gser(gamser,a,x,gln)
    gammq=1.-gamser
  else
    call gcf(gammcf,a,x,gln)
    gammq=gammcf
  endif
end function gammq


subroutine gser(gamser,a,x,gln)
  implicit none
  integer :: n
  real*8 :: gamser, a, x, gln, ap, summ, del
  real*8, external :: gammln
  integer,parameter :: itmax=100
  real,parameter    :: eps=3.e-7
  gln=gammln(a)
  if(x.le.0.)then
    if(x.lt.0.) then
       print*, 'gser'
       stop
    end if
    gamser=0.
    return
  endif
  ap=a
  summ=1./a
  del=summ
  do n=1,itmax
    ap=ap+1.
    del=del*x/ap
    summ=summ+del
    if(abs(del).lt.abs(summ)*eps)go to 1
  end do
  print*, 'gser: a too large, itmax too small'
  stop
1   gamser=summ*exp(-x+a*log(x)-gln)
end subroutine gser



subroutine gcf(gammcf,a,x,gln)
  implicit none
  integer :: n
  real*8 :: gammcf, a, x, gln, gold, a0, a1, b0, b1, fac, an, anf, ana, g
  real*8, external :: gammln
  integer,parameter :: itmax=100
  real,parameter    :: eps=3.e-7
  gln=gammln(a)
  gold=0.
  a0=1.
  a1=x
  b0=0.
  b1=1.
  fac=1.
  do n=1,itmax
    an=real(n)
    ana=an-a
    a0=(a1+a0*ana)*fac
    b0=(b1+b0*ana)*fac
    anf=an*fac
    a1=x*a0+anf*a1
    b1=x*b0+anf*b1
    if(a1.ne.0.)then
      fac=1./a1
      g=b1*fac
      if(abs((g-gold)/g).lt.eps)go to 1
      gold=g
    endif
  end do
  print*, 'gcf: a too large, itmax too small'
  stop
1   gammcf=exp(-x+a*dlog(x)-gln)*g
end subroutine gcf





