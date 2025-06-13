program b_Bayesian_binned

  implicit none
  real , EXTERNAL    ::   ran3, gasdev
  !            ----------------------------------
  !              # BEGINING OF THE USER-DEFINED PARAMETERS
  !            ----------------------------------
!-----------------------------------------
! # Data 
!----------------------------------------- 
  character(LEN=*), parameter :: Magnitude_file = '../inputs/input.txt'	!path of the time-magnitude catalog
  
  integer, parameter :: N_data_max = 800000	! Number maximum of data
  integer, parameter :: T_bins = 200		! Size of the Discrete Time Vector
  real*8, parameter  :: delta_bin = 0.1		! Binning of the catalog			
  real*8, parameter  :: min_bin = 0		! Min bin < min magnitudes for integration lower bound
  real*8, parameter  :: max_bin = 9		! Max bin > max magnitudes for integration upper bound

!-----------------------------------------
! # Parameters of the Markov chain
!-----------------------------------------
  integer, parameter :: it_max = 2000		! Total number of McMC iterations
  integer, parameter :: it_burnin = 1000		! Burn-in iterations
  integer, parameter :: it_thin = 2		! Thinning for independency of accepted models
  integer, parameter :: it_std = 2000		! Number of iterations between two perturbation standard-deviation changes 
!-----------------------------------------
! # Definition of priors
!-----------------------------------------
  integer, parameter :: max_layers = 50 	! Threshold on maximum layers
  
  integer, parameter :: n_draws = 500000 			! Numbers of random draws using Monte-Carlo sampling over model priors
  
  integer, parameter :: Sigma_bins = 100			! Numbers of values explored in the prior uniform distribution if gridsearch
  real*8, parameter :: Sigma_min = 0.01, Sigma_max = 0.8 	! Bounds of the Sigma uniform prior

  integer, parameter :: Mu_bins = 100				! Numbers of values explored in the prior uniform distribution if gridsearch
  real*8, parameter :: Mu_min = 0, Mu_max = 5			! Bounds of the Mu uniform prior
  
  integer, parameter :: Beta_bins = 100				! Numbers of values explored in the prior uniform distribution if gridsearch
  real*8, parameter ::  Beta_min = 0.4*2.3, Beta_max = 1.7*2.3	! Bounds of the uniform prior of Beta parameter (b_value*ln(10))
  !-----------------------------------------
  ! # Definition of initialization
  !-----------------------------------------
  integer, parameter :: Nc_start_min = 3 	! Initial number of layers min
  integer, parameter :: Nc_start_max = 4 	! Initial number of layers max   
  integer, parameter :: Nc_stdi = 1		! Initial standard deviation for layer perturbation 
  !            ----------------------------------
  !              # END OF THE USER-DEFINED PARAMETERS
  !            ----------------------------------
!-----------------------------------------
! # Functions declaration
!-----------------------------------------		
interface
! A function for sorting a list of integers
  subroutine sort_shell(arr)						
  	integer, DIMENSION(:), INTENT(INOUT) :: arr
	integer :: i,j,inc,n
	REAL*8 :: v
  end subroutine sort_shell

! A function for sorting a list of reals   
  subroutine sort_shell_real(arr)					
  	real*8, DIMENSION(:), INTENT(INOUT) :: arr
	integer :: i,j,inc,n
	REAL*8 :: v
  end subroutine sort_shell_real
   
! A function that computes the posterior distribution of the Frequency-Magnitude distribution
! Using a random exploration of the model priors
  subroutine Likelihood_FMD_random(mag_bins, delta_bin, min_bin, n_draws, Beta_prior , Mu_prior, Sigma_prior, &
					Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
	  integer, intent(in) :: mag_bins(:)
	  real*8, intent(in) :: delta_bin
	  real*8, intent(in) :: min_bin
	  integer, intent(in):: n_draws
	  real*8, intent(in) :: Beta_prior(:), Mu_prior(:), Sigma_prior(:)
	  real*8, intent(out) :: Likelihood_FMD_MC_3D(n_draws)
	  real*8, intent(out) :: Beta_posterior(size(Beta_prior)), Mu_posterior(size(Mu_prior)), Sigma_posterior(size(Sigma_prior))
	  real*8, intent(out) :: Maxloglike
  end subroutine Likelihood_FMD_random

! A function that computes the posterior distribution of the Frequency-Magnitude distribution
! Using a gridsearch exploration of the model priors
  subroutine Likelihood_FMD_gridsearch(mag_bins, delta_bin, min_bin, Beta_prior, Mu_prior, Sigma_prior, &
					Likelihood_FMD_3D, Maxloglike) 
	  integer, intent(in) :: mag_bins(:)
	  real*8, intent(in) :: delta_bin
	  real*8, intent(in) :: min_bin
	  real*8, intent(in) :: Mu_prior(:)
	  real*8, intent(in) :: Beta_prior(:)  
	  real*8, intent(in) :: Sigma_prior(:)
	  real*8, intent(out) :: Likelihood_FMD_3D(size(Mu_prior), size(Beta_prior),size(Sigma_prior))
	  real*8, intent(out) :: Maxloglike
  end subroutine Likelihood_FMD_gridsearch

! Function that recomputes the likelihood only on changed segments in case of a death case
  subroutine Likelihood_FMD_random_DEATH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
					Beta_prior, Mu_prior, Sigma_prior, cpert, New_Nc, T_bins, &
                                        New_model_idata, New_model_iprop, & 
                                        Acc_Likelihood, Acc_Maxloglike, &
                                        Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        New_Likelihood, New_Maxloglike, &
                                        New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)                                  
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
  end subroutine Likelihood_FMD_random_DEATH

! Function that recomputes the likelihood only on changed segments in case of a birth case    
  subroutine Likelihood_FMD_random_BIRTH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
					Beta_prior, Mu_prior, Sigma_prior, cpert, New_Nc, T_bins, &
                                        New_model_idata, New_model_iprop, & 
                                        Acc_Likelihood, Acc_Maxloglike, &
                                        Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        New_Likelihood, New_Maxloglike, &
                                        New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)
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
  end subroutine Likelihood_FMD_random_BIRTH
end interface  
   
!-----------------------------------------
! # DECLARATIONS OF VARIABLES 
!-----------------------------------------
integer :: ra							! seed of our run (change if we want a different run)
!-----------------------------------------
! # DECLARATIONS OF VARIABLES DATA
!-----------------------------------------
real*8 :: mag_obs0(N_data_max), times0(N_data_max)  		! max data
real*8, allocatable :: mag_obs(:), times(:) 			! allocated data
real*8, allocatable ::  mag_obs_subset(:)			! allocated data in a subset 
! Histogram of magnitudes 
integer, allocatable::  mag_bins(:), mag_bins_subset(:)		! magnitude binning in the full dataset / by data subsets
integer :: N_data_mag, N_data_subset				! number of data in the full dataset / by subsets
integer :: bin_index, N_bins					! index, number of bins
real*8 :: time_max, time_min, time_lim				! min/max times, bounds for random draw

!-----------------------------------------
! # DECLARATIONS OF VARIABLES PRIORS
!-----------------------------------------
! FMD model parameters : beta/mu/sigma
real*8 :: Beta_prior(Beta_bins), Mu_prior(Mu_bins),  Sigma_prior(Sigma_bins) ! Vector of the prior distributions for Beta, Mu, Sigma (gridsearch)

!----------------------------------------
! # DECLARATIONS OF VARIABLES CONDITIONAL POSTERIOR P(W|D,T)
!-----------------------------------------
real*8 :: Likelihood_all, Maxloglike_all	! Likelihood/Maxloglikelihood of the full dataset (no discontinuities)

! Posterior FMD
real*8 :: Likelihood_FMD_3D(size(Mu_prior), size(Beta_prior),size(Sigma_prior))	  ! output of Likelihood_FMD_gridsearch
real*8 :: Likelihood_FMD_MC_3D(n_draws)						  ! output of Likelihood_FMD_random
real*8 :: Maxloglike

real*8 :: Init_Likelihood(max_layers+1), Init_Maxloglike(max_layers+1) 		! [Init] Likelihood density of FMD for each data subset 
real*8 :: Acc_Likelihood(max_layers+1), Acc_Maxloglike(max_layers+1) 		! [Current] Likelihood density of FMD for each data subset 
real*8 :: New_Likelihood(max_layers+1) , New_Maxloglike(max_layers+1)		! [Proposed] Likelihood density of FMD for each data subset 	
real*8 :: Likelihood_bis(max_layers+1), Maxloglike_bis(max_layers+1) 		! [Intermediary] for move cases Likelihood density of FMD for each data subset 		

! 2D Marginal FMD distributions obtained from Likelihood_gridsearch
real*8 :: Mu_Beta_posterior(size(Mu_prior), size(Beta_prior)) 
real*8 :: Mu_Sigma_posterior(size(Mu_prior), size(Sigma_prior))
real*8 :: Beta_Sigma_posterior(size(Beta_prior), size(Sigma_prior)) 

! 1D Marginals FMD
real*8 :: Beta_posterior(Beta_bins), Mu_posterior(Mu_bins), Sigma_posterior(Sigma_bins)	! Marginals estimated by random exploration: Likelihood_FMD_random output

real*8 :: Init_Beta_T_posterior(T_bins, Beta_bins),  Beta_T_posterior_all(T_bins,Beta_bins) 		! Beta Marginal posterior (at each time bin) : [init]/[all]
real*8 :: New_Beta_T_posterior(T_bins, Beta_bins),  Acc_Beta_T_posterior(T_bins, Beta_bins)	! Beta Marginal posterior (at each time bin) : [Proposed]/[Current]
real*8 :: Beta_T_posterior_bis(T_bins, Beta_bins), Sum_Beta_T_posterior(T_bins, Sigma_bins)	! Beta Marginal posterior (at each time bin) : [Intermediary] for move cases/ [Sum] = output marginal

real*8 :: Init_Mu_T_posterior(T_bins, Mu_bins), Mu_T_posterior_all(T_bins,Mu_bins) 			! Mu Marginal posterior (at each time bin) : [init]/[all]
real*8 :: New_Mu_T_posterior(T_bins, Mu_bins), Acc_Mu_T_posterior(T_bins, Mu_bins)		! Mu Marginal posterior (at each time bin) : [Proposed]/[Current]
real*8 :: Mu_T_posterior_bis(T_bins, Mu_bins), Sum_Mu_T_posterior(T_bins, Sigma_bins)	! Mu Marginal posterior (at each time bin) :[Intermediary] for move cases/ [Sum] = output marginal

real*8 :: Init_Sigma_T_posterior(T_bins, Sigma_bins), Sigma_T_posterior_all(T_bins,Sigma_bins) 	! Sigma Marginal posterior (at each time bin) : [init]/[all]
real*8 :: New_Sigma_T_posterior(T_bins, Sigma_bins),  Acc_Sigma_T_posterior(T_bins, Sigma_bins)	! Sigma Marginal posterior (at each time bin) : [Proposed]/[Current]
real*8 :: Sigma_T_posterior_bis(T_bins, Sigma_bins), Sum_Sigma_T_posterior(T_bins, Sigma_bins)	! Sigma Marginal posterior (at each time bin) : [Intermediary] for move cases/ [Sum] = output marginal

!-----------------------------------------
! # DECLARATIONS OF VARIABLES MARGINAL POSTERIOR P(T|D)
!-----------------------------------------
real*8 :: Nw 											! Weight for parcimony 

real*8 :: Log_Posterior_T_gridsearch(it_max)							! Marginal posterior [Init]/[Gridsearch]											
real*8 :: New_Log_Posterior_T, Acc_Log_Posterior_T(it_max)					! Marginal posterior [Proposed]/[Current]
real*8 :: Log_Posterior_T_all 									! Marginal posterior (2.2.3 b-Bayesian) for 0 layers
!-----------------------------------------
! # DECLARATIONS OF VARIABLES METROPOLIS-HASTINGS
!-----------------------------------------
integer ::  flag_unaccepted, flag_ndata(max_layers)	! some models are ill constructed - reject directly
real*8 :: alpha 					! condition birth/death/pert reversible jump
real*8 :: u, condition					! condition acceptance Metropolis Hastings
real*8:: New_std					!  Methode T.Santos: std adaptation according to AR_moves
integer :: N_thin					! Pour independance des parametres : garder les modeles acceptes tous les dix iterations

!-----------------------------------------
! # DECLARATIONS OF VARIABLES TEMPORAL MODEL 
!-----------------------------------------
! Layers 
integer :: Nc, New_Nc				! Number of Layers for temporal models [Current]/[Proposed]
integer :: Nc_total(it_max)			! Store the number of layers at each iteration of the McMC

integer :: model_T(T_bins), Sum_accepted_T(T_bins)	! Array of accepted discontinuities/ Stack of acc discontinuities (0=no disc/1=disc)
integer :: New_model_T(T_bins)				! [Proposed] Array of proposed discontinuities

real*8 :: model_times(max_layers+1)		! [Current] Array of the times of discontinuities (different of catalog times)
integer :: model_ibins(max_layers+1)		! [Current] Array of indexed of these times in the T_bins vector
integer :: model_idata(max_layers+1)		! [Current] Array of indexed of these times in the data

real*8 :: New_model_times(max_layers+1)		! [Proposed] Array of the times of discontinuities (different of catalog times)
integer :: New_model_ibins(max_layers+1)	! [Proposed] Array of indexed of these times in the T_bins vector
integer :: New_model_idata(max_layers+1)	! [Proposed] Array of indexed of these times in the data 

integer ::  New_model_idata_bis(max_layers+1) 	! [Proposed] Array of indexed of these times in the T_bins vector (intermediary for a move case)
integer ::  New_model_ibins_bis(max_layers+1)	! [Proposed] Array of indexed of these times in the data (intermediary for a move case)

integer :: a_data, b_data, a_ibins, b_ibins	! Indexes for cutting catalogs between two time discontinuities
integer ::  index_data, index_ibins, index_old  ! Indexes for New_model_idata/ New_model_ibins and index_ibins of the previous discontinuity case move
integer :: cpert				! Index of the discontinuity to be disturbed by a move/a death
real*8 :: t_old, t_new, t_init, t_std		! Discontinuity time once disturbed according to t_std, minimum disturbance time, standard time deviation

!-----------------------------------------
! # DECLARATIONS OF VARIABLES Acceptance Rates
!-----------------------------------------
integer :: init_full						! Only to compare gridsearch and random exploration during initialization
integer :: N_move, N_birth, N_death , store_death		! Count number of move/birth /death 
integer :: N_move_accept, N_birth_accept, N_death_accept	! Count number of move/birth /death accepted
real*8 ::  AR_move, AR_birth, AR_death			! Acceptance Rate move/birth/death


!-----------------------------------------
! # Rest
!-----------------------------------------
! iteratives
integer :: i, j, k, it
integer :: io
! Computation time
real :: t1, t2,t3,t4, mean_t_it, t_it(it_max), comp_time, progress					
integer :: line,column ! for outputs writting

!-----------------------------------------
! Variables initialisation
!-----------------------------------------
!open(unit=90,file='seed.txt',status='old')
!read(90,*) ra
close(90)
!WRITE(6,*) 'Seed value=',ra
ra=-1	! seed of our run (change if we want a different run)

init_full = 0 ! For test if full = random during 1st iteration (init)

! Init temporal model (9999999 default)
model_times = 9999999
model_ibins = 9999999
model_idata = 9999999

New_model_times = 9999999
New_model_ibins = 9999999
New_model_idata = 9999999




Sum_accepted_T = 0
model_T = 0 

N_move=0
N_death=0
N_birth=0

N_move_accept = 0
N_death_accept = 0
N_birth_accept = 0


index_data =0
flag_ndata = 0

Log_Posterior_T_gridsearch = 0.0
New_Log_Posterior_T = 0.0
Acc_Log_Posterior_T = 0.0


N_thin = 0 

!-------------------------------------------------
! Etape 0.A : Affectation des vecteurs Muvec , Bvec, Sigma_prior
! Equivalent as np.linspace 
!-------------------------------------------------
do i = 1, Beta_bins
  Beta_prior(i) =  (i-1) * (Beta_max-Beta_min)/(Beta_bins-1) + Beta_min
end do
do i = 1, Mu_bins
  Mu_prior(i) = (i-1) * (Mu_max-Mu_min)/(Mu_bins-1) + Mu_min
end do

do i = 1, Sigma_bins
  Sigma_prior(i) =  (i-1) * (Sigma_max-Sigma_min)/(Sigma_bins-1) + Sigma_min
end do

!-----------------------------------------
! Etape 0.B : Read the magnitude data
! The magnitude catalog must only have two colums
! C1 : Time between 0 and TMAX in minuts or hour or day 
! C2 : Event magnitude
!-----------------------------------------

open(50,file= Magnitude_file ,status='old')
    N_data_mag = 0
    do i =1, N_data_max
	      read(50,*,IOSTAT=io)times0(i),mag_obs0(i)	
	      if (io > 0) then				
			    stop "Check input.  Something was wrong"
	      elseif (io < 0) then
			    exit
	      else
	      	N_data_mag = N_data_mag + 1
	      endif
	      if (i == N_data_max) stop "number of magnitudes data >= ndata_magmax"
    enddo
 close(50)! close the file

!-----------------------------------------
! allocate de tous les tableaux qui dependent de ndata 
allocate(mag_obs(N_data_mag))
allocate(times(N_data_mag))
write(*,*) 'Size of the catalog :', N_data_mag

!-----------------------------------------
mag_obs = mag_obs0(1:N_data_mag)
times =  times0(1:N_data_mag)
write(*,*) 'Minimum magnitude :', minval(mag_obs)
write(*,*) 'Maximum magnitude :', maxval(mag_obs)
write(*,*)

! Re-write input in the output folder for post-processing simplifications
open(unit=61, file="../outputs/input.txt",  status="replace", action="write")
	do i = 1, N_data_mag
		write(61, *) times(i),mag_obs(i)
	end do
 close(61)
 
!-----------------------------------------
! INITIALISATION BINNING
  N_bins = ((max_bin - min_bin)/delta_bin) +1
  allocate(mag_bins(N_bins))
  allocate(mag_bins_subset(N_bins))
  mag_bins = 0
  do i = 1, N_data_mag
    bin_index = int((mag_obs(i) - min_bin) / delta_bin) + 1
    mag_bins(bin_index) = mag_bins(bin_index) + 1   
  end do
  
!-----------------------------------------
! INITIALISATION TIME
time_max = maxval(times)
time_min = minval(times)
! Just for handling some weird cases (1 event in the input catalog)
if (time_max .EQ. time_min ) then 
	if (time_min .NE. 0) then
		time_max = 2 * time_min 
	else 
		time_max = 1
	end if
else if (time_max == 0) then
	time_max = 1
end if

time_lim = (time_max)/100 * 5			! on accepte pas de discontinuite sur les bords std% du catalogue 
t_std = (time_max)/100 * Nc_stdi		! on perturbe la discontinuite de (std%) du catalogue








!------------------------------------------------------------------------
! 1- INITIALIZATION OF THE TRANSDIMENSIONAL INVERSION
!-----------------------------------------------------------------------

  CALL cpu_time(t1)  !Tic. start counting time 
	
  !!1.1 Random draw of the number and location of temporal discontinuities !!
  !-----------------------------------------
  write(*,*) ">> Begin __init__... "
  write(*,*) "------------------- "
  ! Random draw of the number of discontinuities
  Nc = Nc_start_min + (ran3(ra)*(Nc_start_max-Nc_start_min))    
  write(*,*) ". Initial number of discontinuities :", Nc
   
  do k = 1, Nc
  	! Random draw of discontinuity time (not necessarily a catalog times = real between time_min time_max)
	t_init = time_lim + (ran3(ra)*(time_max-(2*time_lim)))	
	
	! a - we find the index associated with this random observation time in the T_bins vector
    	index_ibins = ceiling(T_bins * t_init/(time_max))-1
    			
    	! b - we find the index associated with this random observation time in the data
    	index_data = 0
    	do i = 1, N_data_mag 
    		if (times(i) .le. t_init) then 
    			index_data = i
    		end if 
    	end do	
					
    	model_times(k) = t_init			! store the discontinuities times					
    	model_ibins(k) = index_ibins		! store its index in the T_bins vector
    	model_idata(k) = index_data		! store its index in the data vector
    	model_T(index_ibins) = 1 		! assigns 1 to the location of a discontinuity [0,0,0,1,0,...] in the vector model_T(T_bins)
    	Nc_total(1) = Nc			! store the number of layers accepted Nc at each iteration 
   end do
   
   call sort_shell_real(model_times)	
   call sort_shell(model_ibins)
   call sort_shell(model_idata)
	
   write(*,*) ". Initial model times :", model_times(1:Nc)
   write(*,*) ". Initial model prop :", model_ibins(1:Nc)
   write(*,*) ". Initial model data :", model_idata(1:Nc)


  !! 1.2 Computing conditional P(W|T,D) for each segment !! 
  !------------------------------------------------------

  write(*,*) ". Now, computing first likelihood..."
  write(*,*)".. Likelihood/segments (random) : id | begin | end | n_data | Marginal_b(sample) | MaxLogLike | Log_PT(1)"
  ! a - We find indexes of each temporal segment
  do k = 1, Nc+1
      	if (k==1) then
    		flag_ndata(1) = model_idata(1) 
    		a_data = 1 
      		b_data = a_data + flag_ndata(1) - 1
      		
      		a_ibins = 1
      		b_ibins = model_ibins(k)	
    	
    	else if (k == Nc+1) then
    		flag_ndata(k) = N_data_mag - model_idata(k-1)  
    		a_data = a_data + flag_ndata(k-1)
	    	b_data = a_data + flag_ndata(k) - 1
    		
    		a_ibins = model_ibins(k-1)
      		b_ibins = T_bins
      		
    	else
    		flag_ndata(k) = model_idata(k) - model_idata(k-1) 
    		a_data = a_data + flag_ndata(k-1)
	    	b_data = a_data + flag_ndata(k) - 1
	    	
	      	a_ibins = model_ibins(k-1)
	      	b_ibins = model_ibins(k)
    	end if

    	! b - for the initialization, possibility to compare conditional posterior using gridsearch or random exploration of FMD model parameters
    	write(*,*) 
    	
    	! GRIDSEARCH EXPLORATION
    	if (init_full == 1) then 
    		write(*,*) "... Computing gridsearch conditional for comparison..."
    		Nw = size(Mu_prior)*size(Beta_prior)*size(Sigma_prior)
	    	if (flag_ndata(k)==0) then 								
	    		write(*,*) '.....No events in partition (f)'
	    		
	    	    	Init_Likelihood(k) = 1
	    		Init_Maxloglike(k) = 0	
	    		
	    		do j = a_ibins, b_ibins
	    			Init_Beta_T_posterior(j, :)=1.0/Mu_bins
				Init_Mu_T_posterior(j, :) = 1.0/Mu_bins
				Init_Sigma_T_posterior(j, :) = 1.0/Sigma_bins
			end do	
	    	else 	
			write(*,*) ".....Events in partition :", flag_ndata(k)							
			N_data_subset = size(mag_obs(a_data:b_data))
			
			mag_bins_subset = 0 
			do i = 1, N_data_subset
				bin_index = int((mag_obs(a_data + i-1)-min_bin ) / delta_bin) + 1
				mag_bins_subset(bin_index) = mag_bins_subset(bin_index) + 1	
			end do
			
			call Likelihood_FMD_gridsearch(mag_bins_subset, delta_bin, min_bin, Beta_prior, Mu_prior, Sigma_prior, &
					Likelihood_FMD_3D, Maxloglike)
		    	Init_Likelihood(k) = sum(Likelihood_FMD_3D)
		    	Init_Maxloglike(k) = Maxloglike

		    	Mu_Beta_posterior = sum(Likelihood_FMD_3D, dim=3) /maxval(sum(Likelihood_FMD_3D,dim=3))
			Mu_Sigma_posterior = sum(Likelihood_FMD_3D, dim=2) /maxval(sum(Likelihood_FMD_3D,dim=2))
			Beta_Sigma_posterior = sum(Likelihood_FMD_3D, dim=1) /maxval(sum(Likelihood_FMD_3D,dim=1))
			do j = a_ibins, b_ibins				
				Init_Beta_T_posterior(j,:)=sum(Mu_Beta_posterior, dim=1)/maxval(sum(Mu_Beta_posterior, dim=1))
				Init_Mu_T_posterior(j,:)=sum(Mu_Sigma_posterior, dim=2)/maxval(sum(Mu_Sigma_posterior, dim=2))
				Init_Beta_T_posterior(j,:)=sum(Beta_Sigma_posterior, dim=1)/&
							maxval(sum(Beta_Sigma_posterior, dim=1))
			end do
		end if
		
		if (New_Likelihood(k) == 1) then 
	    		Log_Posterior_T_gridsearch(1) = Log_Posterior_T_gridsearch(1) 
	    	else 
	    		Log_Posterior_T_gridsearch(1) = Log_Posterior_T_gridsearch(1) + log(Init_Likelihood(k))+ Init_Maxloglike(k)- log(Nw)
	    	end if
		write(*,*) ".. Likelihood/segments (gridsearch)",k,a_data,b_data,flag_ndata(k),Init_Likelihood(k),&
							Init_Maxloglike(k), Log_Posterior_T_gridsearch(1)
	
	end if  
	
	! RANDOM EXPLORATION
	Nw = n_draws 
	if (flag_ndata(k)==0) then 
		
	    	write(*,*) '.....No events in partition (r)'
    	    	Init_Likelihood(k) = 1
    		Init_Maxloglike(k) = 0	
    		do j = a_ibins, b_ibins
			Init_Beta_T_posterior(j,  :) = 1.0/Beta_bins
			Init_Mu_T_posterior(j, :) = 1.0/Mu_bins
			Init_Sigma_T_posterior(j, :) = 1.0/Sigma_bins
		end do	
	else
		N_data_subset = flag_ndata(k)
		mag_bins_subset = 0 
		do i = 1, N_data_subset
			bin_index = int((mag_obs(a_data + i-1)-min_bin) / delta_bin) + 1
			mag_bins_subset(bin_index) = mag_bins_subset(bin_index) + 1	
		end do

		call Likelihood_FMD_random(mag_bins_subset, delta_bin, min_bin, n_draws, Beta_prior , Mu_prior, Sigma_prior, &
					Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
		Init_Likelihood(k) = sum(Likelihood_FMD_MC_3D)
	    	Init_Maxloglike(k) = Maxloglike
		do j = a_ibins, b_ibins				
			Init_Beta_T_posterior(j, :) =  Beta_posterior
			Init_Mu_T_posterior(j, :) =  Mu_posterior
			Init_Sigma_T_posterior(j, :) =  Sigma_posterior
		end do
	end if
	
	if (Init_Likelihood(k)  == 1) then 
		
    		Acc_Log_Posterior_T(1) = Acc_Log_Posterior_T(1)
    	else
    		
    		Acc_Log_Posterior_T(1) = Acc_Log_Posterior_T(1) + log(Init_Likelihood(k))+ Init_Maxloglike(k)- log(Nw)
    	end if
	write(*,*) ".. Likelihood/segments (random)", k, a_ibins, b_ibins, flag_ndata(k),Init_Beta_T_posterior(a_ibins,2), &
							Init_Likelihood(k),Init_Maxloglike(k),Acc_Log_Posterior_T(1)

  end do

  !! 1.3 First draw of the type of perturbation of the temporal model 
  !--------------------------------------------------------------------
  write(*,*) ". First random draw Move/Birth/Death :..."
  write(*,*)
  write(*,*)
  flag_unaccepted = 0 
  index_ibins = 0
  New_Nc = Nc
  New_model_T = model_T
  New_model_times =  model_times
  New_model_ibins =  model_ibins
  New_model_idata =  model_idata
  New_model_idata_bis =  model_idata
  New_model_ibins_bis =  model_ibins

  alpha = ran3(ra) 
  if (alpha.lt.1./3) then ! Move
  	write(*,*) ".. Case Move : selecting discontinuity... "
  	N_move = N_move + 1
  	New_Nc = Nc 				! the number of discontinuities does not change
  	cpert = ceiling(ran3(ra)*(Nc))		! we randomly draw one of the Nc discontinuities
  	t_old = model_times(cpert)		! we get the previous time of the discontinuity
  	t_new = t_old + gasdev(ra) * t_std	! we perturb according to a gaussian (gasdev)
	
  	if ((t_new.ge.time_lim) .and. (t_new.le.(time_max-time_lim))) then ! check if within bounds
  		New_model_times(cpert) = t_new
  		write(*,*) "value of t_new",  t_new
  		! we find the index associated with this random observation time in the T_bins vector		
  		index_old = ceiling(T_bins * t_old/time_max)-1	
  		index_ibins = ceiling(T_bins * t_new/time_max)-1	
  		
  		! update proposal with new temporal model 
  		New_model_T(index_old) = 0
  		New_model_T(index_ibins) = 1
  		New_model_ibins(cpert) = index_ibins
  		New_model_idata_bis(cpert) = 9999999
  		New_model_ibins_bis(cpert) = 9999999
  		! we find the index associated with this random observation time in the data
  		index_data = 0
    		do i = 1, N_data_mag 
    			if (times(i) .le. t_new) then 
    				index_data = i
    			end if 
    		end do
    		New_model_idata(cpert) = index_data
  	else 
  		write(*,*) ".. Move rejected : bounds condition"
  		flag_unaccepted = 0
  	end if 
  	
  	call sort_shell_real(New_model_times)
  	call sort_shell(New_model_ibins)
  	call sort_shell(New_model_idata)
  	

	write(*,*) ".. Proposed model Perturbation ", New_model_ibins(1:New_Nc)
			
  else if ((alpha.ge.1./3) .and. (alpha.lt.2./3))  then   ! Death 
  	write(*,*) ".. Case Death : selecting discontinuity..."
  	N_death = N_death +1
	if (Nc/=0) then
		New_Nc = Nc - 1 		! the number of discontinuities is reduced
  		cpert = ceiling(ran3(ra)*(Nc))	! we randomly draw one of the Nc discontinuities		
  		t_new = model_times(cpert)	! we get the previous time of the discontinuity
  		write(*,*) "New discontinuity :", t_new, cpert
  		! update proposal with new temporal model
  		New_model_times(cpert) = 9999999
  		New_model_idata(cpert) = 9999999
  		New_model_ibins(cpert) = 9999999
  		! we find the index associated with this random observation time in the T_bins vector
  		index_ibins = ceiling(T_bins * t_new/time_max)-1   
  		New_model_T(index_ibins) = 0	! we remove the discontinuity in model_T
  		call sort_shell_real(New_model_times)
  		call sort_shell(New_model_ibins)
  		call sort_shell(New_model_idata)
  	else 
  		New_Nc = 0
  		write(*,*) ".. Death rejected : already no layers "
  		flag_unaccepted = 0
  		
  	end if   		
  	write(*,*) ".. Proposed model Death (new Timegrid)", New_model_ibins(1:New_Nc)

  	
  else ! Birth
  	write(*,*) ".. Case Birth : selecting a new discontinuity... "
  	N_birth = N_birth+1	! the number of discontinuities is reduced
	cpert = Nc+1 		! we first put it at the end of all discontinuities
	if (cpert .le. max_layers) then ! check if the number of discontinuities is inf to max_layers
  		t_new = time_lim + (ran3(ra)*(time_max-(2*time_lim)))	! we randomly draw one of the Nc discontinuities
  		
  		! we find the index associated with this random observation time in the T_bins vector				
		index_ibins = ceiling(T_bins * t_new/time_max)-1		
		! update proposal with new temporal model
  		New_model_T(index_ibins) = 1
  		New_model_ibins(cpert) = index_ibins
  		New_model_times(cpert) = t_new
  		
  		index_data = 0
    		do i = 1, N_data_mag 
    			if (times(i) .le. t_new) then 
    				index_data = i
    			end if 
    		end do	
    		New_model_idata(cpert) = index_data
    		
  		call sort_shell_real(New_model_times)
  		call sort_shell(New_model_ibins)
  		call sort_shell(New_model_idata)
		New_Nc = Nc+1
  		do k = 1, Nc+1					! we find the index of the perturbed discontinuity
     			if (New_model_times(k) == t_new) then
        			cpert = k
        				
     			end if
     		end do
  	else 
  		write(*,*) ".. No Birth : max_layers "
  		flag_unaccepted = 0
  	end if 
  	write(*,*) ".. Proposed model Birth (New TimeGrid)", New_model_ibins(1:New_Nc)
  end if
  
  !! 1.4  Computing once the case with no temporal discontinuities Nc=0
  !--------------------------------------------------------------------
  write(*,*)">> Computing once the total Likelihood (case of no discontinuities)..."
  write(*,*)
  call Likelihood_FMD_gridsearch(mag_bins, delta_bin, min_bin, Beta_prior, Mu_prior, Sigma_prior, &
					Likelihood_FMD_3D, Maxloglike)
  Likelihood_all = sum(Likelihood_FMD_3D)
  Maxloglike_all = Maxloglike	
  write(*,*) ".. Total Likelihood (Nc=0) full gridsearch (long) :",  Likelihood_all,  Maxloglike_all
  write(*,*)
  
  call Likelihood_FMD_random(mag_bins, delta_bin, min_bin, n_draws, Beta_prior , Mu_prior, Sigma_prior, &
					Likelihood_FMD_MC_3D, Maxloglike, Beta_posterior, Mu_posterior, Sigma_posterior)
  Likelihood_all = sum(Likelihood_FMD_MC_3D)
  Maxloglike_all = Maxloglike	
  
  Log_Posterior_T_all = dlog(Likelihood_all)  + Maxloglike_all - log(Nw)
  write(*,*) ".. Total Likelihood (Nc=0) (random exploration) :  ", Likelihood_all, Maxloglike_all
  write(*,*) ".. Condition MH no discontinuities Nc=0 :  ", Log_Posterior_T_all
  do j = 1, T_bins				
	Beta_T_posterior_all(j, :) =  Beta_posterior
	Mu_T_posterior_all(j,  :) = Mu_posterior
	Sigma_T_posterior_all(j, :) =  Sigma_posterior
  end do
  write(*,*) ".. If Likelihood full gridsearch different of likelihood random exploration !! Change exploration parameters!!"
  write(*,*) "<< End of Initialization >>"
  write(*,*)
  write(*,*)
   
   
   
















   
!----------------------------------------------------------------------------------------------------
! 2- BEGIN TRANSDIMENSIONAL ITERATIONS
!------------------------------------------------------------------------------------------------------
  write(*,*) ">> Begin reversible-jump..."
  write(*,*) "--------------------------"
  write(*,*) ". Total number of iterations : ", it_max
  write(*,*) ". Thinning : ", it_thin
  write(*,*)
  write(*,*) ">> First iteration : 1"
  write(*,*) ". Computing the likelihood of the new model (for comparison)..."
  
  
  ! Set current model = initial models
  Acc_Likelihood = Init_Likelihood
  Acc_Maxloglike = Init_Maxloglike
  Acc_Beta_T_posterior = Init_Beta_T_posterior
  Acc_Mu_T_posterior = Init_Mu_T_posterior
  Acc_Sigma_T_posterior = Init_Sigma_T_posterior


  Nw = n_draws
  New_std = Nc_stdi

  
  do it = 2, it_max
  	
  	CALL cpu_time(t4)
  	if (mod(it ,it_std)==0) then ! talk
  		  
 		  AR_move = real(N_move_accept)/real(N_move)
		  AR_death = real(N_death_accept)/real(N_death)
		  AR_birth = real(N_birth_accept)/real(N_birth)
		  
		  write(*,*) ". Old perturbation size",New_std
  		  write(*,*) ". Current Acceptance Rates Move/Death/Birth :", it, AR_move, AR_death, AR_birth
  		  New_std = max(0.01, New_std *  AR_move/0.3)	!STD in purcentage
  		  write(*,*) ". New perturbation size", New_std
  		  t_std = (time_max-time_lim)/100  * New_std
  	end if
 
  	!! 2.1 On calcule likelihood methode random 
  	!-------------------------------------------------------------------------
  	New_Likelihood = 0.0
  	New_Maxloglike = 0.0
  	New_Log_Posterior_T = 0.0
  	
  	if (New_Nc==0) then
  		New_Likelihood(1) = Likelihood_all
		New_Maxloglike(1) = Maxloglike_all
		New_Log_Posterior_T = Log_Posterior_T_all
  	else ! On ne calcule la likelihood que la ou elle a change
		if (alpha.lt.1./3) then								!! MOVE CASE
			 ! For a Move Case we consider Move= Death + Birth
			 New_Nc = New_Nc-1  ! Death
			 if (New_Nc==0) then
		  		New_Likelihood(1) = Likelihood_all
				New_Maxloglike(1) = Maxloglike_all

			 else				
				call Likelihood_FMD_random_DEATH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
								Beta_prior, Mu_prior, Sigma_prior, &
								cpert, New_Nc, T_bins, &
                                        			New_model_idata_bis, New_model_ibins_bis, & 
                                        			Acc_Likelihood, Acc_Maxloglike, &
                                       				Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        			Likelihood_bis, Maxloglike_bis, &
                                        			Beta_T_posterior_bis,  Mu_T_posterior_bis, Sigma_T_posterior_bis)                            

  			end if
			
			New_Nc = New_Nc + 1   ! Birth
			 call Likelihood_FMD_random_BIRTH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
							Beta_prior, Mu_prior, Sigma_prior,&
							cpert, New_Nc, T_bins, &
                                        		New_model_idata, New_model_ibins, & 
                                        		Likelihood_bis, Maxloglike_bis, &
                                        		Beta_T_posterior_bis,  Mu_T_posterior_bis, Sigma_T_posterior_bis, &
                                        		New_Likelihood, New_Maxloglike, &
                                        		New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)

			
		else if ((alpha.ge.1./3) .and. (alpha.lt.2./3)) then				!! DEATH CASE

							
			call Likelihood_FMD_random_DEATH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
							Beta_prior, Mu_prior, Sigma_prior, &
							cpert, New_Nc, T_bins, &
                                        		New_model_idata, New_model_ibins, & 
                                        		Acc_Likelihood, Acc_Maxloglike, &
                                       			Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        		New_Likelihood, New_Maxloglike, &
                                        		New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)      

			
		else										!! BIRTH CASE
			 call Likelihood_FMD_random_BIRTH(mag_obs, mag_bins, delta_bin, min_bin, n_draws, & 
							Beta_prior, Mu_prior, Sigma_prior,&
							cpert, New_Nc, T_bins, &
                                        		New_model_idata, New_model_ibins, & 
                                        		Acc_Likelihood, Acc_Maxloglike, &
                                       			Acc_Beta_T_posterior,  Acc_Mu_T_posterior, Acc_Sigma_T_posterior, &
                                        		New_Likelihood, New_Maxloglike, &
                                        		New_Beta_T_posterior,  New_Mu_T_posterior, New_Sigma_T_posterior)
		end if 
		
		
		

		do k = 1, (New_Nc +1)
			if (New_Likelihood(k)==1) then 
				New_Log_Posterior_T = New_Log_Posterior_T
			else 
				New_Log_Posterior_T = New_Log_Posterior_T + dlog(New_Likelihood(k)) + New_Maxloglike(k) - log(Nw) 
			end if
		end do 
    	end if
    	





	
	
  	!! 2.2 TEST CONDITION MCMC - Metropolis Hastings
  	!-------------------------------------------------------------------------
        write(*,*) "-------------------------------------------------------------------------"
	write(*,*) ". Metropolis-Hastings test : new / old max loglikelihood"
	write(*,*) ". Metropolis-Hastings test : ", New_Log_Posterior_T, Acc_Log_Posterior_T(it-1)
    	u = log(ran3(ra))
    	condition = New_Log_Posterior_T - Acc_Log_Posterior_T(it-1)
    	
    	if ((u.gt.condition) .or. (flag_unaccepted ==1) ) then   	! 2.2.A TEST CONDITION MH : REJECTION
    	  								!--------------------------------------
    		write(*,*) ".. Model discarded (old>new) ! "  
    		Acc_Log_Posterior_T(it) = Acc_Log_Posterior_T(it-1) ! on garde le modele precedent
    		
    		flag_unaccepted =0
    		Nc_total(it) = Nc
    		do j = 1, T_bins	
    			Sum_accepted_T(j) = Sum_accepted_T(j) + model_T(j)
   		end do
	
	
	
    	else 								! 2.2.B TEST CONDITION MH : ACCEPT
    									!-------------------------------------
		Nc = New_Nc
		Nc_total(it) = New_Nc
		Acc_Log_Posterior_T(it) = New_Log_Posterior_T
		model_T = New_model_T		
    		model_times = New_model_times
    		model_idata = New_model_idata
    		model_ibins = New_model_ibins
		
		Acc_Likelihood = New_Likelihood
		Acc_Maxloglike = New_Maxloglike
		
		Acc_Beta_T_posterior = New_Beta_T_posterior
		Acc_Mu_T_posterior = New_Mu_T_posterior
		Acc_Sigma_T_posterior = New_Sigma_T_posterior

		if (it .gt. it_burnin) then					
			if (N_thin .lt. it_thin) then 
				N_thin = N_thin+1


			else 
				N_thin = 0
				! Sum all accepted marginals p(B |T,d), p(Mu|T,d), p(Sigma|T,d) and p(T|d)
	    			Sum_Beta_T_posterior = Sum_Beta_T_posterior + Acc_Beta_T_posterior
	    			Sum_Mu_T_posterior = Sum_Mu_T_posterior + Acc_Mu_T_posterior
	    			Sum_Sigma_T_posterior = Sum_Sigma_T_posterior + Acc_Sigma_T_posterior
				do j = 1, T_bins	
	    				Sum_accepted_T(j) = Sum_accepted_T(j) + New_model_T(j)
	   			end do

			end if 
    		end if
    		
    		if (alpha.lt.1./3) then
    			N_move_accept = N_move_accept +1
    		else if ((alpha.ge.1./3) .and. (alpha.lt.2./3)) then
    			N_death_accept = N_death_accept + 1
    		else
    			N_birth_accept = N_birth_accept + 1 
    		end if 
    		
    		
    		
    		write(*,*) ".. Model accepted (old<new) : ", model_times(1:Nc)
    		
    	end if  
    	write(*,*) "-------------------------------------------------------------------------"
    	
    	
    	
    	
    	
    	!! 2.3 ON TIRE UN NOUVEAU MODELE DE DISCONTINUITES
    	!-------------------------------------------------------------------------
	CALL cpu_time(t3)  !Toc. end counting time
	progress = (t3-t1)/3600
	t_it(it) = (t3-t4)
	mean_t_it = sum(t_it(1:it))/it * it_max /3600
	write(*,*) "Time (hours) : ", progress ,'/', mean_t_it
	write(*,*) "<< End of iteration >>"
  	write(*,*)
    	write(*,*)
	if (it.EQ.it_max) then
		write(*,*) ">> Finalisation :", it
	else
    		write(*,*) ">> New iteration : ", it 
	endif
	
	
	
	!---------------------------------------------------------------------------
	  flag_unaccepted = 0 
	  index_ibins = 0
	  New_Nc = Nc
	  New_model_T = model_T
	  New_model_times =  model_times
	  New_model_ibins =  model_ibins
	  New_model_idata =  model_idata
	  New_model_idata_bis =  model_idata
	  New_model_ibins_bis =  model_ibins

	  alpha = ran3(ra) 
	  if (alpha.lt.1./3) then ! Move
	  	write(*,*) ".. Case Move : selecting discontinuity... "
	  	N_move = N_move + 1
	  	if (Nc/=0) then 
		  	New_Nc = Nc 				! the number of discontinuities does not change
		  	cpert = ceiling(ran3(ra)*(Nc))		! we randomly draw one of the Nc discontinuities
		  	t_old = model_times(cpert)		! we get the previous time of the discontinuity
		  	t_new = t_old + gasdev(ra) * t_std	! we perturb according to a gaussian (gasdev)
			
		  	if ((t_new.ge.time_lim) .and. (t_new.le.(time_max-time_lim))) then ! check if within bounds
		  		New_model_times(cpert) = t_new
		  		! we find the index associated with this random observation time in the T_bins vector		
		  		index_old = ceiling(T_bins * t_old/time_max)-1	
		  		index_ibins = ceiling(T_bins * t_new/time_max)-1
		  		
  				New_model_idata_bis(cpert) = 9999999
  				New_model_ibins_bis(cpert) = 9999999
		  		
		  		! update proposal with new temporal model 
		  		New_model_T(index_old) = 0
		  		New_model_T(index_ibins) = 1
		  		New_model_ibins(cpert) = index_ibins
		  		! we find the index associated with this random observation time in the data
		  		index_data = 0
		    		do i = 1, N_data_mag 
		    			if (times(i) .le. t_new) then 
		    				index_data = i
		    			end if 
		    		end do
		    		New_model_idata(cpert) = index_data
		  	else 
		  		write(*,*) ".. Move rejected : bounds condition"
		  		flag_unaccepted = 0
		  	end if 
		  else
		  	flag_unaccepted = 1
	  		New_Nc=0
	  	end if
	  	call sort_shell_real(New_model_times)
		call sort_shell(New_model_ibins)
		call sort_shell(New_model_idata)
		call sort_shell(New_model_ibins_bis)
	  	call sort_shell(New_model_idata_bis)
		write(*,*) ".. Proposed model Move (new Timegrid)", New_model_ibins(1:New_Nc)  
				
	  else if ((alpha.ge.1./3) .and. (alpha.lt.2./3))  then   ! Death 
	  	write(*,*) ".. Case Death : selecting discontinuity..."
	  	N_death = N_death +1
		if (Nc/=0) then
			New_Nc = Nc - 1 		! the number of discontinuities is reduced
	  		cpert = ceiling(ran3(ra)*(Nc))	! we randomly draw one of the Nc discontinuities		
	  		t_new = model_times(cpert)	! we get the previous time of the discontinuity
	  		! update proposal with new temporal model
	  		New_model_times(cpert) = 9999999
	  		New_model_idata(cpert) = 9999999
	  		New_model_ibins(cpert) = 9999999
	  		! we find the index associated with this random observation time in the T_bins vector
	  		index_ibins = ceiling(T_bins * t_new/time_max)-1   
	  		New_model_T(index_ibins) = 0	! we remove the discontinuity in model_T
	  		call sort_shell_real(New_model_times)
	  		call sort_shell(New_model_ibins)
	  		call sort_shell(New_model_idata)
	  	else 
	  		New_Nc = 0
	  		write(*,*) ".. Death rejected : already no layers "
	  		flag_unaccepted = 0
	  		
	  	end if   		
	  	write(*,*) ".. Proposed model Death (new Timegrid)", New_model_ibins(1:New_Nc)

	  	
	  else ! Birth
	  	write(*,*) ".. Case Birth : selecting a new discontinuity... "
	  	N_birth = N_birth+1	! the number of discontinuities is reduced
		cpert = Nc+1 		! we first put it at the end of all discontinuities
		if (cpert .le. max_layers) then ! check if the number of discontinuities is inf to max_layers
	  		t_new = time_lim + (ran3(ra)*(time_max-(2*time_lim)))	! we randomly draw one of the Nc discontinuities
	  		
	  		! we find the index associated with this random observation time in the T_bins vector				
			index_ibins = ceiling(T_bins * t_new/time_max)-1		
			! update proposal with new temporal model
	  		New_model_T(index_ibins) = 1
	  		New_model_ibins(cpert) = index_ibins
	  		New_model_times(cpert) = t_new
	  		
	  		index_data = 0
	    		do i = 1, N_data_mag 
	    			if (times(i) .le. t_new) then ! on retrouve l'index associe a ce temps d'observation random dans les donnees
	    				index_data = i
	    			end if 
	    		end do	
	    		New_model_idata(cpert) = index_data
	    		
	  		call sort_shell_real(New_model_times)
	  		call sort_shell(New_model_ibins)
	  		call sort_shell(New_model_idata)
			New_Nc = Nc+1
	  		do k = 1, Nc+1					! we find the index of the perturbed discontinuity
     				if (New_model_times(k) == t_new) then
        				cpert = k
        				
     				end if
     			end do
	  	else 
	  		write(*,*) ".. No Birth : max_layers "
	  		flag_unaccepted = 0
	  	end if 
	  	write(*,*) ".. Proposed model Birth (New TimeGrid)", New_model_ibins(1:New_Nc)
	  end if
	write(*,*) ".. New number of discontinuities :", New_Nc	
	write(*,*) ".. --------------------------"
	write(*,*) ".. Recomputing likelihood..."  
  end do  
!----------------------------------------------------------------------------------------------------
! 2- FIN TRANSDIMENSIONAL ITERATIONS
!------------------------------------------------------------------------------------------------------ 
  CALL cpu_time(t2)  !Toc. end counting time
  comp_time = nint((t2-t1)*100.0)
  comp_time = comp_time/6000
  write(*,*) "<< End of reversible-jump >>"
  write(*,*) " -------------------------"
  write(*,*) ">> End of inversion (resume)" 
  write(*,*) ". Total time of inversion in minuts : " , comp_time
  write(*,*)
  
  ! Petit bilan des taux d'acceptance 
  AR_move = real(N_move_accept)/real(N_move)
  AR_death = real(N_death_accept)/real(N_death)
  AR_birth = real(N_birth_accept)/real(N_birth)
  write(*,*) "Final acceptance rates : "
  write(*,*) "Name : number_accepted  number_proposed    rate(pct)"
  write(*,*) "Moves  : ", N_move_accept, N_move, AR_move
  write(*,*) "Deaths : ", N_death_accept, N_death, AR_death
  write(*,*) "Births : ", N_birth_accept, N_birth, AR_birth

  
  
!----------------------------------------------------------------------------------------------------
! #3- END OF PARALLELISATION OPENMPI
!------------------------------------------------------------------------------------------------------
  
  
!----------------------------------------------------------------------------------------------------
! #4- WRITING RESULTS - 
!------------------------------------------------------------------------------------------------------
  
  write(*,*) "..Writting outputs files... "
  
  open(unit=62, file="../outputs/Acceptances_Rates.txt",  status="replace", action="write")
  write(62,*) N_birth_accept,N_death_accept, N_move_accept  
  write(62,*) N_birth, N_death, N_move
  write(62,*) AR_birth,  AR_death, AR_move
  close(62)
  
  
  open(unit=63, file="../outputs/Number_of_discontinuities.txt",  status="replace", action="write")
  do line = 1, it_max
  	write(63, *) Nc_total(line)
  end do
  close(63)
  
  open(unit=64, file="../outputs/Probability_discontinuities.txt",  status="replace", action="write")
  do line = 1, T_bins
  	write(64, *) Sum_accepted_T(line)
  end do
  close(64)
  

   open(unit=60, file="../outputs/prior.txt",  status="replace", action="write")
	write(60,*) Mu_min, Mu_max, Mu_bins 
	write(60,*) Beta_min, Beta_max, Beta_bins
	write(60,*) Sigma_min, Sigma_max, Sigma_bins
	write(60,*) it_max, n_draws, T_bins, comp_time
 close(60)

  open(unit=56, file="../outputs/Marginal_b.txt",  status="replace", action="write")
  do line = 1, Beta_bins
  	do column = 1, T_bins
  		write(56, '(F25.5,2X)', advance='no') Sum_Beta_T_posterior(column,line)
  	end do
  	write(56,*)
  end do
  close(56)
  
  open(unit=55, file="../outputs/Marginal_Mu.txt",  status="replace", action="write")
  do line = 1, Mu_bins
  	do column = 1, T_bins
  		write(55,'(F25.5,2X)', advance='no') Sum_Mu_T_posterior(column,line)
  	end do
  	write(55,*)
  end do
  close(55)
  
  

  
  open(unit=59, file="../outputs/Marginal_Sigma.txt",  status="replace", action="write")
  do line = 1, Sigma_bins
  	do column = 1, T_bins
  		write(59, '(F25.5,2X)', advance='no') Sum_Sigma_T_posterior(column,line)
  	end do
  	write(59,*)
  end do
  close(59)
 
  write(*,*) ".. Completed"
  write(*,*) "< End of b-Bayesian program >"
  !------------------------------------------------------------------------------------------------------
end program b_Bayesian_binned










!-----------------------------------------
! UTILS
!-----------------------------------------


 
 FUNCTION ran3(idum)
INTEGER idum
INTEGER MBIG,MSEED,MZin
!     REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=800000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
! write(*,*)' idum ',idum
if(idum.lt.0.or.iff.eq.0)then
	iff=1
	mj=MSEED-iabs(idum)
	mj=mod(mj,MBIG)
	ma(55)=mj
	mk=1
	do 11 i=1,54
	ii=mod(21*i,55)
	ma(ii)=mk
	mk=mj-mk
	if(mk.lt.MZ)mk=mk+MBIG
	mj=ma(ii)
!  write(*,*)' idum av',idum
11      continue
	do 13 k=1,4
	do 12 i=1,55
	ma(i)=ma(i)-ma(1+mod(i+30,55))
	if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
! write(*,*)' idum ap',idum
	inext=0
	inextp=31
	idum=1
endif
! write(*,*)' idum app ',idum
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
!  write(*,*)' idum ',idum
	
return
END


FUNCTION GASDEV(idum)

!     ..Arguments..
integer          idum
real GASDEV

!     ..Local..
real v1,v2,r,fac
real ran3

if (idum.lt.0) iset=0
10 v1=2*ran3(idum)-1
v2=2*ran3(idum)-1
r=v1**2+v2**2
if(r.ge.1.or.r.eq.0) GOTO 10
fac=sqrt(-2*log(r)/r)
GASDEV=v2*fac

RETURN
END












































		
	
	

