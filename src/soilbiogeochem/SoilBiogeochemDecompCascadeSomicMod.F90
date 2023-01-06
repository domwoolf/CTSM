module SoilBiogeochemDecompCascadeSOMicMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Sets the coefficients used in the decomposition submodel.
  ! This uses the SOMic model
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_const_mod                      , only : SHR_CONST_TKFRZ
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools_max
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_met_lit, i_cwd, i_cwdl2
  use clm_varctl                         , only : iulog, spinup_state, anoxia, use_lch4, use_fates
  use clm_varcon                         , only : zsoi
  use decompMod                          , only : bounds_type
  use spmdMod                            , only : masterproc
  use abortutils                         , only : endrun
  use CNSharedParamsMod                  , only : CNParamsShareInst, nlev_soildecomp_standard
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, InitSoilTransfer, use_soil_matrixcn
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilStateType                      , only : soilstate_type
  use TemperatureType                    , only : temperature_type
  use ch4Mod                             , only : ch4_type
  use ColumnType                         , only : col
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term

  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams                      ! Read in parameters from params file
  public :: init_decompcascade_somic        ! Initialization
  public :: decomp_rate_constants_somic     ! Figure out decomposition rates
  !
  ! !PUBLIC DATA MEMBERS
  logical , public :: normalize_q10_to_century_tfunc = .true.! do we normalize the century decomp. rates so that they match the CLM Q10 at a given temp?
  logical , public :: use_century_tfunc = .false.            ! do we use the daycent temperature scalar?
  real(r8), public :: normalization_tref = 15._r8            ! reference temperature for normalizaion (degrees C)
  !

  ! !PRIVATE DATA MEMBERS

  ! soil pool indices
  integer, private :: i_mac_som  ! index of mineral-associated Soil Organic Matter (SOM)
  integer, private :: i_mic_som  ! index of active microbial biomass SOM
  integer, private :: i_doc_som  ! index of dissolved SOM
  integer, private :: i_cel_lit  ! index of cellulose litter pool
  integer, private :: i_lig_lit  ! index of lignin litter pool

  ! cwd variables
  real(r8), private :: cwd_fcel  ! cellulose fraction in coarse woody debris
  real(r8), private :: rf_cwdl2
  real(r8), private :: rf_cwdl3

  ! respiration fractions by transition - these will all be zero except s1s2 (doc uptake by microbes), and cwd decomposition.  The zero rf's are placeholders, that may be removed later.
  ! real(r8), private :: rf_l1s1
  ! real(r8), private :: rf_l2s1
  ! real(r8), private :: rf_l3s1
  ! real(r8), private :: rf_s1s2
  ! real(r8), private :: rf_s1s3
  ! real(r8), private :: rf_s2s1
  ! real(r8), private :: rf_s3s1

  ! indices of transitions
  integer, private :: i_l1s1
  integer, private :: i_l2s1
  integer, private :: i_l3s1
  integer, private :: i_s1s2
  integer, private :: i_s1s3
  integer, private :: i_s2s1
  integer, private :: i_s3s1
  integer, private :: i_cwdl2
  integer, private :: i_cwdl3

  ! define type to pass input parameters for SOMIc model (then define variable "params_inst" of this type)
  type, private :: params_type
     real(r8) :: cn_dom        ! C:N for dissolved organic matter
     real(r8) :: cn_mic        ! C:N for microbial biomass
     real(r8) :: cn_mac        ! C:N for mineral-associated organic matter
     real(r8) :: cue_0         ! default carbon use efficiency for soil microbes
     real(r8) :: mcue          ! temperature-dependence slope of micobial CUE
     real(r8) :: mic_vmax      ! Michaelis-Menten maximum rate coefficient
     real(r8) :: mic_km        ! Michaelis-Menten half saturation coefficient
     real(r8) :: k_l1s1        ! rate constant for litter1 dissolution to DOC (/s)
     real(r8) :: k_l2s1        ! rate constant for litter2 depolymerization to DOC (/s)
     real(r8) :: k_l3s1        ! rate constant for litter2 depolymerization to DOC (/s)
     real(r8) :: k_s1s2        ! rate constant for microbial uptake of DOC
     real(r8) :: k_s1s3        ! rate constant for sorption of DOC to minerals
     real(r8) :: k_s2s1        ! rate constant for microbial turnover (both death and exudates)
     real(r8) :: k_s3s1        ! rate constant for desorption of MAC to DOC
     real(r8) :: mclay         ! slope of texture dependence scalar as linear function of clay content
     real(r8) :: clay0         ! clay content at which texture-dependence scalar equals 1.0
     real(r8) :: cwd_fcel      ! cellulose fraction for CWD
     real(r8) :: rf_cwdl2      ! respiration fraction for CWD to L2
     real(r8) :: rf_cwdl3      ! respiration fraction for CWD to L3
     real(r8), allocatable :: bgc_initial_Cstocks(:)  ! Initial Carbon stocks for a cold-start
     real(r8) :: bgc_initial_Cstocks_depth  ! Soil depth for initial Carbon stocks for a cold-start
  end type params_type
  !
  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio , only: file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompSOMicParamsType'
    character(len=100) :: errCode = 'Error reading in SOMic const file '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! Read SOMic parameters from netcdf file
    tString = 'somic_cwd_fcel'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cwd_fcel = tempr

    tString='somic_rf_cwdl2'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl2 = tempr

    tString='somic_rf_cwdl3'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl3 = tempr

    tString = 'somic_cn_dom'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_dom = tempr

    tString = 'somic_cn_mic'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_mic = tempr

    tString = 'somic_cn_mac'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_mac = tempr

    tString = 'somic_cue_0'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cue_0 = tempr

    tString = 'somic_mcue'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mcue = tempr

    tString = 'somic_mic_vmax'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mic_vmax = tempr

    tString = 'somic_mic_km'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mic_km = tempr

    tString = 'somic_k_l1s1'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_l1s1 = tempr

    tString = 'somic_k_l2s1'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_l2s1 = tempr

    tString = 'somic_k_l3s1'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_l3s1 = tempr

    tString = 'somic_k_s1s2'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_s1s2 = tempr

    tString = 'somic_k_s1s3'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_s1s3 = tempr

    tString = 'somic_k_s2s1'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_s2s1 = tempr

    tString = 'somic_k_s3s1'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_s3s1 = tempr

    tString = 'somic_mclay'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mclay = tempr

    tString = 'somic_clay0'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%clay0 = tempr

    allocate(params_inst%bgc_initial_Cstocks(ndecomp_pools_max))
    tString = 'somic_initial_Cstocks'
    call ncd_io(trim(tString), params_inst%bgc_initial_Cstocks(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    tString = 'somic_initial_Cstocks_depth'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%bgc_initial_Cstocks_depth = tempr

    write(iulog,*) 'Reading SOMIC params:'
    write(iulog,*) 'SOMIC params (somic_cwd_fcel, somic_rf_cwdl2, somic_rf_cwdl3)', &
                    params_inst%cwd_fcel, params_inst%rf_cwdl2, params_inst%rf_cwdl3
    write(iulog,*) 'SOMIC params (somic_cn_dom, somic_cn_mic, somic_cn_mac)', &
                    params_inst%cn_dom, params_inst%cn_mic, params_inst%cn_mac
    write(iulog,*) 'SOMIC params (somic_cue_0, somic_mcue, somic_mic_vmax, somic_mic_km)', &
                    params_inst%cue_0, params_inst%mcue, params_inst%mic_vmax, params_inst%mic_km
    write(iulog,*) 'SOMIC params (somic_k_l1s1, somic_k_l2s1, somic_k_l3s1', &
                    params_inst%k_l1s1, params_inst%k_l2s1, params_inst%k_l3s1
    write(iulog,*) 'SOMIC params (somic_k_s1s2, somic_k_s1s3, somic_k_s2s1, somic_k_s3s1', &
                    params_inst%k_s1s2, params_inst%k_s1s3, params_inst%k_s2s1, params_inst%k_s3s1
    write(iulog,*) 'SOMIC params (somic_mclay, somic_clay0, somic_initial_Cstocks_depth)', &
                    params_inst%mclay, params_inst%clay0, params_inst%bgc_initial_Cstocks_depth
    write(iulog,*) 'SOMIC params (somic_initial_Cstocks)', params_inst%bgc_initial_Cstocks(:)

  end subroutine readParams



  !-----------------------------------------------------------------------
  subroutine init_decompcascade_somic(bounds, soilbiogeochem_state_inst) ! , soilstate_inst) !, soilbiogeochem_carbonstate_inst) ! TODO parameters were deleted, uncomment them once sure
    !
    ! !DESCRIPTION:
    !  Initialize rate constants and decomposition pathways following the decomposition cascade of the SOMic model.
    !  In this subroutine, we only set the time-independent coefficient values
    !  Written by D. Woolf
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds
    type(soilbiogeochem_state_type)      , intent(inout) :: soilbiogeochem_state_inst
    ! type(soilstate_type)                 , intent(in)    :: soilstate_inst
    ! type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_carbonstate_inst

    !
    ! !LOCAL VARIABLES
    integer  :: c, j                         ! indices
    real(r8) :: t                            ! temporary variable
    real(r8) :: cn_s1                        ! C/N ratio of dissolved OM
    real(r8) :: cn_s2                        ! C/N ratio of microbial biomass
    real(r8) :: cn_s3                        ! C/N ratio of mineral-associated OM
    real(r8) :: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.
    ! real(r8) :: clayfact                     ! rate modyfying coefficient due to soil clay content
    ! real(r8) :: ksorb_altered                ! rate constant for sorption, once rate modyfying sclars have been applied
    ! real(r8) :: kmicrobial_uptake_altered    ! rate constant for microbial uptake, once rate modyfying sclars have been applied
    ! real(r8) :: fsorb                        ! fraction of doc removals sorbed to mineral surfaces
    ! real(r8) :: fmic                         ! fraction of doc removals taken up by microbes
    ! real(r8) :: kdoc                         ! rate constant for removal of doc into microbial biomass and mineral-sorption combined
    ! real(r8) :: cue                          ! microbial carbon use efficiency (ratio of growth to uptake)
    !-----------------------------------------------------------------------

    associate(                                                                                     &
         ! Inputs
         ! clay                           => soilstate_inst%cellclay_col                           , & ! Input:  [real(r8)          (:, :)   ]  column 3D clay
         ! Outputs
         donor_pool                     => decomp_cascade_con%cascade_donor_pool                 , & ! Output: [integer           (:)     ]  which pool is C taken from for a given decomposition step
         receiver_pool                  => decomp_cascade_con%cascade_receiver_pool              , & ! Output: [integer           (:)     ]  which pool is C added to for a given decomposition step
         cn_ratio_is_fixed              => decomp_cascade_con%floating_cn_ratio_decomp_pools     , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio
         is_litter                      => decomp_cascade_con%is_litter                          , & ! Output: [logical           (:)     ]  TRUE => pool is a litter pool
         is_soil                        => decomp_cascade_con%is_soil                            , & ! Output: [logical           (:)     ]  TRUE => pool is a soil pool
         is_cwd                         => decomp_cascade_con%is_cwd                             , & ! Output: [logical           (:)     ]  TRUE => pool is a cwd pool
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio                   , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools
         initial_stock                  => decomp_cascade_con%initial_stock                      , & ! Output: [real(r8)          (:)     ]  initial concentration for seeding at spinup
         initial_stock_soildepth        => decomp_cascade_con%initial_stock_soildepth            , & ! Output: [real(r8)          (:)     ]  soil depth for initial concentration for seeding at spinup
         is_metabolic                   => decomp_cascade_con%is_metabolic                       , & ! Output: [logical           (:)     ]  TRUE => pool is metabolic material
         is_cellulose                   => decomp_cascade_con%is_cellulose                       , & ! Output: [logical           (:)     ]  TRUE => pool is cellulose
         is_lignin                      => decomp_cascade_con%is_lignin                          , & ! Output: [logical           (:)     ]  TRUE => pool is lignin
         spinup_factor                  => decomp_cascade_con%spinup_factor                      , & ! Output: [real(r8)          (:)     ]  factor for AD spinup associated with each pool
         ! decomp_cpools_vr               => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col  , &  ! Input: [real(r8) (:,:,:) ] (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
         begc                           => bounds%begc                                           , &
         endc                           => bounds%endc                                             &
         )

      ! allocate(rf_s1s2(begc:endc, 1:nlevdecomp))
      ! allocate( f_s1s2(begc:endc, 1:nlevdecomp))
      ! allocate( f_s1s3(begc:endc, 1:nlevdecomp))

      !------- time-constant coefficients ---------- !
      ! set soil organic matter compartment C:N ratios
      cn_s1 = params_inst%cn_dom
      cn_s2 = params_inst%cn_mic
      cn_s3 = params_inst%cn_mac

      ! set respiration fractions for fluxes between compartments.  In Somic model, CO2 is only evolved from microbial biomass thorugh respiration. rf for other transitions is zero.
      ! rf_l1s1 = 0.0_r8
      ! rf_l2s1 = 0.0_r8
      ! rf_l3s1 = 0.0_r8
      ! rf_s2s1 = 0.0_r8
      ! rf_s3s1 = 0.0_r8
      rf_cwdl2 = params_inst%rf_cwdl2
      rf_cwdl3 = params_inst%rf_cwdl3

      ! set the cellulose and lignin fractions for coarse woody debris
      cwd_fcel = params_inst%cwd_fcel

      ! set path fractions
      ! f_s2s1 = 1.0_r8
      ! f_s3s1 = 1.0_r8
      ! Note that some of the rf and f coefficents vary with time, hence not set here.  These include:
      !   rf_s1s2, f_s1s2, f_s1s3

      initial_stock_soildepth = params_inst%bgc_initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      !----------------  Litter pools
      i_litr_min = 1
      i_met_lit = i_litr_min
      cn_ratio_is_fixed(i_met_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_met_lit) = 'litr1'
      decomp_cascade_con%decomp_pool_name_history(i_met_lit) = 'MET_LIT'
      decomp_cascade_con%decomp_pool_name_long(i_met_lit) = 'metabolic litter'
      decomp_cascade_con%decomp_pool_name_short(i_met_lit) = 'L1'
      is_litter(i_met_lit) = .true.
      is_soil(i_met_lit) = .false.
      is_cwd(i_met_lit) = .false.
      initial_cn_ratio(i_met_lit) = 90._r8
      initial_stock(i_met_lit) = params_inst%bgc_initial_Cstocks(i_met_lit)
      is_metabolic(i_met_lit) = .true.
      is_cellulose(i_met_lit) = .false.
      is_lignin(i_met_lit) = .false.

      i_cel_lit = i_met_lit + 1
      cn_ratio_is_fixed(i_cel_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_cel_lit) = 'litr2'
      decomp_cascade_con%decomp_pool_name_history(i_cel_lit) = 'CEL_LIT'
      decomp_cascade_con%decomp_pool_name_long(i_cel_lit) = 'cellulosic litter'
      decomp_cascade_con%decomp_pool_name_short(i_cel_lit) = 'L2'
      is_litter(i_cel_lit) = .true.
      is_soil(i_cel_lit) = .false.
      is_cwd(i_cel_lit) = .false.
      initial_cn_ratio(i_cel_lit) = 90._r8
      initial_stock(i_cel_lit) = params_inst%bgc_initial_Cstocks(i_cel_lit)
      is_metabolic(i_cel_lit) = .false.
      is_cellulose(i_cel_lit) = .true.
      is_lignin(i_cel_lit) = .false.

      i_lig_lit = i_cel_lit + 1
      cn_ratio_is_fixed(i_lig_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_lig_lit) = 'litr3'
      decomp_cascade_con%decomp_pool_name_history(i_lig_lit) = 'LIG_LIT'
      decomp_cascade_con%decomp_pool_name_long(i_lig_lit) = 'lignin litter'
      decomp_cascade_con%decomp_pool_name_short(i_lig_lit) = 'L3'
      is_litter(i_lig_lit) = .true.
      is_soil(i_lig_lit) = .false.
      is_cwd(i_lig_lit) = .false.
      initial_cn_ratio(i_lig_lit) = 90._r8
      initial_stock(i_lig_lit) = params_inst%bgc_initial_Cstocks(i_lig_lit)
      is_metabolic(i_lig_lit) = .false.
      is_cellulose(i_lig_lit) = .false.
      is_lignin(i_lig_lit) = .true.

      i_litr_max = i_lig_lit
      if (i_litr_min /= 1 .or. i_litr_max < 2 .or. i_litr_max > 3) then
         write(iulog, *) 'Expecting i_litr_min = 1 and i_litr_max = 2 or 3.'
         write(iulog, *) 'See pftconMod, SoilBiogeochemCarbonFluxType, and'
         write(iulog, *) 'clmfates_interfaceMod for ramifications of changing'
         write(iulog, *) 'this assumption.'
         call endrun(msg='ERROR: i_litr_min and/or i_litr_max out of range '// &
              errMsg(sourcefile, __LINE__))
      end if

      !----------------  SOM pools
      i_doc_som = i_lig_lit + 1
      cn_ratio_is_fixed(i_doc_som) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_doc_som) = 'soil1'
      decomp_cascade_con%decomp_pool_name_history(i_doc_som) = 'DOC_SOM'
      decomp_cascade_con%decomp_pool_name_long(i_doc_som) = 'dissolved soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_doc_som) = 'S1'
      is_litter(i_doc_som) = .false.
      is_soil(i_doc_som) = .true.
      is_cwd(i_doc_som) = .false.
      initial_cn_ratio(i_doc_som) = cn_s1
      initial_stock(i_doc_som) = params_inst%bgc_initial_Cstocks(i_doc_som)
      is_metabolic(i_doc_som) = .false.
      is_cellulose(i_doc_som) = .false.
      is_lignin(i_doc_som) = .false.

      i_mic_som = i_doc_som + 1
      cn_ratio_is_fixed(i_mic_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_mic_som) = 'soil2'
      decomp_cascade_con%decomp_pool_name_history(i_mic_som) = 'MIC_SOM'
      decomp_cascade_con%decomp_pool_name_long(i_mic_som) = 'microbial soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_mic_som) = 'S2'
      is_litter(i_mic_som) = .false.
      is_soil(i_mic_som) = .true.
      is_cwd(i_mic_som) = .false.
      initial_cn_ratio(i_mic_som) = cn_s2
      initial_stock(i_mic_som) = params_inst%bgc_initial_Cstocks(i_mic_som)
      is_metabolic(i_mic_som) = .false.
      is_cellulose(i_mic_som) = .false.
      is_lignin(i_mic_som) = .false.

      i_mac_som = i_mic_som + 1
      cn_ratio_is_fixed(i_mac_som) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_mac_som) = 'soil3'
      decomp_cascade_con%decomp_pool_name_history(i_mac_som) = 'MAC_SOM'
      decomp_cascade_con%decomp_pool_name_long(i_mac_som) = 'mineral-associated soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_mac_som) = 'S3'
      is_litter(i_mac_som) = .false.
      is_soil(i_mac_som) = .true.
      is_cwd(i_mac_som) = .false.
      initial_cn_ratio(i_mac_som) = cn_s3
      initial_stock(i_mac_som) = params_inst%bgc_initial_Cstocks(i_mac_som)
      is_metabolic(i_mac_som) = .false.
      is_cellulose(i_mac_som) = .false.
      is_lignin(i_mac_som) = .false.

      !----------------  Coarse woody debris pool
      if (.not. use_fates) then
         ! CWD
         i_cwd = i_mac_som + 1
         cn_ratio_is_fixed(i_cwd) = .true.
         decomp_cascade_con%decomp_pool_name_restart(i_cwd) = 'cwd'
         decomp_cascade_con%decomp_pool_name_history(i_cwd) = 'CWD'
         decomp_cascade_con%decomp_pool_name_long(i_cwd) = 'coarse woody debris'
         decomp_cascade_con%decomp_pool_name_short(i_cwd) = 'CWD'
         is_litter(i_cwd) = .false.
         is_soil(i_cwd) = .false.
         is_cwd(i_cwd) = .true.
         initial_cn_ratio(i_cwd) = 90._r8
         initial_stock(i_cwd) = params_inst%bgc_initial_Cstocks(i_cwd)
         is_metabolic(i_cwd) = .false.
         is_cellulose(i_cwd) = .false.
         is_lignin(i_cwd) = .false.
      endif

      !------------------  Set factor scalars for accelerated spin-up  ---------------!
      ! for now we set all spinup cators equal to one (i.e. no accelrated spinup).  will revisit this decision during model validation and testing.
      speedup_fac = 1._r8
      spinup_factor(i_met_lit) = 1._r8
      spinup_factor(i_cel_lit) = 1._r8
      spinup_factor(i_lig_lit) = 1._r8
      if (.not. use_fates) then
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * CNParamsShareInst%tau_cwd / 2._r8 ))
      end if
      spinup_factor(i_doc_som) = 1._r8
      spinup_factor(i_mic_som) = 1._r8 ! max(1._r8, (speedup_fac * params_inst%tau_s2_bgc))
      spinup_factor(i_mac_som) = 1._r8 ! max(1._r8, (speedup_fac * params_inst%tau_s3_bgc))
      if ( masterproc ) then
         write(iulog, *) 'Spinup_state ', spinup_state
         write(iulog, *) 'Spinup factors ', spinup_factor
      end if

      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1s1 = 1
      decomp_cascade_con%cascade_step_name(i_l1s1) = 'L1S1'
      donor_pool(i_l1s1) = i_met_lit
      receiver_pool(i_l1s1) = i_doc_som

      i_l2s1 = 2
      decomp_cascade_con%cascade_step_name(i_l2s1) = 'L2S1'
      donor_pool(i_l2s1) = i_cel_lit
      receiver_pool(i_l2s1) = i_doc_som

      i_l3s1 = 3
      decomp_cascade_con%cascade_step_name(i_l3s1) = 'L3S2'
      donor_pool(i_l3s1) = i_lig_lit
      receiver_pool(i_l3s1) = i_doc_som

      i_s1s2 = 4
      decomp_cascade_con%cascade_step_name(i_s1s2) = 'S1S2'
      donor_pool(i_s1s2) = i_doc_som
      receiver_pool(i_s1s2) = i_mic_som

      i_s1s3 = 5
      decomp_cascade_con%cascade_step_name(i_s1s3) = 'S1S3'
      donor_pool(i_s1s3) = i_doc_som
      receiver_pool(i_s1s3) = i_mac_som

      i_s2s1 = 6
      decomp_cascade_con%cascade_step_name(i_s2s1) = 'S2S1'
      donor_pool(i_s2s1) = i_mic_som
      receiver_pool(i_s2s1) = i_doc_som

      i_s3s1 = 7
      decomp_cascade_con%cascade_step_name(i_s3s1) = 'S3S1'
      donor_pool(i_s3s1) = i_mac_som
      receiver_pool(i_s3s1) = i_doc_som

      if (.not. use_fates) then
         i_cwdl2 = 8
         decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
         donor_pool(i_cwdl2) = i_cwd
         receiver_pool(i_cwdl2) = i_cel_lit

         i_cwdl3 = 9
         decomp_cascade_con%cascade_step_name(i_cwdl3) = 'CWDL3'
         donor_pool(i_cwdl3) = i_cwd
         receiver_pool(i_cwdl3) = i_lig_lit
      end if

      deallocate(params_inst%bgc_initial_Cstocks)
    end associate
  end subroutine init_decompcascade_somic




  !----- Use the CENTURY T response function
  real(r8) function ft_somic(t1)
    use shr_const_mod , only : SHR_CONST_PI
    real(r8), intent(in) :: t1
    ft_somic = 11.75_r8 + (29.7_r8 / SHR_CONST_PI) * atan(SHR_CONST_PI * 0.031_r8  * (t1 - 15.4_r8))
  end function ft_somic


  real(r8) function cue (t1)
    real(r8), intent(in) :: t1
    real(r8), parameter :: ref_temp = 25._r8   ! named constant for reference temperature for normalization of cue
    real(r8), parameter :: cue_min = 1.e-6_r8  ! named constant for minimum CUE
    real(r8), parameter :: cue_max = 0.7_r8    ! named constant for maximum CUE

    ! Temperature dependence of microbial carbon use efficiency derived from the following sources:
    ! http://onlinelibrary.wiley.com.proxy.library.cornell.edu/doi/10.1890/15-2110.1/full
    ! http://onlinelibrary.wiley.com/doi/10.1111/ele.12113/full
    ! http://link.springer.com/article/10.1007/s10533-013-9948-8
    cue = params_inst%cue_0 - (t1 - ref_temp) * params_inst%mcue
    if (cue < cue_min) then
       cue = cue_min
    end if
    if (cue > cue_max) then
       cue = cue_max
    end if
  end function cue

  !-----------------------------------------------------------------------
  subroutine decomp_rate_constants_somic(bounds, num_soilc, filter_soilc, soilstate_inst, temperature_inst, &
                                         ch4_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    !  calculate rate constants and decomposition pathways for the SOMic decomposition model
    !  written by D. Woolf, based on CLM5 bgc decomposition cascade
    !
    ! !USES:
    use clm_time_manager , only : get_average_days_per_year
    use shr_const_mod    , only : SHR_CONST_PI
    use clm_varcon       , only : secspday
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds
    integer                              , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(temperature_type)               , intent(in)    :: temperature_inst
    type(ch4_type)                       , intent(in)    :: ch4_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_carbonstate_inst

    !
    ! !LOCAL VARIABLES:
    integer :: c, fc, j, k, l                                ! indices
    real(r8), parameter :: eps = 1.e-6_r8                    ! named constant for floating point logic
    real(r8), parameter :: ref_clay = 23._r8                 ! named constant for reference clay content for normalization of sorption rates
    real(r8) :: frw(bounds%begc:bounds%endc)                 ! rooting fraction weight
    real(r8), allocatable:: fr(:,:)                          ! column-level rooting fraction by soil depth
    real(r8) :: m_scalar(bounds%begc:bounds%endc, 1:nlevdecomp) ! Michaelis-Menten rate scalar
    real(r8) :: psi                                          ! temporary soilpsi for water scalar
    real(r8) :: k_l1s1                                       ! base rate constant for decomposition of l1
    real(r8) :: k_l2s1                                       ! base rate constant for decomposition of l2
    real(r8) :: k_l3s1                                       ! base rate constant for decomposition of l3
    real(r8) :: k_sorb                                       ! temporary sorption rate constant
    real(r8) :: k_mic_up                                     ! temporary microbial uptake rate constant
    real(r8) :: k_s1s2                                       ! base rate constant for microbial uptake
    real(r8) :: k_s1s3                                       ! base rate constant for sorption of doc
    real(r8) :: k_s2s1                                       ! base rate constant for loss of microbial biomass (including both death and exudates)
    real(r8) :: k_s3s1                                       ! base rate constant for desorption of mineral-associated OM
    real(r8) :: clay_scalar                                  ! scalar to modify sorption rate depending on clay content
    real(r8) :: f_sorb                                       ! temporary sorption partition fraction
    real(r8) :: f_mic_up                                     ! temporary microbial uptake partition fraction
    real(r8) :: f_resp                                       ! temporary respiration fraction
    real(r8) :: k_frag                                       ! fragmentation rate constant CWD (1/sec)
    real(r8) :: f_growth                                     ! fraction of microbial uptake allocated to growth
    real(r8) :: Q10                                          ! temperature dependence
    real(r8) :: froz_q10                                     ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
    real(r8) :: decomp_depth_efolding                        ! (meters) e-folding depth for reduction in decomposition [
    real(r8) :: ft_somic_30                                  ! reference rate at 30C
    real(r8) :: normalization_factor                         ! factor by which to offset the decomposition rates frm century to a q10 formulation
    real(r8) :: days_per_year                                ! days per year
    real(r8) :: mino2lim                                     ! minimum anaerobic decomposition rate
    real(r8) :: spinup_geogterm_l1(bounds%begc:bounds%endc)  ! geographically-varying spinup term for l1
    real(r8) :: spinup_geogterm_l2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for l2
    real(r8) :: spinup_geogterm_l3(bounds%begc:bounds%endc)  ! geographically-varying spinup term for l3
    real(r8) :: spinup_geogterm_cwd(bounds%begc:bounds%endc) ! geographically-varying spinup term for cwd
    real(r8) :: spinup_geogterm_s1(bounds%begc:bounds%endc)  ! geographically-varying spinup term for s1
    real(r8) :: spinup_geogterm_s2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for s2
    real(r8) :: spinup_geogterm_s3(bounds%begc:bounds%endc)  ! geographically-varying spinup term for s3

    !-----------------------------------------------------------------------

    associate(                                                             &
         cwd_flig         => CNParamsShareInst%cwd_flig                  , &                ! Input:  [real(r8)         ]  lignin fraction of coarse woody debris (frac)
         rf_cwdl2         => CNParamsShareInst%rf_cwdl2                  , &                ! Input:  [real(r8)         ]  respiration fraction in CWD to litter2 transition (frac)
         minpsi           => CNParamsShareInst%minpsi                    , &                ! Input:  [real(r8)         ]  minimum soil suction (mm)
         maxpsi           => CNParamsShareInst%maxpsi                    , &                ! Input:  [real(r8)         ]  maximum soil suction (mm)
         soilpsi          => soilstate_inst%soilpsi_col                  , &                ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)
         t_soisno         => temperature_inst%t_soisno_col               , &                ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         o2stress_sat     => ch4_inst%o2stress_sat_col                   , &                ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         o2stress_unsat   => ch4_inst%o2stress_unsat_col                 , &                ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         finundated       => ch4_inst%finundated_col                     , &                ! Input:  [real(r8) (:)     ]  fractional inundated area
         rf               => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col       , & ! Output: [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)
         pathfrac         => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col , & ! Output: [real(r8) (:,:,:) ]  what fraction of C passes from donor to receiver pool through a given transition (frac)
         t_scalar         => soilbiogeochem_carbonflux_inst%t_scalar_col , &                ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp
         w_scalar         => soilbiogeochem_carbonflux_inst%w_scalar_col , &                ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp
         o_scalar         => soilbiogeochem_carbonflux_inst%o_scalar_col , &                ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
         decomp_k         => soilbiogeochem_carbonflux_inst%decomp_k_col , &                ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         spinup_factor    => decomp_cascade_con%spinup_factor            , &                ! Input:  [real(r8) (:)     ]  factor for AD spinup associated with each pool
         decomp_cpools_vr => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col , &       ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
         clay             => soilstate_inst%cellclay_col                 , &                ! Input:  [real(r8) (:,:)   ]  column 3D clay (%)
         begc             => bounds%begc                                 , &
         endc             => bounds%endc                                   &
         )

      mino2lim = CNParamsShareInst%mino2lim

      if ( use_century_tfunc .and. normalize_q10_to_century_tfunc ) then
         call endrun(msg='ERROR: cannot have both use_century_tfunc and normalize_q10_to_century_tfunc set as true'//&
              errMsg(sourcefile, __LINE__))
      endif

      Q10 = CNParamsShareInst%Q10                                             ! set "Q10" parameter
      froz_q10  = CNParamsShareInst%froz_q10                                  ! set "froz_q10" parameter
      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding         ! Set "decomp_depth_efolding" parameter

      ! base rate constants for transitions (must be in per-second, otherwise we should insert a conversion here)
      k_l1s1 = params_inst%k_l1s1
      k_l2s1 = params_inst%k_l2s1
      k_l3s1 = params_inst%k_l3s1
      k_s1s2 = params_inst%k_s1s2
      k_s1s3 = params_inst%k_s1s3
      k_s2s1 = params_inst%k_s2s1
      k_s3s1 = params_inst%k_s3s1
      k_frag = 1._r8  / (secspday * days_per_year * CNParamsShareInst%tau_cwd)

     ! calculate reference temperature rate scalar
      ft_somic_30 = ft_somic(30._r8)

      ! Set all accelerated spin-up factors to unity.
      ! TO DO: If later we are certain that acceleration is not required, then these factors and their code blocks should be deleted.
      if ( spinup_state >= 1 ) then
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            !
            if ( abs(spinup_factor(i_met_lit) - 1._r8) .gt. eps) then
               spinup_geogterm_l1(c) = spinup_factor(i_met_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_cel_lit) - 1._r8) .gt. eps) then
               spinup_geogterm_l2(c) = spinup_factor(i_cel_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l2(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_lig_lit) - 1._r8) .gt. eps) then
               spinup_geogterm_l3(c) = spinup_factor(i_lig_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l3(c) = 1._r8
            endif
            !
            if ( .not. use_fates ) then
               if ( abs(spinup_factor(i_cwd) - 1._r8) .gt. eps) then
                  spinup_geogterm_cwd(c) = spinup_factor(i_cwd) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
               else
                  spinup_geogterm_cwd(c) = 1._r8
               endif
            endif
            !
            if ( abs(spinup_factor(i_doc_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s1(c) = spinup_factor(i_doc_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_mic_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s2(c) = spinup_factor(i_mic_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s2(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_mac_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s3(c) = spinup_factor(i_mac_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s3(c) = 1._r8
            endif
            !
         end do
      else
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            spinup_geogterm_l1(c)  = 1._r8
            spinup_geogterm_l2(c)  = 1._r8
            spinup_geogterm_l3(c)  = 1._r8
            spinup_geogterm_s1(c)  = 1._r8
            spinup_geogterm_s2(c)  = 1._r8
            spinup_geogterm_s3(c)  = 1._r8
            spinup_geogterm_cwd(c) = 1._r8
         end do
      endif

      !--- time dependent coefficients-----!
      if ( nlevdecomp .eq. 1 ) then

         ! calculate function to weight the temperature and water potential scalars
         ! for decomposition control.

         ! the following normalizes values in fr so that they
         ! sum to 1.0 across top nlevdecomp levels on a column
         frw(begc:endc) = 0._r8
         nlev_soildecomp_standard = 5
         allocate(fr(begc:endc, nlev_soildecomp_standard))
         do j=1, nlev_soildecomp_standard
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               frw(c) = frw(c) + col%dz(c, j)
            end do
         end do
         do j = 1, nlev_soildecomp_standard
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               if (frw(c) /= 0._r8) then
                  fr(c, j) = col%dz(c, j) / frw(c)
               else
                  fr(c, j) = 0._r8
               end if
            end do
         end do

         if ( .not. use_century_tfunc ) then
            ! calculate rate constant scalar for soil temperature
            ! assuming that the base rate constants are assigned for non-moisture
            ! limiting conditions at 25 C.
            do j = 1, nlev_soildecomp_standard
               do fc = 1, num_soilc
                  c = filter_soilc(fc)
                  if (j==1) t_scalar(c, :) = 0._r8
                  if (t_soisno(c, j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c, 1) = t_scalar(c, 1) + (Q10**((t_soisno(c, j) - (SHR_CONST_TKFRZ + 25._r8))/10._r8)) * fr(c, j)
                  else
                     t_scalar(c, 1) = t_scalar(c, 1) + (Q10**(-25._r8/10._r8)) * (froz_q10**((t_soisno(c, j) - SHR_CONST_TKFRZ)/10._r8)) * fr(c, j)
                  endif
               end do
            end do

         else
            ! original daycent uses an arctangent function to calculate the temperature dependence of decomposition
            do j = 1, nlev_soildecomp_standard
               do fc = 1, num_soilc
                  c = filter_soilc(fc)
                  if (j==1) t_scalar(c, :) = 0._r8
                  t_scalar(c, 1) = t_scalar(c, 1) + max(ft_somic(t_soisno(c, j) - SHR_CONST_TKFRZ) / ft_somic_30 * fr(c, j), 0.01_r8)
               end do
            end do

         endif

         ! calculate the rate constant scalar for soil water content.
         ! Uses the log relationship with water potential given in
         ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
         ! a comparison of models. Ecology, 68(5):1190-1200.
         ! and supported by data in
         ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
         ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.
         do j = 1, nlev_soildecomp_standard
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               if (j==1) w_scalar(c,:) = 0._r8
               psi = min(soilpsi(c, j), maxpsi)
               ! decomp only if soilpsi is higher than minpsi
               if (psi > minpsi) then
                  w_scalar(c, 1) = w_scalar(c, 1) + (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c, j)
               end if
            end do
         end do

         if (use_lch4) then
            ! Calculate ANOXIA
            if (anoxia) then
               ! Check for anoxia w/o LCH4 now done in controlMod.

               do j = 1, nlev_soildecomp_standard
                  do fc = 1, num_soilc
                     c = filter_soilc(fc)

                     if (j==1) o_scalar(c, :) = 0._r8

                     o_scalar(c, 1) = o_scalar(c, 1) + fr(c, j) * max(o2stress_unsat(c, j), mino2lim)
                  end do
               end do
            else
               o_scalar(begc:endc, 1:nlevdecomp) = 1._r8
            end if
         else
            o_scalar(begc:endc, 1:nlevdecomp) = 1._r8
         end if
         deallocate(fr)

      else ! nlevdecomp is greater than 1

         if ( .not. use_century_tfunc ) then
            ! calculate rate constant scalar for soil temperature
            ! assuming that the base rate constants are assigned for non-moisture
            ! limiting conditions at 25 C.
            ! Peter Thornton: 3/13/09
            ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
            ! as part of the modifications made to improve the seasonal cycle of
            ! atmospheric CO2 concentration in global simulations. This does not impact
            ! the base rates at 25 C, which are calibrated from microcosm studies.
            do j = 1, nlevdecomp
               do fc = 1, num_soilc
                  c = filter_soilc(fc)
                  if (t_soisno(c, j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c, j) = (Q10**((t_soisno(c, j) - (SHR_CONST_TKFRZ + 25._r8))/10._r8))
                  else
                     ! separate q10 calculation for frozen soils
                     t_scalar(c, j) = (Q10**(-25._r8/10._r8)) * (froz_q10**((t_soisno(c, j) - SHR_CONST_TKFRZ) / 10._r8))
                  endif
               end do
            end do

         else

            do j = 1, nlevdecomp
               do fc = 1, num_soilc
                  c = filter_soilc(fc)
                  ! t scalar is normalized to 1 at 30C, and cannot be less than 0.01
                  t_scalar(c, j) = max(ft_somic(t_soisno(c, j) - SHR_CONST_TKFRZ) / ft_somic_30, 0.01_r8)
               end do
            end do

         endif

         ! calculate the rate constant scalar for soil water content.
         ! Uses the log relationship with water potential given in
         ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
         ! a comparison of models. Ecology, 68(5):1190-1200.
         ! and supported by data in
         ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
         ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.
         do j = 1, nlevdecomp
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               psi = min(soilpsi(c, j), maxpsi)
               ! decomp only if soilpsi is higher than minpsi
               if (psi > minpsi) then
                  w_scalar(c, j) = (log(minpsi / psi) / log(minpsi / maxpsi))
               else
                  w_scalar(c, j) = 0._r8
               end if
            end do
         end do

         if (use_lch4) then
            ! Calculate ANOXIA
            ! Check for anoxia w/o LCH4 now done in controlMod.
            if (anoxia) then
               do j = 1, nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc(fc)

                     o_scalar(c, j) = max(o2stress_unsat(c, j), mino2lim)
                  end do
               end do
            else
               o_scalar(begc:endc, 1:nlevdecomp) = 1._r8
            end if
         else
            o_scalar(begc:endc, 1:nlevdecomp) = 1._r8
         end if
      end if

      if ( normalize_q10_to_century_tfunc ) then
         ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
         normalization_factor = (ft_somic(normalization_tref) / ft_somic_30) / (q10**((normalization_tref - 25._r8) / 10._r8))
         do j = 1, nlevdecomp
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               t_scalar(c, j) = t_scalar(c, j) * normalization_factor
            end do
         end do
      endif

      ! calculate rate constants and path/respiration fractions all litter and som pools and their transitions
      write(iulog,*) 'SOMIC: calculate rate constants and path/respiration fractions...'
      do j = 1, nlevdecomp
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! calculate the rate constant scalar for Michaelis-Menten dynamics to relate microbial biomass to rate
            m_scalar(c, j) = params_inst%mic_vmax * decomp_cpools_vr(c, j, i_mic_som) / &
                            (params_inst%mic_km   + decomp_cpools_vr(c, j, i_mic_som))
            ! ensure that m scalar does not fall below minimum threshold, in case of soil columns starting without any microbial biomass
            if (m_scalar(c, j) < eps) then
               m_scalar(c, j) = eps
            end if

            ! partitioning of doc between competing processes of sorption and microbial uptake
            clay_scalar = 1.0_r8 + params_inst%mclay * (clay(c, j) - ref_clay)
            k_sorb   = params_inst%k_s1s3 * clay_scalar    * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j)    ! modified sorption rate
            ! ensure that k_sorb does not fall below minimum threshold, in case any scalars are zero
            if (k_sorb < eps) then
               k_sorb = eps
            end if
            k_mic_up = params_inst%k_s1s2 * m_scalar(c, j) * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j)    ! modified microbial uptake rate
            ! ensure that k_mic_up does not fall below minimum threshold, in case any scalars are zero
            if (k_mic_up < eps) then
               k_mic_up = eps
            end if
            f_sorb   = k_sorb / (k_sorb + k_mic_up)                                                              ! fraction of doc turnover sorbed
            f_mic_up = 1.0_r8 - f_sorb                                                                           ! fraction of doc turnover taken up by microbes
            f_growth = f_mic_up * cue(t_soisno(c, j) - SHR_CONST_TKFRZ)                                          ! fraction of doc turnover to anabolic growth
            f_resp   = f_mic_up - f_growth                                                                       ! fraction of doc turnover that is respired to CO2

            ! rate constants for decomposition of pools
            decomp_k(c, j, i_met_lit) = k_l1s1 * m_scalar(c, j) * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j) * spinup_geogterm_l1(c)
            decomp_k(c, j, i_cel_lit) = k_l2s1 * m_scalar(c, j) * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j) * spinup_geogterm_l2(c)
            decomp_k(c, j, i_lig_lit) = k_l3s1 * m_scalar(c, j) * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j) * spinup_geogterm_l3(c)
            decomp_k(c, j, i_doc_som) = (f_sorb * k_sorb + f_mic_up * k_mic_up) * spinup_geogterm_s1(c)           ! rate constant for doc loss is weighted mean of sorption and microbial uptake
            decomp_k(c, j, i_mic_som) = k_s2s1 * m_scalar(c, j) * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j) * spinup_geogterm_s2(c)
            decomp_k(c, j, i_mac_som) = k_s3s1 * m_scalar(c, j) * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j) * spinup_geogterm_s3(c)
            ! same for cwd but only if fates is not enabled; fates handles CWD on its own structure
            if (.not. use_fates) then
               decomp_k(c, j, i_cwd) = k_frag * t_scalar(c, j) * w_scalar(c, j) * o_scalar(c, j) * spinup_geogterm_cwd(c)
            end if

            ! Other than doc, which is partitioned between sorption and uptake, all other transitions have a single recipient pool
            ! TO DO: since these never change, we should really just set them once during initialisation
            pathfrac(c, j, i_l1s1) = 1.0_r8
            pathfrac(c, j, i_l2s1) = 1.0_r8
            ! pathfrac(c, j, i_l3s2) = 1.0_r8
            pathfrac(c, j, i_s1s2) = f_growth
            pathfrac(c, j, i_s1s3) = f_sorb
            pathfrac(c, j, i_s2s1) = 1.0_r8
            ! pathfrac(c, j, i_s2s3) = 1.0_r8
            pathfrac(c, j, i_s3s1) = 1.0_r8

            ! Here we set respiration fractions. Note that only microbes respire
            rf(c, j, i_l1s1) = 0.0_r8
            rf(c, j, i_l2s1) = 0.0_r8
            ! rf(c, j, i_l3s2) = 0.0_r8
            rf(c, j, i_s1s2) = f_resp
            rf(c, j, i_s1s3) = 0.0_r8
            rf(c, j, i_s2s1) = 0.0_r8
            ! rf(c, j, i_s2s3) = 0.0_r8
            rf(c, j, i_s3s1) = 0.0_r8

           write(iulog,*) '------------- c, j:', c, j
           write(iulog,*) 'decomp_cpools_vr', decomp_cpools_vr(c, j, :)
           write(iulog,*) 'scalars: m, clay, t, w, o', &
                           m_scalar(c, j), clay_scalar, t_scalar(c, j), w_scalar(c, j), o_scalar(c, j)
           write(iulog,*) 'k_sorb, k_mic_up', k_sorb, k_mic_up
           write(iulog,*) 'f_sorb, f_mic_up, f_growth, f_resp', &
                           f_sorb, f_mic_up, f_growth, f_resp
           write(iulog,*) 'decomp_k', decomp_k(c, j, :)
           write(iulog,*) '-------------'
         end do
      end do

      ! calculate path fractions for cwd
      if (.not. use_fates) then
         pathfrac(begc:endc, 1:nlevdecomp, i_cwdl2) = cwd_fcel
         pathfrac(begc:endc, 1:nlevdecomp, i_cwdl3) = cwd_flig
         rf(begc:endc, 1:nlevdecomp, i_cwdl2) = rf_cwdl2
         rf(begc:endc, 1:nlevdecomp, i_cwdl3) = rf_cwdl2
      end if
    end associate
 end subroutine decomp_rate_constants_somic
end module SoilBiogeochemDecompCascadeSOMicMod
