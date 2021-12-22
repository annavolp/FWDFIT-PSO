
; NAME: 
;   vis_fwdfit_pso_nov2021
;
; PURPOSE:
;   forward fitting method from visibility based on Particle Swarm Optimization
;
; CALLING SEQUENCE:
;   vis_fwdfit_pso_nov2021, type, vis, SwarmSize, uncertainty, param_opt, seedstart
;   
; CALLS:
;   cmreplicate                     [replicates an array or scalar into a larger array, as REPLICATE does.]
;   vis_fwdfit_func_pso             [to calculate model visibilities for a given set of source parameters]
;   pso_func_makealoop_nov2021      [to calculate loop Fourier trasform]
;   swarmintelligence               [to optimize]
;   vis_fwdfit_src_structure        [to create the source structure]
;   vis_fwdfit_src_bifurcate        [to create a modified source structure based on bifurcation of input source structure]
;   vis_fwdfit_pso_source2map       [to create the map]
;   vis_fwdfit_pso_vis_pred         [to create the fit map]
;   
;   
; INPUTS:
;   type: parametric shape to use for the forward fitting method
;         - 'circle' : Gaussian circular source
;         - 'ellipse': Gaussian elliptical source
;         - 'multi'  : double Gaussian circular source
;         - 'loop'   : single curved elliptical gaussian
;
;   vis         : struct containing  the observed visibility
;     vis.obsvis: array containing the values of the observed visibility amplitudes
;     vis.sigamp: array containing the values of the errors on the observed visibility amplitudes
;     vis.u     : u coordinates of the sampling frequencies
;     vis.v     : v coordinates of the sampling frequencies
;
;
; KEYWORDS:
;   lb: array containing the lower bound values of the variables to optimize
;   ub: array containing the upper bound values of the variables to optimize
;
;   For different shapes we have:
;
;       - 'circle'  : lb,ub = [flux, x location, y location, FWHM]
;       - 'ellipse' : lb,ub = [flux, FWHM, ecc * cos(alpha), ecc * sin(alpha), x location, y location] 
;                     'ecc' is the eccentricity of the ellipse and 'alpha' is the orientation angle
;       - 'multi'   : lb,ub = [FWHM1, flux1, FWHM2, flux2, x location, y location, x2 location, y2 location] 
;       - 'loop'    : lb,ub = [flux, FWHM, ecc * cos(alpha), ecc * sin(alpha), x location, y location, loop_angle]
;       
;   param_opt: array containing the values of the parameters to keep fixed during the optimization. 
;              If an entry of 'param_opt' is set equal to 'fit', then the corresponding variable is optimized. 
;              Otherwise, its value is kept fixed equal to the entry of 'param_opt'
;
;   For different shapes we have:
;
;       - 'circle'  : param_opt = [flux, x location, y location, FWHM]
;       - 'ellipse' : param_opt = [flux, FWHM max, FWHM min, alpha, x location, y location] 
;                     'alpha' is the orientation angle of the source
;       - 'multi'   : param_opt = [FWHM1, flux1, FWHM2, flux2, x1 location, y1 location, x2 location, y2 location] 
;       - 'loop'    : param_opt = [flux, FWHM max, FWHM min, alpha, x location, y location, loop_angle]
;
;
;   SwarmSize   : number of particles used in PSO (default is 100)
;   TolFun      : tolerance for the stopping criterion (default is 1e-6)
;   maxiter     : maximum number of iterations of PSO 
;                 (defult is the product between of the numbers of parameters and the number of particles)
;   uncertainty : set to 1 for the computation of the parameters uncertainty (confidence strip approach)
;   silent      : set to 1 for avoiding the print of the retrieved parameters
;   
;   
;   SRCSTR names a source structure array to receive the fitted source parameters.
;   FITSIGMAS returns sigma in fitted quantities in SRCSTR.
;
; HISTORY: November 2021, Massa P., Volpara A. created
;
; CONTACT:
;   massa.p [at] dima.unige.it
;   volpara [at] dima.unige.it

function vis_fwdfit_pso_nov2021, type, vis, $
  lb = lb, ub = ub, $
  param_opt = param_opt, $
  SwarmSize = SwarmSize, TolFun = TolFun, maxiter = maxiter, $
  uncertainty = uncertainty, $
  imsize=imsize, pixel=pixel, $
  silent = silent, $
  srcstr = srcstrout_pso, $
  fitsigmas =fitsigmasout_pso, $
  seedstart = seedstart, $
  no_plot_fit = no_plot_fit


  default, SwarmSize, 100.
  default, TolFun, 1e-06
  default, silent, 0
  default, imsize, [128,128]
  default, pixel, [1.,1.]
  default, seedstart, 0  
  default, no_plot_fit, 0

  phi= max(abs(vis.obsvis)) ;estimate_flux

  case type of

    'circle': begin
      default, param_opt, ['fit', 'fit', 'fit', 'fit'] ; flux, x location, y location, fwhm
      default, lb, [0.1*phi, -100., -100., 1.]
      default, ub, [1.5*phi, 100., 100., 100.]
    end

    'ellipse': begin
      default, param_opt, ['fit', 'fit', 'fit', 'fit', 'fit', 'fit'] ; flux, FWHM max, FWHM min, alpha, x location, y location
      default, lb, [0.1*phi,  1., -5., 0., -100., -100.]
      default, ub, [1.5*phi, 100., 5., 1., 100., 100.]
    end

    'multi': begin
      default, param_opt, ['fit' , 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'] ;FWHM1, flux1, FWHM2, flux2, x1 location, y1 location, x2 location, y2 location
      default, lb, [0.,  0.1*phi, 0.,  0.1*phi, -100., -100., -100., -100.]
      default, ub, [100., 1.5*phi, 100., 1.5*phi, 100., 100., 100., 100.]
    end

    'loop': begin
      default, param_opt, ['fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'] ; flux, FWHM max, FWHM min, alpha, x location, y location, curvature
      default, lb, [0.1*phi,  1., -5., 0., -100., -100., -180.]
      default, ub, [1.5*phi, 100., 5., 1., 100., 100., 180.]
    end

  endcase


  fun_name = 'vis_fwdfit_func_pso'
  Nvars = n_elements(lb)

  visobs = [real_part(vis.obsvis), imaginary(vis.obsvis)]
  nvis   = N_ELEMENTS(visobs)
  vvisobs = transpose(cmreplicate(visobs, SwarmSize))
  sigamp  = [vis.sigamp,vis.sigamp]
  ssigamp = transpose(cmreplicate(sigamp, SwarmSize))

  extra = {type: type, $
    visobs: vvisobs, $
    sigamp: ssigamp, $
    u: vis.u, $
    v: vis.v, $
    n_free: nvis - Nvars, $    ;n_free: degrees of freedom (difference between the number of visibility amplitudes 
                               ;and the number of parameters of the source shape)
    param_opt: param_opt, $
    mapcenter : vis.xyoffset }

  if ~keyword_set(maxiter) then maxiter = Nvars*SwarmSize

  if type eq 'circle' then begin
    
    n_free = nvis - 4.

    if (n_elements(param_opt) ne 4) or (n_elements(lb) ne 4) or (n_elements(ub) ne 4) then begin
        UNDEFINE, lb, ub, param_opt
        message, 'Wrong number of elements of lower bound, upper bound or parameter mask'
    endif
    
    optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
    xopt = optim_f.xopt

    srcstr = {vis_fwdfit_src_structure}
    srcstr.srctype ='circle'

    fitsigmas = {vis_fwdfit_src_structure}
    fitsigmas.srctype ='std.dev'

    srcstr.srcflux         = xopt[0]
    srcstr.srcx            = -xopt[2]+vis[0].xyoffset[0]
    srcstr.srcy            = xopt[1]+vis[0].xyoffset[1]
    srcstr.srcfwhm_max     = xopt[3]
    srcstr.srcfwhm_min     = xopt[3]

    if keyword_set(uncertainty) then begin

      print, ' '
      print, 'Uncertainty: '
      print, '

      ntry = 20

      trial_results = fltarr(Nvars, ntry)
      iseed = findgen(ntry)+seedstart

      for n=0,ntry-1 do begin

        testerror = RANDOMN(iseed[n], nvis)
        vistest   = visobs + testerror * sigamp
        vistest   = transpose(cmreplicate(vistest, SwarmSize))

        extra = {type: type, $
          visobs: vistest, $
          sigamp: ssigamp, $
          u: vis.u, $
          v: vis.v, $
          n_free: nvis - Nvars, $
          param_opt: param_opt, $
          mapcenter : vis.xyoffset}

        optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
        xopt = optim_f.xopt

        trial_results[*,n]  = xopt

      endfor

      std_dev_par = stddev(trial_results, dimension=2)

      fitsigmas.srcflux      = std_dev_par[0]
      fitsigmas.srcx         = std_dev_par[2]
      fitsigmas.srcy         = std_dev_par[1]
      fitsigmas.srcfwhm_max  = std_dev_par[3]
      fitsigmas.srcfwhm_min  = std_dev_par[3]

    endif

  endif

  if type eq 'ellipse' then begin
    
    n_free = nvis-6.

    if (n_elements(param_opt) ne 6) or (n_elements(lb) ne 6) or (n_elements(ub) ne 6) then begin
      UNDEFINE, lb, ub, param_opt
      message, 'Wrong number of elements of lower bound, upper bound or parameter mask'
    endif
      

    optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
    xopt = optim_f.xopt

    srcstr = {vis_fwdfit_src_structure}
    srcstr.srctype    ='ellipse'

    fitsigmas = {vis_fwdfit_src_structure}
    fitsigmas.srctype ='std.dev'

    srcstr.srcflux = xopt[0]

    ecmsr = REFORM(SQRT(xopt[2]^2 + xopt[3]^2))
    eccen = SQRT(1 - EXP(-2*ecmsr))
    
    srcstr.eccen = eccen

    IF ecmsr GT 0 THEN srcstr.srcpa = reform(ATAN(xopt[3], xopt[2]) * !RADEG) + 90.
    IF srcstr.srcpa lt 0. then srcstr.srcpa += 180.

    srcstr.srcfwhm_min = xopt[1] * (1-eccen^2)^0.25
    srcstr.srcfwhm_max = xopt[1] / (1-eccen^2)^0.25

    srcstr.srcx        = -xopt[5]+vis[0].xyoffset[0]
    srcstr.srcy        = xopt[4]+vis[0].xyoffset[1]


    if keyword_set(uncertainty) then begin

      print, ' '
      print, 'Uncertainty: '
      print, '

      ntry = 20

      trial_results = fltarr(Nvars, ntry)
      iseed=findgen(ntry)+seedstart

      for n=0,ntry-1 do begin
        testerror = RANDOMN(iseed[n], nvis)          ; nvis element vector normally distributed with sigma = 1
        vistest   = visobs + testerror * sigamp
        vistest   = transpose(cmreplicate(vistest, SwarmSize))

        extra = {type: type, $
          visobs: vistest, $
          sigamp: ssigamp, $
          u: vis.u, $
          v: vis.v, $
          n_free: nvis - Nvars, $
          param_opt: param_opt, $
          mapcenter : vis.xyoffset}

        optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
        xopt = optim_f.xopt


        ecmsr = REFORM(SQRT(xopt[2]^2 + xopt[3]^2))
        eccen = SQRT(1 - EXP(-2*ecmsr))

        IF ecmsr GT 0 THEN trial_results[3,n] = reform(ATAN(xopt[3], xopt[2]) * !RADEG) + 180.

        trial_results[0,n]  = xopt[0]
        trial_results[1,n]  = xopt[1] / (1-eccen^2)^0.25
        trial_results[2,n]  = xopt[1] * (1-eccen^2)^0.25
        trial_results[4,n]  = xopt[4]
        trial_results[5,n]  = xopt[5]

      endfor

      fitsigmas.srcflux        = stddev(trial_results[0, *])
      fitsigmas.srcfwhm_max    = stddev(trial_results[1, *])
      fitsigmas.srcfwhm_min    = stddev(trial_results[2, *])
      avsrcpa                  = ATAN(TOTAL(SIN(trial_results[3, *] * !DTOR)), $
                                      TOTAL(COS(trial_results[3, *] * !DTOR))) * !RADEG
      groupedpa                = (810 + avsrcpa - trial_results[3, *]) MOD 180.
      fitsigmas.srcpa          = STDDEV(groupedpa)
      fitsigmas.srcx           = stddev(trial_results[5,*])
      fitsigmas.srcy           = stddev(trial_results[4,*])
    endif

  endif

  if type EQ 'multi' then begin
    
    n_free = nvis-8.
    
    if (n_elements(param_opt) ne 8) or (n_elements(lb) ne 8) or (n_elements(ub) ne 8) then begin
      UNDEFINE, lb, ub, param_opt
      message, 'Wrong number of elements of lower bound, upper bound or parameter mask'
    endif

    Nruns = 20
    xx_opt = []
    f = fltarr(Nruns)

    for i = 0,Nruns-1 do begin
      optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
      f[i]    = optim_f.fopt
      xx_opt  = [[xx_opt],optim_f.xopt]
    endfor

    dummy = min(f,location)
    xopt  = xx_opt(location,*)

    srcstr = {vis_fwdfit_src_structure}
    srcstr.srctype ='ellipse'
    srcstr = VIS_FWDFIT_SRC_BIFURCATE(srcstr)

    fitsigmas = {vis_fwdfit_src_structure}
    fitsigmas.srctype ='std.dev'
    fitsigmas = VIS_FWDFIT_SRC_BIFURCATE(fitsigmas)

    srcstr[0].srcflux       = xopt[1]
    srcstr[0].srcfwhm_max   = xopt[0]
    srcstr[0].srcfwhm_min   = xopt[0]
    srcstr[0].srcx          = -xopt[5]+vis[0].xyoffset[0]
    srcstr[0].srcy          = xopt[4]+vis[0].xyoffset[1]

    srcstr[1].srcflux       = xopt[3]
    srcstr[1].srcfwhm_max   = xopt[2]
    srcstr[1].srcfwhm_min   = xopt[2]
    srcstr[1].srcx          = -xopt[7]+vis[0].xyoffset[0]
    srcstr[1].srcy          = xopt[6]+vis[0].xyoffset[1]


    if keyword_set(uncertainty) then begin

      print, ' '
      print, 'Uncertainty: '
      print, '

      ntry = 20
      trial_results = fltarr(Nvars, ntry)
      ;iseed=findgen(ntry)+seedstart

      for n=0,ntry-1 do begin
        nn = n
        testerror  = RANDOMN(nn+seedstart, nvis)
        vistest    = visobs + testerror * sigamp
        vistest    = transpose(cmreplicate(vistest, SwarmSize))

        extra = {type: type, $
          visobs: vistest, $
          sigamp: ssigamp, $
          u: vis.u, $
          v: vis.v, $
          n_free: nvis - Nvars, $
          param_opt: param_opt, $
          mapcenter : vis.xyoffset}

        xx_opt = []
        f = fltarr(Nruns)

        for i = 0,Nruns-1 do begin
          optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
          f[i]    = optim_f.fopt
          xx_opt  = [[xx_opt],optim_f.xopt]
        endfor

        dummy = min(f,location)
        xopt  = xx_opt(location,*)


        trial_results[0, n] = xopt[1]                   ;flux1
        trial_results[1, n] = xopt[0]                   ;FWHM1
        trial_results[2, n] = xopt[3]                   ;flux2
        trial_results[3, n] = xopt[2]                   ;FWHM2

        trial_results[4, n] = xopt[4]                   ;x1
        trial_results[5, n] = xopt[5]                   ;y1

        trial_results[6, n] = xopt[6]                   ;x2
        trial_results[7, n] = xopt[7]                   ;y2
     
     endfor


      fitsigmas[0].srcflux     = stddev(trial_results[0,*])
      fitsigmas[0].srcfwhm_max = stddev(trial_results[1,*])
      fitsigmas[0].srcfwhm_min = stddev(trial_results[1,*])
      fitsigmas[0].srcx        = stddev(trial_results[5,*])
      fitsigmas[0].srcy        = stddev(trial_results[4,*])

      fitsigmas[1].srcflux     = stddev(trial_results[2,*])
      fitsigmas[1].srcfwhm_max = stddev(trial_results[3,*])
      fitsigmas[1].srcfwhm_min = stddev(trial_results[3,*])
      fitsigmas[1].srcx        = stddev(trial_results[7,*])
      fitsigmas[1].srcy        = stddev(trial_results[6,*])

    endif

  endif



  if type eq 'loop' then begin
    
    n_free = nvis-7.

    if (n_elements(param_opt) ne 7) or (n_elements(lb) ne 7) or (n_elements(ub) ne 7) then begin
      UNDEFINE, lb, ub, param_opt
      message, 'Wrong number of elements of lower bound, upper bound or parameter mask'
    endif

    Nruns = 10
    xx_opt = []
    f = fltarr(Nruns)

    for i = 0,Nruns-1 do begin
      optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
      f[i]    = optim_f.fopt
      xx_opt  = [[xx_opt],optim_f.xopt]
    endfor

    dummy = min(f,location)
    xopt  = xx_opt(location,*)

    srcstr = {vis_fwdfit_src_structure}
    srcstr.srctype = 'loop'

    fitsigmas = {vis_fwdfit_src_structure}
    fitsigmas.srctype = 'std.dev'

    srcstr.srcflux = xopt[0]

    ecmsr = REFORM(SQRT(xopt[2]^2 + xopt[3]^2))
    eccen = SQRT(1 - EXP(-2*ecmsr))

    srcstr.eccen = eccen

    IF eccen GT 0.001 THEN srcstr.srcpa = reform(ATAN(xopt[3], xopt[2]) * !RADEG) + 90.
    IF srcstr.srcpa lt 0. then srcstr.srcpa += 180.

    srcstr.srcfwhm_min = xopt[1] * (1-eccen^2)^0.25
    srcstr.srcfwhm_max = xopt[1] / (1-eccen^2)^0.25

    srcstr.srcx        = -xopt[5]+vis[0].xyoffset[0]
    srcstr.srcy        = xopt[4]+vis[0].xyoffset[1]

    srcstr.loop_angle  = xopt[6]

    if keyword_set(uncertainty) then begin

      print, ' '
      print, 'Uncertainty: '
      print, '

      ntry = 20
      trial_results = fltarr(Nvars, ntry)

      for n=0,ntry-1 do begin
        nn=n
        testerror = RANDOMN(nn+seedstart, nvis)          ; nvis element vector normally distributed with sigma = 1
        vistest   = visobs + testerror * sigamp
        vistest   = transpose(cmreplicate(vistest, SwarmSize))

        extra = {type: type, $
          visobs: vistest, $
          sigamp: ssigamp, $
          u: vis.u, $
          v: vis.v, $
          n_free: nvis - Nvars, $
          param_opt: param_opt, $
          mapcenter : vis.xyoffset}

        xx_opt = []
        f = fltarr(Nruns)

        for i = 0,Nruns-1 do begin
          optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
          f[i]    = optim_f.fopt
          xx_opt  = [[xx_opt],optim_f.xopt]
        endfor

        dummy = min(f,location)
        xopt  = xx_opt(location,*)
      
        ecmsr = REFORM(SQRT(xopt[2]^2 + xopt[3]^2))
        eccen = SQRT(1 - EXP(-2*ecmsr))

        IF ecmsr GT 0 THEN trial_results[3,n] = reform(ATAN(xopt[2], -xopt[3]) * !RADEG) + 180.

        trial_results[0,n]  = xopt[0]
        trial_results[1,n]  = xopt[1] / (1-eccen^2)^0.25
        trial_results[2,n]  = xopt[1] * (1-eccen^2)^0.25
        trial_results[4,n]  = xopt[4]
        trial_results[5,n]  = xopt[5]
        trial_results[6,n]  = xopt[6]

      endfor

      fitsigmas.srcflux        = stddev(trial_results[0, *])
      fitsigmas.srcfwhm_max    = stddev(trial_results[1, *])
      fitsigmas.srcfwhm_min    = stddev(trial_results[2, *])
      avsrcpa                  = ATAN(TOTAL(SIN(trial_results[3, *] * !DTOR)), $
                                      TOTAL(COS(trial_results[3, *] * !DTOR))) * !RADEG
      groupedpa                = (810 + avsrcpa - trial_results[3, *]) MOD 180.
      fitsigmas.srcpa          = STDDEV(groupedpa)
      fitsigmas.srcx           = stddev(trial_results[5,*])
      fitsigmas.srcy           = stddev(trial_results[4,*])
      fitsigmas.loop_angle     = stddev(trial_results[6,*])
    endif


  endif

  if ~keyword_set(silent) then begin

    PRINT
    PRINT, 'COMPONENT    TYPE          FLUX       FWHM MAX    FWHM MIN      Angle     X loc      Y loc      Loop FWHM'
    PRINT, '                         cts/s/keV     arcsec      arcsec        deg      arcsec     arcsec        deg   '
    PRINT
    nsrc = N_ELEMENTS(srcstr)
    FOR n = 0, nsrc-1 DO BEGIN
      temp        = [ srcstr[n].srcflux,srcstr[n].srcfwhm_max,  srcstr[n].srcfwhm_min, $ 
                      srcstr[n].srcpa, srcstr[n].srcx, srcstr[n].srcy, srcstr[n].loop_angle]
      PRINT, n+1, srcstr[n].srctype, temp, FORMAT="(I5, A13, F13.2, 1F13.1, F12.1, 2F11.1, F11.1, 2F12.1)"
      
      temp        = [ fitsigmas[n].srcflux,fitsigmas[n].srcfwhm_max, fitsigmas[n].srcfwhm_min, $
                      fitsigmas[n].srcpa, fitsigmas[n].srcx, fitsigmas[n].srcy, fitsigmas[n].loop_angle]
      PRINT, ' ', '(std)', temp, FORMAT="(A7, A11, F13.2, 1F13.1, F12.1, 2F11.1, F11.1, 2F12.1)"
      PRINT, ' '
    endfor

  endif

  UNDEFINE, lb, ub, param_opt
    ;return, {srcstr: srcstr, fitsigmas: fitsigmas}
  
;  param_out = { srcout: srcstr, sigma: fitsigmas, $
;  ;niter: niter, $
;    ;redchi2: redchisq, nfree: nfree, qflag: qflag, 
;    vf_vis_window: fcheck( vf_vis_window, -1) }

  fwdfit_pso_map = vis_FWDFIT_PSO_SOURCE2MAP(srcstr, type=type, pixel=pixel, imsize=imsize, xyoffset=vis[0].xyoffset)

  this_estring=strtrim(fix(vis[0].energy_range[0]),2)+'-'+strtrim(fix(vis[0].energy_range[1]),2)+' keV'
  fwdfit_pso_map.ID = 'STIX VIS_PSO '+this_estring+': '
  fwdfit_pso_map.dx = pixel[0]
  fwdfit_pso_map.dy = pixel[1]
  fwdfit_pso_map.xc = vis[0].xyoffset[0]
  fwdfit_pso_map.yc = vis[0].xyoffset[1]
  this_time_range   = stx_time2any(vis[0].time_range,/vms)
  fwdfit_pso_map.time = anytim((anytim(this_time_range[1])+anytim(this_time_range[0]))/2.,/vms)
  fwdfit_pso_map.DUR  = anytim(this_time_range[1])-anytim(this_time_range[0])
  ;eventually fill in radial distance etc
  add_prop,fwdfit_pso_map,rsun=0.
  add_prop,fwdfit_pso_map,B0=0.
  add_prop,fwdfit_pso_map,L0=0.


srcstrout_pso    = srcstr
fitsigmasout_pso = fitsigmas


if ~no_plot_fit then begin

  visobs   = vis.obsvis
  phaseobs = atan(imaginary(visobs), float(visobs)) * !radeg
  visobs2  = [float(visobs),imaginary(visobs)]

  srcstrout0 = srcstrout_pso
  visobsmap  = vis_fwdfit_pso_vis_pred(srcstrout0, vis, type)  

  phaseobsmap =  atan(visobsmap[n_elements(vis):2*n_elements(vis)-1], visobsmap[0:n_elements(vis)-1]) * !radeg
  ampobsmap   = sqrt(visobsmap[0:n_elements(vis)-1]^2 + visobsmap[n_elements(vis):2*n_elements(vis)-1]^2)

  sigamp   = vis.sigamp
  sigamp2  = [vis.sigamp, vis.sigamp]
  sigphase = sigamp / abs(visobs)  * !radeg


  xx = (findgen(30))/3. + 1.2
  xx = xx[6:29]
  
  nfree = n_elements(vis)-nvars
  chi2  = total(abs(visobsmap - visobs2)^2./sigamp2^2.)/nfree

  charsize = 1.5
  leg_size = 1.5
  thick = 1.8
  symsize = 1.8

  units_phase = 'degrees'
  units_amp   = 'counts s!U-1!n keV!U-1!n'
  xtitle      = 'Detector label'

  window, 1, xsize=1200, ysize=500

  loadct, 5
  linecolors, /quiet

  set_viewport,0.1, 0.45, 0.1, 0.85

  plot, xx, phaseobs, /nodata, xrange=[1.,11.], /xst, xtickinterval=1, xminor=-1, $
    title='VISIBILITY PHASE FIT PSO', $
    xtitle=xtitle, ytitle=units_phase, yrange=yrange, charsize=charsize, thick=thick, /noe

  ; draw vertical dotted lines at each detector boundary
  for i=1,10 do oplot, i+[0,0], !y.crange, linestyle=1


  errplot, xx, (phaseobs - sigphase > !y.crange[0]), (phaseobs + sigphase < !y.crange[1]), $
    width=0, thick=thick, COLOR=7
  oplot, xx, phaseobs, psym=7, thick=thick, symsize=symsize
  oplot, xx, phaseobsmap, psym=4, col=2, thick=thick, symsize=symsize


  leg_text = ['Observed', 'Error on Observed', 'From Image']
  leg_color = [255, 7,2]
  leg_style = [0, 0, 0]
  leg_sym = [7, -3, 4]
  ssw_legend, leg_text, psym=leg_sym, color=leg_color, linest=leg_style, box=0, charsize=leg_size, thick=thick, /left


  set_viewport,0.5, 0.95, 0.1, 0.85

  plot, xx, abs(visobs), /nodata, xrange=[1.,11.], /xst, xtickinterval=1, xminor=-1, $
    title='VISIBILITY AMPLITUDE FIT PSO - CHI2: ' + trim(chi2, '(f12.2)'), $
    xtitle=xtitle, ytitle=units, yrange=yrange, charsize=charsize, thick=thick, /noe

  ; draw vertical dotted lines at each detector boundary
  for i=1,10 do oplot, i+[0,0], !y.crange, linestyle=1

  errplot, xx, (abs(visobs) - sigamp > !y.crange[0]), (abs(visobs) + sigamp < !y.crange[1]), $
    width=0, thick=thick, COLOR=7
  oplot, xx, abs(visobs), psym=7, thick=thick, symsize=symsize
  oplot, xx, ampobsmap, psym=4, col=2, thick=thick, symsize=symsize


  leg_text = ['Observed', 'Error on Observed', 'From Image']
  leg_color = [255, 7,2]
  leg_style = [0, 0, 0]
  leg_sym = [7, -3, 4]
  ssw_legend, leg_text, psym=leg_sym, color=leg_color, linest=leg_style, box=0, charsize=leg_size, thick=thick, /left

  !p.position = [0, 0, 0, 0]

endif


return, fwdfit_pso_map  

end