;+
;
; NAME:
;   stx_vis_fwdfit_pso
;
; PURPOSE:
;   Wrapper around vis_fwdfit_pso
;
; INPUTS:
;   configuration: array containing parametric shapes chosen for the forward fitting method (one for each source component)
;                    - 'circle' : Gaussian circular source
;                    - 'ellipse': Gaussian elliptical source
;                    - 'loop'   : single curved elliptical gaussian
;
;   vis         : struct containing  the observed visibility values
;     vis.obsvis: array containing the values of the observed visibilities
;     vis.sigamp: array containing the values of the errors on the observed visibility amplitudes
;     vis.u     : u coordinates of the sampling frequencies
;     vis.v     : v coordinates of the sampling frequencies
;
; KEYWORDS:
;   SRCIN       : struct containing for each source the parameters to optimize and those fixed, upper and lower bound of the variables.
;                 to create the structure srcin:
;                     srcin = VIS_FWDFIT_PSO_MULTIPLE_SRC_CREATE(vis, configuration) 
;                                  
;                 If not entered, default values are used (see:
;                                                             - vis_fwdfit_pso_circle_struct_define for the circle                                                         
;                                                             - vis_fwdfit_pso_ellipse_struct_define for the ellipse
;                                                             - vis_fwdfit_pso_loop_struct_define for the loop)
;                                                             
;   N_BIRDS     : number of particles used in PSO 
;                 (default is 100)
;   TOLERANCE   : tolerance for the stopping criterion 
;                 (default is 1e-6)
;   MAXITER     : maximum number of iterations of PSO
;                 (default is the product between of the numbers of parameters and the number of particles)
;   UNCERTAINTY : set to 1 for the computation of the parameters uncertainty (confidence strip approach)
;                 (default is 0) 
;   IMSIZE      : array containing the size (number of pixels) of the image to reconstruct
;                 (default is [128., 128.])
;   PIXEL       : array containing the pixel size (in arcsec) of the image to reconstruct
;                 (default is [1., 1.])
;   SILENT      : set to 1 for avoiding the print of the retrieved parameters
;   SEEDSTART   : costant used to initialize the random perturbation of visibilities when uncertainty is computed
;                 (default is fix(randomu(seed) * 100))
;
;
; OUTPUTS:
;   fit parameters and uncertaintly  (srcstr and fitsigams)
;   reduced chi^2 (redchisq)
;   image map (output map structure has north up)


function stx_vis_fwdfit_pso, configuration, vis, $
                             srcin=srcin, $
                             n_birds = n_birds, tolerance = tolerance, maxiter = maxiter, $
                             uncertainty = uncertainty, $
                             imsize=imsize, pixel=pixel, $
                             silent = silent, $
                             srcstr = srcstr,fitsigmas =fitsigmas, redchisq = redchisq, $
                             seedstart = seedstart


default, imsize, [128,128]
default, pixel, [1.,1.]
default, n_birds, 100

; stix view
this_vis = vis
this_vis.xyoffset[0] = vis.xyoffset[1]
this_vis.xyoffset[1] = - vis.xyoffset[0]

phi= max(abs(vis.obsvis)) ;estimate_flux

loc_circle  = where(configuration eq 'circle', n_circle)>0 
loc_ellipse = where(configuration eq 'ellipse', n_ellipse)>0
loc_loop    = where(configuration eq 'loop', n_loop)>0

n_sources = n_elements(configuration)

param_opt=[]
lower_bound=[]
upper_bound=[]

; create param_opt, lower_bound and upper_bound
if ~keyword_set(SRCIN) then begin

    if n_circle gt 0. then begin
      param_opt_circle   = ['fit', 'fit','fit','fit']
      lower_bound_circle = [0.1*phi, -100., -100., 1.]
      upper_bound_circle = [1.5*phi, 100., 100., 100.]

      param_opt_circle   = cmreplicate(param_opt_circle, n_circle)
      lower_bound_circle = cmreplicate(lower_bound_circle, n_circle)
      upper_bound_circle = cmreplicate(upper_bound_circle, n_circle)

      param_opt_circle   = reform(param_opt_circle, [4. *n_circle, 1])
      lower_bound_circle = reform(lower_bound_circle, [4. *n_circle, 1])
      upper_bound_circle = reform(upper_bound_circle, [4. *n_circle, 1])
      
      param_opt   = [param_opt, [param_opt_circle]]
      lower_bound = [lower_bound, [lower_bound_circle]]
      upper_bound = [upper_bound, [upper_bound_circle]]

    endif
    
    if n_ellipse gt 0. then begin
      
      param_opt_ellipse   = ['fit', 'fit', 'fit', 'fit', 'fit', 'fit']
      lower_bound_ellipse = [0.1*phi,  -100, -100., 1., -5., 0.]
      upper_bound_ellipse = [1.5*phi, 100., 100., 100., 5., 1.]

      param_opt_ellipse   = cmreplicate(param_opt_ellipse, n_ellipse)
      lower_bound_ellipse = cmreplicate(lower_bound_ellipse, n_ellipse)
      upper_bound_ellipse = cmreplicate(upper_bound_ellipse, n_ellipse)
  
      param_opt_ellipse   = reform(param_opt_ellipse, [6. *n_ellipse, 1])
      lower_bound_ellipse = reform(lower_bound_ellipse, [6. *n_ellipse, 1])
      upper_bound_ellipse = reform(upper_bound_ellipse, [6. *n_ellipse, 1])
      
      param_opt   = [param_opt, [param_opt_ellipse]]
      lower_bound = [lower_bound, [lower_bound_ellipse]]
      upper_bound = [upper_bound, [upper_bound_ellipse]]

    endif
    
    if n_loop gt 0. then begin
      
      param_opt_loop   = ['fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit']
      lower_bound_loop = [0.1*phi,  -100, -100., 1., -5., 0., -180.]
      upper_bound_loop = [1.5*phi, 100., 100., 100., 5., 1., 180.]

      param_opt_loop  = cmreplicate(param_opt_loop, n_loop)
      lower_bound_loop = cmreplicate(lower_bound_loop, n_loop)
      upper_bound_loop = cmreplicate(upper_bound_loop, n_loop)

      param_opt_loop   = reform(param_opt_loop, [7. *n_loop, 1])
      lower_bound_loop = reform(lower_bound_loop, [7. *n_loop, 1])
      upper_bound_loop = reform(upper_bound_loop, [7. *n_loop, 1])
      
      param_opt   = [param_opt, [param_opt_loop]]        
      lower_bound = [lower_bound, [lower_bound_loop]]
      upper_bound = [upper_bound, [upper_bound_loop]]

    endif

endif else begin
  
  if n_circle gt 0 then begin
    
    for j=0, n_circle-1 do begin
      
      flag=1
      Catch, theError
      IF theError NE 0 THEN BEGIN
        Catch, /Cancel
        x_stix_c:
        flag=0
        this_x_c = 'fit'
      ENDIF
      ; Set up file I/O error handling.
      ON_IOError, x_stix_c
      ; Cause type conversion error.
      if flag then this_x_c = string(double(srcin.circle[j].param_opt.param_y)-58.2)
      
      flag=1
      Catch, theError
      IF theError NE 0 THEN BEGIN
        Catch, /Cancel
        y_stix_c:
        flag=0
        this_y_c = 'fit'
      ENDIF
      ; Set up file I/O error handling.
      ON_IOError, y_stix_c
      ; Cause type conversion error.
      if flag then this_y_c = string( - double(srcin.circle[j].param_opt.param_x) + 26.1)
      
      param_opt=[param_opt, string(srcin.circle[j].param_opt.param_flux), $
                            this_x_c, $
                            this_y_c, $
                            string(srcin.circle[j].param_opt.param_fwhm)]
                            
      lower_bound = [lower_bound, srcin.circle[j].lower_bound.l_b_flux, $
                                  srcin.circle[j].lower_bound.l_b_y, $
                                  - srcin.circle[j].upper_bound.u_b_x, $
                                  srcin.circle[j].lower_bound.l_b_fwhm]
                                  
      upper_bound = [upper_bound, srcin.circle[j].upper_bound.u_b_flux, $
                                  srcin.circle[j].upper_bound.u_b_y, $
                                  - srcin.circle[j].lower_bound.l_b_x, $
                                  srcin.circle[j].upper_bound.u_b_fwhm]

    endfor
  endif
  
  if n_ellipse gt 0 then begin
    for j=0, n_ellipse-1 do begin
           
      flag=1
      Catch, theError
      IF theError NE 0 THEN BEGIN
        Catch, /Cancel
        x_stix_e:
        flag=0
        this_x_e = 'fit'
      ENDIF
      ; Set up file I/O error handling.
      ON_IOError, x_stix_e
      ; Cause type conversion error.
      if flag then this_x_e = string(double(srcin.ellipse[j].param_opt.param_y) - 58.2)

      flag=1
      Catch, theError
      IF theError NE 0 THEN BEGIN
        Catch, /Cancel
        y_stix_e:
        flag=0
        this_y_e = 'fit'
      ENDIF
      ; Set up file I/O error handling.
      ON_IOError, y_stix_e
      ; Cause type conversion error.
      if flag then this_y_e = string( - double(srcin.ellipse[j].param_opt.param_x) + 26.1)

      param_opt   = [param_opt, string(srcin.ellipse[j].param_opt.param_flux), $
                                this_x_e, $
                                this_y_e, $
                                string(srcin.ellipse[j].param_opt.param_fwhm_max), $
                                string(srcin.ellipse[j].param_opt.param_fwhm_min), $
                                string(srcin.ellipse[j].param_opt.param_alpha)]
                                
      lower_bound = [lower_bound, srcin.ellipse[j].lower_bound.l_b_flux, $
                                  srcin.ellipse[j].lower_bound.l_b_y, $
                                  - srcin.ellipse[j].upper_bound.u_b_x, $
                                  srcin.ellipse[j].lower_bound.l_b_fwhm, $
                                  srcin.ellipse[j].lower_bound.l_b_eccos, $
                                  srcin.ellipse[j].lower_bound.l_b_ecsin]
                                  
      upper_bound = [upper_bound, srcin.ellipse[j].upper_bound.u_b_flux, $
                                  srcin.ellipse[j].upper_bound.u_b_y, $
                                  - srcin.ellipse[j].lower_bound.l_b_x, $
                                  srcin.ellipse[j].upper_bound.u_b_fwhm, $
                                  srcin.ellipse[j].upper_bound.u_b_eccos, $
                                  srcin.ellipse[j].upper_bound.u_b_ecsin]
              
    endfor
  endif
  
  if n_loop gt 0 then begin
    for j=0, n_loop-1 do begin
      
      flag=1
      Catch, theError
      IF theError NE 0 THEN BEGIN
        Catch, /Cancel
        x_stix_l:
        flag=0
        this_x_l = 'fit'
      ENDIF
      ; Set up file I/O error handling.
      ON_IOError, x_stix_l
      ; Cause type conversion error.
      if flag then this_x_l = string(double(srcin.loop[j].param_opt.param_y) - 58.2)

      flag=1
      Catch, theError
      IF theError NE 0 THEN BEGIN
        Catch, /Cancel
        y_stix_l:
        flag=0
        this_y_l = 'fit'
      ENDIF
      ; Set up file I/O error handling.
      ON_IOError, y_stix_l
      ; Cause type conversion error.
      if flag then this_y_l = string( - double(srcin.loop[j].param_opt.param_x) + 26.1)


      param_opt   = [param_opt, string(srcin.loop[j].param_opt.param_flux), $
                                this_x_l, $
                                this_y_l, $
                                string(srcin.loop[j].param_opt.param_fwhm_max), $
                                string(srcin.loop[j].param_opt.param_fwhm_min), $
                                string(srcin.loop[j].param_opt.param_alpha), $
                                string(srcin.loop[j].param_opt.param_loopangle)]

      lower_bound = [lower_bound, srcin.loop[j].lower_bound.l_b_flux, $
                                  srcin.loop[j].lower_bound.l_b_y, $
                                  - srcin.loop[j].upper_bound.u_b_x, $
                                  srcin.loop[j].lower_bound.l_b_fwhm, $
                                  srcin.loop[j].lower_bound.l_b_eccos, $
                                  srcin.loop[j].lower_bound.l_b_ecsin, $
                                  srcin.loop[j].lower_bound.l_b_loopangle]

      upper_bound = [upper_bound, srcin.loop[j].upper_bound.u_b_flux, $
                                  srcin.loop[j].upper_bound.u_b_y, $
                                  - srcin.loop[j].lower_bound.l_b_x, $
                                  srcin.loop[j].upper_bound.u_b_fwhm, $
                                  srcin.loop[j].upper_bound.u_b_eccos, $
                                  srcin.loop[j].upper_bound.u_b_ecsin, $
                                  srcin.loop[j].upper_bound.u_b_loopangle]

    endfor
  endif
endelse

  
param_out = vis_fwdfit_pso(configuration, this_vis, $
                            lower_bound = lower_bound, upper_bound = upper_bound, param_opt = param_opt, $
                            n_birds = n_birds, tolerance = tolerance, maxiter = maxiter, $
                            uncertainty = uncertainty, $
                            imsize=imsize, pixel=pixel, $
                            silent = silent, $
                            seedstart = seedstart)

srcstr = param_out.srcstr
fitsigmas = param_out.fitsigmas
redchisq = param_out.redchisq
fwdfit_pso_map = make_map(param_out.data)
this_estring=strtrim(fix(vis[0].energy_range[0]),2)+'-'+strtrim(fix(vis[0].energy_range[1]),2)+' keV'
fwdfit_pso_map.ID = 'STIX VIS_PSO '+this_estring+': '
fwdfit_pso_map.dx = pixel[0]
fwdfit_pso_map.dy = pixel[1]

this_time_range = stx_time2any(vis[0].time_range,/vms)

fwdfit_pso_map.xc = vis[0].xyoffset[0]
fwdfit_pso_map.yc = vis[0].xyoffset[1]

fwdfit_pso_map.time = anytim((anytim(this_time_range[1])+anytim(this_time_range[0]))/2.,/vms)
fwdfit_pso_map.DUR  = anytim(this_time_range[1])-anytim(this_time_range[0])

;rotate map to heliocentric view
fwdfit_pso__map = fwdfit_pso_map
fwdfit_pso__map.data = rotate(fwdfit_pso_map.data,1)

;; Mapcenter corrected for Frederic's mean shift values
fwdfit_pso__map.xc = vis[0].xyoffset[0] + 26.1
fwdfit_pso__map.yc = vis[0].xyoffset[1] + 58.2

data = stx_get_l0_b0_rsun_roll_temp(this_time_range[0])

fwdfit_pso__map.roll_angle    = data.ROLL_ANGLE
add_prop,fwdfit_pso__map,rsun = data.RSUN
add_prop,fwdfit_pso__map,B0   = data.B0
add_prop,fwdfit_pso__map,L0   = data.L0

if ~keyword_set(silent) then begin

  PRINT
  PRINT, 'COMPONENT    TYPE          FLUX       FWHM MAX    FWHM MIN      Angle     X loc      Y loc      Loop FWHM'
  PRINT, '                         cts/s/keV     arcsec      arcsec        deg      arcsec     arcsec        deg   '
  PRINT

endif

nsrc = N_ELEMENTS(srcstr)
FOR n = 0, nsrc-1 DO BEGIN

  ; heliocentric view
  x_new = srcstr[n].srcy - this_vis[0].xyoffset[1]
  y_new = srcstr[n].srcx - this_vis[0].xyoffset[0]

  ;; Center of the sources corrected for Frederic's mean shift values
  srcstr[n].srcx        = - x_new + vis[0].xyoffset[0] + 26.1
  srcstr[n].srcy        = y_new + vis[0].xyoffset[1]  + 58.2

  srcstr[n].SRCPA += 90.

  if ~keyword_set(silent) then begin

    temp        = [ srcstr[n].srcflux,srcstr[n].srcfwhm_max,  srcstr[n].srcfwhm_min, $
      srcstr[n].srcpa, $
      ;x_roll, y_roll, $
      srcstr[n].srcx, srcstr[n].srcy, srcstr[n].loop_angle]
    PRINT, n+1, srcstr[n].srctype, temp, FORMAT="(I5, A13, F13.2, 1F13.1, F12.1, 2F11.1, F11.1, 2F12.1)"

    temp        = [ fitsigmas[n].srcflux,fitsigmas[n].srcfwhm_max, fitsigmas[n].srcfwhm_min, $
      fitsigmas[n].srcpa, fitsigmas[n].srcy, fitsigmas[n].srcx, fitsigmas[n].loop_angle]
    PRINT, ' ', '(std)', temp, FORMAT="(A7, A11, F13.2, 1F13.1, F12.1, 2F11.1, F11.1, 2F12.1)"
    PRINT, ' '

  endif

endfor

undefine, param_opt
undefine, upper_bound
undefine, lower_bound
undefine, srcin

return, fwdfit_pso__map

end