FUNCTION vis_fwdfit_nov_2021,vis0,imsize,pixel,shape=shape,SRCOUT=srcstrout,FITSTDDEV=fitstddev,QFLAG=qflag,REDCHISQ=redchisq,NITER=niter, NFREE=nfree, $
  vf_vis_window=vf_vis_window,xyoffset=xyoffset,image_out = image_out,silent=silent,_REF_EXTRA=extra, no_plot_fit=no_plot_fit

  ; wrapper around VIS_FWDFIT
  ; output map structure has north up
  ;
  ; 8-Nov-2021: Paolo, Anna, Emma: first version

  ;Look for stix visibility structure and extract the vis bag if necessary in hsi-like format
  default, shape, 'circle'
  default, no_plot_fit, 0
  
  vis = vis0
  
  ind = where(vis.v lt 0.)
  vis[ind].u = -vis[ind].u
  vis[ind].v = -vis[ind].v
  vis[ind].obsvis = conj(vis[ind].obsvis)
  
  vis.sigamp = sqrt(vis.sigamp^2 - 0.05^2 * abs(vis.obsvis)^2) ; VIS_FWDFIT add 5% systematic error inside the procedure
  vis.xyoffset *= 0.
  
  case shape of
    
    'circle': begin
              n_free = n_elements(vis)-4
              vis_fwdfit, vis, SRCOUT=srcstrout, FITSTDDEV=fitstddev, QFLAG=qflag, REDCHISQ=redchisq, $
                          NITER=niter, NFREE=nfree, vf_vis_window=vf_vis_window, circle=1, /noplotfit, /quiet
              end
    'ellipse': begin
               n_free = n_elements(vis)-6
               vis_fwdfit, vis, SRCOUT=srcstrout, FITSTDDEV=fitstddev, QFLAG=qflag, REDCHISQ=redchisq, $
                          NITER=niter, NFREE=nfree, vf_vis_window=vf_vis_window, circle=0,/noplotfit, /quiet
              end
     'multi': begin
              n_free = n_elements(vis)-8
              vis_fwdfit, vis, SRCOUT=srcstrout, FITSTDDEV=fitstddev, QFLAG=qflag, REDCHISQ=redchisq, $
                          NITER=niter, NFREE=nfree, vf_vis_window=vf_vis_window, multi=1, /noplotfit, /quiet
              end
      'loop': begin
              n_free = n_elements(vis)-7
              vis_fwdfit, vis, SRCOUT=srcstrout, FITSTDDEV=fitstddev, QFLAG=qflag, REDCHISQ=redchisq, $
                          NITER=niter, NFREE=nfree, vf_vis_window=vf_vis_window, loop=1, /noplotfit, /quiet
              end
  endcase
  
  
  
  
  if ~keyword_set(silent) then begin
  PRINT
  PRINT, 'COMPONENT  PROFILE      FLUX       X(+ve W)  Y(+ve N)   FWHM MAX  FWHM MIN  PosnAngle   Loop_FWHM '
  PRINT, '                      counts/s      arcsec    arcsec     arcsec    arcsec    degrees     degrees  '
  PRINT
  nsrc = N_ELEMENTS(srcstrout)
  FOR n = 0, nsrc-1 DO BEGIN
    
    xx = srcstrout[n].srcx
    srcstrout[n].srcx = -srcstrout[n].srcy + vis0[0].xyoffset[0]
    srcstrout[n].srcy = xx + vis0[0].xyoffset[1]
   ;srcstrout[n].srcpa -= 90.
    srcstrout[n].srcpa += 90.     ;anna
    
    xx = fitstddev[n].srcx
    fitstddev[n].srcx = fitstddev[n].srcy 
    fitstddev[n].srcy = xx 
    
    eccen = srcstrout[n].eccen
    fwhm = srcstrout[n].srcfwhm
    fwhmminor   = fwhm * (1-eccen^2)^0.25
    fwhmmajor   = fwhm / (1-eccen^2)^0.25
    
    std_fwhmmajor = sqrt((fitstddev[n].srcfwhm / (1-eccen^2)^0.25)^2 + (1./2. * fwhm * (1-eccen^2)^(-5./4.) * eccen * fitstddev[n].eccen)^2)
    std_fwhminor  = sqrt((fitstddev[n].srcfwhm * (1-eccen^2)^0.25)^2 + (1./2. * fwhm * (1-eccen^2)^(-3./4.) * eccen * fitstddev[n].eccen)^2)
    
    temp        = [ srcstrout[n].srcflux,  srcstrout[n].srcx,  srcstrout[n].srcy,  $
      fwhmmajor,  fwhmminor, srcstrout[n].srcpa, srcstrout[n].loop_angle]
      
    temp_std    = [ fitstddev[n].srcflux,  fitstddev[n].srcx,  fitstddev[n].srcy,  $
      std_fwhmmajor,  std_fwhminor, fitstddev[n].srcpa, fitstddev[n].loop_angle]
    PRINT, n+1, srcstrout[n].srctype, temp, FORMAT="(I5, A13, F12.2, 4F10.2, 2F12.1, F12.3, F12.1)"
    PRINT, n+1, 'std', temp_std, FORMAT="(I5, A13, F12.2, 4F10.2, 2F12.1, F12.3, F12.1)"
 ENDFOR 
 
 endif
  
  
  param_out = { srcout: srcstrout, sigma: fitstddev, niter: niter, $
    redchi2: redchisq, nfree: nfree, qflag: qflag, vf_vis_window: fcheck( vf_vis_window, -1) }

  vis_source2map, srcstrout, vis0.xyoffset, image_out, pixel = pixel[0], mapsize = imsize[0]
  
  vis_fwdfit_map = make_map(image_out)


  this_estring=strtrim(fix(vis0[0].energy_range[0]),2)+'-'+strtrim(fix(vis0[0].energy_range[1]),2)+' keV'
  vis_fwdfit_map.ID = 'STIX VIS_FWDFIT '+this_estring+': '
  vis_fwdfit_map.dx = pixel[0]
  vis_fwdfit_map.dy = pixel[1]
  vis_fwdfit_map.xc = vis0[0].xyoffset[0]
  vis_fwdfit_map.yc = vis0[0].xyoffset[1]
  this_time_range=stx_time2any(vis0[0].time_range,/vms)
  vis_fwdfit_map.time = anytim((anytim(this_time_range[1])+anytim(this_time_range[0]))/2.,/vms)
  vis_fwdfit_map.DUR = anytim(this_time_range[1])-anytim(this_time_range[0])
  ;eventually fill in radial distance etc
  add_prop,vis_fwdfit_map,rsun=0.
  add_prop,vis_fwdfit_map,B0=0.
  add_prop,vis_fwdfit_map,L0=0.
  
  if ~no_plot_fit then begin
    
  visobs = vis0.obsvis
  phaseobs = atan(imaginary(visobs), float(visobs)) * !radeg
  visobs2=[float(visobs),imaginary(visobs)]
  
  srcstrout0 = srcstrout
  visobsmap = vis_fwdfit_vis_pred(srcstrout0, vis0)
  phaseobsmap =  atan(visobsmap[n_elements(vis0):2*n_elements(vis0)-1], visobsmap[0:n_elements(vis0)-1]) * !radeg
  ampobsmap = sqrt(visobsmap[0:n_elements(vis0)-1]^2 + visobsmap[n_elements(vis0):2*n_elements(vis0)-1]^2)
  
  sigamp = vis0.sigamp
  sigamp2= [sigamp, sigamp]
  sigphase = sigamp / abs(visobs)  * !radeg


  xx = (findgen(30))/3. + 1.2
  xx = xx[6:29]

  chi2 = total(abs(visobsmap - visobs2)^2./sigamp2^2.)/n_free


  charsize = 1.5
  leg_size = 1.5
  thick = 1.8
  symsize = 1.8

  units_phase = 'degrees'
  units_amp   = 'counts s!U-1!n keV!U-1!n'
  xtitle      = 'Detector label'

  window, 0, xsize=1200, ysize=500

  loadct, 5
  linecolors, /quiet

  set_viewport,0.1, 0.45, 0.1, 0.85

  plot, xx, phaseobs, /nodata, xrange=[1.,11.], /xst, xtickinterval=1, xminor=-1, $
    title='VISIBILITY PHASE FIT FWDFIT', $
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
    title='VISIBILITY AMPLITUDE FIT FWDFIT - CHI2: ' + trim(chi2, '(f12.2)'), $
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

  return,vis_fwdfit_map

END