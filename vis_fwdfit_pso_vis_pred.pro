function vis_fwdfit_pso_vis_pred, srcstr, vis, type

  COMMON uvdata, u, v, pa, mapcenter

  u = vis.u
  v = vis.v

;;  TWOPI =  2. * !PI
;;  pa = ((ATAN(v, u) + TWOPI) MOD TWOPI) * !RADEG

  nsrc = N_ELEMENTS(srcstr)
  xh = fltarr(nsrc)
  yh = fltarr(nsrc)
  FOR n = 0, nsrc-1 DO BEGIN
    xx = srcstr[n].srcx - vis[0].xyoffset[0]
    yy = srcstr[n].srcy - vis[0].xyoffset[1]
    xh[n]   = yy
    yh[n]   = -xx
;    xh = srcstr[n].srcy - vis[0].xyoffset[0]
;    yh = srcstr[n].srcx + vis[0].xyoffset[1]
    ;pa = srcstr[n].srcpa + 90.
    pa = srcstr[n].srcpa - 90.
  endfor
  
  vispred_re = fltarr(n_elements(u))
  vispred_im = fltarr(n_elements(u))
  
 case type of
   'circle': begin
;              xh = srcstr.srcy - vis[0].xyoffset[0]
;              yh = srcstr.srcx + vis[0].xyoffset[1]
              vispred_re = srcstr.srcflux * exp(-(!pi^2. * srcstr.srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((xh#u)+(yh#v)))
              vispred_im = srcstr.srcflux * exp(-(!pi^2. * srcstr.srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((xh#u)+(yh#v)))
              vis_pred=[transpose(vispred_re), transpose(vispred_im)]
             end
   
   'ellipse': begin
;                xh = srcstr.srcy - vis[0].xyoffset[0]
;                yh = srcstr.srcx + vis[0].xyoffset[1]
;                pa = srcstr.srcpa + 90.
                u1 = cos(pa * !dtor) * u + sin(pa * !dtor) * v
                v1 = -sin(pa * !dtor)* u + cos(pa * !dtor) * v
                vispred_re = srcstr.srcflux * exp(-(!pi^2. / (4.*alog(2.)))#((u1 * srcstr.srcfwhm_max)^2. + (v1 * srcstr.srcfwhm_min)^2.))*cos(2*!pi*((xh#u)+(yh#v)))
                vispred_im = srcstr.srcflux * exp(-(!pi^2. / (4.*alog(2.)))#((u1 * srcstr.srcfwhm_max)^2. + (v1 * srcstr.srcfwhm_min)^2.))*sin(2*!pi*((xh#u)+(yh#v)))
                vis_pred=[transpose(vispred_re), transpose(vispred_im)]
              end
  
  
  'multi': begin
;            x1h = srcstr[0].srcy - vis[0].xyoffset[0]
;            y1h = srcstr[0].srcx + vis[0].xyoffset[1]
;            x2h = srcstr[1].srcy - vis[0].xyoffset[0]
;            y2h = srcstr[1].srcx + vis[0].xyoffset[1]
            
            re_obs1 = srcstr[0].srcflux * exp(-(!pi^2. * srcstr[0].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((xh[0]#u)+(yh[0]#v)))
            im_obs1 = srcstr[0].srcflux * exp(-(!pi^2. * srcstr[0].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((xh[0]#u)+(yh[0]#v)))
    
            re_obs2 = srcstr[1].srcflux * exp(-(!pi^2. * srcstr[1].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((xh[1]#u)+(yh[1]#v)))
            im_obs2 = srcstr[1].srcflux * exp(-(!pi^2. * srcstr[1].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((xh[1]#u)+(yh[1]#v)))
    
            vispred_re = transpose(re_obs1 + re_obs2)
            vispred_im = transpose(im_obs1 + im_obs2)
    
            vis_pred=[vispred_re, vispred_im]
           end
           
  'loop': begin
;            xh = srcstr.srcy - vis[0].xyoffset[0]
;            yh = srcstr.srcx + vis[0].xyoffset[1]
;            pa = srcstr.srcpa + 90.
            fwhm_loop = srcstr.srcfwhm_max*(1.-srcstr.eccen^2)^0.25
            vis_pred = pso_func_makealoop_nov2021( srcstr.srcflux, fwhm_loop , srcstr.eccen, xh, yh, pa, srcstr.loop_angle, vis.u, vis.v)
            ;vis_pred = [[transpose(vis_pred[0,0:n_elements(u)-1])], [transpose(vis_pred[0,n_elements(u):2*n_elements(u)-1])]]
            vis_pred=transpose(vis_pred)
          end         
           
           
  
  
 endcase
  
  
 
;  if srcstr[0].srctype eq 'circle' then begin
;    nsrc = N_ELEMENTS(srcstr)
;    fwhm=fltarr(nsrc)
;    FOR n = 0, nsrc-1 DO BEGIN
;      vispred_re += srcstr[n].srcflux * exp(-(!pi^2. * srcstr[n].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((srcstr[n].srcx#u)+(srcstr[n].srcy#v)))
;      vispred_im += srcstr[n].srcflux  * exp(-(!pi^2. * srcstr[n].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((srcstr[n].srcx#u)+(srcstr[n].srcy#v)))
;    endfor
;    ;re_obs2 = srcstr[1].srcflux * exp(-(!pi^2. * fwhm2^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((srcstr[1].srcx#u)+(srcstr[1].srcy#v)))
;    ;im_obs2 = srcstr[1].srcflux  * exp(-(!pi^2. * fwhm2^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((srcstr[1].srcx#u)+(srcstr[1].srcy#v)))
;
;    ;vispred_re = re_obs1 + re_obs2
;    ;vispred_im = im_obs1 + im_obs2
;    
;    
;    ;vispred_re = srcstr.srcflux * exp(-(!pi^2. * fwhm^2. / (4.*alog(2.)))*(u^2. + v^2.))*cos(2*!pi*((srcstr.srcx*u)+(srcstr.srcy*v)))
;    ;vispred_im = srcstr.srcflux * exp(-(!pi^2. * fwhm^2. / (4.*alog(2.)))*(u^2. + v^2.))*sin(2*!pi*((srcstr.srcx*u)+(srcstr.srcx*v)))
;    vis_pred=[[vispred_re], [vispred_im]]
;  endif

;  if srcstr[0].srctype eq 'ellipse' then begin
;    
;    u1 = cos(srcstr.srcpa * !dtor) * u + sin(srcstr.srcpa * !dtor) * v
;    v1 = -sin(srcstr.srcpa * !dtor) * u + cos(srcstr.srcpa * !dtor) * v
;    vispred_re = srcstr.srcflux * exp(-(!pi^2. / (4.*alog(2.)))*((u1 * srcstr.srcfwhm_max)^2. + (v1 * srcstr.srcfwhm_min)^2.))*cos(2*!pi*((srcstr.srcx*u)+(srcstr.srcy*v)))
;    vispred_im = srcstr.srcflux * exp(-(!pi^2. / (4.*alog(2.)))*((u1 * srcstr.srcfwhm_max)^2. + (v1 * srcstr.srcfwhm_min)^2.))*sin(2*!pi*((srcstr.srcx*u)+(srcstr.srcy*v)))
;    vis_pred=[[vispred_re], [vispred_im]]
;    
; endif
 
;   if srcstr[0].srctype eq 'multi' then begin
;    fwhm1 = srcstr[0].srcfwhm_max
;    fwhm2 = srcstr[1].srcfwhm_max
;    re_obs1 = srcstr[0].srcflux * exp(-(!pi^2. * fwhm1^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((srcstr[0].srcx#u)+(srcstr[0].srcy#v)))
;    im_obs1 = srcstr[0].srcflux  * exp(-(!pi^2. * fwhm1^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((srcstr[0].srcx#u)+(srcstr[0].srcy#v)))
;
;    re_obs2 = srcstr[1].srcflux * exp(-(!pi^2. * fwhm2^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((srcstr[1].srcx#u)+(srcstr[1].srcy#v)))
;    im_obs2 = srcstr[1].srcflux  * exp(-(!pi^2. * fwhm2^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((srcstr[1].srcx#u)+(srcstr[1].srcy#v)))
;
;    vispred_re = re_obs1 + re_obs2
;    vispred_im = im_obs1 + im_obs2
;    
;    vis_pred=[[vispred_re], [vispred_im]]
;   endif
   
;   if srcstr[0].srctype eq 'loop' then begin
;    fwhm_loop = srcstr.fwhm_max*(1.-strcstr.eccen^2)^0.25
;    vis_pred=pso_func_makealoop_nov2021( srcstr.srcflux, fwhm_loop , srcstr.eccen, srcstr.srcx, srcstr.srcy, srcstr.srcpa, srcstr.loop_angle, vis.u, vis.v)
;   endif
      

  return, vis_pred

end