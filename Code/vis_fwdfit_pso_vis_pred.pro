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
    pa = srcstr[n].srcpa - 90.
  endfor
  
  vispred_re = fltarr(n_elements(u))
  vispred_im = fltarr(n_elements(u))
  
 case type of
   'circle': begin
              vispred_re = srcstr.srcflux * exp(-(!pi^2. * srcstr.srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((xh#u)+(yh#v)))
              vispred_im = srcstr.srcflux * exp(-(!pi^2. * srcstr.srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((xh#u)+(yh#v)))
              vis_pred = [transpose(vispred_re), transpose(vispred_im)]
             end
   
   'ellipse': begin
                u1 = cos(pa * !dtor) * u + sin(pa * !dtor) * v
                v1 = -sin(pa * !dtor)* u + cos(pa * !dtor) * v
                vispred_re = srcstr.srcflux * exp(-(!pi^2. / (4.*alog(2.)))#((u1 * srcstr.srcfwhm_max)^2. + (v1 * srcstr.srcfwhm_min)^2.))*cos(2*!pi*((xh#u)+(yh#v)))
                vispred_im = srcstr.srcflux * exp(-(!pi^2. / (4.*alog(2.)))#((u1 * srcstr.srcfwhm_max)^2. + (v1 * srcstr.srcfwhm_min)^2.))*sin(2*!pi*((xh#u)+(yh#v)))
                vis_pred = [transpose(vispred_re), transpose(vispred_im)]
              end
  
  
  'multi': begin
            re_obs1 = srcstr[0].srcflux * exp(-(!pi^2. * srcstr[0].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((xh[0]#u)+(yh[0]#v)))
            im_obs1 = srcstr[0].srcflux * exp(-(!pi^2. * srcstr[0].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((xh[0]#u)+(yh[0]#v)))
    
            re_obs2 = srcstr[1].srcflux * exp(-(!pi^2. * srcstr[1].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*cos(2*!pi*((xh[1]#u)+(yh[1]#v)))
            im_obs2 = srcstr[1].srcflux * exp(-(!pi^2. * srcstr[1].srcfwhm_max^2. / (4.*alog(2.)))#(u^2. + v^2.))*sin(2*!pi*((xh[1]#u)+(yh[1]#v)))
    
            vispred_re = transpose(re_obs1 + re_obs2)
            vispred_im = transpose(im_obs1 + im_obs2)
    
            vis_pred = [vispred_re, vispred_im]
           end
           
  'loop': begin
            fwhm_loop = srcstr.srcfwhm_max*(1.-srcstr.eccen^2)^0.25
            vis_pred = pso_func_makealoop_nov2021( srcstr.srcflux, fwhm_loop , srcstr.eccen, xh, yh, pa, srcstr.loop_angle, vis.u, vis.v)
            vis_pred = transpose(vis_pred)
          end         

endcase
 

return, vis_pred

end