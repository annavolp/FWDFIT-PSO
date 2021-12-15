function vis_fwdfit_vis_pred, srcstr, vis

COMMON uvdata, u, v, pa, mapcenter

mapcenter = [0., 0.]
u = vis.u
v = vis.v 

TWOPI =  2. * !PI
pa = ((ATAN(v, u) + TWOPI) MOD TWOPI) * !RADEG

xx = srcstr.srcx - vis[0].xyoffset[0]
yy = srcstr.srcy - vis[0].xyoffset[1]

srcstr.srcx   = yy
srcstr.srcy   = -xx
srcstr.srcpa -= 90.

srcparm  = vis_fwdfit_structure2array(srcstr, mapcenter)
jdum     = FINDGEN(2 * n_elements(vis))
vis_pred = vis_fwdfit_func(jdum, srcparm)

return, vis_pred

end