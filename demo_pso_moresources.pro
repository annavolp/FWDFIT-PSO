;************************************************************************************************************************************************;
;                                                                                                                                                ;
;                                                        PSO DEMO                                                              ;
;                                                                                                                                                ;
;************************************************************************************************************************************************;

;; ADD PATH OF THE FOLDER CONTAINING THE PSO MULTIPLE SOURCES CODE
add_path, '/Users/admin/Documents/GitHub/FWDFIT-PSO/PSO_multiple_sources_code' 


; forward fitting method from visibility based on Particle Swarm Optimization

;***************** SET PARAMETERS
imsize = [257,257] ;array containing the size (number of pixels) of the image to reconstruct
pixel  = [1.,1.]   ;array containing the pixel size (in arcsec) of the image to reconstruct

;****************************************** 23-Sep-2021 15:20
;;;***** Insert UID of the science L1 fits file (see STIX website)
uid = '2109230031'
website_url = 'https://datacenter.stix.i4ds.net/download/fits/bsd/'
sock_copy, website_url + uid, out_name, status = status, out_dir = data_folder, local_file=path_sci_file, clobber=0

time_range    = ['23-Sep-2021 15:20:30', '23-Sep-2021 15:22:30']
energy_range  = [18,22]
mapcenter     = [675.,-650.] ; Coordinates of the center of the map to reconstruct (quasi-heliocentric)
xy_flare      = mapcenter    ; Location of the map (quasi-heliocentric). Needed for the visibility phase calibration

;; Compute the visibilities
vis=stix2vis_sep2021(path_sci_file, time_range, energy_range, mapcenter, xy_flare=xy_flare, /silent)

;*************************************** VIS_FWDFIT_PSO ************************************************
; configuration: string array containing parametric shapes chosen for the forward fitting method
;     - Gaussian circular source       : 'circle'
;     - Gaussian elliptical source     : 'ellipse'
;     - Loop source                    : 'loop'

configuration = ['circle', 'circle', 'circle', 'circle']

; 'srcstrout' is a structure containing the values of the optimized parameters

vis_fwdfit_pso_map = stx_vis_fwdfit_pso(configuration, vis, imsize = imsize, pixel = pixel, srcstr = srcstrout)

loadct,5 
window, 0
cleanplot
plot_map, vis_fwdfit_pso_map, /cbar, title='VIS_FWDFIT_PSO', /limb, grid_spacing=5


print, " "
print, "Press SPACE to continue"
print, " "
pause

;********************** UNCERTAINTY
; set /uncertainty for computing an estimate of the uncertainty on the parameters
; 'fitsigmasout_unc' is a structure containing the uncertainty on the optimized parameters

vis_fwdfit_pso_map_unc = stx_vis_fwdfit_pso(configuration, vis, imsize=imsize, pixel=pixel, $
                                            srcstr = srcstrout_unc, $
                                            fitsigmas=fitsigmasout_unc, /uncertainty, $
                                            redchisq = redchisq)

loadct, 5
window, 1
cleanplot
plot_map, vis_fwdfit_pso_map_unc, /cbar, title='VIS_FWDFIT_PSO uncertainty', /limb, grid_spacing=5

print, " "
print, "Press SPACE to continue"
print, " "
pause

;********************** FIX SOME PARAMETERS or CHANGE LOWER and UPPER BOUND
; SRCIN       : struct containing for each source the parameters to be optimized and those fixed, upper and lower bound of the variables.
;               For different shapes:
;                 - 'circle'  : param_opt, lower_bound, upper_bound = [flux, x location, y location, FWHM]
;                 - 'ellipse' : param_opt = [flux, x location, y location, FWHM max, FWHM min, alpha]
;                               lower_bound, upper_bound = [flux, x location, y location, FWHM, ecc * cos(alpha), ecc * sin(alpha)]
;                                 'ecc' is the eccentricity of the ellipse and 'alpha' is the orientation angle
;                 - 'loop'    : param_opt = [flux, x location, y location, FWHM max, FWHM min, alpha, loop_angle]
;                               lower_bound, upper_bound = [flux, x location, y location, FWHM, ecc * cos(alpha), ecc * sin(alpha), loop_angle]

configuration2 = ['circle', 'circle', 'ellipse', 'circle']
srcin = VIS_FWDFIT_PSO_MULTIPLE_SRC_CREATE(vis, configuration2)
srcin.circle[0].param_opt.param_x = '644.'
srcin.circle[0].param_opt.param_y = '-687.'
srcin.circle[1].param_opt.param_x = '664.'
srcin.circle[1].param_opt.param_y = '-658.'
srcin.circle[2].param_opt.param_x = '685.'
srcin.circle[2].param_opt.param_y = '-648.'
srcin.ellipse[0].param_opt.param_x = '719.'
srcin.ellipse[0].param_opt.param_y = '-621.'

; x and y position of each souces fixed, to fit the remaining parameters.

vis_fwdfit_pso_map_fix = stx_vis_fwdfit_pso(configuration2, vis,  srcin = srcin, $
                                            imsize=imsize, pixel=pixel, $
                                            srcstr = srcstrout_fix, $
                                            fitsigmas=fitsigmasout_fix, /uncertainty)
                                             
loadct, 5
window, 2
cleanplot
plot_map, vis_fwdfit_pso_map_fix, /cbar,title='VIS_FWDFIT_PSO x,y fix', /limb, grid_spacing=5


end
