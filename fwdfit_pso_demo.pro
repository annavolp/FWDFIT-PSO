
;; Add the path of the folder containing the code
;add_path, '/Users/admin/Documents/GitHub/Amplitudes-PSO/PSO/'
add_path, '/Users/admin/Documents/GitHub/FWDFIT-PSO/Code'

;;;;;;;;;;;;

data_folder = getenv('SSW_STIX') + '/idl/processing/imaging/data/'

;;;;;;;;;;;; LOAD DATA

;;;;; June 7 2021 - 21:41

path_sci_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178428688_20200607T213708-20200607T215208_V01.fits' ; Path of the science L1 fits file
path_bkg_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178451984_20200607T225959-20200607T235900_V01.fits' ; Path of the background L1 fits file
time_range    = ['7-Jun-2020 21:39:00', '7-Jun-2020 21:42:49'] ; Time range to consider
energy_range  = [6,10]       ; Energy range to consider (keV)
xy_flare      = [-1600., -800.]  ; CFL solution (heliocentric, north up). Needed for the visibility phase calibration
mapcenter     = [-1650., -750.] ; Coordinates of the center of the map to reconstruct (heliocentric, north up)

;;;;;;;;;; CONSTRUCT VISIBILITY STRUCTURE

subc_index=stix_label2ind(['3a','3b','3c','4a','4b','4c','5a','5b','5c','6a','6b','6c','7a','7b','7c',$
  '8a','8b','8c','9a','9b','9c','10a','10b','10c'])


; Create the visibility structure filled with the measured data
vis=stix2vis_sep2021(path_sci_file, path_bkg_file, time_range, energy_range, mapcenter, $
  subc_index=subc_index, xy_flare=xy_flare, pixels=pixels,/silent)
  
;;;;;;;;;; SET PARAMETERS FOR IMAGING

imsize    = [257, 257]    ; number of pixels of the map to recinstruct
pixel     = [1.,1.]       ; pixel size in arcsec

type   = 'multi'
fwdfit_map = vis_fwdfit_nov_2021(vis,imsize,pixel,shape=type,SRCOUT=srcstrout,/no_plot_fit)
vis_pred_analytical = vis_fwdfit_vis_pred(srcstrout, vis)
re_vis_pred_analytical = vis_pred_analytical[0:23]
im_vis_pred_analytical = vis_pred_analytical[24:47]


F = vis_map2vis_matrix(vis.u, vis.v, imsize, pixel)
im_map = rotate(fwdfit_map.data,3)
vis_pred_image = F # im_map[*]
re_vis_pred_image = real_part(vis_pred_image)
im_vis_pred_image = imaginary(vis_pred_image)


fwdfit_pso_map = vis_fwdfit_pso_nov2021(type, vis, imsize=imsize, pixel=pixel,srcstr = srcstrout,/no_plot_fit)
vis_pred_analytical_pso = vis_fwdfit_pso_vis_pred(srcstrout, vis, type)
re_vis_pred_analytical_pso = vis_pred_analytical_pso[0:23]
im_vis_pred_analytical_pso = vis_pred_analytical_pso[24:47]

F = vis_map2vis_matrix(vis.u, vis.v, imsize, pixel)
im_map = rotate(fwdfit_pso_map.data,3)
vis_pred_image_pso = F # im_map[*]
re_vis_pred_image_pso = real_part(vis_pred_image_pso)
im_vis_pred_image_pso = imaginary(vis_pred_image_pso)

;fwdfit_pso_map = vis_fwdfit_pso_nov2021(type, vis, imsize=imsize, pixel=pixel)
;
;
;loadct,5
;window, 2
;plot_map, fwdfit_pso_map, /cbar
;
;stix_plot_fit, fwdfit_pso_map, vis, imsize, pixel, wwindow=3, title=title

end