folder = 'C:\Users\volpa\Desktop\Dottorato\Codice CFL\BP test\'
data_folder = 'C:\Users\volpa\Desktop\Dottorato\PSO 8_12\'

add_path, 'C:\Users\volpa\Desktop\Dottorato\Codice CFL\CFL'
add_path, 'C:\Users\volpa\Desktop\Dottorato\Codice CFL\New codes'
add_path, 'C:\Users\volpa\Desktop\Dottorato\PSO 8_12'

restore, folder+'File sav\flare_files.sav', /v

; Path of the science L1 fits file
path_sci_file = data_folder+ 'solo_L1A_stix-sci-xray-l1-2110280004_20211028T152248-20211028T152656_016712_V01.fits'
; Path of the background L1 fits file
path_bkg_file = data_folder+ 'solo_L1A_stix-sci-xray-l1-2110220013_20211022T003003-20211022T021003_016179_V01.fits' 
; Time range to consider
time_range = ['28-Oct-2021 15:25:00', '28-Oct-2021 15:26:40']
; Energy range to consider (keV)
energy_range = [32,70]    
; CFL solution (heliocentric, north up). Needed for the visibility phase calibration
xy_flare = [260., -600.] 
; Coordinates of the center of the map to reconstruct (heliocentric, north up) 
mapcenter = [260., -600.] 

;**********************************************************************************************************
subc_index = stix_label2ind(['3a','3b','3c','4a','4b','4c','5a','5b','5c','6a','6b','6c',$
  '7a','7b','7c','8a','8b','8c','9a','9b','9c','10a','10b','10c'])
  
; Create the visibility structure filled with the measured data
vis=stix2vis_sep2021(path_sci_file, path_bkg_file, time_range, energy_range, mapcenter, $
  subc_index=subc_index, xy_flare=xy_flare, pixels=pixels,/silent)

;;;;;;;;;; SET PARAMETERS FOR IMAGING
imsize    = [129, 129]    ; number of pixels of the map to recinstruct
pixel     = [2.,2.]       ; pixel size in arcsec

;***************************************** PSO *************************************************
type='ellipse' 
;param_opt=['fit', 'fit', 'fit', 'fit']

fwdfit_pso_map = vis_fwdfit_pso_nov2021(type, vis, SwarmSize=100,  param_opt=param_opt, seedstart=seedstart, $
      imsize=imsize, pixel=pixel, srcstr = srcstrout_pso, fitsigmas =fitsigmasout_pso, no_plot_fit=1)
loadct, 5
window, 0
cleanplot
plot_map, fwdfit_pso_map, /cbar

;chi2
stix_plot_fit, fwdfit_pso_map, vis, imsize, pixel, wwindow=1, title='PSO'

end