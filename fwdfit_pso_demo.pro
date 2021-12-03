
;; Add the path of the folder containing the code
;add_path, '/Users/admin/Documents/GitHub/Amplitudes-PSO/PSO/'
add_path, '/Users/admin/Documents/GitHub/FWDFIT-PSO/Code'

;;;;;;;;;;;;

data_folder = getenv('SSW_STIX') + '/idl/processing/imaging/data/'

;;;;;;;;;;;; LOAD DATA

;;;;; June 7 2021 - 21:41

path_sci_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178428688_20200607T213708-20200607T215208_V01.fits' ; Path of the science L1 fits file
path_bkg_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178451984_20200607T225959-20200607T235900_V01.fits' ; Path of the background L1 fits file
time_range = ['7-Jun-2020 21:39:00', '7-Jun-2020 21:42:49'] ; Time range to consider
energy_range = [6,10]       ; Energy range to consider (keV)
xy_flare = [-1600., -800.]  ; CFL solution (heliocentric, north up). Needed for the visibility phase calibration
mapcenter = [-1650., -750.] ; Coordinates of the center of the map to reconstruct (heliocentric, north up)

;;;;;;;;;; CONSTRUCT VISIBILITY STRUCTURE

; Create the visibility structure filled with the measured data
vis=stix2vis_sep2021(path_sci_file, path_bkg_file, time_range, energy_range, mapcenter, $
  subc_index=subc_index, xy_flare=xy_flare, pixels=pixels)
  
;;;;;;;;;; SET PARAMETERS FOR IMAGING

imsize    = [129, 129]    ; number of pixels of the map to recinstruct
pixel     = [2.,2.]       ; pixel size in arcsec

type   = 'ellipse'
n_free = n_elements(vis) - 4
results = vis_fwdfit_pso_nov2021(type, vis, n_free)
fwdfit_pso_map = vis_FWDFIT_SOURCE2MAP(results.srcstr, type=type, pixel=pixel, imsize=imsize, xyoffset=vis[0].xyoffset)


loadct,5
plot_map, fwdfit_pso_map, /cbar

end