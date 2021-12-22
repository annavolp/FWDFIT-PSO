folder = '/Users/admin/Documents/GitHub/FWDFIT-PSO/'
add_path, folder + 'Code'


data_folder = getenv('SSW_STIX') + '/idl/processing/imaging/data/'

;;;;;;;;;;;; LOAD DATA

;;;;; June 7 2021 - 21:41

path_sci_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178428688_20200607T213708-20200607T215208_V01.fits' ; Path of the science L1 fits file
path_bkg_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178451984_20200607T225959-20200607T235900_V01.fits' ; Path of the background L1 fits file
time_range    = ['7-Jun-2020 21:39:00', '7-Jun-2020 21:42:49'] ; Time range to consider
energy_range  = [6,10]       ; Energy range to consider (keV)
xy_flare      = [-1600., -800.]  ; CFL solution (heliocentric, north up). Needed for the visibility phase calibration
mapcenter     = [-1650., -750.] ; Coordinates of the center of the map to reconstruct (heliocentric, north up)

;**********************************************************************************************************

; Sort the collimators from finest to coarsest: needed just for plotting the visibility amplitude and phase fit
subc_index = stix_label2ind(['3a','3b','3c','4a','4b','4c','5a','5b','5c','6a','6b','6c',$
                             '7a','7b','7c','8a','8b','8c','9a','9b','9c','10a','10b','10c'])
  
; Create the visibility structure filled with the measured data
vis=stix2vis_sep2021(path_sci_file, path_bkg_file, time_range, energy_range, mapcenter, $
  subc_index=subc_index, xy_flare=xy_flare, pixels=pixels,/silent)

;;;;;;;;;; SET PARAMETERS FOR IMAGING
imsize    = [257, 257]    ; number of pixels of the map to recinstruct
pixel     = [1.,1.]       ; pixel size in arcsec

;***************************************** PSO *************************************************

; Select the shape:
; - 'circle': circular Gaussian source
; - 'ellipse': elliptical Gaussian source
; - 'multi': double circular Gaussian source

type='ellipse'

; Comments:
; 1) to avoid the plot of the phase and amplitude fit, set /no_plot_fit;
; 2) use the keyword 'param_opt' to fix some of the parameters (and fit the remaining ones);
;    Please, see the header of the FWDFIT-PSO procedure for details
; 3) 'srcstrout_pso' is a structure containing the values of the optimized parameters
; 4) 'fitsigmasout_pso' is a structure containing the uncertainty on the optimized parameters
; 5) set /uncertainty for computing an estimate of the uncertainty on the parameters

fwdfit_pso_map = vis_fwdfit_pso_nov2021(type, vis,  param_opt=param_opt, imsize=imsize, pixel=pixel, srcstr = srcstrout_pso, fitsigmas=fitsigmasout_pso, /uncertainty)

loadct, 5
window, 0
cleanplot
plot_map, fwdfit_pso_map, /cbar

end