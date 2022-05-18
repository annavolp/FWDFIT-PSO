add_path, 'C:\Users\volpa\Documents\GitHub\FWDFIT-PSO\Code'


;data_folder = getenv('SSW_STIX') + '/idl/processing/imaging/data/'
data_folder_stix = 'C:\Users\volpa\Desktop\Dottorato\PSO 8_12\'

;;;;;;;;;;;; LOAD DATA

;**************************************************************************************
;**************************************** STIX ****************************************
;**************************************************************************************

;;;;; ******************************* May 7 2021 - 18:52 *********************************
path_sci_file = data_folder_stix + 'solo_L1A_stix-sci-xray-l1-2105070034_20210507T184238-20210507T190900_011203_V01.fits'
path_bkg_file = data_folder_stix + 'solo_L1A_stix-sci-xray-l1-2105080012_20210508T020005-20210508T034005_010936_V01.fits'
time_range = ['7-May-2021 18:51:00', '7-May-2021 18:53:40']
xy_flare   = [325., 261.] ; Adapted from imaging
mapcenter  = xy_flare
;;;;; Thermal range
energy_range = [6,10]
;;;;;; Non - Thermal range
;energy_range = [22,50]


;;;;; ******************************* August 26 2021 - 23:19  ****************************
;path_sci_file = data_folder_stix+ 'solo_L1A_stix-sci-xray-l1-2108260030_20210826T231549-20210826T232115_013330_V01.fits'; Path of the science L1 fits file
;path_bkg_file = data_folder_stix + 'solo_L1A_stix-sci-xray-l1-2108200018_20210820T012017-20210820T025617_012845_V01.fits'
;time_range = ['26-Aug-2021 23:18:00', '26-Aug-2021 23:20:00'] ; Time range to consider
;xy_flare = [670., -1185.]
;mapcenter = xy_flare  ; Coordinates of the center of the map to reconstruct (heliocentric, north up)
;;;;;;;;; Thermal range
;;energy_range = [6,10]
;;;;;;;;; Non-thermal range
;energy_range = [15,25]    ; Energy range to consider (keV)

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
; - 'circle' : circular Gaussian source
; - 'ellipse': elliptical Gaussian source
; - 'multi'  : double circular Gaussian source
; - 'loop'   : single curved elliptical gaussian

type='ellipse'

; Comments:
; 1) use the keyword 'param_opt' to fix some of the parameters (and fit the remaining ones);
;    Please, see the header of the FWDFIT-PSO procedure for details
; 2) 'srcstrout_pso' is a structure containing the values of the optimized parameters
; 3) 'fitsigmasout_pso' is a structure containing the uncertainty on the optimized parameters
; 4) set /uncertainty for computing an estimate of the uncertainty on the parameters

vis_fwdfit_pso_map = stx_vis_fwdfit_pso(type, vis,  param_opt=param_opt, imsize=imsize, pixel=pixel, $
  srcstr = srcstrout_pso, fitsigmas=fitsigmasout_pso, /uncertainty)                                        
                                   
loadct, 5
window, 0
cleanplot
plot_map, vis_fwdfit_pso_map, /cbar

stop

;**************************************************************************************
;*************************************** RHESSI ***************************************
;**************************************************************************************
search_network, /enable
imsize = [129 , 129]
pixel = [1.0000, 1.0000]
obj = hsi_image()

obj-> set, cbe_normalize= 1
obj-> set, snr_chk= 1
obj-> set, snr_thresh= 4.00000
obj-> set, mpat_coord= 'CART'
obj-> set, image_dim= imsize
obj-> set, pixel_size= pixel
obj-> set, im_time_interval= ['13-Feb-2002 12:29:40.200', '13-Feb-2002 12:31:22.800']
obj-> set, use_flux_var= 0L
obj-> set, use_phz_stacker= 1L
obj-> set, xyoffset= [-902.999, -374.331]
obj-> set, im_energy_binning= [6.00000, 12.0000]
obj-> set, time_bin_def= [1.00000, 2.00000, 4.00000, 8.00000, 8.00000, 16.0000, 32.0000, $
  64.0000, 128.000]
obj-> set, time_bin_min= 256L
obj-> set, det_index_mask= [0B, 0B, 1B, 1B, 1B, 1B, 1B, 1B, 1B]
obj-> set, vf_loop= 1

vo = obj->get(/obj, class='hsi_visibility') ;salvo vis
vis = vo->getdata()

type='loop'
param_out = vis_fwdfit_pso(type, vis)
srcstr = param_out.srcstr
vis_fwdfit_pso_map = vis_FWDFIT_PSO_SOURCE2MAP(srcstr, type=type, pixel=pixel, imsize=imsize, xyoffset=vis[0].xyoffset)
vis_fwdfit_pso_map=make_map(vis_fwdfit_pso_map)

loadct,5
window,1
plot_map, vis_fwdfit_pso_map, /limb, grid_spacing=5.







end