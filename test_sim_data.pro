;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; DEMO MEM_GE ALGORITHM ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

add_path, '/Users/admin/Documents/GitHub/FWDFIT-PSO/Code'

; Simulated data of a circle gaussian source with:
; - X location = 10 arcsec
; - Y location = 10 arcsec
; - FWHM = 10 arcsec
; - Flux = 10000. counts s^-1 cm^-2 arcsec^-2
ph_src = stx_sim_flare(pixel_data=pixel_data, $
  src_shape = 'gaussian', $
  src_xcen = 10., $
  src_ycen = 10., $
  src_flux = 100000., $
  src_fwhm_wd = 20., $
  src_fwhm_ht = 10., $
  src_phi = 45)

; Sums the pixels of a detector to 4 remaining virtual pixels
sumcase = 1
pixel_data = stx_pixel_sums(pixel_data, sumcase)
subc_str = stx_construct_subcollimator()
vis = stx_visgen(pixel_data, subc_str)

fact = 1.
case sumcase of
  0: begin
    indices = [0, 1] ;two big pixels used
  end
  1: begin
    indices = indgen(3) ;two big pixels and small pixel used
  end
  2: begin
    indices = [0] ;upper row pixels used
  end
  3: begin
    indices = [1] ;'lower row pixels used
  end
  4: begin
    indices = [2] ;small pixels used
    fact=2.
  end
endcase

;Computation of the effective area of the pixels used
tmp = subc_str.det.pixel.area
tmp = reform(tmp[*, 0], 4, 3)
tmp = reform(tmp[0, *])
effective_area = total(tmp[indices])

; Computation of the constant factors M1
M1 = effective_area *fact* 4./(!pi^3.)*sin(!pi/(4.*fact))
vis.obsvis = vis.obsvis / (4. * M1)
vis.sigamp = vis.sigamp / (4. * M1)

vis.sigamp = sqrt(vis.sigamp^2 + 0.05^2 * abs(vis.obsvis)^2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Parameters of the reconstructed image
imsize=[129, 129] ; number of pixels
pixel=[1., 1.] ; pixel size in arcsec

type   = 'ellipse'
fwdfit_map = vis_fwdfit_nov_2021(vis,imsize,pixel,shape=type,SRCOUT=srcstrout,/no_plot_fit)

fwdfit_pso_map = vis_fwdfit_pso_nov2021(type, vis, imsize=imsize, pixel=pixel,srcstr = srcstrout,/no_plot_fit)

loadct,5
window, 0
plot_map, fwdfit_map, /cbar

loadct,5
window, 1
plot_map, fwdfit_pso_map, /cbar


;;MEM_GE
;total_flux = vis_estimate_flux(vis, imsize[0]*pixel[0], silent=0) ;estimate of the total flux of the image
;mem_ge_map = mem_ge(vis, total_flux, percent_lambda = 0.02, imsize = imsize, pixel = pixel, silent = 0, makemap = 1)
;
;loadct, 5
;window, 0
;plot_map, mem_ge_map, /cbar


end