;+
; NAME:
;   swarmintelligence
;
; PURPOSE:
;   Particle Swarm Optimization for constrained optimization
;
; CALLING SEQUENCE:
;   ARRAY = CMREPLICATE(VALUE, DIMS)
;
; DESCRIPTION:
;
;   The CMREPLICATE function constructs an array, which is filled with
;   the specified VALUE template.  CMREPLICATE is very similar to the
;   built-in IDL function REPLICATE.  However there are two
;   differences:
;
;      * the VALUE can be either scalar or an ARRAY.
;
;      * the dimensions are specified as a single vector rather than
;        individual function arguments.
;
;   For example, if VALUE is a 2x2 array, and DIMS is [3,4], then the
;   resulting array will be 2x2x3x4.
;
; INPUTS:
;
;   VALUE - a scalar or array template of any type, to be replicated.
;           NOTE: These two calls do not produce the same result:
;                  ARRAY = CMREPLICATE( 1,  DIMS)
;                  ARRAY = CMREPLICATE([1], DIMS)
;           In the first case the output dimensions will be DIMS and
;           in the second case the output dimensions will be 1xDIMS
;           (except for structures).  That is, a vector of length 1 is
;           considered to be different from a scalar.
;
;   DIMS - Dimensions of output array (which are combined with the
;          dimensions of the input VALUE template).  If DIMS is not
;          specified then VALUE is returned unchanged.
;
; RETURNS:
;   The resulting replicated array.
;
; EXAMPLE:
;   x = [0,1,2]
;   help, cmreplicate(x, [2,2])
;     <Expression>    INT       = Array[3, 2, 2]
;   Explanation: The 3-vector x is replicated 2x2 times.
;
;   x = 5L
;   help, cmreplicate(x, [2,2])
;     <Expression>    LONG      = Array[2, 2]
;   Explanation: The scalar x is replicated 2x2 times.
;
; SEE ALSO:
;
;   REPLICATE
;
; MODIFICATION HISTORY:
;   Written, CM, 11 Feb 2000
;   Fixed case when ARRAY is array of structs, CM, 23 Feb 2000
;   Apparently IDL 5.3 can't return from execute().  Fixed, CM, 24 Feb
;     2000
;   Corrected small typos in documentation, CM, 22 Jun 2000
;   Removed EXECUTE() call by using feature of REBIN() new in IDL 5.6,
;     (thanks to Dick Jackson) CM, 24 Apr 2009
;   Remove some legacy code no longer needed after above change
;     (RETVAL variable no longer defined; thanks to A. van Engelen),
;     CM, 08 Jul 2009
;   Change to square bracket array index notation; there were reports
;     of funny business with parenthesis indexing (thanks Jenny Lovell),
;     CM, 2012-08-16
;
;-
; Copyright (C) 2000, 2009, 2012, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-


function swarmintelligence, obj_fun_name, lower_bound_row, upper_bound_row, $
                            n_birds = n_birds, tolerance = tolerance, maxiter = maxiter, extra = extra, silent = silent 
             
  default, n_birds, 100.                  
  default, tolerance, 1e-6       
  default, silent, 0
  
  n_variable = n_elements(lower_bound_row)
  
  default, maxiter, float(n_variable) * float(n_birds)
  
  if n_elements(lower_bound_row) ne n_elements(upper_bound_row) then message, "Lower and upper bounds must have the same dimension."
  
  lower_bound = transpose(cmreplicate(lower_bound_row, n_birds))
  upper_bound = transpose(cmreplicate(upper_bound_row, n_birds))
  
  birds_pos = lower_bound+randomu(seed,[n_birds,n_variable])*(upper_bound-lower_bound)

  ones = fltarr(n_birds,1)

  birds_vel_aux = ones # reform(min([[upper_bound_row-lower_bound_row],[transpose([n_birds*(fltarr(1,n_variable)+1)])]], dimension=2), [1,n_variable])
  birds_vel_aux1 = 2.*birds_vel_aux*randomu(seed,[n_variable,n_birds])
  birds_vel = -birds_vel_aux+birds_vel_aux1
  
  obj_fun_eval = call_function(obj_fun_name, birds_pos, extra = extra)

  birds_best_fun_eval = obj_fun_eval
  birds_best_pos      = birds_pos

  birds_global_best_fun_eval = min(obj_fun_eval)
  birds_best_fun_eval_iter = (fltarr(20,1)+1)*!values.f_nan
  
  inertia_min = 0.1
  inertia_max = 1.1
  intertia_adapt_count = 0.
  inertia_adapt = inertia_max
  
  
  n_close_birds = max([2.,floor(n_birds/4.)])
  n_close_birds_adapt = n_close_birds
  
  self_behav_const   = 1.49
  social_behav_const = 1.49
  n_iteration_limit = 20
  
  for iter=0,maxiter-1 do begin

    close_birds_ind = fltarr(n_birds, n_close_birds_adapt)
    close_birds_ind[*,0] = indgen(n_birds)
    
    for jj = 0,n_birds-1 do begin
      x_p = LINDGEN(n_birds-1)
      y_p = RANDOMU(dseed, n_birds-1)
      z = x_p[SORT(y_p)]
      close_birds = z[0:n_close_birds_adapt-2]
      close_birds[where([close_birds GE jj])] += 1.
      close_birds_ind[jj,1:-1] = close_birds
    endfor
    
   aux_vect = birds_best_fun_eval[close_birds_ind]
   dummy = min(aux_vect,best_birds_fun_eval_row, dimension=2)
   best_birds_fun_eval_row = ARRAY_INDICES(aux_vect, best_birds_fun_eval_row)
   best_birds_fun_eval_row = best_birds_fun_eval_row[1, *]

 
   randSelf = randomu(seed,[n_birds,n_variable])
   randSocial = randomu(seed,[n_birds,n_variable])
        
   birds_new_vel = inertia_adapt * birds_vel + $
                  self_behav_const * randSelf * (birds_best_pos - birds_pos) + $
                  social_behav_const * randSocial * (birds_best_pos[close_birds_ind[(best_birds_fun_eval_row) * n_birds + findgen(n_birds)],*] - birds_pos)

   ind = where(min(finite(birds_new_vel), dimension=2) eq 1)
   birds_vel[ind,*] = birds_new_vel[ind,*]
    
   swarm_new = birds_pos + birds_vel
    
   lower_bound_ind = where(swarm_new LE lower_bound)
   upper_bound_ind = where(swarm_new GE upper_bound)
   swarm_new[lower_bound_ind] = lower_bound[lower_bound_ind]
   birds_vel[lower_bound_ind] = 0. 
   swarm_new[upper_bound_ind]= upper_bound[upper_bound_ind]
   birds_vel[upper_bound_ind] = 0.

    
    birds_pos = swarm_new
    
    obj_fun_eval = call_function(obj_fun_name, birds_pos, extra = extra)
    
    
    
    ind_best = where(obj_fun_eval LT transpose(birds_best_fun_eval))
    birds_best_fun_eval[ind_best] = obj_fun_eval[ind_best]
    birds_best_pos[ind_best,*]    = birds_pos[ind_best,*]

    birds_best_fun_eval_iter[(iter mod n_iteration_limit)] = min(birds_best_fun_eval)
     
    swarm_new_best = min(birds_best_fun_eval)
   
   if (finite(swarm_new_best)) and (swarm_new_best LT birds_global_best_fun_eval) then begin
      birds_global_best_fun_eval = swarm_new_best 
      intertia_adapt_count       = max([0, intertia_adapt_count-1])
      n_close_birds_adapt        = n_close_birds
    endif else begin
      intertia_adapt_count += 1
      n_close_birds_adapt  = min([n_birds, n_close_birds_adapt+n_close_birds])
    endelse  
    
    if intertia_adapt_count LT 2. then begin
    inertia_adapt = max([inertia_min, min([inertia_max, 2 * inertia_adapt])])
    endif else begin
      if intertia_adapt_count GT 5. then begin
        inertia_adapt = max([inertia_min, min([inertia_max, inertia_adapt / 2.])])
      endif
    endelse
     
    iter_ind = (iter mod n_iteration_limit) + 1

    if iter GT n_iteration_limit then begin
       if iter_ind EQ n_iteration_limit then begin
          birds_max_best_fun_eval_iter = birds_best_fun_eval_iter[0]
       endif else begin
        birds_max_best_fun_eval_iter = birds_best_fun_eval_iter[iter_ind]
       endelse
       obj_fun_update = abs(birds_max_best_fun_eval_iter - birds_best_fun_eval_iter[iter_ind-1]) / max([1,abs(birds_best_fun_eval_iter[iter_ind-1])])
    endif else begin
      obj_fun_update = !VALUES.F_INFINITY  
    endelse

    if obj_fun_update LE tolerance then break
  
    fopt = min(birds_best_fun_eval, location)
    xopt = birds_best_pos[location,*]
    
  endfor
  
stop_crit = 'OK'
if iter ge maxiter then stop_crit = 'Maximum iteration reached'
 
optim = {xopt:xopt, fopt:fopt, stop_crit:stop_crit}

if ~silent then print, 'Obj fun: ', fopt, '  Exit: ', stop_crit

return, optim


end