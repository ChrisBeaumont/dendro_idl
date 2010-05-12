function cnb_decimate_kernels, k_in,  cube, ALL_NEIGHBORS = all_neighbors, $
                               sigma = sigma, delta = delta
  
;+
; NAME:
;   DECIMATE_KERNELS
;
; PURPOSE:
;   To eliminate kernels based on their contrast with the remainder of
;   the emission. This is the first of a 2-pass rejection strategy
;   (the second being decimate_merger). The first pass eliminates the
;   obviously insignificant kernels, which can be identified cheaply.
;   This speeds up future calculations.
;
; CALLING SEQUENCE:
;   new_kernels = DECIMATE_KERNELS(kernels, cube [,/ALL_NEIGHBORS,
;   DELTA = delta, SIGMA = sigma]) 
;
; INPUTS:
;   KERNELS -- The potential kernels as indices in the array CUBE
;   CUBE -- The data cube in question.
;
; KEYWORD PARAMETERS:
;   ALL_NEIGHBORS -- Set this keyword to make pixels have 26 neighbors
;                    and not 6.
;
;   DELTA -- The height that a maximum must be (in units of SIGMA)
;            above a saddle point to be considered a unique kernel.
;
;   SIGMA -- The units of error in the data cube.  Defaults to the
;            median, absolute deviation.
;
; OUTPUTS:
;   NEW_KERNELS -- Kernels decimated by requiring a local maximum to
;                  be N-sigma above the valley where it merges with
;                  another cloud.  
; REQUIRES: 
;   MERGEFIND
;
; MODIFICATION HISTORY:
;   Feb 2010: Adapted from decimate_kernels by Chris Beaumont
;-

; COPY KERNELS TO ALLOW EDITTING/AVOID UNFORTUNATE ACCIDENTS
  kernels = k_in

; CONTRAST CRITERIA
  if n_elements(delta) eq 0 then $
    delta = 2.
  area_rejects = 0
  delta_rejects = 0
  inds = array_indices(cube, kernels)

; MEAN ABSOLUTE DEVIATION (CONVERTED TO SIGMA) OF DATA IN THE CUBE
  if n_elements(sigma) eq 0 then $
    sigma = er_mad(cube)
  
; Establish order of decimation: LOWEST -> HIGHEST
  kernel_value = cube[kernels]
  order = sort(kernel_value)

  ; Set a vector which tracks good kernels to all kernels good.
  valid_kernel = intarr(n_elements(kernels))+1
  dbl_cube = double(cube)
  bad = where(~finite(dbl_cube), badct)
  if badct ne 0 then dbl_cube[bad] = min(dbl_cube,/nan)
  s = sort(cube)

  ;- show progress
  pbar, /new, name='Decimate Kernels'
  for i = 0, n_elements(order)-1 do begin
     if (i mod 5) eq 0 then pbar, 1D * i / (n_elements(order)-1)
     valid_kernel_ind = where(valid_kernel, valid_ct)
     newkernels = kernels[valid_kernel_ind]
     
     seed = array_indices(cube, kernels[order[i]])
     val = cube[kernels[order[i]]] - sigma * delta
     frac = 1D * interpol(s, cube[s], val) / n_elements(s)
     label_seed, dbl_cube, val, seed, result, $
                 all_neighbors = all_neighbors, $
                 external = (frac lt .5)
     assert, total(result[newkernels]) ge 1
     
     ;- test for rejection - region overlaps another seed
     if total(result[newkernels]) gt 1 then begin
        delta_rejects++
        valid_kernel[order[i]] = 0
     endif     
  endfor
  pbar, /close

  message, 'Kernels rejected for contrast: '+string(delta_rejects), /con
  valid_kernel_ind = where(valid_kernel, valid_ct)
  
  newkernels = kernels[valid_kernel_ind]
  return, newkernels
end

