;+
; PURPOSE:
;  This function creates a mask of an image where the "1s" in the
;  output mask correspond to regions with intensities greater than
;  "grow" * sigma, with at least one pixel greater than "sig" * sigma
;
; INPUTS:
;  data: A 2- or 3-d image to mask
;  sig:
;  grow:
;
; KEYWORD PARAPMETERS:
;  nonuniform: If set and "error" is not supplied, then the noise will
;              be estimated assuming it varies over the data
;  error: The 1-sigma noise in the data. Can be a scalar or an array,
;         specifying the noise at each point. This will be calculated
;         if not provided
;  width: Keyword is passed to sigma_cube when calculating irregular
;         noise - gives the sampling size.
;
; MODIFICATION HISTORY:
;  pre-2010: Written by Erik Rosolowsky
;  May 2010: Generalized to handle 3D input arrays. Chris Beaumont
;-
function gmask, data, sig, grow, nonuniform = nonuniform, $
                error = err, width = width

  ;- check inputs
  ndim = size(data, /n_dim)
  if ndim ne 2 && ndim ne 3 then message, 'Input data must be a 2- or 3-D array'

  if n_elements(sig) eq 0 then sig = 3
  if n_elements(grow) eq 0 then grow = 2
  
  if n_elements(err) eq 0 then begin 
     if keyword_set(nonuniform) then err = sigma_cube(data, width = width) $
     else err = er_mad(data)
  endif
  
  ;- the origional "seed" mask of bright pixels
  mask = data gt sig*err
  
  ;- for 3D cubes, set the edges along the 3rd dimension to zero
  if ndim eq 3 then mask = (mask*(shift(mask, 0, 0, -1)+shift(mask, 0, 0, 1)) gt 0)

  ;- expand mask out to all points greater than "grow" * err
  growmask = data gt grow*err
  if ndim eq 3 then $
     growmask = (growmask*(shift(growmask, 0, 0, -1)+shift(growmask, 0, 0, 1)) gt 0)
  mask = dilate_mask(mask, constraint = growmask)

  return, mask
end
