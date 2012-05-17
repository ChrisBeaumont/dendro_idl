;+
; PURPOSE
;  This function uses a dendrogram pointer to recreate the data
;  cube. This cube is aligned such that cube[x,y,v]=t, where x,y,v,t
;  are the arrays in the dendrogram structure.
;
; INPUTS:
;  Ptr: A pointer to a dendrogram structure
;
; OUTPUTS:
;  A data cube
;
; MODIFICATION HISTORY:
;  November 9 2010: Written by Chris Beaumont
;  March 2011: Output cube now has same size as original cube
;  Oct 2011: Realized that the change in March 2011 is not always
;  possible, as x,y,v can be out-of-bounds relative to origional
;  cube. 
;-
function dendro2cube, ptr
  if n_params() eq 0 then begin
     print, 'calling sequence'
     print, ' cube = dendro2cube(ptr)'
     return, !values.f_nan
  endif

  is2D = range((*ptr).v) eq 0
  sz = [0, max((*ptr).x)+1, max((*ptr).y) + 1, max((*ptr).v) + 1]
  if is2D then begin
     cube = fltarr(sz[1], sz[2])
     cube[(*ptr).x, (*ptr).y] = (*ptr).t
  endif else begin
     cube = fltarr(sz[1], sz[2], sz[3])
     cube[(*ptr).x, (*ptr).y, (*ptr).v] = (*ptr).t
  endelse

  return, cube
end
