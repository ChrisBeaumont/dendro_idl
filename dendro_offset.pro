;+
; PURPOSE:
;  This function calculates the offset between the indices used by TOPOLOGIZE
;  et al., and those of the original data cube. 
; 
;  Topologize makes a copy of the original cube, and transforms it into a
;  new cube ("minicube") that may be padded with zeros along the edges and/or
;  cropped. This procedure determines the translation that aligns the internal
;  "minicube" with the original data.
;
; INPUTS:
;  ptr: The pointer output from topologize
;  cube: The original data cube
;
; OUTPUTS:
;  A 3 element vector which gives the location, in the original data, of the
;  first element in the internal topologize cube. Note that these indices may
;  be negative, since the topologize cube is occasionally padded with
;  zeros along the edges.
;
;  In other words, cube[x+offset[0], y+offset[1],
;  v+offset[2]] will equal t, where cube is the original data, and x,y,v,t are
;  tags in the structure created by topologize.
;
; MODIFICATION HISTORY:
;  June 21 2010: Written by Chris Beaumont
;-
function dendro_offset, ptr, cube
  ind = 0
  maxtry = 100
  for ii = 0, maxtry - 1, 1 do begin
     hit = where((*ptr).t[ii] eq cube, ct)
     if ct ne 1 then continue
     ind = array_indices(cube, hit)
     if n_elements(ind) eq 2 then i = [ind, 0]
     offset = ind - [(*ptr).x[ii], (*ptr).y[ii], (*ptr).v[ii]]
     return, offset
  endfor
  message, "Error! Couldn't calculate offset"
end
