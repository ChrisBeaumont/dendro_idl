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
function dendro_offset, ptr

  x0 = (*ptr).x[0]
  y0 = (*ptr).y[0]
  v0 = (*ptr).v[0]
  i0 = (*ptr).cubeindex[0]
  sz = (*ptr).szdata

  dx = i0 mod sz[1] - x0
  dy = (i0 / sz[1]) mod sz[2] - y0
  dz = (i0 / sz[1] / sz[2]) - v0
  return, [dx, dy, dz]
  ind = 0
end
