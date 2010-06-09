;+
; PURPOSE:
;  This function creates masks of individiual dendrogram
;  substructures.
;
; INPUTS:
;  index: The index of the substructure to extract
;  ptr: The ptr variable returned by TOPOLOGIZE
;  labels: The cluster-labeled data (i.e. the contents of
;          data.cll.fits)
;  
; KEYWORD PARAMETERS:
;  single: If set, then only pixels belonging to the index'th
;          structure are considered. Otherwise, all of the structures
;          leafward of index are also included. This is the difference
;          between selecting an annulus and an entire area enclosed by
;          a contour.
;  xlo: The result is a cropped version of the original
;       image/cube. xlo will contain the x coordinate, in the original
;       data, of the left edge of the result.
;  ylo: Same as xlo, but for the bottom edge
;  zlo: Same as zlo, but for the front face (if the data is a cube)
;
; RESULT:
;  A mask, where values of 1 denote pixels belonging to the requested
;  structure. This mask is cropped, since most structures are
;  substantially smaller than the original datacube. The xlo,ylo,zlo
;  keywords contain the offset between this cropped cube and the
;  original data.
;
; MODIFICATION HISTORY:
;  June 2010: Written by Chris Beaumont
;-
function substruct, index, ptr, labels, single = single, $
                    xlo = xlo, ylo = ylo, zlo = zlo
  
  compile_opt idl2

  ndim = size(labels, /n_dim)
  if keyword_set(single) then leaves = index else begin
     leaves = leafward_mergers(index, (*ptr).clusters)
     if leaves[0] ne index then leaves = [index, leaves]
  endelse

  mask = byte(data * 0)

  for i = 0, n_elements(leaves) - 1, 1 do $
     mask = mask or labels eq leaves[i]
  
  hit = where(mask, ct)
  if ct eq 0 then return, -1

  ;- extract a sub-image/cube containing the mask
  ind = array_indices(labels, where(mask))
  xlo = min(ind[0,*], max=xhi)
  ylo = min(ind[1,*], max=yhi)
  if ndim eq 3 then begin
     zlo = min(ind[2,*], max=zhi)
     result = mask[xlo:xhi, ylo:yhi, zlo:zhi]
  endif else begin
     result = mask[xlo:xhi, ylo:yhi]
  endelse
  return, result
end
