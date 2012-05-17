;+
; PURPOSE:
;  This function extracts individual substructures from a dendrogram
;  hierarchy. 
;
; INPUTS:
;  index: The substructure to select. This can be one of several
;  values: 1) A scalar corresponding to a single ID and (if /single is
;  not set), all of the substructures of that Id; 2) An array of ids
;  (in which case only these ids, and not their substructures, will be
;  returned)
;  ptr: The ptr variable returned by TOPOLOGIZE
;  
; KEYWORD PARAMETERS:
;  single: If set, then only pixels belonging to the index'th
;          structure are considered. Otherwise, all of the structures
;          leafward of index are also included. This is the difference
;          between selecting an annulus and an entire area enclosed by
;          a contour. Single is ignored (and treated as un-set) if
;          index is an array.
;
;  count: Set to a variable to hold the number of indices in structure index,
;  or 0 if there are none
;
; RESULT:
;  The indices of the pixels in the (x, y, v, t) tags which contain the
;  requested structure. CAUTION: The x,y,v values are offset from the original
;  data cube. the cubeindex structure will map from the cube described by xyvt
;  to the original data cube.
;
; MODIFICATION HISTORY:
;  June 2010: Written by Chris Beaumont
;  Jan 2011: Added parameter checking. cnb.
;  Feb 2011: Removed support for specifying leaves via index=-1. cnb.
;  July 2011: Removed use of stack for speed.
;-
function substruct, index, ptr, single = single, count = count
                    
  compile_opt idl2

  if n_params() ne 2 then begin
     print, 'calling sequence'
     print, 'ind = substruct(index, ptr, [/single, count = count])'
     return, !values.f_nan
  endif
  count = 0

  if keyword_set(single) || n_elements(index) gt 1 then indices = index else begin
     indices = leafward_mergers(index, (*ptr).clusters)
  endelse
  if min(indices) lt 0 then return, -1

  count = total( (*ptr).cluster_label_h[indices], /preserve)
  if count eq 0 then return, -1
  result = lonarr(count)
  offset = 0
  for i = 0, n_elements(indices) - 1, 1 do begin
     x = indices[i]
     if x lt 0 || x ge n_elements((*ptr).cluster_label_h) then continue     
     if (*ptr).cluster_label_h[x] eq 0 then continue

     ind = (*ptr).cluster_label_ri[(*ptr).cluster_label_ri[x] : $
                                   (*ptr).cluster_label_ri[x+1]-1]
     nind = n_elements(ind)
     result[offset : offset + nind - 1] = ind
     offset += nind
  endfor

  return, result

end
