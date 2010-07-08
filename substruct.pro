;+
; PURPOSE:
;  This function extracts individual substructures from a dendrogram
;  hierarchy. 
;
; INPUTS:
;  index: The index of the substructure to extract
;  ptr: The ptr variable returned by TOPOLOGIZE
;  
; KEYWORD PARAMETERS:
;  single: If set, then only pixels belonging to the index'th
;          structure are considered. Otherwise, all of the structures
;          leafward of index are also included. This is the difference
;          between selecting an annulus and an entire area enclosed by
;          a contour.
;
; RESULT:
;  The indices of the pixels in the (x, y, v, t) tags which contain the
;  requested structure. CAUTION: The x,y,v values are offset from the original
;  data cube by an unspecified amount. Use the dendro_offset function to
;  calculate this offset.
;
; MODIFICATION HISTORY:
;  June 2010: Written by Chris Beaumont
;-
function substruct, index, ptr, single = single
                    
  compile_opt idl2

  if keyword_set(single) then indices = index else $
     indices = leafward_mergers(index, (*ptr).clusters)
  
  s = obj_new('stack')
  for i = 0, n_elements(indices) - 1, 1 do begin
     x = indices[i]
     if (*ptr).cluster_label_h[x] eq 0 then continue
     ind = (*ptr).cluster_label_ri[(*ptr).cluster_label_ri[x] : $
                                   (*ptr).cluster_label_ri[x+1]-1]
     s->push, ind
  endfor

  bad = (s->getSize() eq 0)
  result = s->toArray()
  obj_destroy, s

  return, bad ? !values.f_nana : result
end
