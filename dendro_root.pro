;+
; PURPOSE:
;  Return the ID of the dendrogram root
;
; CALLING SEQUENCE:
;  root = dendro_root(ptr)
;
; INPUTS:
;  ptr: The dendrogram pointer
;
; OUTPUTS:
;  root: The ID of the root
;-
function dendro_root, ptr
  return, n_elements( (*ptr).height) - 1
end
