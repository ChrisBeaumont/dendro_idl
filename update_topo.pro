;+
; PURPOSE:
;  This function updates the pointer returned by old versions of
;  topologize. In particular, it adds a histogram and index of the
;  cluster labels, which are required by several analysis routines.
;
; INPUTS:
;  ptr: A pointer returned by topologize.
;
; OUTPUTS:
;  An updated pointer
;
; MODIFICATION HISTORY:
;  November 2010: Written by Chris Beaumont
;-
function update_topo, ptr
  if n_params() eq 0 then begin
     print, 'calling sequence'
     print, 'newptr = update_topo(ptr)'
     return, ptr_new()
  endif
  
  if contains_tag(*ptr, 'cluster_label_h') then begin
     message, /con, 'Provided pointer is already up-to-date'
     return, ptr
  endif

  cluster_label_h = histogram( (*ptr).cluster_label, $
                               min = 0, max = n_elements((*ptr).height), $
                               rev = cluster_label_ri)
  npix = -1 & fast = 0
  result = create_struct((*ptr), 'cluster_label_h', cluster_label_h, $
                         'cluster_label_ri', cluster_label_ri, $
                         'npix', npix, 'fast', fast)
  return, ptr_new(result, /no_copy)

end
  
  
