;+
; PURPOSE:
;  This function determines the merger partner of a given
;  node. A node's merger partner is the structure that merges
;  with node in the dendrogram.
;
; INPUTS:
;  node: The ID of the node to find
;  clusters: The cluster array
;
; OUTPUTS:
;  The ID of the merger partner, or -1 if no such partner exists
;  (e.g., node is the root of the dendrogram)
;
; KEYWORD PARAMETERS:
;  merge: Set to hold the id of the merger of node with its partner
;
; MODIFICATION HISTORY:
;  Sep 10, 2010: Written by Chris Beaumont
;-
function merger_partner, node, clusters, merge = merged

  if n_params() ne 2 then begin
     print, 'Calling sequence'
     print, 'result = merger_partner(node, clusters, [merged = merged])'
     return, !values.f_nan
  endif

  merged = -1
  nleaf = n_elements(clusters[0,*])+1
  hit = where(clusters eq node, ct)
  if ct eq 0 then return, -1

  ind = array_indices(clusters, hit)
  merged = ind[1] + nleaf
  return, clusters[~ind[0], ind[1]]
end
