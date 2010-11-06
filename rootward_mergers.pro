;+
; PURPOSE:
;  This function returns the the ids of all the nodes in a dendrogram
;  that are visited during a traversal from a given id to the root.
;
; INPUTS:
;  id: The id of the dendrogram for which to find rootward mergers
;  clusters: The cluster matrix describing the dendrogram
;
; OUTPUTS:
;  A list of the rootward mergers of id, including id
;
; MODIFICATION HISTORY:
;  November 2010: Written by Chris Beaumont
;-
function rootward_mergers, id, clusters

  nst = max(clusters)+1
  result = intarr(nst)
  nleaf = n_elements(clusters[0,*]) + 1
  current = id
  i = 0
  while 1 do begin
     result[i++] = current
     hit = where(clusters eq current, ct)
     current = hit/2 + nleaf
     if ct eq 0 then break
  endwhile
  return, result[0:i-1]
end

pro test
  clusters =[ [0,1],[2,3],[4,5]]
  assert, min(rootward_mergers(0, clusters) eq [0, 4])
end
