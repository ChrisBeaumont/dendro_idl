;+
; PURPOSE:
;  This procedure determines the minimum and maximum depth of each
;  node in a dendrogram. These are defined as the minimum / maximum
;  number of nodes visited from some leaf to this node. Leaves have
;  depths of zero.
;
; INPUTS:
;  cluster: The cluster array describing the dendrogram
;
; OUTPUTS:
;  min: The min depth. A (2 x nleaf - 1) array
;  max: The max depth.
;
; MODIFICATION HISTORY:
;  July 21 2010: Written by Chris Beaumont
;-
pro cluster_depth, cluster, min, max
  nleaf = n_elements(cluster)/2+1
  min = intarr(max(cluster+2)) + max(cluster)+5
  max = intarr(max(cluster+2)) + 0
  for i = 0, nleaf - 1, 1 do begin
     mem = cluster_member(i, cluster)
     ind = indgen(n_elements(mem))
     min[mem] <= ind & max[mem] >= ind
  endfor
end

pro test

  c = [[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13]]
  cluster_depth, c, min, max
  print, min, max

  c = [[0,1],[4,2],[5,3]]
  cluster_depth, c, min, max
  print, min, max
end
