;+
; PURPOSE:
;  This function generates an integer array which is meant to distill
;  a dendrogram into its "essential structure". This array lists the
;  traversal of the tree in prefix form (Polish order), with the
;  following modifications:
;   -Instead of labeling non-leaf nodes by their numbers, they are all
;    labeled as -1.
;   -The two components of each structure are sorted by value,
;    where a structure's value is the minimum id of the
;    substructures it contains.
;  these modifications ensure that two trivially-reordered versions
;  of the same dendrogram generate the same key
;
; INPUTS:
;  cluster: The cluster array describing the dendrogram
;
; KEYWORD PARAMETERS:
;  index: If supplied, then only consider the tree rootward of this
;         structure.
;  minval: On output, will contain the "value" of this structure. This
;          is mainly used internally to sort the dendrogram
;
; OUTPUTS:
;  The dendrogram key, as an array of integers
;
; MODIFICATION HISTORY:
;  July 20 2010: Written by Chris Beaumont
;-
function dendrokey, cluster, index = index, minval = minval
  if n_elements(index) eq 0 then return, dendrokey(cluster, index=max(cluster)+1, minval = minval)
  nleaf = n_elements(cluster)/2+1
  if index lt nleaf then begin
     minval = index
     result = index
  endif else begin
     row = index - nleaf
     left = dendrokey(cluster, index = cluster[0, row], minval = min1)
     right = dendrokey(cluster, index = cluster[1, row], minval = min2)
     result = min1 lt min2 ? [-1, left, right] : [-1, right, left]
     minval = min1 < min2
  endelse
  return, result
end

pro test

  clusters = [[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13]]
  
  print, dendrokey(clusters, minval = minval) & print, minval
  tmp = clusters[1,*]
  clusters[1,*] = clusters[0,*]
  clusters[0,*] = tmp
  print, dendrokey(clusters, minval = minval) & print, minval
  clusters=[[2,3],[0,1],[4,5],[6,7],[8,9],[10,11],[12,13]]
  print, dendrokey(clusters, minval=minval) & print, minval
  print, cluster_member(0, clusters)
end
