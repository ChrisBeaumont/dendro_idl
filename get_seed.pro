;+
; PURPOSE:
;  This function finds a the lowest-index leaf node which is contained
;  in each half of each merger.
;
; INPUTS:
;  cluster: The cluster array describing the dendrogram. A (2,
;  nleaf-1) array.
;
; OUTPUTS:
;  an array with the same shape as cluster. result[i,j] returns the
;  lowest-index substructure contained within the structure listed in
;  cluster[i,j]
;
; MODIFICATION HISTORY:
;  July 2010: Written by Chris Beaumont
;-
function get_seed, cluster
  compile_opt idl2

  if n_params() ne 1 then begin
     print, 'calling sequence:'
     print, ' result = get_seed(cluster)'
     print, '   result = (2, nleaf-1) array'
     return, !values.f_nan
  endif

  nleaf = n_elements(cluster)/2+1
  nst = max(cluster+1)
  result = cluster

  h = histogram(cluster, min = 0, reverse = ri)
  assert, range(h) eq 0 && h[0] eq 1

  for i = 0, nleaf - 1, 1 do begin
     ;- what clusters does leaf i show up in?
     m = cluster_member(i, cluster)
     m = m[0:n_elements(m)-2] ;- ignore the "root" node
     ;- at which indices of cluster do _these_ numbers appear?
     ind = ri[ri[m]]
     result[ind] <= i
  endfor
  return, result
end

pro test
  print, 'test1'
  cluster = [[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13]]
  ;- should be 01 23 45 67 02 46 04
  print, get_seed(cluster)
  print, 'test2'
  cluster = [[0,1], [4,2], [5,3]]
  ;- should be 01 02 03
  print, get_seed(cluster)
  
end
