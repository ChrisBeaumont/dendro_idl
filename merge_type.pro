;+
; PURPOSE:
;  This function computes the merge type of each structure in a
;  dendrogram. The merge type is defined to be:
;    0 if the structure is a leaf
;    1 if the structure is a merger of 2 leaves
;    2 if the structure is a leaf/branch merger
;    3 if the structure is a branch/branch merger
;
; INPUTS:
;  clusters: The cluster array describing the dendrogram. A [2,
;  nleaf-1] array.
;
; OUTPUTS:
;  A [2 * nleaf - 1] array, giving the merger type for each structure
;
; MODIFICATION HISTORY:
;  Nov 5 2010: Written by Chris Beaumont
;-
function merge_type, clusters
  compile_opt idl2

  nst = max(clusters) + 2
  nleaf = (size(clusters))[2] + 1
  result = intarr(nst)
  
  for i = nleaf, nst - 1 do begin
     id = clusters[*, i - nleaf]
     isLeaf = (id lt nleaf)
     result[i] = 3 - total(isLeaf)
  endfor
  return, result
end

pro test

  clusters = [[0,1],[2,3],[4,5]]
  assert, array_equal(merge_type(clusters), [0, 0, 0, 0, 1, 1, 3])

  clusters = [[0,1], [5,2],[6,3],[7,4]]
  assert, array_equal(merge_type(clusters), [0,0,0,0,0, 1, 2, 2, 2])
  print, 'All tests passed'
end
  
