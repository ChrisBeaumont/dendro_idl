;+
; PRUPOSE:
;  This function returns the number of each  merger types for all the
;  substructures contained within each node of a dendrogram. See
;  merge_type.pro for how the merge type is defined. 
;
; INPUTS:
;  clusters: The cluster matrix describing the dendrogram
; 
; KEYWORD PARAMETERS:
;  type: Set to a variable to hold The merge type for each
;  _individual_ node (i.e. the output from merge_type).
;
; OUTPUTS:
;  A (4, nstructure) array. The [i,j] entry in this matrix lists how
;  many substructures of node j are of merge type i.
;
; MODIFICATION HISTORY:
;  November 2010: Written by Chris Beaumont
;-

pro _recurse_mtc, clusters, type, id, singletype

  nleaf = (size(clusters))[2] + 1
  isLeaf = id lt nleaf

  if isLeaf then begin
     type[*,id] = [1,0,0,0]
     return
  endif

  sub = clusters[*,id-nleaf]
  _recurse_mtc, clusters, type, sub[0], singletype
  _recurse_mtc, clusters, type, sub[1], singletype
  type[*,id] = type[*,sub[0]] + type[*,sub[1]]
  type[singletype[id], id] += 1
end

function merge_type_cumul, clusters, type = type
  compile_opt idl2

  nst = max(clusters) + 2
  nleaf = (size(clusters))[2] + 1
  result = intarr(4, nst)

  type = merge_type(clusters)
  
  _recurse_mtc, clusters, result, nst - 1, type
  return, result
end

pro test
  clusters = [[0,1],[2,3],[4,5]]
  assert, array_equal(merge_type_cumul(clusters), $
                      [ [1,0,0,0], [1,0,0,0], $
                        [1,0,0,0], [1,0,0,0], $
                        [2,1,0,0], [2,1,0,0], $
                        [4, 2, 0, 1]])
  
  clusters = [[0,1], [5,2],[6,3],[7,4]]
  assert, array_equal(merge_type_cumul(clusters), $
                      [ [1,0,0,0], [1,0,0,0], $
                        [1,0,0,0], [1,0,0,0], $
                        [1,0,0,0], [2,1,0,0], $
                        [3,1,1,0], [4,1,2,0], $
                        [5,1,3,0]])
  print, 'all tests passed'
end
