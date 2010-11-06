;+
; PURPOSE:
;  This function computes the "leaf balance" for every node in a
;  dendrogram. The leaf balance is defined as the absolute difference
;  in the number of leaves between a node's left and right branches. 
;
; INPUTS:
;  clusters: The cluster matrix describing the dendrogram. A
;  [2,nleaf-1] array.
;
; KEYWORD PARAMETERS:
;  leafcount: On output, this variable will hold the number of leaves
;             attached to each node.
;
; OUTPUTS:
;  A [2, 2*nleaf-1] array holding the leaf balance for every structure
;
; MODIFICATION HISTORY:
;  November 5 2010: Written by Chris Beaumont
;-

;+
; Recursive driver for dendro_balance
;-
pro _balance_recurse, result, leaves, clusters, id
  nleaf = (size(clusters))[2] + 1
  isLeaf = id lt nleaf
  if isLeaf then begin
     leaves[id] = 1
     return
  endif
  
  leafwards = clusters[*, id - nleaf]
  _balance_recurse, result, leaves, clusters, leafwards[0]
  _balance_recurse, result, leaves, clusters, leafwards[1]
  
  leaves[id] = total(leaves[leafwards])
  result[id] = abs(leaves[leafwards[0]] - leaves[leafwards[1]])
end
  
function dendro_balance, clusters, leafcount = leafcount

  nst = max(clusters) + 2
  nleaf = n_elements(clusters[*,0]) + 1
  
  leaves = intarr(nst)
  result = leaves
  
  _balance_recurse, result, leaves, clusters, nst-1
  leafcount = leaves
  return, result
end

pro test
  
  ;- perfectly balanced tree with 4 leaves
  clusters=[[0,1],[2,3],[4,5]]
  assert, array_equal(dendro_balance(clusters, leafcount = l), [0,0,0,0,0,0,0])
  assert, array_equal(l, [1,1,1,1,2,2,4])


  ;-fully lopsided tree with 5 leaves
  clusters=[[0,1],[5,2],[6,3],[7,4]]
  assert, array_equal(dendro_balance(clusters, l = l), [0,0,0,0,0,0,1,2,3])
  assert, array_equal(l, [1,1,1,1,1,2,3,4,5])
end

  
