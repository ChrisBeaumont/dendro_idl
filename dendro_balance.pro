;+
; PURPOSE:
;  This function computes the "balance" for every node in a
;  dendrogram. The balance is a property of a quantity defined at
;  every node of the dendrogram. At each node, it is defined as the
;  absolute difference between the quantity of that
;  node's children. Leaves have balances of zero by definition.
;
;  If no quantity is provided via the property keyword, then
;  the "leaf balance" is computed. The leaf balance is the difference
;  between the number of leaves contained within each child node.
;
; INPUTS:
;  clusters: The cluster matrix describing the dendrogram. A
;  [2,nleaf-1] array.
;
; KEYWORD PARAMETERS:
;  property: The property for which to compute the balance. This is an
;  n-structure array. 
;  normalize: Set this keyword to normalize each balance value
;             (balance = |child1 - child2|/(child1+child2)
;
; OUTPUTS:
;  An [n-structure] array holding the balance for every structure
;
; SEE ALSO:
;  accumulate_dendrogram
;
; MODIFICATION HISTORY:
;  November 5 2010: Written by Chris Beaumont
;  November 8 2010: Generalized beyond leaf balances. Added property
;                   keyword. cnb.
;-

;+
; Recursive driver for dendro_balance
;-
pro _balance_recurse, balance, property, clusters, id, normalize = normalize
  nleaf = (size(clusters))[2] + 1
  isLeaf = id lt nleaf
  if isLeaf then return
  leafwards = clusters[*, id - nleaf]
  _balance_recurse, balance, property, clusters, leafwards[0], normalize = normalize
  _balance_recurse, balance, property, clusters, leafwards[1], normalize = normalize
  balance[id] = range(property[leafwards])
  if keyword_set(normalize) then balance[id] /= 1. * total(property[leafwards])
  if ~finite(balance[id]) then balance[id] = 0
end
  
function dendro_balance, clusters, property = property, normalize = normalize

  nst = max(clusters) + 2
  nleaf = n_elements(clusters[*,0]) + 1
  
  if ~keyword_set(property) then begin
     property = intarr(nst)
     property[get_leaves(clusters)] = 1
     property = accumulate_dendro(property, clusters)
  endif
  balance = property * 0.
  _balance_recurse, balance, property, clusters, nst-1, normalize = normalize
  return, balance
end

pro test
  
  ;- perfectly balanced tree with 4 leaves
  clusters=[[0,1],[2,3],[4,5]]
  assert, array_equal(dendro_balance(clusters), [0,0,0,0,0,0,0])

  ;-fully lopsided tree with 5 leaves
  clusters=[[0,1],[5,2],[6,3],[7,4]]
  assert, array_equal(dendro_balance(clusters), [0,0,0,0,0,0,1,2,3])

  ;- property value
  clusters=[[0,1],[2,3],[4,5]]
  prop = [10,5,1,1,0,0,0]
  assert, array_equal(dendro_balance(clusters, p=prop), [0,0,0,0,5,0,0])


  clusters=[[0,1],[2,3],[4,5]]
  prop = [10,5,1,1,0,0,0]
  b = dendro_balance(clusters, p=prop,/norm)
  guess = [0,0,0,0,5./15., 0,0]
  assert, min(abs(guess - b) lt 1e-3)


  
end

  
