;+
; PURPOSE:
;  This function calculates which structures (leaves and nodes) in a
;  dendrogram lie above (i.e., leafward) of a given node, along the
;  same branch. 
;
; INPUTS:
;  node: The node for which to find the leafward mergers (can also be
;        a leaf)
;  clusters: The cluster array returned by cluster_tree
;
; OUTPUTS:
;  A vector of the leafward mergers for node. If node is a leaf, it
;  has no leafward mergers and node is returned by itself. Otherwise,
;  node is not included in the output list.
;
; MODIFICATION HISTORY:
;  April 2010: Written by Chris Beaumont
;-
function leafward_mergers, node, clusters

  compile_opt idl2

  ;- clusters is a [2, n_leaf-1] array,
  ;- the i'th row of which lists the index
  ;- of the 2 structures (leaves or nodes) that
  ;- merge into the n_leaf+i'th node

  ;- so, we look up which 2 structures merge
  ;- to make the requested node - these are
  ;- two of the leafward mergers. We then
  ;- recursively find each of those structures'
  ;- leafward mergers, and return the whole set.

  sz = size(clusters)
  offset = sz[2] + 1

  ;- if node is a leaf, we are done
  if node lt offset then return, node

  children = clusters[*,node - offset]
  s = obj_new('stack')
  s->push, children
  while ~s->isEmpty() do begin
     newnode = s->pop()
     if n_elements(result) eq 0 then result = newnode $
     else result=[result, newnode]
     if newnode lt offset then continue

     ;- sanity check - indices should always be decreasing
     assert, max(clusters[*,newnode-offset]) lt newnode

     s->push, clusters[*, newnode - offset]
  endwhile
  obj_destroy, s
  return, result
end
