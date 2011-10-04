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
; KEYWORD PARAMETERS:
;  parents: If set, include only the 2 structures which merge to form
;  node. If node is a leaf, then return -1.
;
; OUTPUTS:
;  A vector of the leafward mergers for node. The input node is included in
;  this list. 
;
; MODIFICATION HISTORY:
;  April 2010: Written by Chris Beaumont
;  June 2010: Input node now included in output by default. cnb.
;  October 2010: Added PARENTS keywords. cnb.
;-
function leafward_mergers, node, clusters, parents = parents
  
  compile_opt idl2
  on_error, 2
  
  if n_params() ne 2 then begin
     print, 'calling sequence'
     print, 'result = leafward_mergers(node, clusters, [/parents])'
     return, !values.f_nan
  endif

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

  if node - offset ge n_elements(clusters[0,*]) then return, -1

  if keyword_set(parents) then begin
     if node lt offset then return, -1
     return, clusters[*, node-offset]
  endif

  ;- if node is a leaf, we are done
  if node lt offset then return, node
  result = node

  children = clusters[*,node - offset]
  s = obj_new('stack')
  s->push, children
  while ~s->isEmpty() do begin
     newnode = s->pop()
     result = [result, newnode]
     if newnode lt offset then continue

     ;- sanity check - indices should always be decreasing
     assert, max(clusters[*,newnode-offset]) lt newnode
     s->push, clusters[*, newnode - offset]
  endwhile
  obj_destroy, s
  return, result
end
