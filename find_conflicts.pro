;+
; PURPOSE:
;  this code is a recursive helper routine for find_conflicts. It does
;  the actual work of conflict detection, for a single node. It calls
;  itself to traverse the dendrogram.
;
; INPUTS:
;  index: The ID of the current dendrogram structure to examine
;  clusters: The clusters array for the dendrogram
;  height: The height array for the dendrogram
;  conflicts: A stack containing the conflicts so far. New conflicts
;  will be pushed onto this stack.
;
; MODIFICATION HISTORY:
;  July 2010: Written by Chris Beaumont
;-
pro recursive_conflict, index, clusters, height, conflicts
  nleaf = n_elements(clusters[0,*]) + 1

  ;- if index is a leaf, we are done
  if index lt nleaf then return

  ;- find the two structures that merger to form cluster_index
  ;- find their heights
  leafward = clusters[*, index - nleaf]
  lev = height[index]
  lev2 = height[leafward]

  ;- leafward level equals current level - a conflict
  if lev2[0] eq lev then begin
     assert, leafward[0] ge nleaf
     leaf2 = clusters[*, leafward[0] - nleaf]
     conflicts->push, [leafward[1], leaf2[0], leaf2[1]]
  endif
  if lev2[1] eq lev then begin
     assert, leafward[1] ge nleaf
     leaf2 = clusters[*, leafward[1] - nleaf]
     conflicts->push, [leafward[0], leaf2[0], leaf2[1]]
  endif

  ;- recurse on two leafwards
  recursive_conflict, leafward[0], clusters, height, conflicts
  recursive_conflict, leafward[1], clusters, height, conflicts
  return
end
  
;+
; PURPOSE:
;  This function locates any non-binary mergers in a merger
;  matrix. Non-binary splits are found in the following way: first, a
;  dendrogram is generated from the merger matrix. The height of each
;  component is compared to the height of each of its child (leafward)
;  substructures. If these two heights are equal, this indicates a
;  three-way split between one "child" and two "grandchildren" of the
;  original node. 
;
; INPUTS:
;  merger: The merger matrix
;
; OUTPUTS:
;  A [3, nconflict] array of conflicts. Each row in this array gives
;  the location of three seed locations (rows/columns in the symmeric
;  merger matrix) which produce a non-binary split in the
;  dendrogram. merger[result[0,i], result[1,i]] is equal to
;  merger[result[1,i], result[2,i]] (and so on for other
;  perumutations). The merger matrix must be refined so that these
;  contours are no longer equal. This will resolve the conflict.
;
; MODIFICATION HISTORY:
;  July 2010: Written by Chris Beaumont
;-
function find_conflicts, merger                             

  ;- make the dendrogram
  generate_dendrogram, merger, $
                       clusters = clusters, height = height, $
                       xlocation = xlocation, leafnodes = leafnodes
  
  ;- calculate the ID of the root
  root = max(clusters)+1

  ;- recursively find conflicts
  conflicts = obj_new('stack')
  recursive_conflict, root, clusters, height, conflicts
  result = conflicts->toArray()
  bad = conflicts->getSize() eq 0
  obj_destroy, conflicts
  
  if bad then return, -1 
  for i = 0, n_elements(result) -1 do result[i] = min(leafward_mergers(result[i], clusters))
  return, reform(result, 3, n_elements(result) / 3)


  ;- sanity check code
  ;- confirm assertions for all conflicts
  for i = 0, n_elements(result[0,*])-1 do begin
     l1 = leafward_mergers(result[0,i], clusters)
     l2 = leafward_mergers(result[1,i], clusters)
     l3 = leafward_mergers(result[2,i], clusters)
     assert, ~finite(intersection(l1, l2))
     assert, ~finite(intersection(l3, l2))
     assert, ~finite(intersection(l1, l3))
     assert, merger[min(l1), min(l2)] eq $
             merger[min(l1), min(l3)]
     assert, merger[min(l1), min(l2)] eq $
             merger[min(l2), min(l3)]
     assert, merger[min(l1), min(l3)] eq $
             merger[min(l2), min(l3)]
  endfor

  print, 'all tests passed'
end
