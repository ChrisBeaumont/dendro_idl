;+
; PURPOSE:
;  This function computes the height for each structure in a
;  dendrogram. The height of a structure is defined as the number of
;  nodes visited when traveling from the dendrogram root to the
;  structure. The root has a height of zero
;
; INPUTS:
;  Cluster: the cluster array describing the dendrogram
;
; OUTPUTS:
;  An n-structure element array giving the height of each node.
;
; MODIFICATION HISTORY:
;  November 2010; Written by Chris Beaumont
;-
function cluster_height, cluster

  nst = max(cluster) + 2
  offset = n_elements(cluster[0,*]) + 1
  todo = obj_new('stack')
  
  result = intarr(nst)
  rec={height_rec, id:max(cluster)+1, depth:0}
  todo->push, rec

  while ~todo->isEmpty() do begin
     rec = todo->pop()
     result[rec.id] = rec.depth
     
     ;- 2 children
     isLeaf = rec.id lt offset
     if isLeaf then continue
     children = cluster[*, rec.id - offset]
     todo->push, {height_rec, id:children[0], depth:rec.depth+1}
     todo->push, {height_rec, id:children[1], depth:rec.depth+1}
  endwhile
  obj_destroy, todo
  return, result
end

pro test
  cluster = [[0,1],[2,3],[4,5]]
  assert, array_equal(cluster_height(cluster), [2,2,2,2,1,1,0])

  cluster = [[0,1],[5,2],[6,3],[7,4]]
  assert, array_equal(cluster_height(cluster), [4,4,3,2,1,3,2,1,0])
end
