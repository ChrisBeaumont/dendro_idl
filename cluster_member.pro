;+
; PURPOSE:
;  This function follows a dendrogram structure down (rootward) the
;  dendrogram, and returns a list of all the structures
;  encountered. In other words, it returns a list of structures for
;  which a given node is a nested substructure.
;
; INPUTS:
;  leafnode: The ID of the structure to examine
;  clusters: The cluster matrix describing the dendrogram
;
; OUTPUTS:
;  A list of clusters for which leafnode is a member / substructure
;
; MODIFICATION HISTORY:
;  <July 2010: Written by Erik Rosolowsky
;  July 2010: Documented by Chris Beaumont.
function cluster_member, leafnode, clusters
  compile_opt idl2

; find the number all clusters that LEAFNODE belongs to.

  roots = [leafnode]
  sz = size(clusters)
  node = leafnode
  repeat begin
    ind = where(clusters eq node)
    clnum = (ind/2)+sz[2]+1
    roots = [roots, clnum]
    node = clnum[0]
  endrep until (clnum eq 2*sz[2])
  
  return, roots
end
