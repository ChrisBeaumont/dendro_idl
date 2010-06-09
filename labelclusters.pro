;+
; PURPOSE:
;  This function labels each pixel in an image/cube according to the
;  highest dendrogram structure to which it belongs.
;
;
; MODIFICATION HISTORY:
;  <June 2010: Written by Erik Rosolowsky
;  June 2010: Got rid of the annoying "32767" label for masked out
;             pixels. They are now -1. Chris Beaumont
;-
function labelclusters, height, clusters, decimkern, levelsin, $
  x, y, v, t, szin, $
  all_neighbors = all_neighbors

  sz = szin
  if sz[0] eq 2 then sz[3] = 1
  cubify, x, y, v, t, sz, c = c
  x0 = decimkern mod sz[1]
  y0 = (decimkern mod (sz[1]*sz[2]))/sz[1]
  v0 = decimkern/(sz[1]*sz[2])
  clusterlabel = intarr(sz[1], sz[2], sz[3]) - 1
  levels = reverse(levelsin[sort(levelsin)])
; We need an array to tell secondary clusters (not leaves) where to
; find a point within them.  This initializes that array
  usekernels = intarr(max(clusters)+1)+max(clusters)+2

  for k = 0, max(clusters) do begin
     print, k
; Find all parent clusters
     roots = cluster_member(k, clusters)
     print, roots
; Tell every parent cluster to use the smallest index kernel within it
; as a seed point.
     usekernels[roots] = (min(roots))+intarr(n_elements(roots)) $
                         < usekernels[roots]
     mrglevels = height[roots[1]]
; For this cluster, find the lowest level in it that's above a merge
     sup = max(where(levels gt float(mrglevels), ct))
     if ct eq 0 then sup = 0
     above_mrg = levels[sup]
; Contour above this level
     mask = label_region(c ge above_mrg, all_neighbors = all_neighbors, /ulong)
; Reject regions that don't contain the kernel
     mask = mask eq mask[x0[usekernels[k]], y0[usekernels[k]], v0[usekernels[k]]] ; Label every point in this with the smallest value present (k or a
; preassigned value)
     clusterlabel[where(mask)] = clusterlabel[where(mask)] < k
  endfor
;  clusterlabel[where(c lt min(height))] = max(clusters)+1
; Send back a vector.
  return, clusterlabel[x, y, v]
end



