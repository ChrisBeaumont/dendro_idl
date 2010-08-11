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
;  July 2010: Input is now the cube, instead of the vectorized
;             version. simpler to understand. cnb.
;-
function labelclusters, height, clusters, decimkern, levelsin, $
                        cube, $
                        all_neighbors = all_neighbors, fast = fast, $
                        contour_res = contour_res
                        

  ;- just making sure that the kernel locations have been 
  ;- calculated correctly, given all the vectorifying/cubifying
  nkern = n_elements(decimkern)
  if keyword_set(fast) then assert, max(abs(cube[decimkern] - height[0:nkern-1])) lt 1e-4
  sz = size(cube)
  if size(cube, /n_dim) eq 2 then sz[3] = 1
  badval = 32767

  clusterlabel = intarr(sz[1], sz[2], sz[3]) + badval
  if ~keyword_set(fast) then levels = reverse(levelsin[sort(levelsin)])
  
  ; We need an array to tell secondary clusters (not leaves) where to
  ; find a point within them.  This initializes that array
  usekernels = intarr(max(clusters)+1)+max(clusters)+2
  merge_array = fltarr(max(clusters)+1)

  for k = 0, max(clusters) do begin
     
     ; Find all parent clusters
     roots = cluster_member(k, clusters)
     
     ; Tell every parent cluster to use the smallest index kernel within it
     ; as a seed point.
     usekernels[roots] <= min(roots)

     hi = max(height[roots], loc)
     assert, roots[0] eq k
     assert, loc eq 0

     mrglevels = height[roots[1]]

     ;- want to find the lowest level that's above the next merger
     if keyword_set(fast) then begin
        ;- cnb_mergefind is meant to return an upper limit to the merge 
        ;- levels, so mrglevels is just above the next merger
        above_mrg = mrglevels
     endif else begin
        ;- find the lowest level above mrglevels
        sup = max(where(levels gt float(mrglevels), ct))
        if ct eq 0 then sup = 0
        above_mrg = levels[sup]
     endelse
     merge_array[k] = above_mrg

     ; Contour above this level
     mask = label_region(cube gt above_mrg, all_neighbors = all_neighbors, /ulong)

     ; Reject regions that don't contain the kernel
     mask = mask eq mask[decimkern[usekernels[k]]] ; Label every point in this with the smallest value present (k or a
; preassigned value)

     clusterlabel[where(mask)] <= k
  endfor

  ;- mask out the unused values to -1
  bad = where(clusterlabel eq badval, ct)
  if ct ne 0 then clusterlabel[bad] = -1

  return, clusterlabel
end



