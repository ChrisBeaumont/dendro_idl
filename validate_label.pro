;+
; PURPOSE:
;  This procedure ensures that the labelcube created by labelclusters is
;  sensible. It enforces a bunch of assertions, and will halt if any fail.
;
; INPUTS:
;  cube: The data cube
;  kernels: Indices of the kernels
;  labels: The label cube
;  clusters: The clusters array
;
; MODIFICATION HISTORY:
;  July 22 2010: Written by Chris Beaumont.
;-
pro validate_label, cube, kernels, labels, clusters, all_neighbors = all_neighbors

  nst = max(clusters)
  nleaf = n_elements(kernels)
  seeds = get_seed(clusters)
  nonneg = (labels ge 0)
  pbar, 'validate_label', /new
  for i = 0, max(labels), 1 do begin
     pbar, 1. * i / max(labels)

     if i lt nleaf then seed = i else seed = seeds[0,i-nleaf]

     ;- this mask contains structure i and all (leafward) substructs
     mask = labels le i and nonneg
     assert, mask[kernels[seed]]
     r = label_region(mask, all_neighbors = all_neighbors, /ulong)
     mask = r eq r[kernels[seed]]

     ;- this region _must_ contain the alread-merged kernels 
     if i lt nleaf then begin
        assert, mask[kernels[seed]]
     endif else begin
        assert, mask[kernels[seeds[0,i-nleaf]]]
        assert, mask[kernels[seeds[1,i-nleaf]]]
     endelse

     ;- this region _cannot_ contain the next merger seed
     next = where(clusters eq i)
     next_row = next / 2 & next_col = next mod 2
     sibling = clusters[~next_col, next_row]
     assert, sibling ne i
     assert, ~mask[kernels[seeds[~next_col, next_row]]]
  endfor
  pbar, /close
end
