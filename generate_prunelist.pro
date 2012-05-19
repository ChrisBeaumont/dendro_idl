;+
; PURPOSE:
;  Generate a list of dendrogram indices that should be pruned, based
;  on 1 or more criteria
;
; INPUTS:
;  ptr: A dendrogram pointer to prune
;  count: On output, will hold the number of pruned indices
;
; KEYWORD PARAMETERS:
;  delta: To avoid pruning, the brightest pixel in a structure must be
;         at least delta above the next merger contour
;  npix: To avoid pruning, a structure must contain at least npix
;        pixels
;  minflux: To avoid pruning, the summed intensity must be >= minflux
;  minpeak: To avoid pruning, the brightest pixel must be >= minpeak
;
;  If any of these are missing, the criteria is not enforced.
;
;  outifle: If present, it gives the name of a file to write a new
;  'seeds' file to. This file, when passed tot he C++ utility via the
;  '-k' argument, will re-generate a new dendrogram with the
;  appropriate structure removed.
;
; OUTPUTS:
;  A list of indices which should be pruned.
;-
function generate_prunelist, ptr, count, delta = delta, npix = npix, $
                             minflux = minflux, minpeak = minpeak, $
                             outfile = outfile

  if n_params() ne 2 then begin
     print, 'calling sequence'
     print, 'result = generate_prunelist(ptr, count, '
     print, '                           (delta = delta, npix = npix)'
     return, !values.f_nan
  endif

  count = 0
  if n_elements(delta) eq 0 && n_elements(npix) eq 0 then $
     return, -1

  if ~keyword_set(delta) then delta = 0
  if ~keyword_set(npix) then npix = 0
  if ~keyword_set(minflux) then minflux = -!values.f_infinity
  if ~keyword_set(minpeak) then minpeak = -!values.f_infinity

  nst = n_elements((*ptr).height)
  kill = bytarr(nst)
  visited = bytarr(nst)

  d_index = obj_new('dendro_index', ptr)

  l = get_leaves((*ptr).clusters)
  ;- for each leaf, travel rootward until prune test fails,
  ;- or we re-visit a node
  for i = 0, n_elements(l) - 1, 1 do begin
     index = l[i]
     h = (*ptr).height[index]
     repeat begin
        if visited[index] then break
        visited[index] = 1B
        lm = d_index->leafward_mergers(index)
        p = merger_partner(index, (*ptr).clusters, merge = m)
        if p eq -1 then break ;- don't prune the root
        np = total( (*ptr).cluster_label_h[lm] )
        de = (h - (*ptr).height[m])
        ind = substruct(index, ptr)
        dokill = np lt npix || de lt delta || $
                 total((*ptr).t[ind]) lt minflux || $
                 max( (*ptr).t[ind] ) lt minpeak
        kill[index] = dokill
        index = m
     endrep until dokill eq 0
     ;- all structures rootward of a non-pruned structure
     ;- are implicitly safe from pruning
     visited[rootward_mergers(l[i], (*ptr).clusters)] = 1B
  endfor
  obj_destroy, d_index

  result = where(kill, count)

  if keyword_set(outfile) then $
     junk = pruned_seeds(ptr, result, outfile=outfile)

  return, result
end
