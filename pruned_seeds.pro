;+
; PURPOSE:
;  This function takes a list of to-be-pruned indices (generated,
;  e.g., by generate_prunelist), and computes a new list of
;  kernels (seeds, or local maxima, or new leaves). Optionally,
;  it will create an output file that can be read by the C++
;  utility.
;
; INPUTS:
;  ptr: The dendrogram pointer to be pruned
;  prunelist: The list of structures to prune
;
; KEYWORD PARAMETERS:
;  count: A variable which, on output, will hold the number of new
;         seeds
;  outfile: The name of a file. If present, the new seed list will be
;           written to file. This can be read by the C++ utility.
;  keep: By default, prunelist contains the indices to prune. If /keep
;        is set, prunelist contains the indices to keep.
;
; OUTPUT:
;  A list of the new seeds
;-
function pruned_seeds, ptr, prunelist, count = count, $
                       outfile = outfile, keep = keep

  ls = get_leaves( (*ptr).clusters)
  nleaf = n_elements(ls)
  nst = n_elements( (*ptr).height)
  dokeep = keyword_set(keep)

  keep = replicate(1B, nst)
  if dokeep then begin
     keep *= 0
     keep[prunelist] = 1
  endif else $
     keep[prunelist] = 0

  keep = where(keep, ct)

  result = replicate(0B, nleaf)

  d_index = obj_new('dendro_index', ptr)
  for i = 0, ct - 1, 1 do begin
     index = keep[i]
     lm = d_index->leafward_mergers(index)
     lm = lm[where(lm lt nleaf)]

     ;- if one of the seeds is already included,
     ;- we don't need to do anything
     if max(result[lm]) eq 1 then continue

     ;- otherwise, add brightest leaf
     hi = max( (*ptr).height[lm], hiloc)
     leaf = lm[hiloc]
     result[leaf] = 1B
  endfor
  obj_destroy, d_index

  result = where(result, count)
  if keyword_set(outfile) then begin
     seeds = (*ptr).seeds[result]
     writecol, outfile, seeds, (*ptr).t[seeds], fmt='(i, ",", e)'
  endif

  return, result
end


