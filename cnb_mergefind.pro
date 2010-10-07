pro mergefind_step, value, lower, upper, npix, cube, kernels, diagonal
  assert, n_params() eq 7
;  assert, max(upper gt cube[kernels]) eq 0
  assert, min(upper - lower) ge 0
  
  
  ;- label cube based on "islands" higher than value
  r = label_region(cube gt value, all_n = all_neighbors, /ulong)
  h = histogram(r, min = 0)
  kernel_ct = n_elements(kernels)
  id = rebin(r[kernels], kernel_ct, kernel_ct)
  tid = transpose(id)
  
  ;- value is a lower limit to kernels on the same island,
  ;- and an upper limit to kernels on different islands
  joined = (id eq tid and (id ne 0) and (tid ne 0) and ~diagonal)
  split = (id ne tid and (id ne 0) and (tid ne 0))
  
  j = where(joined, jct)
  s = where(split, sct)
  assert, jct + sct ne 0
  sanity = upper * 0
  if jct ne 0 then begin
     x = j mod kernel_ct & y = j / kernel_ct
     assert, min(id[x] eq id[y]) eq 1
     npix[x, y] = npix[x, y] < h[id[x]]
     npix[y,x] = npix[x,y]
     sanity[x,y] += lower[x,y] lt value
     lower[x,y] = lower[x,y] > value
     lower[y,x] = lower[x,y]
  endif
  if sct ne 0 then begin
     x = s mod kernel_ct & y = s / kernel_ct
     sanity[x,y] += upper[x,y] gt value
     upper[x,y] = upper[x,y] < value
     upper[y,x] = upper[x,y]
  endif
  assert, max(sanity) eq 1
end


;+
; PURPOSE:
;  This function computes the contour levels at which
;  different seed locations in a data cube merge with one another.
;
; INPUTS:
;  cube: The data. A 2- or 3D array
;  kernels: The (1D) indices of the kernel locations (dendrogram
;           leaves)
;
; KEYWORD PARAMETERS:
;  ALL_NEIGHBORS: Set to include 26, rather than 6, neighbors for each
;                 pixel. 
;  contour_res: The contour level will be determined to within this much of the true
;               value. Note that this is just a first guess, and the
;               precision is increased if necessary (see blow).
;  npix: on output, will return an (nkernel x nkernel) array, listing
;        the number of pixels in the region defined by the contour
;        which merges kernels i and j. The diagonals of this array are
;        filled with bogus, large numbers. Used by decimate_merger.
;
;
; OUTPUTS:
;  An n_kernel x n_kernel array, whose (i,j) element denotes the
;  contour level at which kernel_i merges with kernel_j. cnb_mergefind returns
;  the upper bound for convenience later. Structures i and j merge somewhere
;  between outputs[i,j] - contour_res and outputs[i,j]. 
;
; PROCEDURE:
;  We repeatedly contour the data, and keep track of the lower and
;  upper limits for the contour values which critically separate each
;  pair of kernels. A pair's 'critical contour' encircles each
;  point in a figure 8 pattern, such that lower contours enclose the
;  pair in a single structure, and higher contours split them apart.
;
;  At each iteration, a new test contour is chosen to refine the
;  merger matrix. This happens during 3 phases:
;   - phase 1 is completed when each critical contour is determined to
;     within a resolution specified by the contour_res keyword. At
;     each iteration, a not-yet-converged pair is chosen, and the
;     test contour level bisects the current lower and upper limits
;     for this pair. 
;   - phase 2 refines the merger matrix such that there are no
;     non-binary splits when the final dendrogram is created. This
;     phase creates a dendrogram from the current merger matrix, finds
;     non-binary splits, identifies which kernel pairs are
;     responsible, and re-determines these contour levels at higher
;     precision. The process repeats until no conflicts remain
;   - phase 3 provides the final refinement. It ensures that the
;     uncertainties in merger levels are sufficiently small that the
;     merger order of each kernel pair (and the structure of the resulting
;     dendrogram) is unambiguously determined.
;
;  Note that, while the test contour level at each iteration is
;  targeted to refine a single merger pair, it often tightens the
;  bounds on many more.
;
; MODIFICATION HISTORY:
;  Feb 2010: Written by Chris Beaumont. Adapted from mergefind.pro,
;            written by Erik Rosolowsky
;  June 2010: Added contour_res keyword. Added code to resolve multi-way
;  mergers. cnb. 
;  July 2010: Now returning the upper limits for each merger. I think
;  this is the most useful quantity. cnb.
;  July 19 2010: Added npix keyword, and code to ensure
;                the dendrogram is uniquely determined. cnb.
;-
function cnb_mergefind, cube, kernels, $
                        all_neighbors = all_neighbors, $
                        contour_res = contour_res, $
                        npix = npix
  compile_opt idl2
  
  ;- pick a default resolution
  if ~keyword_set(contour_res) then tol = range(cube[kernels]) * 1d-3 $
  else tol = contour_res

  ; CHECK THAT WE HAVE SOME KERNELS TO TRY MERGING.
  kernel_ct = n_elements(kernels)
  if kernel_ct eq 0 then begin
     message, 'No kernels to merge!'
  endif

  ;- initialize the merger matrix
;  merger = dblarr(kernel_ct, kernel_ct)+!values.f_nan
;  merger[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]

  ;- initialize npix to large numers. off-diagonal elements 
  ;- will be updated.
  npix = lonarr(kernel_ct, kernel_ct) + n_elements(cube)
      
  ;- arrays to hold the bounds of possible contour levels for each
  ;- merger pair. Each merger lies somewhere in between 0 and min(kern1, kern2)
  lower = dblarr(kernel_ct, kernel_ct)
  lower[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]
  upper = rebin(cube[kernels], kernel_ct, kernel_ct)
  upper = upper < transpose(upper)
  upper = double(upper)

  ;- array to track convergence
  converged = bytarr(kernel_ct, kernel_ct)
  converged[indgen(kernel_ct), indgen(kernel_ct)] = 1
  pbar, name='Mergefind 1', /new    
  diagonal = converged          ;-1s along the diagonal

  ;- first step: find mergers to a precision of at least contour_res
  repeat begin
     ;- pick a new contour value - find an unmerged pair, 
     ;- and bisect the possible range of contour values for that pair
     todo = where(converged eq 0, todoct)
     if todoct eq 0 then break
     i = todo[0] mod kernel_ct & j = todo[0] / kernel_ct

     ;- some sanity checks, to make sure we're tracking everything correctly
     assert, max(upper gt cube[kernels]) eq 0
     assert, min(upper - lower) ge 0

     ;- bisect a new search value
     value = lower[i,j] + (upper[i,j] - lower[i,j]) / 2
     assert, value gt lower[i,j] && value lt upper[i,j]
     
     ;- tighten the bounds, based on this test value
     mergefind_step, value, lower, upper, npix, cube, kernels, diagonal

     ;- check for convergence
     converged = (upper - lower) lt contour_res

     ;- guess at how much longer we have
     nc = total(converged)
     pbar,  1D * nc / n_elements(converged)
     
  endrep until 0
  pbar, /close

  ;- second step -- refine merger to resolve non-binary mergers
  do2 = 0
  if do2 then conflicts = find_conflicts(upper, count = con_ct)
  while do2 && con_ct ne 0 do begin
     nconflict = con_ct
     pbar, name='Mergefind 2 -'+strtrim(con_ct,2)+' conflicts', /new    
     for i = 0, nconflict - 1, 1 do begin
        pbar, 1. * i / nconflict
        s1 = conflicts[0,i]
        s2 = conflicts[1,i]
        s3 = conflicts[2,i]
        ;- refine until conflict is resolved
        ;- resolved when: merger is not at the same level,
        ;- or the merger level is known exactly
        while upper[s1, s2] eq upper[s1, s3] && $
           upper[s1, s2] eq upper[s2, s3] && $
           ((upper - lower)[s1,s2] gt 0 || $
            (upper - lower)[s1, s3] gt 0 || $
            (upper - lower)[s2,s3] gt 0) do begin
           
           ;- find the lest-well determined pair, and bisect
           range = (upper - lower)[[s1, s1, s2], [s2, s3, s3]]
           value = ((upper + lower)/2.)[[s1, s1, s2], [s2, s3, s3]]
           assert, n_elements(range) eq 3
           top = max(range, loc)
           value = value[loc]
           mergefind_step, value, lower, upper, npix, cube, kernels, diagonal
        endwhile
     endfor
     pbar, /close
     ;- recalculate conflicts
     conflicts = find_conflicts(upper, count = con_ct)
  endwhile
  
  ;- third level of refinement (which, now that I think about it, 
  ;- would automatically fix the 2nd level. o well). Find any 
  ;- ambiguous mergers, where the order of merging is not determined
  ;- given the uncertainties between lower and upper. Refine these 
  ;- levels, so that the dendrogram structure is uniquely determined
  nelem = n_elements(lower)
  conflicts = find_ambiguities(lower, upper, count = con_ct)
  while con_ct ne 0 do begin
     nconflict = con_ct
     print, con_ct
     pbar, name='Mergefind 3 -'+strtrim(con_ct,2)+' conflicts', /new    
     for i = 0, nconflict - 1, 1 do begin
        pbar, 1. * i / nconflict
        s1 = conflicts[0,i]
        s2 = conflicts[1,i]
        s3 = conflicts[2,i]
        ;- resolve the conflict
        while ambiguous_triplet(conflicts[*,i], lower, upper) do begin
           ;- find the least-constrained merger, and bisect it
           range = (upper - lower)[[s1, s1, s2], [s2, s3, s3]]
           value = ((upper + lower)/2.)[[s1, s1, s2], [s2, s3, s3]]
           top = max(range, loc)
           value = value[loc]
           mergefind_step, value, lower, upper, npix, cube, kernels, diagonal
        endwhile
     endfor
     pbar, /close
     ;- recalculate conflicts, if we've introduced any more 
     ;- (note that I don't think this is possible)
     conflicts = find_ambiguities(lower, upper, count = con_ct)
  endwhile
    
  ;- the upper limit is the most convenient. Once the dendrogram
  ;- is created, these contour levels are guaranteed to split
  ;- apart the two components of each merger. This makes
  ;- cluster labeling easier.
  return, upper
end
