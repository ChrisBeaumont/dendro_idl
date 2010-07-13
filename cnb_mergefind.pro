pro mergefind_step, value, lower, upper, cube, kernels, diagonal
;  assert, max(upper gt cube[kernels]) eq 0
  assert, min(upper - lower) ge 0
  
  
  ;- label cube based on "islands" higher than value
  r = label_region(cube gt value, all_n = all_neighbors, /ulong)
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
;               value.
;
; OUTPUTS:
;  An n_kernel x n_kernel array, whose (i,j) element denotes the
;  contour level at which kernel_i merges with kernel_j. cnb_mergefind returns
;  the upper bound for convenience later. Structures i and j merge somewhere
;  between outputs[i,j] - contour_res and outputs[i,j]. 
;
; PROCEDURE:
;  We use a binary search to identify contour levels. While there are
;  (n_kernel x (n_kernel - 1))/2 contour levels to find, there are far
;  fewer unique values. cnb_mergefind takes advantage of this, and
;  searches for multiple contour values simultaneously.
;
; MODIFICATION HISTORY:
;  Feb 2010: Written by Chris Beaumont. Adapted from mergefind.pro,
;            written by Erik Rosolowsky
;  June 2010: Added contour_res keyword. Added code to resolve multi-way
;  mergers. cnb. 
;  July 2010: Now returning the upper limits for each merger. I think
;  this is the most useful quantity. cnb.
;-
function cnb_mergefind, cube, kernels, $
                        all_neighbors = all_neighbors, $
                        contour_res = contour_res
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
  merger = dblarr(kernel_ct, kernel_ct)+!values.f_nan
  merger[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]
      
  ;- arrays to hold the bounds of possible contour levels for each
  ;- merger pair. Each merger lies somewhere in between 0 and min(kern1, kern2)
  lower = dblarr(kernel_ct, kernel_ct)
  lower[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]
  upper = rebin(cube[kernels], kernel_ct, kernel_ct)
  upper = upper < transpose(upper)

  ;- array to track convergence
  converged = bytarr(kernel_ct, kernel_ct)
  converged[indgen(kernel_ct), indgen(kernel_ct)] = 1
  pbar, name='Mergefind', /new    
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
     mergefind_step, value, lower, upper, cube, kernels, diagonal

     ;- check for convergence
     converged = (upper - lower) lt contour_res

     ;- guess at how much longer we have
     nc = total(converged)
     pbar,  1D * nc / n_elements(converged)
     
  endrep until 0
  pbar, /close

  ;- second step -- refine merger to resolve non-binary mergers
  conflicts = find_conflicts(upper)
  while conflicts[0] ne -1 do begin
     nconflict = n_elements(conflicts[0,*])
     for i = 0, nconflict - 1, 1 do begin
        print, i, nconflict
        s1 = conflicts[0,i]
        s2 = conflicts[1,i]
        s3 = conflicts[2,i]
        ;- resolve the conflict
        while upper[s1, s2] eq upper[s1, s3] && $
           upper[s1, s2] eq upper[s2, s3] do begin
           
           range = (upper - lower)[[s1, s1, s2], [s2, s3, s3]]
           value = ((upper + lower)/2.)[[s1, s1, s2], [s2, s3, s3]]
           assert, n_elements(range) eq 3
           top = max(range, loc)
           value = value[loc]
           mergefind_step, value, lower, upper, cube, kernels, diagonal
        endwhile
     endfor
     ;- recalculate conflicts
     conflicts = find_conflicts(upper)
  endwhile
  
  ;- maybe do one more level of refinement: generate dendros 
  ;- on lower and upper, assert that they are the same?
  ;- then, the structure of the dendrogram is unique.
  ;- how long would that take???


  ;- XXX eliminate seeds which are too small at this point?     


  ;- the upper limit is the most convenient. Once the dendrogram
  ;- is created, these contour levels are guaranteed to split
  ;- apart the two components of each merger. This makes
  ;- cluster labeling easier.
  return, upper
end
