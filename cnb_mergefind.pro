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
;  tol: The contour level will be determined to within tol of the true
;       value.
;
; OUTPUTS:
;  An n_kernel x n_kernel array, whose (i,j) element denotes the
;  contour level at which kernel_i merges with kernel_j. 
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
;-
function cnb_mergefind, cube, kernels, $
                        all_neighbors = all_neighbors, $
                        tol = tol

  ;- pick a default tolerance level
  if n_elements(tol) eq 0 then tol = range(cube[kernels]) * 1d-3

  ; CHECK THAT WE HAVE SOME KERNELS TO TRY MERGING.
  kernel_ct = n_elements(kernels)
  if kernel_ct eq 0 then begin
     message, 'WARNING: No kernels to merge!', /con
     return, !values.f_nan
  endif

  ;- initialize the merger matrix
  merger = fltarr(kernel_ct, kernel_ct)+!values.f_nan
  merger[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]
  
  kval = cube[kernels]
  maxvalue = max(kval, /nan)
  minvalue = min(kval, /nan)
  
  ;- discretize the cube data values, to a granularity set by tol
  vals = cube - (cube mod tol)
  levs = vals[uniq(vals, sort(vals))]
  levs = levs[where(finite(levs))]
  ;- get the location, in levs, of the kernels
  kinds = value_locate(levs, vals[kernels])
  
  ;- arrays to hold the bounds of possible contour levels for each
  ;- merger pair.
  ind_lower = lonarr(kernel_ct, kernel_ct)
  ind_lower[indgen(kernel_ct), indgen(kernel_ct)] = kinds
  ind_upper = rebin(kinds, kernel_ct, kernel_ct)
  ind_upper = ind_upper < transpose(ind_upper)

  ;- array to track convergence
  converged = bytarr(kernel_ct, kernel_ct)
  converged[indgen(kernel_ct), indgen(kernel_ct)] = 1
  pbar, name='Mergefind', /new    
  diagonal = converged          ;-1s along the diagonal
  
  repeat begin
     ;- pick a new contour value - find an unmerged pair, 
     ;- and bisect the possible range of contour values for that pair
     todo = where(converged eq 0, todoct)
     if todoct eq 0 then break
     i = todo[0] mod kernel_ct & j = todo[0] / kernel_ct

     ;- some sanity checks, to make sure we're tracking everything correctly
     assert, max(levs[ind_upper] gt cube[kernels]) eq 0
     assert, min(ind_upper - ind_lower) ge 0

     ;- bisect a new search value
     ind = ind_lower[i,j] + (ind_upper[i,j] - ind_lower[i,j]) / 2
     assert, ind gt ind_lower[i,j] && ind lt ind_upper[i,j]
     lev = levs[ind]

     ;- label cube based on "islands" higher than lev
     r = label_region(cube gt lev, all_n = all_neighbors, /ulong)
     id = rebin(r[kernels], kernel_ct, kernel_ct)
     tid = transpose(id)

     ;- lev is an lower limit to kernels on the same island,
     ;- and an upper limit to kernels on different islands
     joined = (id eq tid and (id ne 0) and (tid ne 0) and ~diagonal)
     split = (id ne tid and (id ne 0) and (tid ne 0))
;     assert, min(joined eq transpose(joined)) eq 1
;     assert, min(split eq transpose(split)) eq 1
     
     j = where(joined, jct)
     s = where(split, sct)
     assert, jct + sct ne 0 ;- otherwise, we should have exited
     sanity = ind_upper * 0
     if jct ne 0 then begin
        x = j mod kernel_ct & y = j / kernel_ct
        assert, min(id[x] eq id[y]) eq 1
        sanity[x,y] += ind_lower[x,y] lt ind
        ind_lower[x,y] = ind_lower[x,y] > ind
        ind_lower[y,x] = ind_lower[x,y]
     endif
     if sct ne 0 then begin
        x = s mod kernel_ct & y = s / kernel_ct
        sanity[x,y] += ind_upper[x,y] gt ind
        ind_upper[x,y] = ind_upper[x,y] < ind
        ind_upper[y,x] = ind_upper[x,y]
     endif
     assert, max(sanity) eq 1
     converged = (ind_upper - ind_lower) lt 2

     ;- guess at how much longer we have
     nc = total(converged)
     pbar,  1D * nc / n_elements(converged)
     
  endrep until 0
  pbar, /close

  ;- XXX eliminate seeds which are too small at this point?     
  return, (levs[ind_lower] + levs[ind_upper])/2.
end
