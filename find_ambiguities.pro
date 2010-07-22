;+
; PURPOSE:
;  Given the lower and upper bounds on a merger matrix, this function
;  returns a list of kernel triplets whose merger order is not
;  determined. When there are no more ambiguous kernel triplets, the
;  structure of the dendrogram is completely determined.
;
; INPUTS:
;  lower: The lower bound for the merger matrix
;  upper: The upper bound for the merger matrix
;
; KEYWORD PARAMETERS:
;  count: The number of ambiguous triplets
;
; OUTPUTS:
;  A (3,n) list of ambiguous kernel triplets. If the entries in one of
;  these rows are i,j,k, then the true order of merger_ij, merger_jk,
;  and merger_ik are not well determined. 
;
; MODIFICATION HISTORY:
;  July 21 2010: Written by Chris Beaumont
;-
function find_ambiguities, lower, upper, count = count
  compile_opt idl2
  
  if n_params() ne 2 then begin
     print, 'calling sequence:'
     print, ' result = find_ambiguities(lower, upper, [count = count])'
     return, !values.f_nan
  endif

  generate_dendrogram, upper, $
                       clusters = cluster, $
                       height = height

  seeds = get_seed(cluster)
  nleaf = n_elements(lower[0,*])
  sz = size(cluster) & nrow = sz[2]
  result = obj_new('stack')

  ;- examine each structure in the dendro
  for i = 0, nrow - 2, 1 do begin
     s1 = seeds[0,i] & s2 = seeds[1,i]
     l1 = lower[s1, s2] & u1 = upper[s1, s2]
     l1 = l1[0] & u1 = u1[0]

     ;-compare this structure to all the remaining 
     ;- structures
     s3 = seeds[0,i+1:*] & s4 = seeds[1, i+1:*]
     l2 = lower[s3, s4] & u2 = upper[s3, s4]
     
     ;- are there any overlaps?
     overlap = (l1 lt u2 and u1 gt l2)
     hit = where(overlap, ct)
     if ct eq 0 then continue

     for j = 0, ct - 1 do begin
        s = [s1, s2, s3[hit[j]], s4[hit[j]]]
        s = s[uniq(s, sort(s))]
        ns = n_elements(s)
        assert, ns ge 3
        ;- there are 3 or 4 seed values between these
        ;- two structures. All possible triplet permutations
        ;- are candidates
        if ns eq 3 then result->push, s $
        else begin 
           result->push, s[[0,1,2]]
           result->push, s[[0,1,3]]
           result->push, s[[0,2,3]]
           result->push, s[[1,2,3]]
        endelse        
     endfor
  endfor 
  array = result->toArray()
  obj_destroy, result

  nelem = n_elements(array)
  array = reform(array, 3, nelem/3)

  ;- find out which candidates are actually ambiguous
  ambiguous = ambiguous_triplet(array, lower, upper)
  hit = where(ambiguous, count)
  if count ne 0 then return, array[*, hit] else return, -1
end

pro test

  restore, 'mergetest2.sav'
  t0 = systime(/seconds)
;  profiler, /reset & profiler, /system & profiler
  x = find_ambiguities(lower, upper)
  print, time2string(systime(/seconds) - t0)
  help, x
;  profiler, /report
end
  
     
