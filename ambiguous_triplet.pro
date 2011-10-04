;+
; PURPOSE:
;  This function determines whether the mergers between a triplet of
;  seed local maxima have been constrained precisely enough. The
;  function returns a 1 if a triplet is ambiguous -- that is, if the
;  order in which the thee seeds merge into a 'mini-dendrogram'
;  is not known given the uncertainties.
;
;  Note: rarely, an unambiguous merger order cannot be obtained,
;        given finite floating-point precision. In such a case, this
;        function will label the situation as unambiguous (since
;        resolving the ambiguity is impossible anyways).
;
; INPUTS:
;  inds: a (3xn) set of n triplets. Each triplet gives the index of a
;  local maximum, and a row/column in the merger matrix.
;  lower: The lower bounds for the merger matrix
;  upper: The upper bounds for the merger matrix
;
; OUTPUTS:
;  An n-element vector, with 1s for the triplets which remain
;  ambiguous, and need tighter bounds.
;
; MODIFICATION HISTORY:
;  July 21 2010: Written by Chris Beaumont
;  Feb 7 2010: Caught case when lower and upper converge to within
;              machine precision, and hence are intrinsically
;              ambiguous. cnb.
;- 
function ambiguous_triplet, inds, lower, upper
  compile_opt idl2

  ;- get the lower and upper limits for each merger
  m1l = lower[inds[0,*], inds[1,*]] & m1u = upper[inds[0,*], inds[1,*]]
  m2l = lower[inds[2,*], inds[1,*]] & m2u = upper[inds[2,*], inds[1,*]]
  m3l = lower[inds[0,*], inds[2,*]] & m3u = upper[inds[0,*], inds[2,*]]

  ;- determine which pairs of mergers potentially overlap
  ;- not that 1 pair _always_ overlaps, since the two seeds
  ;- which merge at the highest contour necessarily merge
  ;- with the third seed at the _same_, lower contour level.
  a12 = (m1l lt m2u) and (m1u gt m2l)
  a23 = (m2l lt m3u) and (m2u gt m3l)
  a13 = (m1l lt m3u) and (m1u gt m3l)

  tot = a12 + a23 + a13

  ;- tot should be 1 (non-ambigous) or 3 (ambiguous)
  h = histogram(tot, min = 0)
  assert, h[0] eq 0 && h[2] eq 0

  b1 = (m1u + m1l) / 2D
  b2 = (m2u + m2l) / 2D
  b3 = (m3u + m3l) / 2D
  eps = (b1 eq m1u) or (b1 eq m1l) or $
        (b2 eq m2u) or (b2 eq m2l) or $
        (b3 eq m3u) or (b3 eq m3l)
        
  return, tot eq 3 and ~eps
end
