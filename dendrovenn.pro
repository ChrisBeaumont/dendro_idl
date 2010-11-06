;+
; PURPOSE:
;  This function computes venn-diagram like comparisons between two
;  dendrogram structures. It determines which of the substrutures that
;  make up structures a and b belong to both sets, one and not the
;  other, etc. This can be useful when making image masks with large
;  datasets.
;
; INPTUS:
;  a: The index of the first substructure
;  b: The index of the second substructure
;  clusters: The dendrogram clusters array
;
; OUTPUTS:
;  A structure with six tags. anotb lists the substructures within a
;  that are not in b. bnota lists the substructures within b not in
;  a. ab lists the substructures common to both structures. If any of
;  these sets are empy, -1 is returned in that field. The tags
;  anotbct, bnotact, and abct list the number of substructures in each category.
;
; MODIFICATION HISTORY:
;  October 2010: Written by Chris Beaumont
;-
function dendrovenn, a, b, clusters

  if a eq -1 then aleaf = get_leaves()
  if a lt -1 then aleaf = -1
  if a gt -1 then aleaf = leafward_mergers(a, clusters)

  if b eq -1 then bleaf = get_leaves()
  if b lt -1 then bleaf = -1
  if b gt -1 then bleaf = leafward_mergers(b, clusters)


  ;- edge cases
  if a lt -1 then begin
     return, {anotbct:0, anotb:-1, bnotact:n_elements(bleaf), $
              bnota:bleaf, abct:0, ab:-1}
  endif
  if b lt -1 then begin
     return, {anotbct:n_elements(aleaf), anotb:aleaf, $
              bnotact:0, bnota:-1, abct:0, ab:-1}
  endif
     
  aleaf = (a eq -1) ? get_leaves() : leafward_mergers(a, clusters)
  bleaf = (b eq -1) ? get_leaves() : leafward_mergers(b, clusters)

  top = max([aleaf, bleaf])
  ah = histogram([aleaf], min=0, max = top, loc = l )
  bh = histogram([bleaf], min=0, max = top)
  assert, max(ah) le 1 && max(bh) le 1

  anotb = where(ah and ~bh, ct1)
  anotb = ct1 eq 0 ? -1 : l[anotb]

  bnota = where(bh and ~ah, ct2)
  bnota = ct2 eq 0 ? -1 : l[bnota]

  ab = where(ah and bh, ct3)
  ab = ct3 eq 0 ? -1 : l[ab]
  return, {anotb: anotb, anotbct: ct1, $
           bnota: bnota, bnotact: ct2, $
           ab:ab, abct: ct3}
end
