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
;  A structure with three tags. anotb lists the substructures within a
;  that are not in b. bnota lists the substructures within b not in
;  a. ab lists the substructures common to both structures. If any of
;  these sets are empy, -1 is returned in that field.
;
; MODIFICATION HISTORY:
;  October 2010: Written by Chris Beaumont
;-
function dendrovenn, a, b, clusters

  aleaf = leafward_mergers(a, clusters)
  bleaf = leafward_mergers(b, clusters)

  top = max([aleaf, bleaf])
  ah = histogram([aleaf], min=0, max = top, loc = l )
  bh = histogram([bleaf], min=0, max = top)
  assert, max(ah) le 1 && max(bh) le 1

  anotb = where(ah and ~bh, ct)
  anotb = ct eq 0 ? -1 : l[anotb]

  bnota = where(bh and ~ah, ct)
  bnota = ct eq 0 ? -1 : l[bnota]

  ab = where(ah and bh, ct)
  ab = ct eq 0 ? -1 : l[ab]
  return, {anotb: anotb, bnota:bnota, ab:ab}
end
