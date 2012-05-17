;+
; PURPOSE:
;  This class creates creates an internal index for a dendrogram that
;  speeds up subsequent calls to functions like leafward_merger.
;
; MODIFICATION HISTORY:
;  June 2011: Written by Chris Beaumont
;-

;+
; PURPOSE:
;  Calculate the leafward mergers of a given leaf/branch
;
; INPUTS:
;  index: The index for which to return leafward mergers
;
; KEYWORD PARAMETERS:
;  Count: The number of leafward nodes
;
; OUTPUT:
;  A list of leafward mergers, including index itself.
;-
function dendro_index::leafward_mergers, index, count = count
  count = (*self.index)[2 * index + 1]
  pos = (*self.index)[2 * index]
  return, (*self.index)[pos:pos+count-1]
end


;+
; PURPOSE:
;  Generate the indices belonging to a given structure
;
; INPUTS:
;  index: The index to consider
;
; KEYWORD PARAMETERS:
;  count: The number of elements returned
;
; OUTPUT:
;  A list of indices, or -1 if the index has none.
;-
function dendro_index::substruct, index, count = count, single = single
  count = 0

  if n_elements(index) ne 1 then $
     message, 'index must be a scalar'
  if ~keyword_set(single) then begin
     indices = self->leafward_mergers(index)
  endif else begin
     indices = index
  endelse

  ptr = self.ptr
  count = total((*ptr).cluster_label_h[indices])
  if count eq 0 then return, -1

  offset = 0
  result = lonarr(count)
  for i = 0L, n_elements(count) - 1, 1 do begin
     x = indices[i]
     if (*ptr).cluster_label_h[x] eq 0 then continue
     ind = (*ptr).cluster_label_ri[(*ptr).cluster_label_ri[x] : $
                                   (*ptr).cluster_label_ri[x+1]-1]
     nind = n_elements(ind)
     result[offset:offset + nind - 1] = ind
     offset += n_elements(ind)
  endfor
  return, result
end


;+
; PURPOSE:
;  An internal method to build the index
;
; METHOD:
;  The index list is similar to histogram's reverse indices.
;  index[2*i] lists the number of children to structure i (including i
;  itself). index[2*i+1] lists an offset such that, if index[2*i] = n
;  and index[2*i+1] = off, index[n:n + off-1] = a list of structure i
;  and its substructures. index[n] = i.
;
; INPUTS:
;  clusters: The cluster array
;  index: The index vector to populate
;  id: The dendrogram id to consider
;  pos: The variable off discussed above
;
; OUTPUTS:
;  The number of substructures for structure id
;-
function dendro_index::_build_index, clusters, index, id, pos
  nleaf = self.nleaf
  index[2L * id] = pos
  index[pos] = id
  if id lt nleaf then begin
     index[2 * id + 1] = 1
     return, 1
  endif

  left = clusters[0, id - nleaf]
  right = clusters[1, id - nleaf]

  lpos = pos + 1
  left = self->_build_index(clusters, index, left, lpos)
  rpos = lpos + left
  right = self->_build_index(clusters, index, right, rpos)

  index[2 * id + 1] = left + right + 1
  return, left + right + 1
end

;+
; PURPOSE:
;  Create a new index
;
; INPUTS:
;  ptr: A dendrogram pointer to index
;
; OUTPUTS:
;  1 for success
;-
function dendro_index::init, ptr
  self.ptr = ptr

  nst = n_elements( (*ptr).height)
  nleaf = (nst + 1) / 2
  self.nleaf =  nleaf

  index = lonarr(3L * nst)
  junk = self->_build_index((*ptr).clusters, index, nst - 1, 2 * nst)
  assert, junk eq nst

  self.index = ptr_new(index, /no_copy)
  return, 1
end

;+
; PURPOSE:
;  Delete an index and free memory
;-
pro dendro_index::cleanup
  ptr_free, self.index
end

;+
; PURPOSE:
;  Defines the index class
;-
pro dendro_index__define
  obj = {dendro_index, $
         ptr:ptr_new(), $ ;- dendrogram pointer
         index:ptr_new(), $ ;- index list
         nleaf: 0L} ;- number of leaves
end

pro test
  file = '~/stella_sims/paper/raw_dendro_ppv_256_noise3.fits'
  p = dendrocpp2idl(file)
  nst = n_elements((*p).height)
  t0 = systime(/seconds)
  in = obj_new('dendro_index', p)
  t1 = systime(/seconds)
  print, 'Time to index: ' + time2string(t1 - t0)
  dt1 = 0
  dt2 = 0
  for i = 0, 100, 1 do begin
     t0 = systime(/seconds)
     l1 = leafward_mergers(nst-i-1, (*p).clusters)
     t1 = systime(/seconds)
     l2 = in->leafward_mergers(nst-i-1)
     t2 = systime(/seconds)
     dt2 += t2 - t1
     dt1 += t1 - t0
     assert, array_equal(l1[sort(l1)], l2[sort(l2)])
  endfor
  print, 'Time in leafward_mergers: '+time2string(dt1)
  print, 'Time in index->leafward: '+time2string(dt2)
  obj_destroy, in
  ptr_free, p
  print, 'all tests passed'
end
