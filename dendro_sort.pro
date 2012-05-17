function inten_key, ptr
  nst = n_elements((*ptr).height)
  result = fltarr(nst)
  for i = 0, nst - 1, 1 do begin
     result[i] = total((*ptr).t[substruct(i, ptr)], /nan)
  endfor
  return, result
end

function height_key, ptr
  nst = n_elements((*ptr).height)
  result = fltarr(nst)
  for i = 0, nst - 1, 1 do begin
     result[i] = max((*ptr).t[substruct(i, ptr)], /nan)
  endfor
  return, result
end

pro _width, ptr, width, i
  nst = n_elements((*ptr).height)
  nleaf = (nst + 1) / 2
  if i lt nleaf then begin
     width[i] = 1
  endif else begin
     p = leafward_mergers(i, (*ptr).clusters, /parents)
     _width, ptr, width, p[0]
     _width, ptr, width, p[1]
     width[i] = width[p[0]] + width[p[1]]
  endelse
end

pro _pos, ptr, i, width, pos, key
  nst = n_elements((*ptr).height)
  nleaf = (nst + 1) / 2
  if i lt nleaf then return
  p = leafward_mergers(i, (*ptr).clusters, /parents)
  if key[p[0]] gt key[p[1]] then $
     p = [p[1], p[0]]
  w = width[p[0]] + width[p[1]]
  gap = width[i] - w
  assert, gap ge 0
  pos[p[0]] = pos[i] - width[i] / 2 + gap / 3. + width[p[0]] / 2.
  pos[p[1]] = pos[i] - width[i] / 2. + 2. * gap / 3 + width[p[0]] + width[p[1]] / 2.
  _pos, ptr, p[0], width, pos, key
  _pos, ptr, p[1], width, pos, key
end

pro dendro_sort, ptr, key, width = width, inten = inten, height = height
  nst = n_elements((*ptr).height)
  pos = fltarr(nst)

  if keyword_set(inten) && keyword_set(height) then $
     message, 'cannot set both inten and height'

  if keyword_set(inten) then key = inten_key(ptr)
  if keyword_set(height) then key = height_key(ptr)

  if ~keyword_set(width) then begin
     width = fltarr(nst)
     _width, ptr, width, nst - 1
  endif


  _pos, ptr, nst - 1, width, pos, key
  (*ptr).xlocation = pos
end
