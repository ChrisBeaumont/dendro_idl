function dendro_catalog_2, arg, levels, $
                         len_scale = len_scale, vel_scale = vel_scale, $
                         flux2mass = flux2mass

  if n_params() ne 2 then begin
     print, 'Calling sequence'
     print, 'result = dendro_catalog_2(arg, levels, [len_scale=len_scale,'
     print, '                  vel_scale=vel_scale, flux2mass=flux2mass])'
     return, !values.f_nan
  endif

  if size(arg, /tname) eq 'STRING' then begin
     file = arg
     if ~file_test(file) then $
        message, 'Cannot find dendrogram file: ', file

     catch, error
     if error ne 0 then begin
        catch, /cancel
        print, 'Could not convert file to pointer: '+file
        message, !error_state.msg
     endif
     ptr = dendrocpp2idl(file)
     catch, /cancel
  endif else if size(arg, /tname) eq 'POINTER' then begin
     ptr = arg
  endif else $
     message, 'Input must be a string or a pointer'

  nlevel = n_elements(levels)
  nst = n_elements( (*ptr).height)

  rec = structure_props(/get_record)
  data = replicate(rec, nst, nlevel)

  for i = 0, nst - 1, 1 do begin
     ind0 = substruct(i, ptr, count = ct0)
     if ct0 eq 0 then continue
     ;ind0 = ind0[sort((*ptr).t[ind0])]
     ;t = (*ptr).t[ind0]

     p = merger_partner(i, (*ptr).clusters, merge = m)
     h_hi = (*ptr).height[i]
     h_lo = m eq -1 ? min((*ptr).height) : (*ptr).height[m]
     lw = leafward_mergers(i, (*ptr).clusters)

     for j = 0, nlevel - 1, 1 do begin
        if levels[j] le h_lo || levels[j] gt h_hi then continue
        hit = where((*ptr).t[ind0] gt levels[j], ct)
        if ct eq 0 then continue
        ind = ind0[hit]

        ;cutoff = value_locate(t, levels[j]) < (ct0 - 1)
        ;if cutoff le 0 then continue
        ;ind = ind0[0:cutoff]

        rec = structure_props(ptr, ind, len_scale=len_scale, $
                              vel_scale=vel_scale, flux2mass=flux2mass)

        ;- copy to leafward mergers
        assert, max(finite(data[lw, j].x)) eq 0
        data[lw, j] = rec
     endfor
  endfor

  if size(arg, /tname) eq 'STRING' then ptr_free, ptr
  return, data
end
