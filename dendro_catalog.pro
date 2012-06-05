;+
; PURPOSE:
;  Measure several properties about each structure in a dendrogram
;
; INPUTS:
;   arg: Either a string naming a C++ generated dendrogram file,
;        or a pointer to an IDL-generated dendrogram
;
; KEYWORD PARAMETERS:
;  len_scale: physical size of each pixel
;  vel_scale: velocity channel width
;  flux2mass: Multiply by this factor to convert intensity to mass in M_sun
;
; OUTPUTS:
;  An array of structures (see structure_prop for fields)
;
; MODIFICATION HISTORY:
;  Jan 2011: Written by Chris Beaumont
;-
function dendro_catalog, arg, $
                         len_scale = len_scale, vel_scale = vel_scale, $
                         flux2mass = flux2mass

  if n_params() ne 1 then begin
     print, 'calling sequence:'
     print, ' result = dendro_catalog(file/ptr, [len_scale = len_scale, '
     print, '                         vel_scale = vel_scale, flux2mass = flux2mass])'
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

  nst = n_elements( (*ptr).height )
  rec = structure_props(/get_record)
  data = replicate(rec, nst)

  for i = 0, nst - 1, 1 do begin
     ind = substruct(i, ptr, count = ct)
     if ct eq 0 then continue
     data[i] = structure_props(ptr, ind, len_scale=len_scale, $
                               vel_scale = vel_scale, flux2mass = flux2mass)
  endfor

  for i = 0, nst - 1, 1 do begin
     partner = merger_partner(i, (*ptr).clusters, merge = m)
     leaf = leafward_mergers(i, (*ptr).clusters, /parents)
     if partner ne -1 then $
        data[i].shoulder_height = (*ptr).height[i] - (*ptr).height[m]
     if leaf[0] eq -1 then continue
     data[i].vol_left = max(data[leaf].vol, min = lo)
     data[i].vol_right = lo
  endfor

  if size(arg, /tname) eq 'STRING' then ptr_free, ptr
  return, data
end

