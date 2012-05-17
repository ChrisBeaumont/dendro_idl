;+
; PURPOSE:
;  Measure several properties about each structure in a dendrogram
;
; INPUTS:
;   arg: Either a string naming a C++ generated dendrogram file,
;        or a pointer to an IDL-generated dendrogram
;
; KEYWORD PARAMETERS:
;  len_scale: physical size of each pixel, in pc
;  vel_scale: velocity channel width, in km/s
;  flux2mass: Multiply by this factor to convert intensity to mass in M_sun
;
; OUTPUTS:
;  A catalog with the following tags:
;    sig_maj -- 2nd moment about the first principal axis, when shape
;               is projected onto 2D sky. Units are pixels or pc
;    sig_min -- 2nd moment about the second principal axis, when shape
;               is projected onto 2D sky. Units are pixels or pc
;    sig_v   -- 2nd moment of velocity. Units are km/s
;    sig_r   -- sqrt(sig_maj * sig_min). Units are pixels or pc
;    flux    -- integrated flux.
;    vol     -- voxel volume of structure
;    virial  -- virial parameter estimated from sig_r, sig_v, flux,
;               and flux2mass
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

  if ~keyword_set(len_scale) then len_scale = 1
  if ~keyword_set(vel_scale) then vel_scale = 1
  if ~keyword_set(flux2mass) then flux2mass = 1

  nst = n_elements( (*ptr).height )

  nan = !values.f_nan
  rec = {x: nan, $
         y: nan, $
         v: nan, $
         sig_maj:nan, $
         sig_min:nan, $
         sig_v:nan, $
         sig_r:nan, $
         area_mask:nan, $
         perimeter_mask:nan, $
         flux:nan, $
         peak_inten:nan, $
         vol:nan, $
         virial:nan, $
         shoulder_height:nan, $
         vol_left:nan, $
         vol_right:nan }

  data = replicate(rec, nst)

  for i = 0, nst - 1, 1 do begin
     ind = substruct(i, ptr, count = ct)
     if ct eq 0 then continue

     x = (*ptr).x[ind]
     y = (*ptr).y[ind]
     v = (*ptr).v[ind]
     t = (*ptr).t[ind]

     amask_ind = long(x) + long(y) * (max(x)+1)
     amask_ind = uniq(amask_ind, sort(amask_ind))
     data[i].area_mask = n_elements(amask_ind)

     if total(t, /double) le 0 then begin
        print, 'Integrated flux is negative. Skipping structure ' + strtrim(i, 2)
        continue
     endif

     tt = total(t, /double)
     data[i].x = total( x * t, /double) / tt
     data[i].y = total( y * t, /double) / tt
     data[i].v = total( v * t, /double) / tt
     data[i].peak_inten = max(t, /nan)

     stamp = dblarr( range(x)+2, range(y)+2, range(v) + 2)
     stamp[ x- min(x), y - min(y), v - min(v) ] = t
     indices, stamp, ix, iy, iz

     ;- get shape statistics for 2D projection
     sz = size(stamp)
     shape_stat3, reform(total(stamp, 3), sz[1], sz[2], 1), mean = mean, $
                  paxis = paxis, obl = obl, $
                  sph = sph

     ;- sanity check -- mean is calculated correctly
     tt = total(stamp,/double)
     mean2 = [total(ix * stamp, /double) / tt, $
              total(iy * stamp, /double) / tt, $
              total(iz * stamp, /double) / tt]
     assert, max(abs(mean2[0:1] - mean[0:1]) / (abs(mean[0:1]) > .005) ) lt (1d-3)
     mean = mean2

     ;- find perimeter of 2d mask
     mask = reform(max(stamp ne 0, dim=3), sz[1], sz[2])
     contour, mask, lev = .5, path_xy = pth, /path_data_coords, $
              path_info = info
     best = max(info.n, maxloc)
     first = info[maxloc].offset
     last = first + info[maxloc].n - 1
     pth = pth[*, first : last]
     data[i].perimeter_mask = perimeter(pth[0, *], pth[1, *])

     ;- find normalized principle axis for 2D projection
     ix -= mean[0] & iy -= mean[1] & iz -= mean[2]
     ax1 = reform(paxis[*,0] * [1,1,0])
     if max(abs(ax1)) eq 0 then ax1 = [1,0]
     ax1 /= sqrt( total(ax1^2) )

     ;- project (xy) onto the major/minor axes in XY plane
     p_maj = ix * ax1[0] + iy * ax1[1]
     p_min = (-1) * ix * ax1[1] + iy * ax1[0]

     assert, abs(tt - total(t,/double)) / tt lt 1d-3
     assert, total(stamp * ix) / tt lt 1d-3
     assert, total(stamp * iy) / tt lt 1d-3
     assert, total(stamp * iz) / tt lt 1d-3

     assert, total(stamp * p_maj) / tt lt 1d-2
     assert, total(stamp * p_min) / tt lt 1d-2
     assert, total(stamp * iz) / tt lt 1d-2

     sig_maj = sqrt(total(stamp * p_maj^2) / tt) * len_scale
     sig_min = sqrt(total(stamp * p_min^2) / tt) * len_scale
     sig_vel = sqrt(total(stamp * iz^2) / tt) * vel_scale
     sig_r = sqrt(sig_maj * sig_min)

     if sig_maj lt sig_min then begin
        erase
        p = [0,0,1,1]
        tvimage, bytscl(total(stamp, 3, /nan)), /keep, /noi, pos = p
        contour, total(stamp, 3), ix, iy, /nodata, /noerase, pos = p
        oplot, [0, ax1[0]], [0, ax1[1]], color = fsc_color('red')
        stop
     endif

     data[i].sig_maj = sig_maj
     data[i].sig_min = sig_min
     data[i].sig_v = sig_vel
     data[i].sig_r = sig_r
     data[i].flux = tt
     data[i].vol = n_elements(t)
  endfor
  eta = 1.91 ;- correct for concentration of R. see rosolowsky 2008
  g = 6.673d-8
  data.virial = 5 * eta * data.sig_r * data.sig_v^2 / (data.flux * flux2mass * g)

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

