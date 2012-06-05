;+
; PURPOSE:
;  Calculates several properties of a region of pixels in a
;  dendrogrammed data set
;
; INPUTS:
;  ptr: A pointer to a dendrogram structure
;  ind: An array of indices, defining which indices in the x/y/v/t
;       arrays belong to the structure of interest. (Calculated, e.g.,
;       by substruct)
;
; KEYWORD PARAMETERS:
;  len_scale: physical size of each pixel
;  vel_scale: velocity channel width
;  flux2mass: Multiply by this factor to convert intensity to mass
;
;  get_record: If set, then return a nan'd version of the
;              output structure (for determining it's format)
;
;
; OUTPUTS:
;  A structure with the following tags:
;    x : intensity-weighted mean x position
;    y : "                     " y position
;    v : "                     " v position
;    sig_maj : Semi-major axis of the projection onto the xy plane
;    sig_min : Semi-minor axis of the projection onto the xy plane
;    sig_v : RMS dispersion of velocity
;    sig_r : sqrt(sig_maj * sig_min)
;    flux : Sum of intensities over all indices
;    vol : Number of indices
;    area_mask : Area (in pixels) of projection onto xy plane
;    perimeter_mask : Perimeter (in pixels) of projection onto xy plane
;    peak_inten : Brightest pixel value
;    virial : virial ratio
;    shoulder_height : placeholder (NOT CALCULATED)
;    vol_left : placeholder (NOT CALCULATED)
;    vol_right : placeholder (NOT CALCULATED)
;
; NOTES:
;  If len_scale is provided, sig_maj, sig_min, and sig_r will be
;     converted from pixel units
;  If vel_scale is provided, sig_v and virial will be converted
;  The virial conversion assumes len_scale and vel_scale convert into
;  cgs units
;-
function structure_props, ptr, ind, len_scale = len_scale, vel_scale=vel_scale, $
                          flux2mass = flux2mass, get_record = get_record
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

  if keyword_set(get_record) then $
     return, rec

  if n_params() ne 2 then begin
     print, 'calling sequence'
     print, 'result = structure_props(ptr, ind, [len_scale=len_scale, '
     print, '                         vel_scale = vel_scale, flux2mass=flux2mass, /get_record])'
     return, nan
  endif

  if ind[0] eq -1 then return, rec

  if ~keyword_set(len_scale) then len_scale = 1
  if ~keyword_set(vel_scale) then vel_scale = 1
  if ~keyword_set(flux2mass) then flux2mass = 1

  x = (*ptr).x[ind]
  y = (*ptr).y[ind]
  v = (*ptr).v[ind]
  t = (*ptr).t[ind]

  amask_ind = long(x) + long(y) * (max(x)+1)
  amask_ind = uniq(amask_ind, sort(amask_ind))
  rec.area_mask = n_elements(amask_ind)

  if total(t, /double) le 0 then begin
     print, 'Integrated flux is negative. Skipping structure '
     return, rec
  endif

  tt = total(t, /double)
  rec.x = total( x * t, /double) / tt
  rec.y = total( y * t, /double) / tt
  rec.v = total( v * t, /double) / tt
  rec.peak_inten = max(t, /nan)

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
  rec.perimeter_mask = perimeter(pth[0, *], pth[1, *])

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

  rec.sig_maj = sig_maj
  rec.sig_min = sig_min
  rec.sig_v = sig_vel
  rec.sig_r = sig_r
  rec.flux = tt
  rec.vol = n_elements(t)

  eta = 1.91 ;- correct for concentration of R. see rosolowsky 2008
  g = 6.673d-8
  rec.virial = 5 * eta * rec.sig_r * rec.sig_v^2 / (rec.flux * flux2mass * g)

  return, rec
end
