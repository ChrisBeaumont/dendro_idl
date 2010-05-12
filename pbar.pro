;+
; PURPOSE:
;  This procedure draws progress bars to the screen.
;
; INPUTS:
;  progress: A number from 0-1, indicating the fractional completion
;  of some process. A partially-filled progress bar will be drawn to
;  the screen.
;
; KEYWORDS:
;  name: Choose a display name for the current monitored process
;  new: Create a new progress bar
;  close: Close the current progress bar
;
; RESTRICTIONS:
;  Only one pbar can be displayed at a time
;
; EXAMPLE:
;  pbar, 0, /new, name = 'Example'
;  for i = 0, 100, 1 do begin
;    pbar, i/100.
;    wait, .1
;  endfor
;  pbar, /close
;
; MODIFICATION HISTORY:
;  Feb 2010: Written by Chris Beaumont
;-
pro pbar, progress, name = name, new = new, close= close
  common pbar, t0, wid, pwid

  if keyword_set(close) then begin
     wdelete, wid
     wdelete, pwid
     wid = -1
     return
  endif

  if keyword_set(new) || n_elements(progress) eq 0 then begin
     ;pbar, /close
     name = keyword_set(name) ? string(name) : 'Progress'
     window, xsize = 500, ysize = 100, title = name, xpos = 0, ypos = 0, /free
     wid = !d.window
     
     window, xsize = 500, ysize = 100, title = name, xpos= 0, ypos = 0, /free,/pix
     pwid = !d.window
     
     t0 = systime(/seconds)
     return
  endif
     
  wold = !d.window
  wset, pwid
  im = bytarr(500, 100)
  im[0:((progress < 1) * 499), *] = 255
  tv, im
  elapsed = systime(/seconds) - t0
  todo=time2string(elapsed / progress - elapsed)+' remaining'
  done=time2string(elapsed)+' elapsed'
  if progress gt 0 then begin
     xyouts, .1, .6, done, $
             /norm, charsize = 2, color = fsc_color('red'), $
             charthick = 2
     xyouts, .1, .4, todo, /norm, charsize = 2, charthick = 2, $
             color = fsc_color('red')
  endif
  wset, wid
  device, copy=[0,0,!D.x_size, !D.y_size,0,0,pwid]
  wset, wold

end

pro test

  pbar, name='Testing...', /new
  
  for i = 0, 100, 1 do begin
     pbar, i/100D
     wait, .1
  endfor
  pbar, /close
end
