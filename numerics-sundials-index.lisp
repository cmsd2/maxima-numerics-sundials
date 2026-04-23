(in-package :cl-info)
(let (
(deffn-defvr-pairs '(
; CONTENT: (<INDEX TOPIC> . (<FILENAME> <BYTE OFFSET> <LENGTH IN CHARACTERS> <NODE NAME>))
("np_cvode" . ("numerics-sundials.info" 2362 46 "Function np_cvode f vars y0 tspan"))
("np_cvode <1>" . ("numerics-sundials.info" 2639 95 "Function np_cvode f vars y0 tspan rtol atol method events max_steps"))
("np_cvode <2>" . ("numerics-sundials.info" 2992 45 "Function np_cvode f_func y0 tspan"))
("np_cvode <3>" . ("numerics-sundials.info" 3219 2038 "Function np_cvode f_func y0 tspan rtol atol method max_steps"))
))
(section-pairs '(
; CONTENT: (<NODE NAME> . (<FILENAME> <BYTE OFFSET> <LENGTH IN CHARACTERS>))
("Definitions for numerics-sundials" . ("numerics-sundials.info" 1873 300))
("Introduction" . ("numerics-sundials.info" 607 626))
)))
(load-info-hashtables (maxima::maxima-load-pathname-directory) deffn-defvr-pairs section-pairs))
