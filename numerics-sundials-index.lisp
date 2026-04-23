(in-package :cl-info)
(let (
(deffn-defvr-pairs '(
; CONTENT: (<INDEX TOPIC> . (<FILENAME> <BYTE OFFSET> <LENGTH IN CHARACTERS> <NODE NAME>))
("numerics_sundials_hello" . ("numerics-sundials.info" 1195 78 "Function numerics_sundials_hello"))
))
(section-pairs '(
; CONTENT: (<NODE NAME> . (<FILENAME> <BYTE OFFSET> <LENGTH IN CHARACTERS>))
("Definitions for numerics-sundials" . ("numerics-sundials.info" 958 125))
("Introduction to numerics-sundials" . ("numerics-sundials.info" 649 161))
)))
(load-info-hashtables (maxima::maxima-load-pathname-directory) deffn-defvr-pairs section-pairs))
