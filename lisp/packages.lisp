;;; packages.lisp — CL package for numerics-sundials

(defpackage #:numerics-sundials
  (:use #:cl)
  (:export
   ;; Context
   #:ensure-sun-context
   ;; N_Vector helpers
   #:make-nvector
   #:nvector-to-list
   #:nvector-to-array
   #:copy-array-to-nvector
   #:copy-nvector-to-array
   ;; Error handling
   #:sundials-check-flag
   #:sundials-flag-message))
