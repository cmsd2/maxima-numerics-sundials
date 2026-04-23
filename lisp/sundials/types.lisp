;;; types.lisp — SUNDIALS type definitions, constants, and error handling

(in-package #:numerics-sundials)

;;; --- Basic types ---
;;; SUNDIALS uses sunrealtype = double, sunindextype = int64_t by default.

(cffi:defctype sunrealtype :double)
(cffi:defctype sunindextype :long)
(cffi:defctype sunbooleantype :int)

;; Opaque pointer types
(cffi:defctype sun-context :pointer)
(cffi:defctype n-vector :pointer)
(cffi:defctype sun-matrix :pointer)
(cffi:defctype sun-linear-solver :pointer)

;;; --- CVODE constants ---

(defconstant +cv-adams+ 1 "CVODE Adams method (non-stiff)")
(defconstant +cv-bdf+ 2 "CVODE BDF method (stiff)")

(defconstant +cv-normal+ 1 "CVode task: integrate to tout")
(defconstant +cv-one-step+ 2 "CVode task: take one internal step")

;; CVODE return flags
(defconstant +cv-success+ 0)
(defconstant +cv-tstop-return+ 1)
(defconstant +cv-root-return+ 2)
(defconstant +cv-warning+ 99)
(defconstant +cv-too-much-work+ -1)
(defconstant +cv-too-much-acc+ -2)
(defconstant +cv-err-failure+ -3)
(defconstant +cv-conv-failure+ -4)
(defconstant +cv-linit-fail+ -5)
(defconstant +cv-lsetup-fail+ -6)
(defconstant +cv-lsolve-fail+ -7)
(defconstant +cv-rhsfunc-fail+ -8)
(defconstant +cv-first-rhsfunc-err+ -9)
(defconstant +cv-reptd-rhsfunc-err+ -10)
(defconstant +cv-unrec-rhsfunc-err+ -11)

;;; --- Error handling ---

(defun sundials-flag-message (flag)
  "Return a human-readable message for a SUNDIALS return flag."
  (case flag
    (-1 "too much work (increase max_steps or try a different method)")
    (-2 "too much accuracy requested (increase tolerances)")
    (-3 "error test failures (problem may be ill-conditioned)")
    (-4 "convergence failures (try BDF method for stiff problems)")
    (-5 "linear solver initialization failed")
    (-6 "linear solver setup failed")
    (-7 "linear solver solve failed")
    (-8 "RHS function evaluation failed")
    (-9 "first RHS evaluation failed (check initial conditions and RHS function)")
    (-10 "repeated RHS evaluation errors")
    (-11 "unrecoverable RHS function error")
    (otherwise (format nil "unknown error code ~A" flag))))

(defun sundials-check-flag (flag context)
  "Check a SUNDIALS return flag and signal a Maxima error if negative."
  (when (< flag 0)
    (maxima::merror "~A: SUNDIALS error (flag=~A): ~A"
                    context flag (sundials-flag-message flag))))
