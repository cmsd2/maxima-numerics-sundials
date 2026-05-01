;;; cvode.lisp — CVODE low-level CFFI bindings

(in-package #:numerics-sundials)

;;; --- Core CVODE functions ---

(cffi:defcfun ("CVodeCreate" %cvode-create) :pointer
  (lmm :int)
  (ctx sun-context))

(cffi:defcfun ("CVodeInit" %cvode-init) :int
  (cvode-mem :pointer)
  (f :pointer)               ; CVRhsFn callback
  (t0 sunrealtype)
  (y0 n-vector))

(cffi:defcfun ("CVodeSStolerances" %cvode-ss-tolerances) :int
  (cvode-mem :pointer)
  (rtol sunrealtype)
  (atol sunrealtype))

(cffi:defcfun ("CVodeSetLinearSolver" %cvode-set-linear-solver) :int
  (cvode-mem :pointer)
  (ls sun-linear-solver)
  (a sun-matrix))

(cffi:defcfun ("CVodeSetUserData" %cvode-set-user-data) :int
  (cvode-mem :pointer)
  (user-data :pointer))

(cffi:defcfun ("CVodeSetMaxNumSteps" %cvode-set-max-num-steps) :int
  (cvode-mem :pointer)
  (mxsteps :long))

(cffi:defcfun ("CVodeSetInitStep" %cvode-set-init-step) :int
  (cvode-mem :pointer)
  (hin sunrealtype))           ; initial step size (0 = let CVODE choose)

;;; --- Rootfinding (event detection) ---

(cffi:defcfun ("CVodeRootInit" %cvode-root-init) :int
  (cvode-mem :pointer)
  (nrtfn :int)
  (g :pointer))              ; CVRootFn callback

(cffi:defcfun ("CVodeGetRootInfo" %cvode-get-root-info) :int
  (cvode-mem :pointer)
  (rootsfound :pointer))     ; int* array — sign indicates direction:
                             ;   +1: rising crossing (g goes - → +)
                             ;   -1: falling crossing (g goes + → -)
                             ;    0: this root didn't fire

(cffi:defcfun ("CVodeSetRootDirection" %cvode-set-root-direction) :int
  (cvode-mem :pointer)
  (rootdir   :pointer))      ; int* of length nrtfn — per-root filter:
                             ;   +1: only detect rising crossings
                             ;   -1: only detect falling crossings
                             ;    0: detect any crossing (default)

;;; --- Reinitialisation (after a discrete state change) ---

(cffi:defcfun ("CVodeReInit" %cvode-reinit) :int
  (cvode-mem :pointer)
  (t0 sunrealtype)
  (y0 n-vector))             ; reuse the same cvode-mem (and its
                             ; linear solver) but with a new initial
                             ; condition.  Multistep history is reset
                             ; — appropriate after an event-triggered
                             ; state reset (the trajectory is no longer
                             ; smoothly continuous through t0).

;;; --- Stepping ---

(cffi:defcfun ("CVode" %cvode-solve) :int
  (cvode-mem :pointer)
  (tout sunrealtype)
  (yout n-vector)
  (tret :pointer)            ; sunrealtype* (output: actual t reached)
  (itask :int))

;;; --- Cleanup ---

(cffi:defcfun ("CVodeFree" %cvode-free) :void
  (cvode-mem-ptr :pointer))  ; void** (pointer to pointer)

;;; --- Diagnostics ---

(cffi:defcfun ("CVodeGetNumSteps" %cvode-get-num-steps) :int
  (cvode-mem :pointer)
  (nsteps :pointer))         ; long*

(cffi:defcfun ("CVodeGetLastOrder" %cvode-get-last-order) :int
  (cvode-mem :pointer)
  (qlast :pointer))           ; int*
