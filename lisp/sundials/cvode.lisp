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

;;; --- Rootfinding (event detection) ---

(cffi:defcfun ("CVodeRootInit" %cvode-root-init) :int
  (cvode-mem :pointer)
  (nrtfn :int)
  (g :pointer))              ; CVRootFn callback

(cffi:defcfun ("CVodeGetRootInfo" %cvode-get-root-info) :int
  (cvode-mem :pointer)
  (rootsfound :pointer))     ; int* array

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
