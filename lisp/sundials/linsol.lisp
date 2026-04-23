;;; linsol.lisp — SUNLinearSolver dense bindings

(in-package #:numerics-sundials)

(cffi:defcfun ("SUNLinSol_Dense" %sun-linsol-dense) sun-linear-solver
  (y n-vector)
  (a sun-matrix)
  (ctx sun-context))

(cffi:defcfun ("SUNLinSolFree" %sun-linsol-free) :int
  (ls sun-linear-solver))
