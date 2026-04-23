;;; matrix.lisp — SUNMatrix dense bindings

(in-package #:numerics-sundials)

(cffi:defcfun ("SUNDenseMatrix" %sun-dense-matrix) sun-matrix
  (m sunindextype)
  (n sunindextype)
  (ctx sun-context))

(cffi:defcfun ("SUNMatDestroy" %sun-mat-destroy) :void
  (a sun-matrix))
