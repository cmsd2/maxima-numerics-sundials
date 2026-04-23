;;; numerics-sundials-loader.lisp — Bootstrap SUNDIALS CFFI bindings via ASDF/Quicklisp

(in-package :maxima)

;; Derive the project directory from this file's location
(let* ((here (make-pathname :directory (pathname-directory *load-truename*))))
  (pushnew here asdf:*central-registry* :test #'equal)
  (pushnew (merge-pathnames "lisp/" here)
           asdf:*central-registry* :test #'equal))

;; Load the system via Quicklisp (resolves cffi + dependencies)
(funcall (intern "QUICKLOAD" :ql) "numerics-sundials" :silent t)
