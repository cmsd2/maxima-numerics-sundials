;;; numerics-sundials.asd — ASDF system definition

(defsystem "numerics-sundials"
  :description "SUNDIALS ODE/DAE solvers for Maxima numerics ndarrays"
  :version "0.1.0"
  :license "MIT"
  :depends-on ("numerics/core" "cffi" "trivial-garbage" "alexandria")
  :serial t
  :components
  ((:file "packages")
   (:module "sundials"
    :serial t
    :components
    ((:file "library")
     (:file "types")
     (:file "context")
     (:file "nvector")
     (:file "matrix")
     (:file "linsol")
     (:file "cvode")
     (:file "cvode-wrapper")))))
