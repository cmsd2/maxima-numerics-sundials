;;; nvector.lisp — N_Vector serial bindings and ndarray bridge
;;;
;;; Uses copy (not zero-copy) between N_Vector and CL arrays.
;;; ODE state vectors are small (tens of elements), so copy cost is negligible.

(in-package #:numerics-sundials)

;;; --- CFFI bindings ---

(cffi:defcfun ("N_VNew_Serial" %nv-new-serial) n-vector
  (n sunindextype)
  (ctx sun-context))

(cffi:defcfun ("N_VDestroy" %nv-destroy) :void
  (v n-vector))

(cffi:defcfun ("N_VGetArrayPointer_Serial" %nv-get-array-pointer) :pointer
  (v n-vector))

(cffi:defcfun ("N_VGetLength_Serial" %nv-get-length) sunindextype
  (v n-vector))

;;; --- Helpers ---

(defun make-nvector (n &optional initial-values)
  "Create an N_Vector_Serial of length N. Optionally initialize from a
   list or simple-array of double-float values."
  (let* ((ctx (ensure-sun-context))
         (v (%nv-new-serial n ctx)))
    (when (cffi:null-pointer-p v)
      (maxima::merror "Failed to allocate N_Vector of length ~A" n))
    (when initial-values
      (let ((ptr (%nv-get-array-pointer v)))
        (etypecase initial-values
          (list
           (loop for val in initial-values
                 for i from 0
                 do (setf (cffi:mem-aref ptr :double i)
                          (coerce val 'double-float))))
          ((simple-array double-float (*))
           (loop for i from 0 below n
                 do (setf (cffi:mem-aref ptr :double i)
                          (aref initial-values i)))))))
    v))

(defun nvector-to-list (v)
  "Extract N_Vector contents as a list of double-floats."
  (let* ((n (%nv-get-length v))
         (ptr (%nv-get-array-pointer v)))
    (loop for i from 0 below n
          collect (cffi:mem-aref ptr :double i))))

(defun nvector-to-array (v)
  "Extract N_Vector contents as a simple-array of double-floats."
  (let* ((n (%nv-get-length v))
         (ptr (%nv-get-array-pointer v))
         (arr (make-array n :element-type 'double-float)))
    (loop for i from 0 below n
          do (setf (aref arr i) (cffi:mem-aref ptr :double i)))
    arr))

(defun copy-array-to-nvector (arr v)
  "Copy values from a CL array into an existing N_Vector."
  (let ((ptr (%nv-get-array-pointer v))
        (n (length arr)))
    (loop for i from 0 below n
          do (setf (cffi:mem-aref ptr :double i)
                   (coerce (aref arr i) 'double-float)))))

(defun copy-nvector-to-array (v arr)
  "Copy N_Vector values into an existing CL array."
  (let ((ptr (%nv-get-array-pointer v))
        (n (length arr)))
    (loop for i from 0 below n
          do (setf (aref arr i) (cffi:mem-aref ptr :double i)))))
