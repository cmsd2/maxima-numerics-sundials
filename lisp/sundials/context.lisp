;;; context.lisp — SUNContext lifecycle management
;;;
;;; SUNDIALS v6.0+ requires a SUNContext for all operations.
;;; We use a single global context (Maxima is single-threaded, no MPI).
;;;
;;; When SUNDIALS is built with MPI support (e.g., Homebrew on macOS),
;;; SUNContext_Create expects MPI_COMM_NULL, not a null pointer.
;;; We detect this at load time and pass the correct value.

(in-package #:numerics-sundials)

(cffi:defcfun ("SUNContext_Create" %sun-context-create) :int
  (comm :pointer)       ; SUNComm: MPI_Comm (if MPI) or int (if no MPI)
  (ctx :pointer))       ; SUNContext* (output)

(cffi:defcfun ("SUNContext_Free" %sun-context-free) :int
  (ctx :pointer))       ; SUNContext* (pointer to context pointer)

;;; --- MPI_COMM_NULL detection ---
;;; When SUNDIALS links against MPI, passing NULL (0) as the communicator
;;; triggers MPI_Comm_dup before MPI_Init. We need to pass MPI_COMM_NULL,
;;; which in Open MPI is the address of the global `ompi_mpi_comm_null`.

(defun find-mpi-comm-null ()
  "Find the value of MPI_COMM_NULL for the linked MPI library.
   Returns a CFFI pointer, or null-pointer if MPI is not linked."
  ;; Open MPI: MPI_COMM_NULL = (MPI_Comm)&ompi_mpi_comm_null
  (let ((sym (cffi:foreign-symbol-pointer "ompi_mpi_comm_null")))
    (when sym (return-from find-mpi-comm-null sym)))
  ;; MPICH: MPI_COMM_NULL = 0x04000000 (a sentinel value)
  (let ((sym (cffi:foreign-symbol-pointer "MPIR_Comm_null")))
    (when sym (return-from find-mpi-comm-null (cffi:make-pointer #x04000000))))
  ;; No MPI or unknown implementation: null pointer (works for non-MPI builds)
  (cffi:null-pointer))

(defvar *sun-comm-null* (find-mpi-comm-null)
  "The correct null communicator for SUNContext_Create.
   MPI_COMM_NULL when SUNDIALS is MPI-enabled, null-pointer otherwise.")

;;; --- Context management ---

(defvar *sun-context* nil
  "Global SUNContext for serial SUNDIALS operations.
   Created lazily on first use, persists for the Maxima session.")

(defun ensure-sun-context ()
  "Return the global SUNContext, creating it if necessary."
  (unless *sun-context*
    (cffi:with-foreign-object (ctx-ptr :pointer)
      (let ((flag (%sun-context-create *sun-comm-null* ctx-ptr)))
        (sundials-check-flag flag "SUNContext_Create")
        (setf *sun-context* (cffi:mem-ref ctx-ptr :pointer)))))
  *sun-context*)
