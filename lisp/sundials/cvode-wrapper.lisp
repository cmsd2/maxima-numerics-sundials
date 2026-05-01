;;; cvode-wrapper.lisp — High-level CVODE Maxima API
;;;
;;; Provides np_cvode for ODE integration via SUNDIALS CVODE.
;;; Supports expression mode (symbolic RHS) and function mode (callable RHS).
;;; Follows the same patterns as np_odeint in numerics/integrate/ode.lisp.

(in-package #:numerics-sundials)

;;; ================================================================
;;; Callback trampoline infrastructure
;;; ================================================================

(defvar *cvode-rhs-table* (make-hash-table)
  "Maps integer keys to Lisp RHS closures for CVODE callbacks.")
(defvar *cvode-rhs-counter* 0)

(defvar *cvode-root-table* (make-hash-table)
  "Maps integer keys to Lisp root-finding closures.")
(defvar *cvode-root-counter* 0)

(defun register-rhs-callback (fn)
  "Register a Lisp closure as a CVODE RHS callback. Returns an integer key."
  (let ((key (incf *cvode-rhs-counter*)))
    (setf (gethash key *cvode-rhs-table*) fn)
    key))

(defun unregister-rhs-callback (key)
  "Remove a callback from the RHS table."
  (remhash key *cvode-rhs-table*))

(defun register-root-callback (fn)
  "Register a Lisp closure for CVODE rootfinding. Returns an integer key."
  (let ((key (incf *cvode-root-counter*)))
    (setf (gethash key *cvode-root-table*) fn)
    key))

(defun unregister-root-callback (key)
  "Remove a callback from the root table."
  (remhash key *cvode-root-table*))

;;; --- CVODE RHS trampoline ---
;;; C signature: int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)

(cffi:defcallback cvode-rhs-trampoline :int
    ((t-val sunrealtype) (y n-vector) (ydot n-vector) (user-data :pointer))
  (let ((fn (gethash (cffi:pointer-address user-data) *cvode-rhs-table*)))
    (if fn
        (handler-case
            (progn (funcall fn t-val y ydot) 0)
          (error () -1))
        -1)))

;;; --- CVODE root trampoline ---
;;; C signature: int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)

(cffi:defcallback cvode-root-trampoline :int
    ((t-val sunrealtype) (y n-vector) (gout :pointer) (user-data :pointer))
  (let ((fn (gethash (cffi:pointer-address user-data) *cvode-root-table*)))
    (if fn
        (handler-case
            (progn (funcall fn t-val y gout) 0)
          (error () -1))
        -1)))

;;; ================================================================
;;; RHS closure builders
;;; ================================================================

(defun make-cvode-rhs-closure-expr (compiled-rhs neq)
  "Create a CVODE RHS closure from a compiled expression function.
   COMPILED-RHS: compiled function from coerce-float-fun, taking (t, y1, ..., yn).
   Returns: closure (t-val y-nvec ydot-nvec) -> void."
  (lambda (t-val y-nvec ydot-nvec)
    (let* ((y-ptr (%nv-get-array-pointer y-nvec))
           (ydot-ptr (%nv-get-array-pointer ydot-nvec))
           (y-list (loop for i from 0 below neq
                         collect (cffi:mem-aref y-ptr :double i)))
           (result (apply compiled-rhs (coerce t-val 'double-float) y-list)))
      (loop for v in (cdr result)
            for i from 0
            do (setf (cffi:mem-aref ydot-ptr :double i)
                     (coerce v 'double-float))))))

(defun make-cvode-rhs-closure-fn (f-func neq)
  "Create a CVODE RHS closure from a Maxima callable.
   F-FUNC: Maxima function f(t, y) -> list of derivatives.
   Pre-allocates a callback ndarray to avoid allocation in the hot loop."
  (let* ((cb-tensor (magicl:empty (list neq) :type 'double-float
                                              :layout :column-major))
         (cb-handle (numerics:make-ndarray cb-tensor))
         (cb-wrapped (maxima::numerics-wrap cb-handle))
         (cb-storage (maxima::numerics-tensor-storage cb-tensor)))
    (lambda (t-val y-nvec ydot-nvec)
      ;; Copy N_Vector y into callback ndarray
      (let ((y-ptr (%nv-get-array-pointer y-nvec)))
        (loop for i from 0 below neq
              do (setf (aref cb-storage i)
                       (cffi:mem-aref y-ptr :double i))))
      ;; Call Maxima function
      (let ((result (if (symbolp f-func)
                        (maxima::mfuncall f-func t-val cb-wrapped)
                        (maxima::mapply f-func
                                        (list t-val cb-wrapped)
                                        f-func))))
        (unless (and (consp result) (eq (caar result) 'maxima::mlist))
          (maxima::merror
            "np_cvode: f(t,y) must return a list of derivatives, got: ~M"
            result))
        ;; Copy result into ydot N_Vector
        (let ((ydot-ptr (%nv-get-array-pointer ydot-nvec)))
          (loop for v in (cdr result) for i from 0
                do (setf (cffi:mem-aref ydot-ptr :double i)
                         (coerce (maxima::$float v) 'double-float))))))))

;;; ================================================================
;;; Event (root-finding) closure builder
;;; ================================================================

(defun make-cvode-event-closure-expr (compiled-events neq nrtfn)
  "Create a combined root-finding closure from compiled event functions.
   COMPILED-EVENTS: list of compiled functions, each (t, y1, ..., yn) -> scalar.
   Returns: closure (t-val y-nvec gout-ptr) -> void."
  (lambda (t-val y-nvec gout-ptr)
    (let* ((y-ptr (%nv-get-array-pointer y-nvec))
           (y-list (loop for i from 0 below neq
                         collect (cffi:mem-aref y-ptr :double i))))
      (loop for fn in compiled-events
            for i from 0
            do (let ((val (apply fn (coerce t-val 'double-float) y-list)))
                 (setf (cffi:mem-aref gout-ptr :double i)
                       (coerce (if (consp val)
                                   ;; coerce-float-fun returns (mlist val)
                                   (cadr val)
                                   val)
                               'double-float)))))))

;;; ================================================================
;;; Helpers
;;; ================================================================

(defun np-cvode-extract-values (arg name)
  "Extract a list of double-float values from an ndarray or Maxima list."
  (cond
    ((maxima::$ndarray_p arg)
     (let* ((handle (maxima::numerics-unwrap arg))
            (tensor (numerics:ndarray-tensor handle)))
       (coerce (maxima::numerics-tensor-storage tensor) 'list)))
    ((and (listp arg) (eq (caar arg) 'maxima::mlist))
     (mapcar (lambda (v) (coerce (maxima::$float v) 'double-float))
             (cdr arg)))
    (t (maxima::merror "np_cvode: ~A must be a list or ndarray" name))))

(defun np-cvode-parse-method (method)
  "Map Maxima method symbol to CVODE linear multistep constant."
  (cond
    ((null method) +cv-adams+)
    ((eq method 'maxima::$adams) +cv-adams+)
    ((eq method 'maxima::$bdf) +cv-bdf+)
    (t (maxima::merror "np_cvode: unknown method '~A'. Use 'adams' or 'bdf'."
                       method))))

;;; ================================================================
;;; Core stepping loop
;;; ================================================================

(defun np-cvode-core (rhs-closure neq y0-list t-list lmm rtol atol
                       &key event-closure (nrtfn 0) (max-steps 5000)
                            rootdir)
  "CVODE integration core.
   RHS-CLOSURE: (lambda (t y ydot) ...).
   EVENT-CLOSURE: single combined closure (t y gout) or nil.
   NRTFN: number of root functions (event expressions).
   ROOTDIR: optional list of NRTFN integers in {-1, 0, +1} — per-root
     direction filter (-1 = falling only, 0 = any, +1 = rising).
     When nil or omitted, all roots use the default 0.
   Returns 2D ndarray [n_times, 1+neq].
   If events are present, returns [trajectory, events-list] where each
     event entry is [t_event, [y...], [event_indices], [directions]];
     directions are +1/-1 indicating the crossing direction per
     event index that fired."
  (let* ((ctx (ensure-sun-context))
         (n-times (length t-list))
         ;; Create SUNDIALS objects
         (cvode-mem (%cvode-create lmm ctx))
         (y-nv (make-nvector neq y0-list))
         (sun-mat (%sun-dense-matrix neq neq ctx))
         (sun-ls (%sun-linsol-dense y-nv sun-mat ctx))
         ;; Register RHS callback
         (rhs-key (register-rhs-callback rhs-closure))
         ;; Result storage
         (ncols (1+ neq))
         (result (magicl:empty (list n-times ncols) :type 'double-float
                                                     :layout :column-major))
         (result-storage (maxima::numerics-tensor-storage result))
         (events-collected nil))
    (unwind-protect
        (progn
          (when (cffi:null-pointer-p cvode-mem)
            (maxima::merror "np_cvode: CVodeCreate failed"))
          ;; Initialize
          (sundials-check-flag
            (%cvode-init cvode-mem
                         (cffi:callback cvode-rhs-trampoline)
                         (coerce (first t-list) 'double-float)
                         y-nv)
            "CVodeInit")
          ;; Set tolerances
          (sundials-check-flag
            (%cvode-ss-tolerances cvode-mem rtol atol)
            "CVodeSStolerances")
          ;; Set linear solver
          (sundials-check-flag
            (%cvode-set-linear-solver cvode-mem sun-ls sun-mat)
            "CVodeSetLinearSolver")
          ;; Set user data (integer key encoded as pointer)
          (%cvode-set-user-data cvode-mem (cffi:make-pointer rhs-key))
          ;; Set max steps
          (%cvode-set-max-num-steps cvode-mem max-steps)
          ;; Set a sensible initial step size.
          ;; SUNDIALS' default init-step heuristic divides by ||f(t0, y0)||
          ;; — so it triggers DIV-BY-ZERO whenever the RHS evaluates to
          ;; zero at the start (e.g. dy/dt = t with t0 = 0, dy/dt = 0
          ;; identically, or any model that starts at an instantaneous
          ;; equilibrium of its drift term).  Picking ~1/1000 of the
          ;; first output interval works for any reasonable problem
          ;; (CVODE adapts away from it within the first step or two).
          (let ((h0 (when (>= (length t-list) 2)
                      (* 1.0d-3 (- (coerce (second t-list) 'double-float)
                                   (coerce (first  t-list) 'double-float))))))
            (when (and h0 (> h0 0.0d0))
              (sundials-check-flag
                (%cvode-set-init-step cvode-mem h0)
                "CVodeSetInitStep")))
          ;; Set up event detection if requested.
          ;; SUNDIALS uses the same user_data pointer for both RHS and root
          ;; callbacks. Our root trampoline looks up in *cvode-root-table*
          ;; using the same key as the RHS, so we register the event closure
          ;; under the rhs-key.
          (when (and event-closure (> nrtfn 0))
            (setf (gethash rhs-key *cvode-root-table*) event-closure)
            (sundials-check-flag
              (%cvode-root-init cvode-mem nrtfn
                                (cffi:callback cvode-root-trampoline))
              "CVodeRootInit")
            ;; Apply the per-root direction filter if any non-zero entry
            ;; was specified.  SUNDIALS ignores the call when every dir
            ;; is 0 (its default) — and skipping it keeps the behaviour
            ;; identical to pre-direction releases for callers that
            ;; don't pass `rootdir'.
            (when (and rootdir (some (lambda (d) (not (zerop d))) rootdir))
              (cffi:with-foreign-object (rd-ptr :int nrtfn)
                (loop for d in rootdir
                      for i from 0
                      do (setf (cffi:mem-aref rd-ptr :int i) (coerce d 'integer)))
                (sundials-check-flag
                  (%cvode-set-root-direction cvode-mem rd-ptr)
                  "CVodeSetRootDirection"))))
          ;; Store initial condition in row 0
          (setf (aref result-storage 0)
                (coerce (first t-list) 'double-float))
          (loop for j from 0 below neq
                do (setf (aref result-storage (* (1+ j) n-times))
                         (nth j y0-list)))
          ;; Step through output times
          (cffi:with-foreign-object (tret :double)
            (loop for idx from 1 below n-times
                  for tout-d = (coerce (nth idx t-list) 'double-float)
                  do (let ((keep-stepping t))
                       (loop while keep-stepping do
                         (let ((flag (%cvode-solve
                                       cvode-mem tout-d y-nv tret
                                       +cv-normal+)))
                           (cond
                             ;; Normal return — reached tout
                             ((= flag +cv-success+)
                              (setf keep-stepping nil))
                             ;; Root found — record event and continue.
                             ;; CVodeGetRootInfo returns +1/-1 per root
                             ;; that fired (sign of the crossing) and 0
                             ;; for roots that didn't.  Capture both the
                             ;; index list (which roots fired) and the
                             ;; direction list (in the same order, only
                             ;; for fired roots) so callers can filter
                             ;; spurious detections in user space.
                             ((= flag +cv-root-return+)
                              (let ((t-event (cffi:mem-ref tret :double))
                                    (y-event (nvector-to-list y-nv)))
                                (cffi:with-foreign-object
                                    (roots-found :int nrtfn)
                                  (%cvode-get-root-info cvode-mem roots-found)
                                  (let ((indices  '())
                                        (dirs     '()))
                                    (loop for i from 0 below nrtfn
                                          for r = (cffi:mem-aref
                                                    roots-found :int i)
                                          unless (zerop r)
                                            do (push i indices)
                                            and do (push r dirs))
                                    (push (list t-event y-event
                                                (nreverse indices)
                                                (nreverse dirs))
                                          events-collected))))
                              ;; Continue stepping toward tout
                              )
                             (t
                              (sundials-check-flag flag "CVode"))))))
                  ;; Store row
                  do (let ((t-actual (cffi:mem-ref tret :double)))
                       (setf (aref result-storage idx) t-actual)
                       (let ((y-ptr (%nv-get-array-pointer y-nv)))
                         (loop for j from 0 below neq
                               do (setf (aref result-storage
                                              (+ idx (* (1+ j) n-times)))
                                        (cffi:mem-aref y-ptr :double j)))))))
          ;; Build return value
          (let ((traj (maxima::numerics-wrap
                        (numerics:make-ndarray result))))
            (if (> nrtfn 0)
                ;; Events were requested: always return [traj, events_list]
                ;; Each event entry is
                ;;   [t_event, [y...], [event_indices], [directions]]
                ;; where directions parallels indices (same length, in
                ;; the same order) and contains +1 (rising) or -1
                ;; (falling) per fired root.
                `((maxima::mlist)
                  ,traj
                  ((maxima::mlist)
                   ,@(mapcar
                       (lambda (ev)
                         `((maxima::mlist)
                           ,(first  ev)
                           ((maxima::mlist) ,@(second ev))
                           ((maxima::mlist) ,@(third  ev))
                           ((maxima::mlist) ,@(fourth ev))))
                       (nreverse events-collected))))
                traj)))
      ;; Cleanup (unwind-protect cleanup forms)
      (unregister-rhs-callback rhs-key)
      (remhash rhs-key *cvode-root-table*)
      (cffi:with-foreign-object (mem-ptr :pointer)
        (setf (cffi:mem-ref mem-ptr :pointer) cvode-mem)
        (%cvode-free mem-ptr))
      (%sun-linsol-free sun-ls)
      (%sun-mat-destroy sun-mat)
      (%nv-destroy y-nv))))

;;; ================================================================
;;; Expression mode
;;; ================================================================

(defun np-cvode-expression-mode (f vars y0 tspan
                                  &key (rtol 1.0d-8) (atol 1.0d-8)
                                       method events (max-steps 5000)
                                       rootdir)
  "Expression mode: np_cvode([f1,...], [t,y1,...], y0, tspan ...).
   ROOTDIR is an optional Maxima list of -1/0/+1 integers (one per
   event) selecting the crossing direction CVODE should detect."
  (unless (and (listp f) (eq (caar f) 'maxima::mlist))
    (maxima::merror "np_cvode: first argument must be a list of RHS expressions"))
  (unless (and (listp vars) (eq (caar vars) 'maxima::mlist))
    (maxima::merror "np_cvode: second argument must be a variable list [t, y1, ...]"))
  (let* ((neq (length (cdr f)))
         (y0-list (np-cvode-extract-values y0 "y0"))
         (t-list (np-cvode-extract-values tspan "tspan"))
         (lmm (np-cvode-parse-method method))
         (rtol-d (coerce (maxima::$float rtol) 'double-float))
         (atol-d (coerce (maxima::$float atol) 'double-float)))
    ;; Validate dimensions
    (unless (= neq (length y0-list))
      (maxima::merror "np_cvode: ~A equations but y0 has ~A elements"
                      neq (length y0-list)))
    ;; Compile RHS
    (let* ((compiled-rhs (compile nil (maxima::coerce-float-fun f vars)))
           (rhs-closure (make-cvode-rhs-closure-expr compiled-rhs neq))
           ;; Compile event functions if provided
           (nrtfn (if (and events (listp events)
                          (eq (caar events) 'maxima::mlist))
                     (length (cdr events))
                     0))
           (event-closure
             (when (> nrtfn 0)
               (let* ((event-exprs (cdr events))
                      (compiled-events
                        (mapcar (lambda (g-expr)
                                  (compile nil
                                    (maxima::coerce-float-fun
                                      `((maxima::mlist) ,g-expr)
                                      vars)))
                                event-exprs)))
                 (make-cvode-event-closure-expr
                   compiled-events neq nrtfn)))))
      ;; rootdir from Maxima — accept either a Maxima list ((mlist) -1 0 1)
      ;; or nil.  Convert to a Lisp list of integers; complain if its
      ;; length doesn't match nrtfn.
      (let ((rootdir-lisp
              (when (and rootdir
                         (consp rootdir)
                         (consp (car rootdir))
                         (eq (caar rootdir) 'maxima::mlist))
                (let ((rd (mapcar (lambda (d)
                                    (truncate (maxima::$float d)))
                                  (cdr rootdir))))
                  (unless (= (length rd) nrtfn)
                    (maxima::merror "np_cvode: rootdir has ~A entries but ~A events were given"
                                    (length rd) nrtfn))
                  (dolist (d rd)
                    (unless (member d '(-1 0 1))
                      (maxima::merror "np_cvode: rootdir entries must be -1, 0, or +1; got ~A"
                                      d)))
                  rd))))
        (np-cvode-core rhs-closure neq y0-list t-list lmm rtol-d atol-d
                        :event-closure event-closure
                        :nrtfn nrtfn
                        :max-steps max-steps
                        :rootdir rootdir-lisp)))))

;;; ================================================================
;;; Function mode
;;; ================================================================

(defun np-cvode-function-mode (f-func y0 tspan
                                &key (rtol 1.0d-8) (atol 1.0d-8)
                                     method (max-steps 5000))
  "Function mode: np_cvode(f_func, y0, tspan ...)."
  (let* ((y0-list (np-cvode-extract-values y0 "y0"))
         (neq (length y0-list))
         (t-list (np-cvode-extract-values tspan "tspan"))
         (lmm (np-cvode-parse-method method))
         (rtol-d (coerce (maxima::$float rtol) 'double-float))
         (atol-d (coerce (maxima::$float atol) 'double-float))
         (rhs-closure (make-cvode-rhs-closure-fn f-func neq)))
    (np-cvode-core rhs-closure neq y0-list t-list lmm rtol-d atol-d
                    :max-steps max-steps)))

;;; ================================================================
;;; Stateful API: persistent CVODE handles with ReInit support
;;; ================================================================
;;;
;;; np_cvode is single-shot — every call builds a fresh cvode-mem,
;;; integrates a tspan, and tears everything down.  For event-driven
;;; loops (bouncing ball, switched dynamics, hybrid controllers) the
;;; consumer pays the setup cost on each segment, even though the
;;; only thing changing is the initial state.  SUNDIALS exposes
;;; CVodeReInit precisely for this: keep cvode-mem alive, just reset
;;; the IC.  We expose it through a small lifecycle API:
;;;
;;;   handle : np_cvode_create(rhs, vars, y0, t0,
;;;                            rtol, atol, method, events, max_steps, rootdir)$
;;;   [t, y, evs] : np_cvode_step(handle, t_target)$
;;;   np_cvode_reinit(handle, t_new, y_new)$
;;;   [t, y, evs] : np_cvode_step(handle, next_target)$
;;;   ...
;;;   np_cvode_close(handle)$
;;;
;;; Handles are integer keys into *cvode-handle-table*; resources are
;;; freed by np_cvode_close.  The API is intentionally explicit about
;;; lifecycle — Maxima's GC doesn't reach into the C heap, so users
;;; must call np_cvode_close (typically inside an unwind-protect-style
;;; pattern) to avoid leaking SUNDIALS objects.

(defstruct cvode-handle
  cvode-mem
  ctx
  y-nv
  sun-mat
  sun-ls
  rhs-key
  neq
  nrtfn
  max-steps
  has-events-p)

;;; Two layers wrapping a CVODE handle:
;;;
;;;   *cvode-handle-table*  — integer key → cvode-handle struct
;;;     (holds the foreign pointers and rhs-key; freed by
;;;     cvode-handle-free).  Keyed by integer because integers are
;;;     cheap to pass around as Maxima atoms.
;;;
;;;   cvode-ref struct      — the user-visible Lisp object that
;;;     anchors GC reachability.  When Maxima drops its last reference
;;;     to the cvode-ref, GC eventually runs a finalizer that calls
;;;     cvode-handle-free with the closed-over key, so resources get
;;;     reclaimed even if the user forgot np_cvode_close.
;;;
;;; The finalizer closes over the integer key, NOT the cvode-ref —
;;; otherwise the closure would hold a strong reference to the wrapper,
;;; preventing GC from ever firing the finalizer.

(defvar *cvode-handle-table* (make-hash-table)
  "Maps integer keys to cvode-handle structs.")
(defvar *cvode-handle-counter* 0)

(defstruct cvode-ref
  "User-visible wrapper around a CVODE handle's integer key.
   Anchors GC reachability so dropping the wrapper triggers cleanup."
  key)

(defun cvode-handle-register (h)
  (let ((key (incf *cvode-handle-counter*)))
    (setf (gethash key *cvode-handle-table*) h)
    key))

(defun cvode-handle-lookup (key)
  (or (gethash key *cvode-handle-table*)
      (maxima::merror "np_cvode handle ~A not found (already closed?)" key)))

(defun cvode-handle-free (key)
  "Free the foreign resources attached to KEY.  Idempotent: if the
   handle was already closed (or never existed), this is a no-op,
   so manual np_cvode_close + later GC-triggered finalizer is safe."
  (let ((h (gethash key *cvode-handle-table*)))
    (when h
      (unregister-rhs-callback (cvode-handle-rhs-key h))
      (remhash (cvode-handle-rhs-key h) *cvode-root-table*)
      (cffi:with-foreign-object (mem-ptr :pointer)
        (setf (cffi:mem-ref mem-ptr :pointer) (cvode-handle-cvode-mem h))
        (%cvode-free mem-ptr))
      (%sun-linsol-free  (cvode-handle-sun-ls  h))
      (%sun-mat-destroy  (cvode-handle-sun-mat h))
      (%nv-destroy       (cvode-handle-y-nv    h))
      (remhash key *cvode-handle-table*))))

(defun cvode-make-ref-with-finalizer (key)
  "Wrap KEY in a cvode-ref and register a finalizer that frees the
   handle on GC.  Returns the wrapper.  The finalizer closes over
   KEY (an integer, harmless), not over the cvode-ref itself."
  (let ((ref (make-cvode-ref :key key)))
    (trivial-garbage:finalize ref (lambda () (cvode-handle-free key)))
    ref))

(defun cvode-form->key (form)
  "Extract the integer key from a Maxima `(($cvode_handle simp) <ref>)'
   form.  Errors on shape mismatch."
  (cond
    ((and (consp form)
          (consp (car form))
          (eq (caar form) 'maxima::$cvode_handle)
          (consp (cdr form))
          (cvode-ref-p (cadr form)))
     (cvode-ref-key (cadr form)))
    (t
     (maxima::merror "np_cvode: not a valid handle: ~M" form))))

(defun np-cvode-create-internal (rhs-closure neq y0-list t0 lmm rtol atol
                                  &key event-closure (nrtfn 0)
                                       (max-steps 5000)
                                       rootdir)
  "Build a persistent CVODE handle.  Mirrors the setup phase of
   np-cvode-core but doesn't step.  Returns an integer key."
  (let* ((ctx (ensure-sun-context))
         (cvode-mem (%cvode-create lmm ctx))
         (y-nv (make-nvector neq y0-list))
         (sun-mat (%sun-dense-matrix neq neq ctx))
         (sun-ls (%sun-linsol-dense y-nv sun-mat ctx))
         (rhs-key (register-rhs-callback rhs-closure))
         (h (make-cvode-handle :cvode-mem cvode-mem
                                :ctx ctx
                                :y-nv y-nv
                                :sun-mat sun-mat
                                :sun-ls sun-ls
                                :rhs-key rhs-key
                                :neq neq
                                :nrtfn nrtfn
                                :max-steps max-steps
                                :has-events-p (and event-closure
                                                   (> nrtfn 0)))))
    (when (cffi:null-pointer-p cvode-mem)
      (maxima::merror "np_cvode_create: CVodeCreate failed"))
    (sundials-check-flag
      (%cvode-init cvode-mem
                   (cffi:callback cvode-rhs-trampoline)
                   (coerce t0 'double-float)
                   y-nv)
      "CVodeInit")
    (sundials-check-flag (%cvode-ss-tolerances cvode-mem rtol atol)
                          "CVodeSStolerances")
    (sundials-check-flag (%cvode-set-linear-solver cvode-mem sun-ls sun-mat)
                          "CVodeSetLinearSolver")
    (%cvode-set-user-data cvode-mem (cffi:make-pointer rhs-key))
    (%cvode-set-max-num-steps cvode-mem max-steps)
    ;; Same init-step nudge as np-cvode-core; it's a one-time setup so
    ;; we can't size it from a `next time' here.  Pick something small
    ;; and let CVODE adapt.
    (sundials-check-flag (%cvode-set-init-step cvode-mem 1.0d-6)
                          "CVodeSetInitStep")
    (when (and event-closure (> nrtfn 0))
      (setf (gethash rhs-key *cvode-root-table*) event-closure)
      (sundials-check-flag
        (%cvode-root-init cvode-mem nrtfn
                           (cffi:callback cvode-root-trampoline))
        "CVodeRootInit")
      (when (and rootdir (some (lambda (d) (not (zerop d))) rootdir))
        (cffi:with-foreign-object (rd-ptr :int nrtfn)
          (loop for d in rootdir
                for i from 0
                do (setf (cffi:mem-aref rd-ptr :int i)
                         (coerce d 'integer)))
          (sundials-check-flag
            (%cvode-set-root-direction cvode-mem rd-ptr)
            "CVodeSetRootDirection"))))
    (cvode-handle-register h)))

(defun np-cvode-step-internal (handle t-target)
  "Advance the integrator to T-TARGET, collecting any events along the
   way.  Returns three values: actual t reached, state vector (list),
   and an events-list (list of [t_event y_event indices directions])."
  (let* ((h (cvode-handle-lookup handle))
         (cvode-mem (cvode-handle-cvode-mem h))
         (y-nv (cvode-handle-y-nv h))
         (neq (cvode-handle-neq h))
         (nrtfn (cvode-handle-nrtfn h))
         (events-collected nil)
         (t-actual 0.0d0))
    (cffi:with-foreign-object (tret :double)
      (let ((keep-stepping t))
        (loop while keep-stepping do
          (let ((flag (%cvode-solve cvode-mem
                                     (coerce t-target 'double-float)
                                     y-nv tret +cv-normal+)))
            (cond
              ((= flag +cv-success+)
               (setf t-actual (cffi:mem-ref tret :double))
               (setf keep-stepping nil))
              ((= flag +cv-root-return+)
               (let ((t-event (cffi:mem-ref tret :double))
                     (y-event (nvector-to-list y-nv)))
                 (cffi:with-foreign-object (roots-found :int nrtfn)
                   (%cvode-get-root-info cvode-mem roots-found)
                   (let ((indices '())
                         (dirs    '()))
                     (loop for i from 0 below nrtfn
                           for r = (cffi:mem-aref roots-found :int i)
                           unless (zerop r)
                             do (push i indices)
                             and do (push r dirs))
                     (push (list t-event y-event
                                 (nreverse indices) (nreverse dirs))
                           events-collected)))))
              (t
               (sundials-check-flag flag "CVode")))))))
    (values t-actual
            (nvector-to-list y-nv)
            (nreverse events-collected))))

(in-package #:maxima)

;;; Stateful API entry points (Maxima-callable).

(defun $np_cvode_create (&rest args)
  "np_cvode_create(f, vars, y0, t0, rtol, atol, method, events, max_steps, rootdir)
   Build a persistent CVODE handle.  All arguments after `t0' are
   optional; defaults match `np_cvode'.

   Returns an opaque handle of the form `(($cvode_handle) <ref>)'.
   Pass it to np_cvode_step / np_cvode_reinit.  For typical use you
   can let it go — handles are reclaimed by GC eventually, or at
   session end.  See np_cvode_close / np_cvode_with_handle for
   explicit-cleanup variants if you're hitting memory pressure in
   long-running scripts."
  (unless (>= (length args) 4)
    (merror "np_cvode_create: expected at least 4 arguments (f, vars, y0, t0)"))
  (let* ((f      (first args))
         (vars   (second args))
         (y0     (third args))
         (t0     (coerce (maxima::$float (fourth args)) 'double-float))
         (rtol-d (coerce (maxima::$float (or (fifth args) 1.0d-8)) 'double-float))
         (atol-d (coerce (maxima::$float (or (sixth args) 1.0d-8)) 'double-float))
         (method (seventh args))
         (events (eighth args))
         (max-steps (or (ninth args) 5000))
         (rootdir (tenth args)))
    (unless (and (listp f) (eq (caar f) 'maxima::mlist))
      (merror "np_cvode_create: first argument must be a list of RHS expressions"))
    (unless (and (listp vars) (eq (caar vars) 'maxima::mlist))
      (merror "np_cvode_create: second argument must be a variable list [t, y1, ...]"))
    (let* ((neq (length (cdr f)))
           (y0-list (numerics-sundials::np-cvode-extract-values y0 "y0"))
           (lmm (numerics-sundials::np-cvode-parse-method method)))
      (unless (= neq (length y0-list))
        (merror "np_cvode_create: ~A equations but y0 has ~A elements"
                neq (length y0-list)))
      (let* ((compiled-rhs (compile nil (coerce-float-fun f vars)))
             (rhs-closure (numerics-sundials::make-cvode-rhs-closure-expr
                            compiled-rhs neq))
             (nrtfn (if (and events (listp events)
                             (eq (caar events) 'maxima::mlist))
                       (length (cdr events))
                       0))
             (event-closure
               (when (> nrtfn 0)
                 (let* ((event-exprs (cdr events))
                        (compiled-events
                          (mapcar (lambda (g-expr)
                                    (compile nil
                                      (coerce-float-fun
                                        `((maxima::mlist) ,g-expr)
                                        vars)))
                                  event-exprs)))
                   (numerics-sundials::make-cvode-event-closure-expr
                     compiled-events neq nrtfn))))
             (rootdir-lisp
               (when (and rootdir
                          (consp rootdir)
                          (consp (car rootdir))
                          (eq (caar rootdir) 'maxima::mlist))
                 (mapcar (lambda (d) (truncate (maxima::$float d)))
                         (cdr rootdir)))))
        (let* ((key (numerics-sundials::np-cvode-create-internal
                      rhs-closure neq y0-list t0 lmm rtol-d atol-d
                      :event-closure event-closure
                      :nrtfn nrtfn
                      :max-steps max-steps
                      :rootdir rootdir-lisp))
               (ref (numerics-sundials::cvode-make-ref-with-finalizer key)))
          `(($cvode_handle simp) ,ref))))))

(defun $np_cvode_step (handle t-target)
  "np_cvode_step(handle, t_target) — integrate to t_target.
   Returns [t_actual, y_list, events_list].
   events_list entries are [t_event, [y...], [event_indices], [directions]]."
  (multiple-value-bind (t-actual y-list events-list)
      (numerics-sundials::np-cvode-step-internal
        (numerics-sundials::cvode-form->key handle) t-target)
    `((mlist) ,t-actual
              ((mlist) ,@y-list)
              ((mlist) ,@(mapcar
                           (lambda (ev)
                             `((mlist)
                               ,(first  ev)
                               ((mlist) ,@(second ev))
                               ((mlist) ,@(third  ev))
                               ((mlist) ,@(fourth ev))))
                           events-list)))))

(defun $np_cvode_reinit (handle t-new y-new)
  "np_cvode_reinit(handle, t_new, y_new) — reset the integrator's
   initial condition.  Keeps the linear solver and event-detection
   setup; multistep history is discarded (which is correct after a
   discrete state change)."
  (let* ((key (numerics-sundials::cvode-form->key handle))
         (h (numerics-sundials::cvode-handle-lookup key))
         (cvode-mem (numerics-sundials::cvode-handle-cvode-mem h))
         (y-nv (numerics-sundials::cvode-handle-y-nv h))
         (neq (numerics-sundials::cvode-handle-neq h))
         (y-list (numerics-sundials::np-cvode-extract-values y-new "y_new")))
    (unless (= neq (length y-list))
      (merror "np_cvode_reinit: y_new has ~A elements but handle expects ~A"
              (length y-list) neq))
    ;; Copy the new state into the persistent N_Vector.
    (let ((y-ptr (numerics-sundials::%nv-get-array-pointer y-nv)))
      (loop for v in y-list
            for i from 0
            do (setf (cffi:mem-aref y-ptr :double i)
                     (coerce v 'double-float))))
    (numerics-sundials::sundials-check-flag
      (numerics-sundials::%cvode-reinit
        cvode-mem
        (coerce (maxima::$float t-new) 'double-float)
        y-nv)
      "CVodeReInit")
    '$done))

(defun $np_cvode_close (handle)
  "np_cvode_close(handle) — free SUNDIALS resources for a handle
   immediately.  After this the handle is invalid; subsequent calls
   error.

   For typical interactive use this is unnecessary — handles are
   reclaimed by GC (eventually) and at session end (for sure).  Call
   close when you want determinism: tight loops creating thousands
   of handles, very large state vectors where each handle is
   non-trivial, or scripts running long enough that bounded GC
   timing matters."
  (numerics-sundials::cvode-handle-free
    (numerics-sundials::cvode-form->key handle))
  '$done)

(defun $np_cvode_with_handle (create-args body-fn)
  "np_cvode_with_handle(create_args, body_fn) — resource-safe wrapper.
   Builds a persistent CVODE handle from create_args, calls body_fn
   with the handle, then closes the handle in an unwind-protect so
   cleanup runs even if body_fn errors.

   Implemented at the Lisp level rather than in Maxima because
   Maxima's `unwind_protect' has historically been unreliable about
   running its cleanup form when the protected expression signals an
   error (versions vary).  Common Lisp's unwind-protect is solid, so
   we delegate the resource discipline here."
  (unless (and (consp create-args) (eq (caar create-args) 'mlist))
    (merror "np_cvode_with_handle: create_args must be a Maxima list"))
  (let ((handle (apply '$np_cvode_create (cdr create-args))))
    (unwind-protect
        (mfuncall body-fn handle)
      (ignore-errors ($np_cvode_close handle)))))

(in-package #:maxima)

(defun $np_cvode (&rest args)
  "Integrate ODE system dy/dt = f(t, y) using SUNDIALS CVODE.

   Expression mode:
     np_cvode([f1,...], [t,y1,...], y0, tspan
              [, rtol, atol, method, events, max_steps, rootdir])

   Function mode:
     np_cvode(f_func, y0, tspan [, rtol, atol, method, max_steps])

   Method: adams (default) or bdf.
   Events: list of g_i(t,y) expressions — solver reports zero crossings.
   rootdir: list of -1, 0, or +1 (one per event) — direction filter.
            -1 = falling crossings only (g goes + → -),
             0 = either direction (default),
            +1 = rising crossings only (g goes - → +).

   Returns: 2D ndarray [n_times, 1+neq], or [trajectory, events_list]
   if events used.  Each event entry is
     [t_event, [y...], [event_indices], [directions]]
   where directions are +1 / -1 indicating which way each fired root
   crossed zero."
  (unless (>= (length args) 3)
    (merror "np_cvode: expected at least 3 arguments"))
  (let ((first-arg (first args)))
    (if (numerics-callable-p first-arg)
        ;; Function mode: np_cvode(f, y0, tspan [, rtol, atol, method, max_steps])
        (numerics-sundials::np-cvode-function-mode
          first-arg (second args) (third args)
          :rtol (or (fourth args) 1.0d-8)
          :atol (or (fifth args) 1.0d-8)
          :method (sixth args)
          :max-steps (or (seventh args) 5000))
        ;; Expression mode: np_cvode(f, vars, y0, tspan
        ;;                           [, rtol, atol, method, events, max_steps, rootdir])
        (progn
          (unless (>= (length args) 4)
            (merror "np_cvode: expression mode requires at least 4 arguments"))
          (numerics-sundials::np-cvode-expression-mode
            first-arg (second args) (third args) (fourth args)
            :rtol (or (fifth args) 1.0d-8)
            :atol (or (sixth args) 1.0d-8)
            :method (seventh args)
            :events (eighth args)
            :max-steps (or (ninth args) 5000)
            :rootdir (tenth args))))))
