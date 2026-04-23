;;; library.lisp — CFFI library loading for SUNDIALS
;;;
;;; Loads SUNDIALS shared libraries at ASDF load time.
;;; The user must install SUNDIALS via their system package manager:
;;;   macOS:   brew install sundials
;;;   Debian:  sudo apt install libsundials-dev
;;;   Fedora:  sudo dnf install sundials-devel
;;;   Windows: vcpkg install sundials
;;;
;;; For non-standard install locations, set LD_LIBRARY_PATH / DYLD_LIBRARY_PATH
;;; or push to cffi:*foreign-library-directories* in maxima-init.mac.

(in-package #:numerics-sundials)

;;; --- Library definitions ---

;; SUNDIALS v7+ has a separate core library; v6 doesn't.
(cffi:define-foreign-library libsundials-core
  (:darwin (:or "libsundials_core.dylib" (:default "libsundials_core")))
  (:unix (:or "libsundials_core.so" (:default "libsundials_core")))
  (:windows (:or "sundials_core.dll" (:default "sundials_core")))
  (t (:default "libsundials_core")))

(cffi:define-foreign-library libsundials-nvecserial
  (:darwin (:or "libsundials_nvecserial.dylib" (:default "libsundials_nvecserial")))
  (:unix (:or "libsundials_nvecserial.so" (:default "libsundials_nvecserial")))
  (:windows (:or "sundials_nvecserial.dll" (:default "sundials_nvecserial")))
  (t (:default "libsundials_nvecserial")))

(cffi:define-foreign-library libsundials-sunmatrixdense
  (:darwin (:or "libsundials_sunmatrixdense.dylib" (:default "libsundials_sunmatrixdense")))
  (:unix (:or "libsundials_sunmatrixdense.so" (:default "libsundials_sunmatrixdense")))
  (:windows (:or "sundials_sunmatrixdense.dll" (:default "sundials_sunmatrixdense")))
  (t (:default "libsundials_sunmatrixdense")))

(cffi:define-foreign-library libsundials-sunlinsoldense
  (:darwin (:or "libsundials_sunlinsoldense.dylib" (:default "libsundials_sunlinsoldense")))
  (:unix (:or "libsundials_sunlinsoldense.so" (:default "libsundials_sunlinsoldense")))
  (:windows (:or "sundials_sunlinsoldense.dll" (:default "sundials_sunlinsoldense")))
  (t (:default "libsundials_sunlinsoldense")))

(cffi:define-foreign-library libsundials-cvode
  (:darwin (:or "libsundials_cvode.dylib" (:default "libsundials_cvode")))
  (:unix (:or "libsundials_cvode.so" (:default "libsundials_cvode")))
  (:windows (:or "sundials_cvode.dll" (:default "sundials_cvode")))
  (t (:default "libsundials_cvode")))

;;; --- Load all libraries ---

(defun load-sundials-libraries ()
  "Load all SUNDIALS shared libraries. Called at ASDF load time."
  (handler-case
      (progn
        ;; v7+ has a separate core library; v6 doesn't — ignore if missing.
        (ignore-errors (cffi:use-foreign-library libsundials-core))
        (cffi:use-foreign-library libsundials-nvecserial)
        (cffi:use-foreign-library libsundials-sunmatrixdense)
        (cffi:use-foreign-library libsundials-sunlinsoldense)
        (cffi:use-foreign-library libsundials-cvode))
    (cffi:load-foreign-library-error (e)
      (maxima::merror
        "numerics-sundials: SUNDIALS library not found.~%~
         Install it:~%  ~
         macOS:   brew install sundials~%  ~
         Debian:  sudo apt install libsundials-dev~%  ~
         Fedora:  sudo dnf install sundials-devel~%  ~
         Windows: vcpkg install sundials~%~
         ~%For non-standard locations, set LD_LIBRARY_PATH or push to~%~
         cffi:*foreign-library-directories* in maxima-init.mac.~%~
         ~%(Error: ~A)" e))))

(eval-when (:load-toplevel :execute)
  (load-sundials-libraries))
