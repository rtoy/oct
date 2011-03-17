;;;; -*- Mode: lisp -*-
;;;;
;;;; Copyright (c) 2007,2011 Raymond Toy
;;;;
;;;; Permission is hereby granted, free of charge, to any person
;;;; obtaining a copy of this software and associated documentation
;;;; files (the "Software"), to deal in the Software without
;;;; restriction, including without limitation the rights to use,
;;;; copy, modify, merge, publish, distribute, sublicense, and/or sell
;;;; copies of the Software, and to permit persons to whom the
;;;; Software is furnished to do so, subject to the following
;;;; conditions:
;;;;
;;;; The above copyright notice and this permission notice shall be
;;;; included in all copies or substantial portions of the Software.
;;;;
;;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
;;;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
;;;; OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
;;;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;;;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
;;;; WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
;;;; FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
;;;; OTHER DEALINGS IN THE SOFTWARE.

(in-package #:oct)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (setf *readtable* *oct-readtable*))

;; For the tests, we need to turn off underflow for clisp.
#+clisp
(ext:without-package-lock ()
  (setq sys::*inhibit-floating-point-underflow* t))

;; Compute how many bits are the same for two numbers EST and TRUE.
;; Return T if they are identical.
(defun bit-accuracy (est true)
  (let* ((diff (abs (- est true)))
	 (err (float (if (zerop true)
			 diff
			 (/ diff (abs true)))
		     1d0)))
    (if (zerop diff)
	t
	(- (log err 2)))))

(defun check-accuracy (limit est true)
  (let ((bits (bit-accuracy est true)))
    (if (numberp bits)
	(if (< bits limit)
	    (list bits limit est true)))))

(defvar *null* (make-broadcast-stream))

;;; Some simple tests from the Yozo Hida's qd package.

;; Pi via Machin's formula
(rt:deftest oct.pi.machin
    (let* ((*standard-output* *null*)
	   (val (make-instance 'qd-real :value (octi::test2 nil)))
	   (true oct:+pi+))
      (check-accuracy 213 val true))
  nil)

;; Pi via Salamin-Brent algorithm
(rt:deftest oct.pi.salamin-brent
    (let* ((*standard-output* *null*)
	   (val (make-instance 'qd-real :value (octi::test3 nil)))
	   (true oct:+pi+))
      (check-accuracy 202 val true))
  nil)

;; Pi via Borweign's Quartic formula
(rt:deftest oct.pi.borweign
    (let* ((*standard-output* *null*)
	   (val (make-instance 'qd-real :value (octi::test4 nil)))
	   (true oct:+pi+))
      (check-accuracy 211 val true))
  nil)

;; e via Taylor series
(rt:deftest oct.e.taylor
    (let* ((*standard-output* *null*)
	   (val (make-instance 'qd-real :value (octi::test5 nil)))
	   (true (make-instance 'qd-real :value octi::+qd-e+)))
      (check-accuracy 212 val true))
  nil)

;; log(2) via Taylor series
(rt:deftest oct.log2.taylor
    (let* ((*standard-output* *null*)
	   (val (make-instance 'qd-real :value (octi::test6 nil)))
	   (true (make-instance 'qd-real :value octi::+qd-log2+)))
      (check-accuracy 212 val true))
  nil)

;;; Tests of atan where we know the analytical result
(rt:deftest oct.atan.1
    (let* ((arg (/ (sqrt #q3)))
	   (y (/ (atan arg) +pi+))
	   (true (/ #q6)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan.2
    (let* ((arg (sqrt #q3))
	 (y (/ (atan arg) +pi+))
	 (true (/ #q3)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan.3
    (let* ((arg #q1)
	 (y (/ (atan arg) +pi+))
	 (true (/ #q4)))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan.4
    (let* ((arg #q1q100)
	   (y (/ (atan arg) +pi+))
	   (true #q.5))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan.5
    (let* ((arg #q-1q100)
	   (y (/ (atan arg) +pi+))
	   (true #q-.5))
    (check-accuracy 212 y true))
  nil)


(defun atan-qd/duplication (arg)
  (make-instance 'qd-real
		 :value (octi::atan-qd/duplication (qd-value arg))))

;;; Tests of atan where we know the analytical result.  Same tests,
;;; but using the atan duplication formula.
(rt:deftest oct.atan/dup.1
    (let* ((arg (/ (sqrt #q3)))
	   (y (/ (atan-qd/duplication arg) +pi+))
	   (true (/ #q6)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/dup.2
    (let* ((arg (sqrt #q3))
	   (y (/ (atan-qd/duplication arg) +pi+))
	   (true (/ #q3)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/dup.3
    (let* ((arg #q1)
	   (y (/ (atan-qd/duplication arg) +pi+))
	   (true (/ #q4)))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/dup.4
    (let* ((arg #q1q100)
	   (y (/ (atan-qd/duplication arg) +pi+))
	   (true #q.5))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/dup.5
    (let* ((arg #q-1q100)
	   (y (/ (atan-qd/duplication arg) +pi+))
	   (true #q-.5))
    (check-accuracy 212 y true))
  nil)

;;; Tests of atan where we know the analytical result.  Same tests,
;;; but using a CORDIC implementation.
(defun atan-qd/cordic (arg)
  (make-instance 'qd-real
		 :value (octi::atan-qd/cordic (qd-value arg))))

(rt:deftest oct.atan/cordic.1
    (let* ((arg (/ (sqrt #q3)))
	   (y (/ (atan-qd/cordic arg) +pi+))
	   (true (/ #q6)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/cordic.2
    (let* ((arg (sqrt #q3))
	   (y (/ (atan-qd/cordic arg) +pi+))
	   (true (/ #q3)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/cordic.3
    (let* ((arg #q1)
	   (y (/ (atan-qd/cordic arg) +pi+))
	   (true (/ #q4)))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/cordic.4
    (let* ((arg #q1q100)
	   (y (/ (atan-qd/cordic arg) +pi+))
	   (true #q.5))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.atan/cordic.5
    (let* ((arg #q-1q100)
	   (y (/ (atan-qd/cordic arg) +pi+))
	   (true #q-.5))
    (check-accuracy 212 y true))
  nil)


;;; Tests of sin where we know the analytical result.
(rt:deftest oct.sin.1
    (let* ((arg (/ +pi+ 6))
	   (y (sin arg))
	   (true #q.5))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.sin.2
    (let* ((arg (/ +pi+ 4))
	   (y (sin arg))
	   (true (sqrt #q.5)))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.sin.3
    (let* ((arg (/ +pi+ 3))
	   (y (sin arg))
	   (true (/ (sqrt #q3) 2)))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.big-sin.1
    (let* ((arg (oct:make-qd (ash 1 120)))
	   (y (sin arg))
	   (true #q3.778201093607520226555484700569229919605866976512306642257987199414885q-1))
      (check-accuracy 205 y true))
  nil)

(rt:deftest oct.big-sin.2
    (let* ((arg (oct:make-qd (ash 1 1023)))
	   (y (sin arg))
	   (true #q5.631277798508840134529434079444683477103854907361251399182750155357133q-1))
      (check-accuracy 205 y true))
  nil)

;;; Tests of tan where we know the analytical result.
(rt:deftest oct.tan.1
    (let* ((arg (/ +pi+ 6))
	 (y (tan arg))
	 (true (/ (sqrt #q3))))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.tan.2
    (let* ((arg (/ +pi+ 4))
	 (y (tan arg))
	 (true #q1))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.tan.3
  (let* ((arg (/ +pi+ 3))
	 (y (tan arg))
	 (true (sqrt #q3)))
    (check-accuracy 212 y true))    
  nil)

;;; Tests of tan where we know the analytical result.  Uses CORDIC
;;; algorithm.

(defun tan/cordic (arg)
  (make-instance 'qd-real
		 :value (octi::tan-qd/cordic (qd-value arg))))

(rt:deftest oct.tan/cordic.1
    (let* ((arg (/ +pi+ 6))
	 (y (tan/cordic arg))
	 (true (/ (sqrt #q3))))
    (check-accuracy 211 y true))
  nil)

(rt:deftest oct.tan/cordic.2
    (let* ((arg (/ +pi+ 4))
	 (y (tan/cordic arg))
	 (true #q1))
    (check-accuracy 211 y true))
  nil)

(rt:deftest oct.tan/cordic.3
  (let* ((arg (/ +pi+ 3))
	 (y (tan/cordic arg))
	 (true (sqrt #q3)))
    (check-accuracy 210 y true))    
  nil)

;;; Tests of asin where we know the analytical result.

(rt:deftest oct.asin.1
    (let* ((arg #q.5)
	 (y (asin arg))
	 (true (/ +pi+ 6)))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.asin.2
    (let* ((arg (sqrt #q.5))
	   (y (asin arg))
	   (true (/ +pi+ 4)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.asin.3
    (let* ((arg (/ (sqrt #q3) 2))
	   (y (asin arg))
	   (true (/ +pi+ 3)))
      (check-accuracy 212 y true))
  nil)

;;; Tests of log.

(rt:deftest oct.log.1
    (let* ((arg #q2)
	   (y (log arg))
	   (true (make-instance 'qd-real :value octi::+qd-log2+)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.log.2
    (let* ((arg #q10)
	 (y (log arg))
	 (true (make-instance 'qd-real :value octi::+qd-log10+)))
    (check-accuracy 207 y true))
  nil)

(rt:deftest oct.log.3
    (let* ((arg (+ 1 (scale-float #q1 -80)))
	   (y (log arg))
	   (true #q8.2718061255302767487140834995607996176476940491239977084112840149578911975528492q-25))
      (check-accuracy 212 y true))
  nil)

;;; Tests of log using Newton iteration.

(defun log/newton (arg)
  (make-instance 'qd-real
		 :value (octi::log-qd/newton (qd-value arg))))

(rt:deftest oct.log/newton.1
    (let* ((arg #q2)
	   (y (log/newton arg))
	   (true (make-instance 'qd-real :value octi::+qd-log2+)))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.log/newton.2
    (let* ((arg #q10)
	 (y (log/newton arg))
	 (true (make-instance 'qd-real :value octi::+qd-log10+)))
    (check-accuracy 207 y true))
  nil)

(rt:deftest oct.log/newton.3
    (let* ((arg (+ 1 (scale-float #q1 -80)))
	   (y (log/newton arg))
	   (true #q8.2718061255302767487140834995607996176476940491239977084112840149578911975528492q-25))
      (check-accuracy 212 y true))
  nil)

;;; Tests of log using AGM.

(defun log/agm (arg)
  (make-instance 'qd-real
		 :value (octi::log-qd/agm (qd-value arg))))

(rt:deftest oct.log/agm.1
    (let* ((arg #q2)
	   (y (log/agm arg))
	   (true (make-instance 'qd-real :value octi::+qd-log2+)))
      (check-accuracy 203 y true))
  nil)

(rt:deftest oct.log/agm.2
    (let* ((arg #q10)
	 (y (log/agm arg))
	 (true (make-instance 'qd-real :value octi::+qd-log10+)))
    (check-accuracy 205 y true))
  nil)

(rt:deftest oct.log/agm.3
    (let* ((arg (+ 1 (scale-float #q1 -80)))
	   (y (log/agm arg))
	   (true #q8.2718061255302767487140834995607996176476940491239977084112840149578911975528492q-25))
      (check-accuracy 123 y true))
  nil)

;;; Tests of log using AGM2, a faster variaton of AGM.

(defun log/agm2 (arg)
  (make-instance 'qd-real
		 :value (octi::log-qd/agm2 (qd-value arg))))

(rt:deftest oct.log/agm2.1
    (let* ((arg #q2)
	   (y (log/agm2 arg))
	   (true (make-instance 'qd-real :value octi::+qd-log2+)))
      (check-accuracy 203 y true))
  nil)

(rt:deftest oct.log/agm2.2
    (let* ((arg #q10)
	 (y (log/agm2 arg))
	 (true (make-instance 'qd-real :value octi::+qd-log10+)))
    (check-accuracy 205 y true))
  nil)

(rt:deftest oct.log/agm2.3
    (let* ((arg (+ 1 (scale-float #q1 -80)))
	   (y (log/agm2 arg))
	   (true #q8.2718061255302767487140834995607996176476940491239977084112840149578911975528492q-25))
      (check-accuracy 123 y true))
  nil)

;;; Tests of log using AGM3, a faster variation of AGM2.
(defun log/agm3 (arg)
  (make-instance 'qd-real
		 :value (octi::log-qd/agm3 (qd-value arg))))

(rt:deftest oct.log/agm3.1
    (let* ((arg #q2)
	   (y (log/agm3 arg))
	   (true (make-instance 'qd-real :value octi::+qd-log2+)))
      (check-accuracy 203 y true))
  nil)

(rt:deftest oct.log/agm3.2
    (let* ((arg #q10)
	 (y (log/agm3 arg))
	 (true (make-instance 'qd-real :value octi::+qd-log10+)))
    (check-accuracy 205 y true))
  nil)

(rt:deftest oct.log/agm3.3
    (let* ((arg (+ 1 (scale-float #q1 -80)))
	   (y (log/agm3 arg))
	   (true #q8.2718061255302767487140834995607996176476940491239977084112840149578911975528492q-25))
      (check-accuracy 123 y true))
  nil)

;;; Tests of sqrt to make sure we overflow or underflow where we
;;; shouldn't.

(rt:deftest oct.sqrt.1
    (let* ((arg #q1q200)
	   (y (sqrt arg))
	   (true #q1q100))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.sqrt.2
    (let* ((arg #q1q200)
	   (y (sqrt arg))
	   (true #q1q100))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.sqrt.3
    (let* ((arg #q1q300)
	   (y (sqrt arg))
	   (true #q1q150))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.sqrt.4
    (let* ((arg #q1q-200)
	   (y (sqrt arg))
	   (true #q1q-100))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.sqrt.5
    (let* ((arg #q1q-250)
	   (y (sqrt arg))
	   (true #q1q-125))
      (check-accuracy 212 y true))
  nil)

;;; Tests of log1p(x) = log(1+x), using the duplication formula.

(defun log1p/dup (arg)
  (make-instance 'qd-real
		 :value (octi::log1p-qd/duplication (qd-value arg))))

(rt:deftest oct.log1p.1
    (let* ((arg #q9)
	   (y (log1p/dup arg))
	   (true #q2.3025850929940456840179914546843642076011014886287729760333279009675726096773525q0))
      (check-accuracy 212 y true))
  nil)

(rt:deftest oct.log1p.2
    (let* ((arg (scale-float #q1 -80))
	   (y (log1p/dup arg))
	   (true #q8.2718061255302767487140834995607996176476940491239977084112840149578911975528492q-25))
      (check-accuracy 212 y true))
  nil)

;;; Tests of expm1(x) = exp(x) - 1, using a Taylor series with
;;; argument reduction.

(defun expm1/series (arg)
  (make-instance 'qd-real
		 :value (octi::expm1-qd/series (qd-value arg))))

(rt:deftest oct.expm1/series.1
  (let* ((arg #q0)
	 (y (expm1/series arg))
	 (true #q0))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.expm1/series.2
  (let* ((arg #q1)
	 (y (expm1/series arg))
	 (true #q1.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952q0))
    (check-accuracy 211 y true))
  nil)

(rt:deftest oct.expm1/series.3
    (let* ((arg (scale-float #q1 -100))
	   (y (expm1/series arg))
	   (true #q7.888609052210118054117285652830973804370994921943802079729680186943164342372119432861876389514693341738324702996270767390039172777809233288470357147q-31))
      (check-accuracy 211 y true))
  nil)

;;; Tests of expm1(x) = exp(x) - 1, using duplication formula.

(defun expm1/dup (arg)
  (make-instance 'qd-real
		 :value (octi::expm1-qd/duplication (qd-value arg))))


(rt:deftest oct.expm1/dup.1
  (let* ((arg #q0)
	 (y (expm1/dup arg))
	 (true #q0))
    (check-accuracy 212 y true))
  nil)

(rt:deftest oct.expm1/dup.2
  (let* ((arg #q1)
	 (y (expm1/dup arg))
	 (true #q1.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952q0))
    (check-accuracy 211 y true))
  nil)

(rt:deftest oct.expm1/dup.3
    (let* ((arg (scale-float #q1 -100))
	   (y (expm1/dup arg))
	   (true #q7.888609052210118054117285652830973804370994921943802079729680186943164342372119432861876389514693341738324702996270767390039172777809233288470357147q-31))
      (check-accuracy 211 y true))
  nil)

;; If we screw up integer-decode-qd, printing is wrong.  Here is one
;; case where integer-decode-qd was screwed up and printing the wrong
;; thing.
(rt:deftest oct.integer-decode.1
    (multiple-value-bind (frac exp s)
	(octi:integer-decode-qd (octi::%make-qd-d -0.03980126756814893d0
						-2.7419792323327893d-18
						0d0 0d0))
      (unless (and (eql frac 103329998279901916046530991816704)
		   (eql exp -111)
		   (eql s -1))
	(list frac exp s)))
  nil)
      
;;;
;;; Add a few tests for the branch cuts.  Many of these tests assume
;;; that Lisp has support for signed zeroes.  If not, these tests are
;;; probably wrong.

(defun check-signs (fun arg expected)
  (let* ((z (funcall fun arg))
	 (x (realpart z))
	 (y (imagpart z)))
    ;; If the Lisp doesn't support signed zeroes, then this test
    ;; should always pass.
    (if (or (eql -0d0 0d0)
	    (and (= (float-sign x) (float-sign (realpart expected)))
		 (= (float-sign y) (float-sign (imagpart expected)))))
	t
	(list z expected fun arg))))
      
;; asin has a branch cut on the real axis |x|>1.  For x < -1, it is
;; continuous with quadrant II; for x > 1, continuous with quadrant
;; IV.
(rt:deftest oct.asin-branch-neg.1
    (let ((true (cl:asin #c(-2d0 1d-20))))
      (check-signs #'asin -2d0 true))
  t)

(rt:deftest oct.asin-branch-neg.2
    (let ((true (cl:asin #c(-2d0 1d-20))))
      (check-signs #'asin #q-2 true))
  t)

(rt:deftest oct.asin-branch-neg.3
    (let ((true (cl:asin #c(-2d0 1d-20))))
      (check-signs #'asin #c(-2d0 0d0) true))
  t)

(rt:deftest oct.asin-branch-neg.4
    (let ((true (cl:asin #c(-2d0 1d-20))))
      (check-signs #'asin #q(-2 0) true))
  t)

(rt:deftest oct.asin-branch-neg.5
    (let ((true (cl:asin #c(-2d0 1d-20))))
      (check-signs #'asin #c(-2d0 -0d0) (conjugate true)))
  t)

(rt:deftest oct.asin-branch-neg.6
    (let ((true (cl:asin #c(-2d0 1d-20))))
      (check-signs #'asin #q(-2d0 -0d0) (conjugate true)))
  t)

(rt:deftest oct.asin-branch-pos.1
    (let ((true (cl:asin #c(2d0 -1d-20))))
      (check-signs #'asin #c(2d0 0d0) (conjugate true)))
  t)

(rt:deftest oct.asin-branch-pos.2
    (let ((true (cl:asin #c(2d0 -1d-20))))
      (check-signs #'asin #q(2 0d0) (conjugate true)))
  t)

(rt:deftest oct.asin-branch-pos.3
    (let ((true (cl:asin #c(2d0 -1d-20))))
      (check-signs #'asin #c(2d0 -0d0) true))
  t)

(rt:deftest oct.asin-branch-pos.4
    (let ((true (cl:asin #c(2d0 -1d-20))))
      (check-signs #'asin #q(2d0 -0d0) true))
  t)

;; acos branch cut is the real axis, |x| > 1.  For x < -1, it is
;; continuous with quadrant II; for x > 1, quadrant IV.

(rt:deftest oct.acos-branch-neg.1
    (let ((true (cl:acos #c(-2d0 1d-20))))
      (check-signs #'acos -2d0 true))
  t)

(rt:deftest oct.acos-branch-neg.2
    (let ((true (cl:acos #c(-2d0 1d-20))))
      (check-signs #'acos #q-2 true))
  t)

(rt:deftest oct.acos-branch-neg.3
    (let ((true (cl:acos #c(-2d0 1d-20))))
      (check-signs #'acos #c(-2d0 0d0) true))
  t)

(rt:deftest oct.acos-branch-neg.4
    (let ((true (cl:acos #c(-2d0 1d-20))))
      (check-signs #'acos #q(-2 0) true))
  t)

(rt:deftest oct.acos-branch-neg.5
    (let ((true (cl:acos #c(-2d0 1d-20))))
      (check-signs #'acos #c(-2d0 -0d0) (conjugate true)))
  t)

(rt:deftest oct.acos-branch-neg.6
    (let ((true (cl:acos #c(-2d0 1d-20))))
      (check-signs #'acos #q(-2d0 -0d0) (conjugate true)))
  t)

(rt:deftest oct.acos-branch-pos.1
    (let ((true (cl:acos #c(2d0 -1d-20))))
      (check-signs #'acos #c(2d0 0d0) (conjugate true)))
  t)

(rt:deftest oct.acos-branch-pos.2
    (let ((true (cl:acos #c(2d0 -1d-20))))
      (check-signs #'acos #q(2 0d0) (conjugate true)))
  t)

(rt:deftest oct.acos-branch-pos.3
    (let ((true (cl:acos #c(2d0 -1d-20))))
      (check-signs #'acos #c(2d0 -0d0) true))
  t)

(rt:deftest oct.acos-branch-pos.4
    (let ((true (cl:acos #c(2d0 -1d-20))))
      (check-signs #'acos #q(2d0 -0d0) true))
  t)

;; atan branch cut is the imaginary axis, |y| > 1.  For y < -1, it is
;; continuous with quadrant IV; for x > 1, quadrant II.

(rt:deftest oct.atan-branch-neg.1
    (let ((true (cl:atan #c(1d-20 -2d0))))
      (check-signs #'atan #c(0d0 -2d0) true))
  t)

(rt:deftest oct.atan-branch-neg.2
    (let ((true (cl:atan #c(1d-20 -2d0))))
      (check-signs #'atan #q(0 -2) true))
  t)

(rt:deftest oct.atan-branch-neg.3
    (let ((true (cl:atan #c(-1d-20 -2d0))))
      (check-signs #'atan #c(-0d0 -2d0) true))
  t)

(rt:deftest oct.atan-branch-neg.4
    (let ((true (cl:atan #c(-1d-20 -2d0))))
      (check-signs #'atan #q(-0d0 -2d0) true))
  t)

(rt:deftest oct.atan-branch-pos.1
    (let ((true (cl:atan #c(1d-20 2d0))))
      (check-signs #'atan #c(0d0 2d0) true))
  t)

(rt:deftest oct.atan-branch-pos.2
    (let ((true (cl:atan #c(1d-20 2d0))))
      (check-signs #'atan #q(0d0 2d0) true))
  t)

(rt:deftest oct.atan-branch-pos.3
    (let ((true (cl:atan #c(-1d-20 2d0))))
      (check-signs #'atan #c(-0d0 2d0) true))
  t)

(rt:deftest oct.atan-branch-pos.4
    (let ((true (cl:atan #c(-1d-20 2d0))))
      (check-signs #'atan #q(-0d0 2d0) true))
  t)

;; Test x < -1
(rt:deftest oct.atanh-branch-neg.1
    (let ((true (cl:atanh #c(-2d0 -1d-20))))
      (check-signs #'atanh -2d0 true))
  t)

(rt:deftest oct.atanh-branch-neg.2
    (let ((true (cl:atanh #c(-2d0 -1d-20))))
      (check-signs #'atanh #q-2 true))
  t)

;; Test x > 1
(rt:deftest oct.atanh-branch-pos.1
    (let ((true (cl:atanh #c(2d0 1d-20))))
      (check-signs #'atanh 2d0 true))
  t)

(rt:deftest oct.atanh-branch-pos.2
    (let ((true (cl:atanh #c(2d0 1d-20))))
      (check-signs #'atanh #q2 true))
  t)

;; elliptic_k(-1) = gamma(1/4)^2/2^(5/2)/sqrt(%pi)
(rt:deftest oct.elliptic-k.1d
    (let* ((val (elliptic-k -1d0))
	   (true #q1.311028777146059905232419794945559706841377475715811581408410851900395293535207125115147766480714547q0))
      (check-accuracy 53 val true))
  nil)

(rt:deftest oct.elliptic-k.1q
    (let* ((val (elliptic-k #q-1q0))
	   (true #q1.311028777146059905232419794945559706841377475715811581408410851900395293535207125115147766480714547q0))
      (check-accuracy 210 val true))
  nil)

;; elliptic_k(1/2) = %pi^(3/2)/2/gamma(3/4)^2
(rt:deftest oct.elliptic-k.2d
    (let* ((val (elliptic-k 0.5d0))
	   (true #q1.854074677301371918433850347195260046217598823521766905585928045056021776838119978357271861650371897q0))
      (check-accuracy 53 val true))
  nil)

(rt:deftest oct.elliptic-k.2q
    (let* ((val (elliptic-k #q.5))
	   (true #q1.854074677301371918433850347195260046217598823521766905585928045056021776838119978357271861650371897q0))
      (check-accuracy 210 val true))
  nil)

;; jacobi_sn(K,1/2) = 1, where K = elliptic_k(1/2)
(rt:deftest oct.jacobi-sn.1d
    (let* ((ek (elliptic-k .5d0))
	   (val (jacobi-sn ek .5d0)))
      (check-accuracy 54 val 1d0))
  nil)

(rt:deftest oct.jacobi-sn.1q
    (let* ((ek (elliptic-k #q.5))
	   (val (jacobi-sn ek #q.5)))
      (check-accuracy 212 val #q1))
  nil)

;; jacobi_cn(K,1/2) = 0
(rt:deftest oct.jacobi-cn.1d
    (let* ((ek (elliptic-k .5d0))
	   (val (jacobi-cn ek .5d0)))
      (check-accuracy 50 val 0d0))
  nil)

(rt:deftest oct.jacobi-sn.1q
    (let* ((ek (elliptic-k #q.5))
	   (val (jacobi-cn ek #q.5)))
      (check-accuracy 210 val #q0))
  nil)

;; jacobi-dn(K, 1/2) = sqrt(1/2)
(rt:deftest oct.jacobi-dn.1d
    (let* ((ek (elliptic-k .5d0))
	   (true (sqrt .5d0))
	   (val (jacobi-dn ek .5d0)))
      (check-accuracy 52 val true))
  nil)

(rt:deftest oct.jacobi-dn.1q
    (let* ((ek (elliptic-k #q.5))
	   (true (sqrt #q.5))
	   (val (jacobi-dn ek #q.5)))
      (check-accuracy 212 val true))
  nil)

(rt:deftest oct.carlson-rf.1d
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0d0 2d0 1d0))
	  (true 1.31102877714605990523241979494d0))
      (check-accuracy 53 rf true))
  nil)

(rt:deftest oct.carlson-rf.1q
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    (let ((rf (carlson-rf #q0 #q2 #q1))
	  (true #q1.311028777146059905232419794945559706841377475715811581408410851900395q0))
      (check-accuracy 212 rf true))
  nil)

(rt:deftest oct.carlson-rd.1d
    ;; Rd(0,2,1) = 3*integrate(s^2/sqrt(1-s^4), s, 0 ,1)
    ;;             = 3*beta(3/4,1/2)/4
    ;;             = 3*sqrt(%pi)*gamma(3/4)/gamma(1/4)
    (let ((rd (carlson-rd 0d0 2d0 1d0))
	  (true 1.7972103521033883d0))
      (check-accuracy 51 rd true))
  nil)

(rt:deftest oct.carlson-rd.1q
    (let ((rd (carlson-rd #q0 #q2 #q1))
	  (true #q1.797210352103388311159883738420485817340818994823477337395512429419599q0))
      (check-accuracy 212 rd true))
  nil)

;; Test some of the contagion stuff.

(rt:deftest oct.carlson-rf.contagion.1
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0 2 1))
	  (true 1.31102877714605990523241979494d0))
      (check-accuracy 23 rf true))
  nil)

(rt:deftest oct.carlson-rf.contagion.1d
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0d0 2 1))
	  (true 1.31102877714605990523241979494d0))
      (check-accuracy 53 rf true))
  nil)

(rt:deftest oct.carlson-rf.contagion.2d
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0 2d0 1))
	  (true 1.31102877714605990523241979494d0))
      (check-accuracy 53 rf true))
  nil)

(rt:deftest oct.carlson-rf.contagion.3d
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0 2 1d0))
	  (true 1.31102877714605990523241979494d0))
      (check-accuracy 53 rf true))
  nil)

(rt:deftest oct.carlson-rf.contagion.1q
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf #q0q0 2 1))
	  (true #q1.311028777146059905232419794945559706841377475715811581408410851900395q0))
      (check-accuracy 212 rf true))
  nil)

(rt:deftest oct.carlson-rf.contagion.2q
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0 #q2q0 1))
	  (true #q1.311028777146059905232419794945559706841377475715811581408410851900395q0))
      (check-accuracy 212 rf true))
  nil)

(rt:deftest oct.carlson-rf.contagion.3q
    ;; Rf(0,2,1) = integrate(1/sqrt(1-s^4), s, 0 ,1)
    ;;           = 1/4*beta(1/2,1/2)
    ;;           = sqrt(%pi)/4*gamma(1/4)/gamma(3/4)
    (let ((rf (carlson-rf 0 2 #q1q0))
	  (true #q1.311028777146059905232419794945559706841377475715811581408410851900395q0))
      (check-accuracy 212 rf true))
  nil)

;; Elliptic integral of the third kind

;; elliptic-pi(0,phi,m) = elliptic-f(phi, m)
(rt:deftest oct.elliptic-pi.1d
    (loop for k from 0 to 100
       for phi = (random (/ pi 2))
       for m = (random 1d0)
       for epi = (elliptic-pi 0 phi m)
       for ef = (elliptic-f phi m)
       for result = (check-accuracy 48 epi ef)
       unless (eq nil result)
       append (list (list phi m) result))
  nil)

(rt:deftest oct.elliptic-pi.1q
    (loop for k from 0 below 100
       for phi = (random (/ +pi+ 2))
       for m = (random #q1)
       for epi = (elliptic-pi 0 phi m)
       for ef = (elliptic-f phi m)
       for result = (check-accuracy 53 epi ef)
       unless (eq nil result)
       append (list (list phi m) result))
  nil)

;; DLMF 19.6.3
;;
;; PI(n; pi/2 | 0) = pi/(2*sqrt(1-n))
(rt:deftest oct.elliptic-pi.19.6.3.d
    (loop for k from 0 below 100
       for n = (random 1d0)
       for epi = (elliptic-pi n (/ pi 2) 0)
       for true = (/ pi (* 2 (sqrt (- 1 n))))
       for result = (check-accuracy 47 epi true)
       unless (eq nil result)
       append (list (list (list k n) result)))
  nil)

(rt:deftest oct.elliptic-pi.19.6.3.q
    (loop for k from 0 below 100
       for n = (random #q1)
       for epi = (elliptic-pi n (/ (float-pi n) 2) 0)
       for true = (/ (float-pi n) (* 2 (sqrt (- 1 n))))
       for result = (check-accuracy 208 epi true)
       unless (eq nil result)
       append (list (list (list k n) result)))
  nil)

;; elliptic-pi(n, phi, 0) = 
;;   atan(sqrt(1-n)*tan(phi))/sqrt(1-n)   n < 1
;;   atanh(sqrt(n-1)*tan(phi))/sqrt(n-1)  n > 1
;;   tan(phi)                             n = 1
;;
;; These are easy to derive if you look at the integral:
;;
;; ellipti-pi(n, phi, 0) = integrate(1/(1-n*sin(t)^2), t, 0, phi)
;;
;; and this can be easily integrated to give the above expressions for
;; the different values of n.
(rt:deftest oct.elliptic-pi.n0.d
    ;; Tests for random values for phi in [0, pi/2] and n in [0, 1]
    (loop for k from 0 below 100
       for phi = (random (/ pi 2))
       for n = (random 1d0)
       for epi = (elliptic-pi n phi 0)
       for true = (/ (atan (* (tan phi) (sqrt (- 1 n))))
		     (sqrt (- 1 n)))
       for result = (check-accuracy 47.5 epi true)
       unless (eq nil result)
       append (list (list (list k n phi) result)))
  nil)

(rt:deftest oct.elliptic-pi.n1.d
    (loop for k from 0 below 100
       for phi = (random (/ pi 2))
       for epi = (elliptic-pi 1 phi 0)
       for true = (tan phi)
       for result = (check-accuracy 36 epi true)
       unless (eq nil result)
       append (list (list (list k phi) result)))
  nil)

(rt:deftest oct.elliptic-pi.n2.d
    (loop for k from 0 below 100
       for phi = (random (/ pi 2))
       for n = (+ 1d0 (random 100d0))
       for epi = (elliptic-pi n phi 0)
       for true = (/ (atanh (* (tan phi) (sqrt (- n 1))))
		     (sqrt (- n 1)))
       for result = (check-accuracy 47 epi true)
       ;; Not sure if this formula holds when atanh gives a complex
       ;; result.  Wolfram doesn't say
       when (and (not (complexp true)) result)
       append (list (list (list k n phi) result)))
  nil)

(rt:deftest oct.elliptic-pi.n0.q
    ;; Tests for random values for phi in [0, pi/2] and n in [0, 1]
    (loop for k from 0 below 100
       for phi = (random (/ +pi+ 2))
       for n = (random #q1)
       for epi = (elliptic-pi n phi 0)
       for true = (/ (atan (* (tan phi) (sqrt (- 1 n))))
		     (sqrt (- 1 n)))
       for result = (check-accuracy 208 epi true)
       unless (eq nil result)
       append (list (list (list k n phi) result)))
  nil)

(rt:deftest oct.elliptic-pi.n1.q
    (loop for k from 0 below 100
       for phi = (random (/ +pi+ 2))
       for epi = (elliptic-pi 1 phi 0)
       for true = (tan phi)
       for result = (check-accuracy 194 epi true)
       unless (eq nil result)
       append (list (list (list k phi) result)))
  nil)

(rt:deftest oct.elliptic-pi.n2.q
    (loop for k from 0 below 100
       for phi = (random (/ +pi+ 2))
       for n = (+ #q1 (random #q1))
       for epi = (elliptic-pi n phi 0)
       for true = (/ (atanh (* (tan phi) (sqrt (- n 1))))
		     (sqrt (- n 1)))
       for result = (check-accuracy 206 epi true)
       ;; Not sure if this formula holds when atanh gives a complex
       ;; result.  Wolfram doesn't say
       when (and (not (complexp true)) result)
       append (list (list (list k n phi) result)))
  nil)

;; Tests for theta functions.

(rt:deftest oct.theta3.1.d
    ;; A&S 16.38.5
    ;; sqrt(2*K/%pi) = theta3(0,q)
    (loop for k from 0 below 100
       for m = (random 1d0)
       for t3 = (elliptic-theta-3 0 (elliptic-nome m))
       for true = (sqrt (/ (* 2 (elliptic-k m)) (float-pi m)))
       for result = (check-accuracy 51 t3 true)
       when result
       append (list (list (list k m) result)))
  nil)

(rt:deftest oct.theta3.1.q
    ;; A&S 16.38.5
    ;; sqrt(2*K/%pi) = theta3(0,q)
    (loop for k from 0 below 100
       for m = (random #q1)
       for t3 = (elliptic-theta-3 0 (elliptic-nome m))
       for true = (sqrt (/ (* 2 (elliptic-k m)) (float-pi m)))
       for result = (check-accuracy 206 t3 true)
       when result
       append (list (list (list k m) result)))
  nil)

(rt:deftest oct.theta2.1.d
    ;; A&S 16.38.7
    ;; sqrt(2*sqrt(m)*K/%pi) = theta2(0,q)
    (loop for k from 0 below 100
       for m = (random 1d0)
       for t3 = (elliptic-theta-2 0 (elliptic-nome m))
       for true = (sqrt (/ (* 2 (sqrt m) (elliptic-k m)) (float-pi m)))
       for result = (check-accuracy 49 t3 true)
       when result
       append (list (list (list k m) result)))
  nil)

(rt:deftest oct.theta2.1.q
    ;; A&S 16.38.7
    ;; sqrt(2*sqrt(m)*K/%pi) = theta2(0,q)
    (loop for k from 0 below 100
       for m = (random #q1)
       for t3 = (elliptic-theta-2 0 (elliptic-nome m))
       for true = (sqrt (/ (* 2 (sqrt m) (elliptic-k m)) (float-pi m)))
       for result = (check-accuracy 205 t3 true)
       when result
       append (list (list (list k m) result)))
  nil)

(rt:deftest oct.theta4.1.d
    ;; A&S 16.38.8
    ;; sqrt(2*sqrt(1-m)*K/%pi) = theta2(0,q)
    (loop for k from 0 below 100
       for m = (random 1d0)
       for t3 = (elliptic-theta-4 0 (elliptic-nome m))
       for true = (sqrt (/ (* 2 (sqrt (- 1 m)) (elliptic-k m))
			   (float-pi m)))
       for result = (check-accuracy 49 t3 true)
       when result
       append (list (list (list k m) result)))
  nil)

(rt:deftest oct.theta4.1.q
    ;; A&S 16.38.8
    ;; sqrt(2*sqrt(1-m)*K/%pi) = theta2(0,q)
    (loop for k from 0 below 100
       for m = (random #q1)
       for t3 = (elliptic-theta-4 0 (elliptic-nome m))
       for true = (sqrt (/ (* 2 (sqrt (- 1 m)) (elliptic-k m))
			   (float-pi m)))
       for result = (check-accuracy 204 t3 true)
       when result
       append (list (list (list k m) result)))
  nil)

(rt:deftest gamma.1.d
    (let ((g (gamma 0.5d0))
	  (true (sqrt pi)))
      ;; This should give full accuracy but doesn't.
      (check-accuracy 51 g true))
  nil)

(rt:deftest gamma.1.q
    (let ((g (gamma #q0.5))
	  (true (sqrt +pi+)))
      ;; This should give full accuracy but doesn't.
      (check-accuracy 197 g true))
  nil)

(rt:deftest gamma.2.d
    (loop for k from 0 below 100
       for y = (+ 1 (random 100d0))
       for g = (abs (gamma (complex 0 y)))
       for true = (sqrt (/ pi y (sinh (* pi y))))
       for result = (check-accuracy 45 g true)
       when result
       append (list (list (list k y) result)))
  nil)

(rt:deftest gamma.2.q
    (loop for k from 0 below 100
       for y = (+ 1 (random #q100))
       for g = (abs (gamma (complex 0 y)))
       for true = (sqrt (/ +pi+ y (sinh (* +pi+ y))))
       for result = (check-accuracy 196 g true)
       when result
       append (list (list (list k y) result)))
  nil)

(rt:deftest gamma.3.d
    (loop for k from 0 below 100
       for y = (+ 1 (random 100d0))
       for g = (abs (gamma (complex 1/2 y)))
       for true = (sqrt (/ pi (cosh (* pi y))))
       for result = (check-accuracy 44 g true)
       when result
       append (list (list (list k y) result)))
  nil)

(rt:deftest gamma.3.q
    (loop for k from 0 below 100
       for y = (+ 1 (random #q100))
       for g = (abs (gamma (complex 1/2 y)))
       for true = (sqrt (/ +pi+ (cosh (* +pi+ y))))
       for result = (check-accuracy 196 g true)
       when result
       append (list (list (list k y) result)))
  nil)

;; gamma_incomplete(2,z) = integrate(t*exp(-t), t, z, inf)
;;                       = (z+1)*exp(-z)
(rt:deftest gamma-incomplete-tail.1.d
    (let* ((z 5d0)
	   (gi (incomplete-gamma-tail 2 z))
	   (true (* (+ z 1) (exp (- z)))))
      (check-accuracy 52 gi true))
  nil)

(rt:deftest gamma-incomplete-tail.2.d
    (let* ((z #c(1 5d0))
	   (gi (incomplete-gamma-tail 2 z))
	   (true (* (+ z 1) (exp (- z)))))
      (check-accuracy 50 gi true))
  nil)

(rt:deftest gamma-incomplete-tail.1.q
    (let* ((z 5d0)
	   (gi (incomplete-gamma-tail 2 z))
	   (true (* (+ z 1) (exp (- z)))))
      (check-accuracy 212 gi true))
  nil)

(rt:deftest gamma-incomplete-tail.1.q
    (let* ((z #q(1 5))
	   (gi (incomplete-gamma-tail 2 z))
	   (true (* (+ z 1) (exp (- z)))))
      (check-accuracy 206 gi true))
  nil)
