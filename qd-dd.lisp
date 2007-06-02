;;; -*- Mode: lisp -*-

(in-package #:qdi)

;;; double-double float routines needed for quad-double.
;;;
;;; Not needed for CMUCL.
;;;
;;; These routines were taken directly from CMUCL.

(declaim (inline quick-two-sum))
(defun quick-two-sum (a b)
  "Computes fl(a+b) and err(a+b), assuming |a| >= |b|"
  (declare (double-float a b))
  (let* ((s (+ a b))
	 (e (- b (- s a))))
    (values s e)))

(declaim (inline two-sum))
(defun two-sum (a b)
  "Computes fl(a+b) and err(a+b)"
  (declare (double-float a b))
  (let* ((s (+ a b))
	 (v (- s a))
	 (e (+ (- a (- s v))
	       (- b v))))
    (values s e)))

(declaim (inline two-prod))
(declaim (inline split))
;; This algorithm is the version given by Yozo Hida.  It has problems
;; with overflow because we multiply by 1+2^27.
;;
;; But be very careful about replacing this with a new algorithm.  The
;; values computed here are very important to get the rounding right.
;; If you change this, the rounding may be different, which will
;; affect other parts of the algorithm.
;;
;; I (rtoy) tried a different algorithm that split the number in two
;; as described, but without overflow.  However, that caused
;; -9.4294948327242751340284975915175w0/1w14 to return a value that
;; wasn't really close to -9.4294948327242751340284975915175w-14.
;;
;; This also means we can't print numbers like 1w308 with the current
;; printing algorithm, or even divide 1w308 by 10.
#+nil
(defun split (a)
  "Split the double-float number a into a-hi and a-lo such that a =
  a-hi + a-lo and a-hi contains the upper 26 significant bits of a and
  a-lo contains the lower 26 bits."
  (declare (double-float a))
  (let* ((tmp (* a (+ 1 (expt 2 27))))
	 (a-hi (- tmp (- tmp a)))
	 (a-lo (- a a-hi)))
    (values a-hi a-lo)))


(defun split (a)
  "Split the double-float number a into a-hi and a-lo such that a =
  a-hi + a-lo and a-hi contains the upper 26 significant bits of a and
  a-lo contains the lower 26 bits."
  (declare (double-float a)
	   (optimize (speed 3)))
  ;; This splits the number a into 2 halves of 26 bits each, but the
  ;; halves are, I think, supposed to be properly rounded in an IEEE
  ;; fashion.
  ;;
  ;; For numbers that are very large, we use a different algorithm.
  ;; For smaller numbers, we can use the original algorithm of Yozo
  ;; Hida.
  (if (> (abs a) (scale-float 1d0 (- 1023 27)))
      ;; I've tested this algorithm against Yozo's method for 1
      ;; billion randomly generated double-floats between 2^(-995) and
      ;; 2^996, and identical results are obtained.  For numbers that
      ;; are very small, this algorithm produces different numbers
      ;; because of underflow.  For very large numbers, we, of course
      ;; produce different results because Yozo's method causes
      ;; overflow.
      (let* ((tmp (* a (+ 1 (scale-float 1d0 -27))))
	     (as (* a (scale-float 1d0 -27)))
	     (a-hi (* (- tmp (- tmp as)) (expt 2 27)))
	     (a-lo (- a a-hi)))
	(values a-hi a-lo))
      ;; Yozo's algorithm.
      (let* ((tmp (* a (+ 1 (expt 2 27))))
	     (a-hi (- tmp (- tmp a)))
	     (a-lo (- a a-hi)))
	(values a-hi a-lo))))


(defun two-prod (a b)
  "Compute fl(a*b) and err(a*b)"
  (declare (double-float a b))
  (let ((p (* a b)))
    (multiple-value-bind (a-hi a-lo)
	(split a)
      (multiple-value-bind (b-hi b-lo)
	  (split b)
	(let ((e (+ (+ (- (* a-hi b-hi) p)
		       (* a-hi b-lo)
		       (* a-lo b-hi))
		    (* a-lo b-lo))))
	  (values p e))))))

(declaim (inline two-sqr))
(defun two-sqr (a)
  "Compute fl(a*a) and err(a*b).  This is a more efficient
  implementation of two-prod"
  (declare (double-float a))
  (let ((q (* a a)))
    (multiple-value-bind (a-hi a-lo)
	(split a)
      (values q (+ (+ (- (* a-hi a-hi) q)
		      (* 2 a-hi a-lo))
		   (* a-lo a-lo))))))
