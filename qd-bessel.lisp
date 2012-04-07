;;;; -*- Mode: lisp -*-
;;;;
;;;; Copyright (c) 2011 Raymond Toy
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

;;; References:
;;;
;;; [1] Borwein, Borwein, Crandall, "Effective Laguerre Asymptotics",
;;; http://people.reed.edu/~crandall/papers/Laguerre-f.pdf
;;;
;;; [2] Borwein, Borwein, Chan, "The Evaluation of Bessel Functions
;;; via Exp-Arc Integrals", http://web.cs.dal.ca/~jborwein/bessel.pdf
;;;

(defvar *debug-exparc* nil)

;; B[k](p) = 1/2^(k+3/2)*integrate(exp(-p*u)*u^(k-1/2),u,0,1)
;;         = 1/2^(k+3/2)/p^(k+1/2)*integrate(t^(k-1/2)*exp(-t),t,0,p)
;;         = 1/2^(k+3/2)/p^(k+1/2) * g(k+1/2, p)
;;
;; where g(a,z) is the lower incomplete gamma function.
;;
;; There is the continued fraction expansion for g(a,z) (see
;; cf-incomplete-gamma in qd-gamma.lisp):
;;
;;  g(a,z) = z^a*exp(-z)/ CF
;;
;; So
;;
;;  B[k](p) = 1/2^(k+3/2)/p^(k+1/2)*p^(k+1/2)*exp(-p)/CF
;;          = exp(-p)/2^(k+3/2)/CF
;;
(defun bk (k p)
  (/ (exp (- p))
     (* (sqrt (float 2 (realpart p))) (ash 1 (+ k 1)))
     (let ((a (float (+ k 1/2) (realpart p))))
       (lentz #'(lambda (n)
		  (+ n a))
	      #'(lambda (n)
		  (if (evenp n)
		      (* (ash n -1) p)
		      (- (* (+ a (ash n -1)) p))))))))

;; exp-arc I function, as given in the Laguerre paper
;;
;; I(p, q) = 4*exp(p) * sum(g[k](-2*%i*q)/(2*k)!*B[k](p), k, 0, inf)
;;
;; where g[k](p) = product(p^2+(2*j-1)^2, j, 1, k) and B[k](p) as above.
;;
;; For computation, note that g[k](p) = g[k-1](p) * (p^2 + (2*k-1)^2)
;; and (2*k)! = (2*k-2)! * (2*k-1) * (2*k).  Then, let
;;
;;  R[k](p) = g[k](p)/(2*k)!
;;
;; Then
;;
;;  R[k](p) = g[k](p)/(2*k)!
;;          = g[k-1](p)/(2*k-2)! * (p^2 + (2*k-1)^2)/((2*k-1)*(2*k)
;;          = R[k-1](p) * (p^2 + (2*k-1)^2)/((2*k-1)*(2*k)
;;
;; In the exp-arc paper, the function is defined (equivalently) as
;; 
;; I(p, q) = 2*%i*exp(p)/q * sum(r[2*k+1](-2*%i*q)/(2*k)!*B[k](p), k, 0, inf)
;;
;; where r[2*k+1](p) = p*product(p^2 + (2*j-1)^2, j, 1, k)
;;
;; Let's note some properties of I(p, q).
;;
;; I(-%i*z, v) = 2*%i*exp(-%i*z)/q * sum(r[2*k+1](-2*%i*v)/(2*k)!*B[k](-%i*z))
;;
;; Note thate B[k](-%i*z) = 1/2^(k+3/2)*integrate(exp(%i*z*u)*u^(k-1/2),u,0,1)
;;                        = conj(B[k](%i*z).
;;
;; Hence I(-%i*z, v) = conj(I(%i*z, v)) when both z and v are real.
(defun exp-arc-i (p q)
  (let* ((sqrt2 (sqrt (float 2 (realpart p))))
	 (exp/p/sqrt2 (/ (exp (- p)) p sqrt2))
	 (v (* #c(0 -2) q))
	 (v2 (expt v 2))
	 (eps (epsilon (realpart p))))
    (when *debug-exparc*
      (format t "sqrt2 = ~S~%" sqrt2)
      (format t "exp/p/sqrt2 = ~S~%" exp/p/sqrt2))
    (do* ((k 0 (1+ k))
	  (bk (/ (incomplete-gamma 1/2 p)
		 2 sqrt2 (sqrt p))
	      (- (/ (* bk (- k 1/2)) 2 p)
		 (/ exp/p/sqrt2 (ash 1 (+ k 1)))))
	  ;; ratio[k] = r[2*k+1](v)/(2*k)!.
	  ;; r[1] = v and r[2*k+1](v) = r[2*k-1](v)*(v^2 + (2*k-1)^2)
	  ;; ratio[0] = v
	  ;; and ratio[k] = r[2*k-1](v)*(v^2+(2*k-1)^2) / ((2*k-2)! * (2*k-1) * 2*k)
	  ;;              = ratio[k]*(v^2+(2*k-1)^2)/((2*k-1) * 2 * k)
	  (ratio v
		 (* ratio (/ (+ v2 (expt (1- (* 2 k)) 2))
			     (* 2 k (1- (* 2 k))))))
	  (term (* ratio bk)
		(* ratio bk))
	  (sum term (+ sum term)))
	 ((< (abs term) (* (abs sum) eps))
	  (* sum #c(0 2) (/ (exp p) q)))
      (when *debug-exparc*
	(format t "k      = ~D~%" k)
	(format t " bk    = ~S~%" bk)
	(format t " ratio = ~S~%" ratio)
	(format t " term  = ~S~%" term)
	(format t " sum   - ~S~%" sum)))))

(defun exp-arc-i-2 (p q)
  (let* ((sqrt2 (sqrt (float 2 (realpart p))))
	 (exp/p/sqrt2 (/ (exp (- p)) p sqrt2))
	 (v (* #c(0 -2) q))
	 (v2 (expt v 2))
	 (eps (epsilon (realpart p))))
    (when *debug-exparc*
      (format t "sqrt2 = ~S~%" sqrt2)
      (format t "exp/p/sqrt2 = ~S~%" exp/p/sqrt2))
    (do* ((k 0 (1+ k))
	  (bk (bk 0 p)
	      (bk k p))
	  (ratio v
		 (* ratio (/ (+ v2 (expt (1- (* 2 k)) 2))
			     (* 2 k (1- (* 2 k))))))
	  (term (* ratio bk)
		(* ratio bk))
	  (sum term (+ sum term)))
	 ((< (abs term) (* (abs sum) eps))
	  (* sum #c(0 2) (/ (exp p) q)))
      (when *debug-exparc*
	(format t "k      = ~D~%" k)
	(format t " bk    = ~S~%" bk)
	(format t " ratio = ~S~%" ratio)
	(format t " term  = ~S~%" term)
	(format t " sum   - ~S~%" sum)))))


;; This currently only works for v an integer.
;;
(defun bessel-j-exp-arc (v z)
  (let* ((iz (* #c(0 1) z))
	 (i+ (exp-arc-i-2 iz v))
	 (i- (exp-arc-i-2 (- iz ) v)))
    (/ (+ (* (cis (* v (float-pi i+) -1/2))
	     i+)
	  (* (cis (* v (float-pi i+) 1/2))
	     i-))
       (float-pi i+)
       2)))

(defun paris-series (v z n)
  (labels ((pochhammer (a k)
	     (/ (gamma (+ a k))
		(gamma a)))
	   (a (v k)
	     (* (/ (pochhammer (+ 1/2 v) k)
		   (gamma (float (1+ k) z)))
		(pochhammer (- 1/2 v) k))))
    (* (loop for k from 0 below n
	     sum (* (/ (a v k)
		       (expt (* 2 z) k))
		    (/ (cf-incomplete-gamma (+ k v 1/2) (* 2 z))
		       (gamma (+ k v 1/2)))))
       (/ (exp z)
	  (sqrt (* 2 (float-pi z) z))))))

