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

(eval-when (:compile-toplevel :load-toplevel :execute)
  (setf *readtable* *oct-readtable*))

;; For log-gamma we use the asymptotic formula
;;
;; log(gamma(z)) ~ (z - 1/2)*log(z) + log(2*%pi)/2
;;                   + sum(bern(2*k)/(2*k)/(2*k-1)/z^(2k-1), k, 1, inf)
;;
;;               = (z - 1/2)*log(z) + log(2*%pi)/2
;;                  + 1/12/z*(1 - 1/30/z^2 + 1/105/z^4 + 1/140/z^6 + ...
;;                              + 174611/10450/z^18 + ...)
;;
;; For double-floats, let's stop the series at the power z^18.  The
;; next term is 77683/483/z^20.  This means that for |z| > 8.09438,
;; the series has double-float precision.
;;
;; For quad-doubles, let's stop the series at the power z^62.  The
;; next term is about 6.364d37/z^64.  So for |z| > 38.71, the series
;; has quad-double precision.
(defparameter *log-gamma-asymp-coef*
  #(-1/30 1/105 -1/140 1/99 -691/30030 1/13 -3617/10200 43867/20349 
    -174611/10450 77683/483 -236364091/125580 657931/25 -3392780147/7830 
    1723168255201/207669 -7709321041217/42160 151628697551/33 
    -26315271553053477373/201514950 154210205991661/37 
    -261082718496449122051/1758900 1520097643918070802691/259161 
    -2530297234481911294093/9890 25932657025822267968607/2115 
    -5609403368997817686249127547/8725080 19802288209643185928499101/539 
    -61628132164268458257532691681/27030 29149963634884862421418123812691/190323 
    -354198989901889536240773677094747/31900 
    2913228046513104891794716413587449/3363 
    -1215233140483755572040304994079820246041491/16752085350 
    396793078518930920708162576045270521/61 
    -106783830147866529886385444979142647942017/171360 
    133872729284212332186510857141084758385627191/2103465
    ))

#+nil
(defun log-gamma-asymp-series (z nterms)
  ;; Sum the asymptotic formula for n terms
  ;;
  ;; 1 + sum(c[k]/z^(2*k+2), k, 0, nterms)
  (let ((z2 (* z z))
	(sum 1)
	(term 1))
    (dotimes (k nterms)
      (setf term (* term z2))
      (incf sum (/ (aref *log-gamma-asymp-coef* k) term)))
    sum))

(defun log-gamma-asymp-series (z nterms)
  (loop with y = (* z z)
     for k from 1 to nterms
     for x = 0 then
       (setf x (/ (+ x (aref *log-gamma-asymp-coef* (- nterms k)))
		  y))
     finally (return (+ 1 x))))
       

(defun log-gamma-asymp-principal (z nterms log2pi/2)
  (+ (- (* (- z 1/2)
	   (log z))
	z)
     log2pi/2))

(defun log-gamma-asymp (z nterms log2pi/2)
  (+ (log-gamma-asymp-principal z nterms log2pi/2)
     (* 1/12 (/ (log-gamma-asymp-series z nterms) z))))

(defun log2pi/2 (precision)
  (ecase precision
    (single-float
     (coerce (/ (log (* 2 pi)) 2) 'single-float))
    (double-float
     (coerce (/ (log (* 2 pi)) 2) 'double-float))
    (qd-real
     (/ (log +2pi+) 2))))

(defun log-gamma-aux (z limit nterms)
  (let ((precision (float-contagion z)))
    (cond ((minusp (realpart z))
	   ;; Use reflection formula if realpart(z) < 0
	   ;;   log(gamma(-z)) = log(pi)-log(-z)-log(sin(pi*z))-log(gamma(z))
	   ;; Or
	   ;;   log(gamma(z)) = log(pi)-log(-z)-log(sin(pi*z))-log(gamma(-z))
	   (- (apply-contagion (log pi) precision)
	      (log (- z))
	      (apply-contagion (log (sin (* pi z))) precision)
	      (log-gamma (- z))))
	  (t
	   (let ((absz (abs z)))
	     (cond ((>= absz limit)
		    ;; Can use the asymptotic formula directly with 9 terms
		    (log-gamma-asymp z nterms (log2pi/2 precision)))
		   (t
		    ;; |z| is too small.  Use the formula
		    ;; log(gamma(z)) = log(gamma(z+1)) - log(z)
		    (- (log-gamma (+ z 1))
		       (log z)))))))))

(defmethod log-gamma ((z cl:number))
  (log-gamma-aux z 9 9))

(defmethod log-gamma ((z qd-real))
  (log-gamma-aux z 26 26))

(defmethod log-gamma ((z qd-complex))
  (log-gamma-aux z 26 26))

(defun gamma-aux (z limit nterms)
  (let ((precision (float-contagion z)))
    (cond ((minusp (realpart z))
	   ;; Use reflection formula if realpart(z) < 0:
	   ;;  gamma(-z) = -pi*csc(pi*z)/gamma(z+1)
	   ;; or
	   ;;  gamma(z) = pi*csc(pi*z)/gamma(1-z)
	   (/ (float-pi z)
	      (sin (* (float-pi z) z))
	      (gamma-aux (- 1 z) limit nterms)))
	  (t
	   (let ((absz (abs z)))
	     (cond ((>= absz limit)
		    ;; Use log gamma directly:
		    ;;  log(gamma(z)) = principal part + 1/12/z*(series part)
		    ;; so
		    ;;  gamma(z) = exp(principal part)*exp(1/12/z*series)
		    (exp (log-gamma z))
		    #+nil
		    (* (exp (log-gamma-asymp-principal z nterms
						       (log2pi/2 precision)))
		       (exp (* 1/12 (/ (log-gamma-asymp-series z nterms) z)))))
		   (t
		    ;; 1 <= |z| <= limit
		    ;; gamma(z) = gamma(z+1)/z
		    (/ (gamma-aux (+ 1 z) limit nterms) z))))))))
		 
(defmethod gamma ((z cl:number))
  (gamma-aux z 9 9))

(defmethod gamma ((z qd-real))
  (gamma-aux z 39 32))

(defmethod gamma ((z qd-complex))
  (gamma-aux z 39 32))

