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


;; Lentz's algorithm for evaluating continued fractions.
;;
;; Let the continued fraction be:
;;
;;      a1    a2    a3
;; b0 + ----  ----  ----
;;      b1 +  b2 +  b3 +
;;
(defun lentz (bf af)
  (flet ((value-or-tiny (v)
	   (if (zerop v)
	       (etypecase v
		 ((or double-float cl:complex)
		  least-positive-normalized-double-float)
		 ((or qd-real qd-complex)
		  (make-qd least-positive-normalized-double-float)))
	       v)))
    (let* ((f (value-or-tiny (funcall bf 0)))
	   (c f)
	   (d 0)
	   (eps (epsilon f)))
      (loop
	 for j from 1
	 for an = (funcall af j)
	 for bn = (funcall bf j)
	 do (progn
	      (setf d (value-or-tiny (+ bn (* an d))))
	      (setf c (value-or-tiny (+ bn (/ an c))))
	      (setf d (/ d))
	      (let ((delta (* c d)))
		(setf f (* f delta))
		(when (<= (abs (- delta 1)) eps)
		  (return)))))
      f)))

;; Continued fraction for erf(b):
;;
;; z[n] = 1+2*n-2*z^2
;; a[n] = 4*n*z^2
;;
;; This works ok, but has problems for z > 3 where sometimes the
;; result is greater than 1.
#+nil
(defun erf (z)
  (let* ((z2 (* z z))
	 (twoz2 (* 2 z2)))
    (* (/ (* 2 z)
	  (sqrt (float-pi z)))
       (exp (- z2))
       (/ (lentz #'(lambda (n)
		     (- (+ 1 (* 2 n))
			twoz2))
		 #'(lambda (n)
		     (* 4 n z2)))))))

;; Tail of the incomplete gamma function:
;; integrate(x^(a-1)*exp(-x), x, z, inf)
;;
;; The continued fraction, valid for all z except the negative real
;; axis:
;;
;; b[n] = 1+2*n+z-a
;; a[n] = n*(a-n)
;;
;; See http://functions.wolfram.com/06.06.10.0003.01
(defun cf-incomplete-gamma-tail (a z)
  (/ (* (expt z a)
	(exp (- z)))
     (let ((z-a (- z a)))
       (lentz #'(lambda (n)
		  (+ n n 1 z-a))
	      #'(lambda (n)
		  (* n (- a n)))))))

;; Incomplete gamma function:
;; integrate(x^(a-1)*exp(-x), x, 0, z)
;;
;; The continued fraction, valid for all z except the negative real
;; axis:
;;
;; b[n] = n - 1 + z + a
;; a[n] = -z*(a + n)
;;
;; See http://functions.wolfram.com/06.06.10.0005.01.  We modified the
;; continued fraction slightly and discarded the first quotient from
;; the fraction.
(defun cf-incomplete-gamma (a z)
  (/ (* (expt z a)
	(exp (- z)))
     (let ((za1 (+ z a 1)))
       (- a (/ (* a z)
	       (lentz #'(lambda (n)
			  (+ n za1))
		      #'(lambda (n)
			  (- (* z (+ a n))))))))))

;; Series expansion for incomplete gamma.  Intended for |a|<1 and
;; |z|<1.  The series is
;;
;; g(a,z) = z^a * sum((-z)^k/k!/(a+k), k, 0, inf)
(defun s-incomplete-gamma (a z)
  (let ((-z (- z))
	(eps (epsilon z)))
    (loop for k from 0
       for term = 1 then (* term (/ -z k))
       for sum = (/ a) then (+ sum (/ term (+ a k)))
       when (< (abs term) (* (abs sum) eps))
       return (* sum (expt z a)))))

  

;; Tail of the incomplete gamma function.
(defun incomplete-gamma-tail (a z)
  "Tail of the incomplete gamma function defined by:

  integrate(t^(a-1)*exp(-t), t, z, inf)"
  (let* ((prec (float-contagion a z))
	 (a (apply-contagion a prec))
	 (z (apply-contagion z prec)))
    (if (and (zerop (imagpart a))
	     (zerop (imagpart z)))
	;; For real values, we split the result to compute either the
	;; tail directly from the continued fraction or from gamma(a)
	;; - incomplete-gamma.  The continued fraction doesn't
	;; converge on the negative real axis, so we can't use that
	;; there.  And accuracy appears to be better if z is "small".
	;; We take this to mean |z| < |a-1|.  Note that |a-1| is the
	;; peak of the integrand.
	(if (and (> (abs z) (abs (- a 1)))
		(not (minusp (realpart z))))
	    (cf-incomplete-gamma-tail a z)
	    (- (gamma a) (incomplete-gamma a z)))
	(cf-incomplete-gamma-tail a z))))

(defun incomplete-gamma (a z)
  "Incomplete gamma function defined by:

  integrate(t^(a-1)*exp(-t), t, 0, z)"
  (let* ((prec (float-contagion a z))
	 (a (apply-contagion a prec))
	 (z (apply-contagion z prec)))
    (if (and (< (abs a) 1) (< (abs z) 1))
	(s-incomplete-gamma a z)
	(if (and (realp a) (realp z))
	    (if (< z (- a 1))
		(cf-incomplete-gamma a z)
		(- (gamma a) (cf-incomplete-gamma-tail a z)))
	    ;; The continued fraction doesn't converge very fast if a
	    ;; and z are small.  In this case, use the series
	    ;; expansion instead, which converges quite rapidly.
	    (if (< (abs z) (abs a))
		(cf-incomplete-gamma a z)
		(- (gamma a) (cf-incomplete-gamma-tail a z)))))))

(defun erf (z)
  "Error function:

    erf(z) = 2/sqrt(%pi)*sum((-1)^k*z^(2*k+1)/k!/(2*k+1), k, 0, inf)

  For real z, this is equivalent to

    erf(z) = 2/sqrt(%pi)*integrate(exp(-t^2), t, 0, z) for real z."
  ;;
  ;; Erf is an odd function: erf(-z) = -erf(z)
  (if (minusp (realpart z))
      (- (erf (- z)))
      (/ (incomplete-gamma 1/2 (* z z))
	 (sqrt (float-pi z)))))

(defun erfc (z)
  "Complementary error function:

    erfc(z) = 1 - erf(z)"
  ;; Compute erfc(z) via 1 - erf(z) is not very accurate if erf(z) is
  ;; near 1.  Wolfram says
  ;;
  ;; erfc(z) = 1 - sqrt(z^2)/z * (1 - 1/sqrt(pi)*gamma_incomplete_tail(1/2, z^2))
  ;;
  ;; For real(z) > 0, sqrt(z^2)/z is 1 so
  ;;
  ;; erfc(z) = 1 - (1 - 1/sqrt(pi)*gamma_incomplete_tail(1/2,z^2))
  ;;         = 1/sqrt(pi)*gamma_incomplete_tail(1/2,z^2)
  ;;
  ;; For real(z) < 0, sqrt(z^2)/z is -1 so
  ;;
  ;; erfc(z) = 1 + (1 - 1/sqrt(pi)*gamma_incomplete_tail(1/2,z^2))
  ;;         = 1 + 1/sqrt(pi)*gamma_incomplete(1/2,z^2)
  (if (>= (realpart z) 0)
      (/ (incomplete-gamma-tail 1/2 (* z z))
	 (sqrt (float-pi z)))
      (+ 1
	 (/ (incomplete-gamma 1/2 (* z z))
	    (sqrt (float-pi z))))))

(defun exp-integral-e (v z)
  "Exponential integral E:

   E(v,z) = integrate(exp(-t)/t^v, t, 1, inf)"
  ;; Wolfram gives E(v,z) = z^(v-1)*gamma_incomplete_tail(1-v,z)
  (* (expt z (- v 1))
     (incomplete-gamma-tail (- 1 v) z)))

;; Series for Fresnel S
;;
;;   S(z) = z^3*sum((%pi/2)^(2*k+1)(-z^4)^k/(2*k+1)!/(4*k+3), k, 0, inf)
;;
;; Compute as
;;
;;   S(z) = z^3*sum(a(k)/(4*k+3), k, 0, inf)
;;
;; where
;;
;;   a(k+1) = -a(k) * (%pi/2)^2 * z^4 / (2*k+2) / (2*k+3)
;;
;;   a(0) = %pi/2.
(defun fresnel-s-series (z)
  (let* ((pi/2 (* 1/2 (float-pi z)))
	 (factor (- (* (expt z 4) pi/2 pi/2)))
	 (eps (epsilon z))
	 (sum 0)
	 (term pi/2))
    (loop for k2 from 0 by 2
       until (< (abs term) (* eps sum))
       do
       (incf sum (/ term (+ 3 k2 k2)))
       (setf term (/ (* term factor)
		     (* (+ k2 2)
			(+ k2 3)))))
    (* sum (expt z 3))))
    
(defun fresnel-s (z)
  "Fresnel S:

   S(z) = integrate(sin(%pi*t^2/2), t, 0, z) "
  (let ((prec (float-contagion z))
	(sqrt-pi (sqrt (float-pi z))))
    (flet ((fs (z)
	     ;; Wolfram gives
	     ;;
	     ;;  S(z) = (1+%i)/4*(erf(c*z) - %i*erf(conjugate(c)*z))
	     ;;
	     ;; where c = sqrt(%pi)/2*(1+%i).
	     ;;
	     ;; But for large z, we should use erfc.  Then
	     ;;  S(z) = 1/2 - (1+%i)/4*(erfc(c*z) - %i*erfc(conjugate(c)*z))
	     (if (and t (> (abs z) 2))
		 (- 1/2
		    (* #c(1/4 1/4)
		       (- (erfc (* #c(1/2 1/2) sqrt-pi z))
			  (* #c(0 1)
			     (erfc (* #c(1/2 -1/2) sqrt-pi z))))))
		 (* #c(1/4 1/4)
		    (- (erf (* #c(1/2 1/2) sqrt-pi z))
		       (* #c(0 1)
			  (erf (* #c(1/2 -1/2) sqrt-pi z)))))))
          (rfs (z)
            ;; When z is real, recall that erf(conjugate(z)) =
            ;; conjugate(erf(z)).  Then
            ;;
            ;;  S(z) = 1/2*(realpart(erf(c*z)) - imagpart(erf(c*z)))
            ;;
            ;; But for large z, we should use erfc.  Then
            ;;
            ;;  S(z) = 1/2 - 1/2*(realpart(erfc(c*z)) - imagpart(erf(c*z)))
            (if (> (abs z) 2)
                (let ((s (erfc (* #c(1/2 1/2) sqrt-pi z))))
                  (- 1/2
                     (* 1/2 (- (realpart s) (imagpart s)))))
                (let ((s (erf (* #c(1/2 1/2) sqrt-pi z))))
                  (* 1/2 (- (realpart s) (imagpart s)))))))
      ;; For small z, the erf terms above suffer from subtractive
      ;; cancellation.  So use the series in this case.  Some simple
      ;; tests were done to determine that for double-floats we want
      ;; to use the series for z < 1 to give max accuracy.  For
      ;; qd-real, the above formula is good enough for z > 1d-5.
      (if (< (abs z) (ecase prec
		       (single-float 1.5f0)
		       (double-float 1d0)
		       (qd-real #q1)))
	  (fresnel-s-series z)
	  (if (realp z)
	      ;; FresnelS is real for a real argument. And it is odd.
	      (if (minusp z)
		  (- (rfs (- z)))
		  (rfs z))
	      (fs z))))))

(defun fresnel-c (z)
  "Fresnel C:

   C(z) = integrate(cos(%pi*t^2/2), t, 0, z) "
  (let ((sqrt-pi (sqrt (float-pi z))))
    (flet ((fs (z)
	     ;; Wolfram gives
	     ;;
	     ;;  C(z) = (1-%i)/4*(erf(c*z) + %i*erf(conjugate(c)*z))
	     ;;
	     ;; where c = sqrt(%pi)/2*(1+%i).
	     (* #c(1/4 -1/4)
		(+ (erf (* #c(1/2 1/2) sqrt-pi z))
		   (* #c(0 1)
		      (erf (* #c(1/2 -1/2) sqrt-pi z)))))))
      (if (realp z)
	  ;; FresnelS is real for a real argument. And it is odd.
	  (if (minusp z)
	      (- (realpart (fs (- z))))
	      (realpart (fs z)))
	  (fs z)))))

(defun sin-integral (z)
  "Sin integral:

  Si(z) = integrate(sin(t)/t, t, 0, z)"
  ;; Wolfram has
  ;;
  ;; Si(z) = %i/2*(gamma_inc_tail(0, -%i*z) - gamma_inc_tail(0, %i*z) + log(-%i*z)-log(%i*z))
  ;;
  (flet ((si (z)
	   (* #c(0 1/2)
	      (let ((iz (* #c(0 1) z))
		    (-iz (* #c(0 -1) z)))
		(+ (- (incomplete-gamma-tail 0 -iz)
		      (incomplete-gamma-tail 0 iz))
		   (- (log -iz)
		      (log iz)))))))
    (if (realp z)
	;; Si is odd and real for real z.  In this case, we have
	;;
	;; Si(x) = %i/2*(gamma_inc_tail(0, -%i*x) - gamma_inc_tail(0, %i*x) - %i*%pi)
	;;       = %pi/2 + %i/2*(gamma_inc_tail(0, -%i*x) - gamma_inc_tail(0, %i*x))
	;; But gamma_inc_tail(0, conjugate(z)) = conjugate(gamma_inc_tail(0, z)), so
	;;
	;; Si(x) = %pi/2 + imagpart(gamma_inc_tail(0, %i*x))
	(cond ((< z 0)
	       (- (sin-integral (- z))))
	      ((= z 0)
	       (* 0 z))
	      (t
	       (+ (* 1/2 (float-pi z))
		  (imagpart (incomplete-gamma-tail 0 (complex 0 z))))))
	(si z))))

(defun cos-integral (z)
  "Cos integral:

    Ci(z) = integrate((cos(t) - 1)/t, t, 0, z) + log(z) + gamma

  where gamma is Euler-Mascheroni constant"
  ;; Wolfram has
  ;;
  ;; Ci(z) = log(z) - 1/2*(gamma_inc_tail(0, -%i*z) + gamma_inc_tail(0, %i*z) + log(-%i*z)+log(%i*z))
  ;;
  (flet ((ci (z)
	   (- (log z)
	      (* 1/2
		 (let ((iz (* #c(0 1) z))
		       (-iz (* #c(0 -1) z)))
		   (+ (+ (incomplete-gamma-tail 0 -iz)
			 (incomplete-gamma-tail 0 iz))
		      (+ (log -iz)
			 (log iz))))))))
    (if (and (realp z) (plusp z))
	(realpart (ci z))
	(ci z))))
  