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

(declaim (inline descending-transform ascending-transform))

;; Determine which of x and y has the higher precision and return the
;; value of the higher precision number.  If both x and y are
;; rationals, just return 1f0, for a single-float value.
(defun float-contagion-2 (x y)
  (etypecase x
    (cl:rational
     (etypecase y
       (cl:rational
	1f0)
       (cl:float
	y)
       (qd-real
	y)))
    (single-float
     (etypecase y
       ((or cl:rational single-float)
	x)
       ((or double-float qd-real)
	y)))
    (double-float
     (etypecase y
       ((or cl:rational single-float double-float)
	x)
       (qd-real
	y)))
    (qd-real
     x)))
    
;; Return a floating point (or complex) type of the highest precision
;; among all of the given arguments.
(defun float-contagion (&rest args)
  ;; It would be easy if we could just add the args together and let
  ;; normal contagion do the work, but we could easily introduce
  ;; overflows or other errors that way.  So look at each argument and
  ;; determine the precision and choose the highest precision.  
  (let ((complexp (some #'complexp args))
	(max-type
	 (etypecase (reduce #'float-contagion-2 (mapcar #'realpart (if (cdr args)
								       args
								       (list (car args) 0))))
	   (single-float 'single-float)
	   (double-float 'double-float)
	   (qd-real 'qd-real))))
    max-type))

(defun apply-contagion (number precision)
  (etypecase number
    ((or cl:real qd-real)
     (coerce number precision))
    ((or cl:complex qd-complex)
     (complex (coerce (realpart number) precision)
	      (coerce (imagpart number) precision)))))

;;; Jacobian elliptic functions

(defun ascending-transform (u m)
  ;; A&S 16.14.1
  ;;
  ;; Take care in computing this transform.  For the case where
  ;; m is complex, we should compute sqrt(mu1) first as
  ;; (1-sqrt(m))/(1+sqrt(m)), and then square this to get mu1.
  ;; If not, we may choose the wrong branch when computing
  ;; sqrt(mu1).
  ;;
  ;; mu = 4*sqrt(m)/(1+sqrt(m))^2
  ;; sqrt(mu1) = (1-sqrt(m))/(1+sqrt(m))
  ;; v = u/(1+sqrt(mu1))
  ;;
  ;; Return v, mu, sqrt(mu1)
  (let* ((root-m (sqrt m))
	 (mu (/ (* 4 root-m)
		(expt (1+ root-m) 2)))
	 (root-mu1 (/ (- 1 root-m) (+ 1 root-m)))
	 (v (/ u (1+ root-mu1))))
    (values v mu root-mu1)))

(defun descending-transform (u m)
  ;; A&S 16.12.1
  ;;
  ;; Note: Don't calculate mu first, as given in 16.12.1.  We
  ;; should calculate sqrt(mu) = (1-sqrt(m1)/(1+sqrt(m1)), and
  ;; then compute mu = sqrt(mu)^2.  If we calculate mu first,
  ;; sqrt(mu) loses information when m or m1 is complex.
  ;;
  ;; sqrt(mu) = (1-sqrt(m1))/(1+sqrt(m1))
  ;; v = u/(1+sqrt(mu))
  ;;
  ;; where m1 = 1-m
  ;;
  ;; Return v, mu, sqrt(mu)
  (let* ((root-m1 (sqrt (- 1 m)))
	 (root-mu (/ (- 1 root-m1) (+ 1 root-m1)))
	 (mu (* root-mu root-mu))
	 (v (/ u (1+ root-mu))))
    (values v mu root-mu)))


;; Could use the descending transform, but some of my tests show
;; that it has problems with roundoff errors.

;; WARNING: This doesn't work very well for u > 1000 or so.  For
;; example (elliptic-dn-ascending 1000b0 .5b0) -> 3.228b324, but dn <= 1.
#+nil
(defun elliptic-dn-ascending (u m)
  (cond ((zerop m)
	 ;; A&S 16.6.3
	 1.0)
	((< (abs (- 1 m)) (* 4 (epsilon u)))
	 ;; A&S 16.6.3
	 (/ (cosh u)))
	(t
	 (multiple-value-bind (v mu root-mu1)
	     (ascending-transform u m)
	   ;; A&S 16.14.4
	   (let* ((new-dn (elliptic-dn-ascending v mu)))
	     (* (/ (- 1 root-mu1) mu)
		(/ (+ root-mu1 (* new-dn new-dn))
		   new-dn)))))))

;; Don't use the descending version because it requires cn, dn, and
;; sn.
;;
;; WARNING: This doesn't work very well for large u.
;; (elliptic-cn-ascending 1000b0 .5b0) -> 4.565b324.  But |cn| <= 1.
#+nil
(defun elliptic-cn-ascending (u m)
  (cond ((zerop m)
	 ;; A&S 16.6.2
	 (cos u))
	((< (abs (- 1 m)) (* 4 (epsilon u)))
	 ;; A&S 16.6.2
	 (/ (cl:cosh u)))
	(t
	 (multiple-value-bind (v mu root-mu1)
	     (ascending-transform u m)
	   ;; A&S 16.14.3
	   (let* ((new-dn (elliptic-dn-ascending v mu)))
	     (* (/ (+ 1 root-mu1) mu)
		(/ (- (* new-dn new-dn) root-mu1)
		   new-dn)))))))

;;
;; This appears to work quite well for both real and complex values
;; of u.
(defun elliptic-sn-descending (u m)
  (cond ((= m 1)
	 ;; A&S 16.6.1
	 (tanh u))
	((< (abs m) (epsilon u))
	 ;; A&S 16.6.1
	 (sin u))
	(t
	 ;; A&S 16.12.2
	 ;;
	 ;; sn(u|m) = (1 + sqrt(mu))*sn(v|u)/(1 + sqrt(mu)*sn(v|mu)^2)
	 (multiple-value-bind (v mu root-mu)
	     (descending-transform u m)
	   (let* ((new-sn (elliptic-sn-descending v mu)))
	     (/ (* (1+ root-mu) new-sn)
		(1+ (* root-mu new-sn new-sn))))))))

;; We don't use the ascending transform here because it requires
;; evaluating sn, cn, and dn.  The ascending transform only needs
;; sn.
#+nil
(defun elliptic-sn-ascending (u m)
  (if (< (abs (- 1 m)) (* 4 flonum-epsilon))
      ;; A&S 16.6.1
      (tanh u)
      (multiple-value-bind (v mu root-mu1)
	  (ascending-transform u m)
	;; A&S 16.14.2
	(let* ((new-cn (elliptic-cn-ascending v mu))
	       (new-dn (elliptic-dn-ascending v mu))
	       (new-sn (elliptic-sn-ascending v mu)))
	  (/ (* (+ 1 root-mu1) new-sn new-cn)
	     new-dn)))))

(defun jacobi-sn (u m)
  "Compute Jacobian sn for argument u and parameter m"
  (let ((s (elliptic-sn-descending u m)))
    (if (and (realp u) (realp m))
	(realpart s)
	s)))

(defun jacobi-dn (u m)
  "Compute Jacobi dn for argument u and parameter m"
  ;; Use the Gauss transformation from
  ;; http://functions.wolfram.com/09.29.16.0013.01:
  ;;
  ;;
  ;; dn((1+sqrt(m))*z, 4*sqrt(m)/(1+sqrt(m))^2)
  ;;   =  (1-sqrt(m)*sn(z, m)^2)/(1+sqrt(m)*sn(z,m)^2)
  ;;
  ;; So
  ;;
  ;; dn(y, mu) = (1-sqrt(m)*sn(z, m)^2)/(1+sqrt(m)*sn(z,m)^2)
  ;;
  ;; where z = y/(1+sqrt(m)) and mu=4*sqrt(m)/(1+sqrt(m))^2.
  ;;
  ;; Solve for m, and we get
  ;;
  ;; sqrt(m) = -(mu+2*sqrt(1-mu)-2)/mu or (-mu+2*sqrt(1-mu)+2)/mu.
  ;;
  ;; I don't think it matters which sqrt we use, so I (rtoy)
  ;; arbitrarily choose the first one above.
  ;;
  ;; Note that (1-sqrt(1-mu))/(1+sqrt(1-mu)) is the same as
  ;; -(mu+2*sqrt(1-mu)-2)/mu.  Also, the former is more
  ;; accurate for small mu.
  (let* ((root (let ((root-1-m (sqrt (- 1 m))))
		 (/ (- 1 root-1-m)
		    (+ 1 root-1-m))))
	 (z (/ u (+ 1 root)))
	 (s (elliptic-sn-descending z (* root root)))
	 (p (* root s s )))
    (/ (- 1 p)
       (+ 1 p))))

(defun jacobi-cn (u m)
  "Compute Jacobi cn for argument u and parameter m"
  ;; Use the ascending Landen transformation, A&S 16.14.3.
  ;;
  ;; cn(u,m) = (1+sqrt(mu1))/mu * (dn(v,mu)^2-sqrt(mu1))/dn(v,mu)
  (multiple-value-bind (v mu root-mu1)
      (ascending-transform u m)
    (let ((d (jacobi-dn v mu)))
      (* (/ (+ 1 root-mu1) mu)
	 (/ (- (* d d) root-mu1)
	    d)))))

;;; Elliptic Integrals
;;;
;; Translation of Jim FitzSimons' bigfloat implementation of elliptic
;; integrals from http://www.getnet.com/~cherry/elliptbf3.mac.
;;
;; The algorithms are based on B.C. Carlson's "Numerical Computation
;; of Real or Complex Elliptic Integrals".  These are updated to the
;; algorithms in Journal of Computational and Applied Mathematics 118
;; (2000) 71-85 "Reduction Theorems for Elliptic Integrands with the
;; Square Root of two quadritic factors"
;;

(defun errtol (&rest args)
  ;; Compute error tolerance as sqrt(<float-precision>).  Not sure
  ;; this is quite right, but it makes the routines more accurate as
  ;; precision increases increases.
  (sqrt (reduce #'min (mapcar #'(lambda (x)
				  (if (rationalp x)
				      single-float-epsilon
				      (epsilon x)))
			      args))))

(defun carlson-rf (x y z)
  "Compute Carlson's Rf function:

  Rf(x, y, z) = 1/2*integrate((t+x)^(-1/2)*(t+y)^(-1/2)*(t+z)^(-1/2), t, 0, inf)"
  (let* ((precision (float-contagion x y z))
	 (xn (apply-contagion x precision))
	 (yn (apply-contagion y precision))
	 (zn (apply-contagion z precision))
	 (a (/ (+ xn yn zn) 3))
	 (epslon (/ (max (abs (- a xn))
			 (abs (- a yn))
			 (abs (- a zn)))
		    (errtol x y z)))
	 (an a)
	 (power4 1)
	 (n 0)
	 xnroot ynroot znroot lam)
    (loop while (> (* power4 epslon) (abs an))
       do
       (setf xnroot (sqrt xn))
       (setf ynroot (sqrt yn))
       (setf znroot (sqrt zn))
       (setf lam (+ (* xnroot ynroot)
		    (* xnroot znroot)
		    (* ynroot znroot)))
       (setf power4 (* power4 1/4))
       (setf xn (* (+ xn lam) 1/4))
       (setf yn (* (+ yn lam) 1/4))
       (setf zn (* (+ zn lam) 1/4))
       (setf an (* (+ an lam) 1/4))
       (incf n))
    ;; c1=-3/14,c2=1/6,c3=9/88,c4=9/22,c5=-3/22,c6=-9/52,c7=3/26
    (let* ((xndev (/ (* (- a x) power4) an))
	   (yndev (/ (* (- a y) power4) an))
	   (zndev (- (+ xndev yndev)))
	   (ee2 (- (* xndev yndev) (* 6 zndev zndev)))
	   (ee3 (* xndev yndev zndev))
	   (s (+ 1
		 (* -1/10 ee2)
		 (* 1/14 ee3)
		 (* 1/24 ee2 ee2)
		 (* -3/44 ee2 ee3))))
      (/ s (sqrt an)))))

;; rd(x,y,z) = integrate(3/2*(t+x)^(-1/2)*(t+y)^(-1/2)*(t+z)^(-3/2), t, 0, inf)
;;
;; E(K) = rf(0, 1-K^2, 1) - (K^2/3)*rd(0,1-K^2,1)
;;
;; B = integrate(s^2/sqrt(1-s^4), s, 0 ,1)
;;   = beta(3/4,1/2)/4
;;   = sqrt(%pi)*gamma(3/4)/gamma(1/4)
;;   = 1/3*rd(0,2,1)
(defun carlson-rd (x y z)
  "Compute Carlson's Rd function:

  Rd(x,y,z) = integrate(3/2*(t+x)^(-1/2)*(t+y)^(-1/2)*(t+z)^(-3/2), t, 0, inf)"
  (let* ((precision (float-contagion x y z))
	 (xn (apply-contagion x precision))
	 (yn (apply-contagion y precision))
	 (zn (apply-contagion z precision))
	 (a (/ (+ xn yn (* 3 zn)) 5))
	 (epslon (/ (max (abs (- a xn))
			 (abs (- a yn))
			 (abs (- a zn)))
		    (errtol x y z)))
	 (an a)
	 (sigma 0)
	 (power4 1)
	 (n 0)
	 xnroot ynroot znroot lam)
    (loop while (> (* power4 epslon) (abs an))
       do
       (setf xnroot (sqrt xn))
       (setf ynroot (sqrt yn))
       (setf znroot (sqrt zn))
       (setf lam (+ (* xnroot ynroot)
		    (* xnroot znroot)
		    (* ynroot znroot)))
       (setf sigma (+ sigma (/ power4
			       (* znroot (+ zn lam)))))
       (setf power4 (* power4 1/4))
       (setf xn (* (+ xn lam) 1/4))
       (setf yn (* (+ yn lam) 1/4))
       (setf zn (* (+ zn lam) 1/4))
       (setf an (* (+ an lam) 1/4))
       (incf n))
    ;; c1=-3/14,c2=1/6,c3=9/88,c4=9/22,c5=-3/22,c6=-9/52,c7=3/26
    (let* ((xndev (/ (* (- a x) power4) an))
	   (yndev (/ (* (- a y) power4) an))
	   (zndev (- (* (+ xndev yndev) 1/3)))
	   (ee2 (- (* xndev yndev) (* 6 zndev zndev)))
	   (ee3 (* (- (* 3 xndev yndev)
		      (* 8 zndev zndev))
		   zndev))
	   (ee4 (* 3 (- (* xndev yndev) (* zndev zndev)) zndev zndev))
	   (ee5 (* xndev yndev zndev zndev zndev))
	   (s (+ 1
		 (* -3/14 ee2)
		 (* 1/6 ee3)
		 (* 9/88 ee2 ee2)
		 (* -3/22 ee4)
		 (* -9/52 ee2 ee3)
		 (* 3/26 ee5)
		 (* -1/16 ee2 ee2 ee2)
		 (* 3/10 ee3 ee3)
		 (* 3/20 ee2 ee4)
		 (* 45/272 ee2 ee2 ee3)
		 (* -9/68 (+ (* ee2 ee5) (* ee3 ee4))))))
      (+ (* 3 sigma)
	 (/ (* power4 s)
	    (expt an 3/2))))))

;; Complete elliptic integral of the first kind.  This can be computed
;; from Carlson's Rf function:
;;
;;   K(m) = Rf(0, 1 - m, 1)
(defun elliptic-k (m)
  "Complete elliptic integral of the first kind K for parameter m

  K(m) = integrate(1/sqrt(1-m*sin(x)^2), x, 0, %pi/2).

  Note: K(m) = F(%pi/2, m), where F is the (incomplete) elliptic
  integral of the first kind."
  (cond ((= m 0)
	 (/ (float +pi+ m) 2))
	(t
	 (let ((precision (float-contagion m)))
	   (carlson-rf 0 (- 1 m) 1)))))

;; Elliptic integral of the first kind.  This is computed using
;; Carlson's Rf function:
;;
;;  F(phi, m) = sin(phi) * Rf(cos(phi)^2, 1 - m*sin(phi)^2, 1)
(defun elliptic-f (x m)
  "Incomplete Elliptic integral of the first kind:

  F(x, m) = integrate(1/sqrt(1-m*sin(phi)^2), phi, 0, x)

  Note for the complete elliptic integral, you can use elliptic-k"
  (let* ((precision (float-contagion x m))
	 (x (apply-contagion x precision))
	 (m (apply-contagion m precision)))
    (cond ((and (realp m) (realp x))
	   (cond ((> m 1)
		  ;; A&S 17.4.15
		  ;;
		  ;; F(phi|m) = 1/sqrt(m)*F(theta|1/m)
		  ;;
		  ;; with sin(theta) = sqrt(m)*sin(phi)
		  (/ (elliptic-f (asin (* (sqrt m) (sin x))) (/ m))
		     (sqrt m)))
		 ((< m 0)
		  ;; A&S 17.4.17
		  (let* ((m (- m))
			 (m+1 (+ 1 m))
			 (root (sqrt m+1))
			 (m/m+1 (/ m m+1)))
		    (- (/ (elliptic-f (/ (float-pi m) 2) m/m+1)
			  root)
		       (/ (elliptic-f (- (/ (float-pi x) 2) x) m/m+1)
			  root))))
		 ((= m 0)
		  ;; A&S 17.4.19
		  x)
		 ((= m 1)
		  ;; A&S 17.4.21
		  ;;
		  ;; F(phi,1) = log(sec(phi)+tan(phi))
		  ;;          = log(tan(pi/4+pi/2))
		  (log (tan (+ (/ x 2) (/ (float-pi x) 4)))))
		 ((minusp x)
		  (- (elliptic-f (- x) m)))
		 ((> x (float-pi x))
		  ;; A&S 17.4.3
		  (multiple-value-bind (s x-rem)
		      (truncate x (float-pi x))
		    (+ (* 2 s (elliptic-k m))
		       (elliptic-f x-rem m))))
		 ((<= x (/ (float-pi x) 2))
		  (let ((sin-x (sin x))
			(cos-x (cos x))
			(k (sqrt m)))
		    (* sin-x
		       (carlson-rf (* cos-x cos-x)
				   (* (- 1 (* k sin-x))
				      (+ 1 (* k sin-x)))
				   1))))
		 ((< x (float-pi x))
		  (+ (* 2 (elliptic-k m))
		     (elliptic-f (- x (float pi x)) m)))))
	  (t
	   (let ((sin-x (sin x))
		 (cos-x (cos x))
		 (k (sqrt m)))
	     (* sin-x
		(carlson-rf (* cos-x cos-x)
			    (* (- 1 (* k sin-x))
			       (+ 1 (* k sin-x)))
			    1)))))))

;; Incomplete elliptic integral of the second kind.
;;
;; E(phi, m) = integrate(sqrt(1-m*sin(x)^2), x, 0, phi)
;;
(defun elliptic-e (phi m)
  "Incomplete elliptic integral of the second kind:

E(phi, m) = integrate(sqrt(1-m*sin(x)^2), x, 0, phi)"
  (let* ((precision (float-contagion phi m))
	 (phi (apply-contagion phi precision))
	 (m (apply-contagion m precision)))
    (cond ((= m 0)
	   ;; A&S 17.4.23
	   phi)
	  ((= m 1)
	   ;; A&S 17.4.25
	   (sin phi))
	  (t
	   (let* ((sin-phi (sin phi))
		  (cos-phi (cos phi))
		  (k (sqrt m))
		  (y (* (- 1 (* k sin-phi))
			(+ 1 (* k sin-phi)))))
	     (- (* sin-phi
		   (carlson-rf (* cos-phi cos-phi) y 1))
		(* (/ m 3)
		   (expt sin-phi 3)
		   (carlson-rd (* cos-phi cos-phi) y 1))))))))

;; Complete elliptic integral of second kind.
;;
;; E(phi) = integrate(sqrt(1-m*sin(x)^2), x, 0, %pi/2)
;;
(defun elliptic-ec (m)
  "Complete elliptic integral of the second kind:

E(m) = integrate(sqrt(1-m*sin(x)^2), x, 0, %pi/2)"
  (cond ((= m 0)
	 ;; A&S 17.4.23
	 (/ (float-pi m) 2))
	((= m 1)
	 ;; A&S 17.4.25
	 (float 1 m))
	(t
	 (let* ((k (sqrt m))
		(y (* (- 1 k)
		      (+ 1 k))))
	   (- (carlson-rf 0 y 1)
	      (* (/ m 3)
		 (carlson-rd 0 y 1)))))))
