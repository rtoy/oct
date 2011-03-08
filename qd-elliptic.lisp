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
  (let ((s (elliptic-sn-descending u m)))
    (if (and (realp u) (realp m))
	(realpart s)
	s)))

(defun jacobi-dn (u m)
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
  ;; Use the ascending Landen transformation, A&S 16.14.3.
  ;;
  ;; cn(u,m) = (1+sqrt(mu1))/mu * (dn(v,mu)^2-sqrt(mu1))/dn(v,mu)
  (multiple-value-bind (v mu root-mu1)
      (ascending-transform u m)
    (let ((d (jacobi-dn v mu)))
      (* (/ (+ 1 root-mu1) mu)
	 (/ (- (* d d) root-mu1)
	    d)))))

(defun errtol (&rest args)
  ;; Compute error tolerance as sqrt(2^(-fpprec)).  Not sure this is
  ;; quite right, but it makes the routines more accurate as fpprec
  ;; increases.
  (sqrt (reduce #'min (mapcar #'(lambda (x)
				  (if (rationalp x)
				      double-float-epsilon
				      (epsilon x)))
			      args))))

(defun carlson-rf (x y z)
  (let* ((xn x)
	 (yn y)
	 (zn z)
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

(defun elliptic-k (m)
  (cond ((= m 0)
	 (/ (float +pi+ m) 2))
	(t
	 (carlson-rf 0 (- 1 m) 1))))
