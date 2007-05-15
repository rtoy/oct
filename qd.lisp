;;; -*- Mode: lisp -*-

(in-package "QD")

(eval-when (:compile-toplevel)
  (setf *inline-expansion-limit* 1600))

;; We need some object that can hold 4 double-float numbers.  A
;; (complex double-double-float) is perfect for that because CMUCL can
;; handle them without consing.

(deftype %quad-double ()
  '(complex double-double-float))

;; All of the following functions should be inline.
(declaim (inline
	  qd-0 qd-1 qd-2 qd-3
	  %make-qd-d
	  qd-parts
	  three-sum
	  three-sum2))


;; QD-0, QD-1, QD-2, and QD-3 extract the various parts of a
;; quad-double.  QD-0 is the most significant part and QD-3 is the
;; least.
(defun qd-0 (q)
  (declare (type %quad-double q)
	   (optimize (speed 3)))
  (kernel:double-double-hi (realpart q)))
(defun qd-1 (q)
  (declare (type %quad-double q)
	   (optimize (speed 3)))
  (kernel:double-double-lo (realpart q)))
(defun qd-2 (q)
  (declare (type %quad-double q)
	   (optimize (speed 3)))
  (kernel:double-double-hi (imagpart q)))
(defun qd-3 (q)
  (declare (type %quad-double q)
	   (optimize (speed 3)))
  (kernel:double-double-lo (imagpart q)))

(eval-when (:compile-toplevel :load-toplevel :execute)
(defun %make-qd-d (a0 a1 a2 a3)
  "Make a %quad-double from 4 double-floats, exactly using the given
  values.  No check is made to see if the values make sense.  A0 is
  the most significant part and A3 is the least.
" 
  (declare (double-float a0 a1
			 a2 a3))
  (complex (kernel:%make-double-double-float a0 a1)
	   (kernel:%make-double-double-float a2 a3)))
)


(defun qd-parts (qd)
  "Extract the four doubles comprising a quad-double and return them
  as multiple values.  The most significant double is the first value."
  (declare (type %quad-double qd))
  (let ((re (realpart qd))
	(im (imagpart qd)))
    (values (kernel:double-double-hi re)
	    (kernel:double-double-lo re)
	    (kernel:double-double-hi im)
	    (kernel:double-double-lo im))))


#+nil
(defun make-qd-z (z)
  "Create a %quad-double from a complex double-double-float"
  (declare (type (complex double-double-float) z))
  (make-qd-dd (realpart z) (imagpart z)))

;; Internal routines for implementing quad-double.
(defun three-sum (a b c)
  (declare (double-float a b c)
	   (optimize (speed 3)))
  (multiple-value-bind (t1 t2)
      (c::two-sum a b)
    (multiple-value-bind (a t3)
	(c::two-sum c t1)
      (multiple-value-bind (b c)
	  (c::two-sum t2 t3)
	(values a b c)))))

(defun three-sum2 (a b c)
  (declare (double-float a b c)
	   (optimize (speed 3)))
  (multiple-value-bind (t1 t2)
      (c::two-sum a b)
    (multiple-value-bind (a t3)
	(c::two-sum c t1)
      (values a (+ t2 t3) c))))

;; Not needed????
#+nil
(declaim (inline quick-three-accum))
#+nil
(defun quick-three-accum (a b c)
  (declare (double-float a b c)
	   (optimize (speed 3) (space 0)))
  (multiple-value-bind (s b)
      (c::two-sum b c)
    (multiple-value-bind (s a)
	(c::two-sum a s)
      (let ((za (/= a 0))
	    (zb (/= b 0)))
	(when (and za zb)
	  (return-from quick-three-accum (values (+ s 0d0) (+ a 0d0) (+ b 0d0))))
	(when (not za)
	  (setf a s)
	  (setf s 0d0))
	(when (not zb)
	  (setf b a)
	  (setf a s))
	(values 0d0 a b)))))


;; These functions are quite short, so we inline them to minimize
;; consing.
(declaim (inline make-qd-d
		 add-d-qd
		 add-dd-qd
		 neg-qd
		 sub-qd
		 sub-qd-dd
		 sub-qd-d
		 sub-d-qd
		 make-qd-dd
		 abs-qd
		 qd->
		 zerop-qd
		 onep-qd
		 plusp-qd
		 minusp-qd
		 qd-=
		 scale-float-qd))

;; Should these functions be inline?  The QD C++ source has these as
;; inline functions, but these are fairly large functions.  However,
;; inlining makes quite a big difference in speed and consing.
(declaim (#+qd-inline inline
	 #-qd-inline maybe-inline
	 renorm-4
	 renorm-5
	 add-qd-d
	 add-qd-dd
	 add-qd
	 mul-qd-d
	 mul-qd-dd
	 mul-qd
	 sqr-qd
	 div-qd
	 div-qd-d
	 div-qd-dd))

#-qd-inline
(declaim (ext:start-block renorm-4 renorm-5
			  make-qd-d
			  add-qd-d add-d-qd add-qd-dd
			  add-dd-qd
			  add-qd
			  neg-qd
			  sub-qd sub-qd-dd sub-qd-d sub-d-qd
			  mul-qd-d mul-qd-dd mul-qd
			  sqr-qd
			  div-qd div-qd-d div-qd-dd
			  make-qd-dd
			  ))

#+(or)
(defun quick-renorm (c0 c1 c2 c3 c4)
  (declare (double-float c0 c1 c2 c3 c4)
	   (optimize (speed 3)))
  (multiple-value-bind (s t3)
      (c::quick-two-sum c3 c4)
    (multiple-value-bind (s t2)
	(c::quick-two-sum c2 s)
      (multiple-value-bind (s t1)
	  (c::quick-two-sum c1 s)
	(multiple-value-bind (c0 t0)
	    (c::quick-two-sum c0 s)
	  (multiple-value-bind (s t2)
	      (c::quick-two-sum t2 t3)
	    (multiple-value-bind (s t1)
		(c::quick-two-sum t1 s)
	      (multiple-value-bind (c1 t0)
		  (c::quick-two-sum t0 s)
		(multiple-value-bind (s t1)
		    (c::quick-two-sum t1 t2)
		  (multiple-value-bind (c2 t0)
		      (c::quick-two-sum t0 s)
		    (values c0 c1 c2 (+ t0 t1))))))))))))

(defun renorm-4 (c0 c1 c2 c3)
  (declare (double-float c0 c1 c2 c3)
	   (optimize (speed 3) (safety 0)))
  (let ((s2 0d0)
	(s3 0d0))
    (multiple-value-bind (s0 c3)
	(c::quick-two-sum c2 c3)
      (multiple-value-bind (s0 c2)
	  (c::quick-two-sum c1 s0)
	(multiple-value-bind (c0 c1)
	    (c::quick-two-sum c0 s0)
	  (let ((s0 c0)
		(s1 c1))
	    (cond ((/= s1 0)
		   (multiple-value-setq (s1 s2)
		     (c::quick-two-sum s1 c2))
		   (if (/= s2 0)
		       (multiple-value-setq (s2 s3)
			 (c::quick-two-sum s2 c3))
		       (multiple-value-setq (s1 s2)
			 (c::quick-two-sum s1 c3))))
		  (t
		   (multiple-value-setq (s0 s1)
		     (c::quick-two-sum s0 c2))
		   (if (/= s1 0)
		       (multiple-value-setq (s1 s2)
			 (c::quick-two-sum s1 c3))
		       (multiple-value-setq (s0 s1)
			 (c::quick-two-sum s0 c3)))))
	    (values s0 s1 s2 s3)))))))

(defun renorm-5 (c0 c1 c2 c3 c4)
  (declare (double-float c0 c1 c2 c3)
	   (optimize (speed 3) (safety 0)))
  (let ((s2 0d0)
	(s3 0d0))
    (declare (double-float s2 s3))
    (multiple-value-bind (s0 c4)
	(c::quick-two-sum c3 c4)
      (multiple-value-bind (s0 c3)
	  (c::quick-two-sum c2 s0)
	(multiple-value-bind (s0 c2)
	    (c::quick-two-sum c1 s0)
	  (multiple-value-bind (c0 c1)
	      (c::quick-two-sum c0 s0)
	    (let ((s0 c0)
		  (s1 c1))
	      (declare (double-float s0 s1))
	      (multiple-value-setq (s0 s1)
		(c::quick-two-sum c0 c1))
	      (cond ((/= s1 0)
		     (multiple-value-setq (s1 s2)
		       (c::quick-two-sum s1 c2))
		     (cond ((/= s2 0)
			    (multiple-value-setq (s2 s3)
			      (c::quick-two-sum s2 c3))
			    (if (/= s3 0)
				(incf s3 c4)
				(incf s2 c4)))
			   (t
			    (multiple-value-setq (s1 s2)
			      (c::quick-two-sum s1 c3))
			    (if (/= s2 0)
				(multiple-value-setq (s2 s3)
				  (c::quick-two-sum s2 c4))
				(multiple-value-setq (s1 s2)
				  (c::quick-two-sum s1 c4))))))
		    (t
		     (multiple-value-setq (s0 s1)
		       (c::quick-two-sum s0 c2))
		     (cond ((/= s1 0)
			    (multiple-value-setq (s1 s2)
			      (c::quick-two-sum s1 c3))
			    (if (/= s2 0)
				(multiple-value-setq (s2 s3)
				  (c::quick-two-sum s2 c4))
				(multiple-value-setq (s1 s2)
				  (c::quick-two-sum s1 c4))))
			   (t
			    (multiple-value-setq (s0 s1)
			      (c::quick-two-sum s0 c3))
			    (if (/= s1 0)
				(multiple-value-setq (s1 s2)
				  (c::quick-two-sum s1 c4))
				(multiple-value-setq (s0 s1)
				  (c::quick-two-sum s0 c4)))))))
	      (values s0 s1 s2 s3))))))))

(defun make-qd-d (a0 &optional (a1 0d0 a1-p) (a2 0d0) (a3 0d0))
  "Create a %quad-double from four double-floats, appropriately
  normalizing the result from the four double-floats.
"
  (declare (double-float a0 a1 a2 a3)
	   (optimize (speed 3)
		     #+cmu
		     (ext:inhibit-warnings 3)))
  (if a1-p
      (multiple-value-bind (s0 s1 s2 s3)
	  (renorm-4 a0 a1 a2 a3)
	(%make-qd-d s0 s1 s2 s3))
      (%make-qd-d a0 0d0 0d0 0d0)))

;;;; Addition

;; Quad-double + double
(defun add-qd-d (a b)
  "Add a quad-double A and a double-float B"
  (declare (type %quad-double a)
	   (double-float b)
	   (optimize (speed 3)
		     (space 0)))
  (multiple-value-bind (c0 e)
      (c::two-sum (qd-0 a) b)
    (multiple-value-bind (c1 e)
	(c::two-sum (qd-1 a) e)
      (multiple-value-bind (c2 e)
	  (c::two-sum (qd-2 a) e)
	(multiple-value-bind (c3 e)
	    (c::two-sum (qd-3 a) e)
	  (multiple-value-call #'%make-qd-d
	    (renorm-5 c0 c1 c2 c3 e)))))))

(defun add-d-qd (a b)
  (declare (double-float a)
	   (type %quad-double b)
	   (optimize (speed 3)))
  (add-qd-d b a))

(defun add-qd-dd (a b)
  "Add a quad-double A and a double-double B"
  (declare (type %quad-double a)
	   (double-double-float b)
	   (optimize (speed 3)
		     (space 0)))
  (multiple-value-bind (s0 t0)
      (c::two-sum (qd-0 a) (kernel:double-double-hi b))
    (multiple-value-bind (s1 t1)
	(c::two-sum (qd-1 a) (kernel:double-double-lo b))
      (multiple-value-bind (s1 t0)
	  (c::two-sum s1 t0)
	(multiple-value-bind (s2 t0 t1)
	    (three-sum (qd-2 a) t0 t1)
	  (multiple-value-bind (s3 t0)
	      (c::two-sum t0 (qd-3 a))
	    (let ((t0 (+ t0 t1)))
	      (multiple-value-call #'%make-qd-d
		(renorm-5 s0 s1 s2 s3 t0)))))))))

(defun add-dd-qd (a b)
  (declare (double-double-float a)
	   (type %quad-double b)
	   (optimize (speed 3)
		     (space 0)))
  (add-qd-dd b a))


#+(or)
(defun add-qd-1 (a b)
  (declare (type %quad-double a b)
	   (optimize (speed 3)))
  (multiple-value-bind (s0 t0)
      (c::two-sum (qd-0 a) (qd-0 b))
    (multiple-value-bind (s1 t1)
	(c::two-sum (qd-1 a) (qd-1 b))
      (multiple-value-bind (s2 t2)
	  (c::two-sum (qd-2 a) (qd-2 b))
	(multiple-value-bind (s3 t3)
	    (c::two-sum (qd-3 a) (qd-3 b))
	  (multiple-value-bind (s1 t0)
	      (c::two-sum s1 t0)
	    (multiple-value-bind (s2 t0 t1)
		(three-sum s2 t0 t1)
	      (multiple-value-bind (s3 t0)
		  (three-sum2 s3 t0 t2)
		(let ((t0 (+ t0 t1 t3)))
		  (multiple-value-call #'%make-qd-d
		    (renorm-5 s0 s1 s2 s3 t0)))))))))))

;; Same as above, except that everything is expanded out for compilers
;; which don't do a very good job with dataflow.  CMUCL is one of
;; those compilers.

(defun add-qd (a b)
  (declare (type %quad-double a b)
	   (optimize (speed 3)
		     (space 0)))
  ;; This is the version that is NOT IEEE.  Should we use the IEEE
  ;; version?  It's quite a bit more complicated.
  ;;
  ;; In addition, this is reorganized to minimize data dependency.
  (multiple-value-bind (a0 a1 a2 a3)
      (qd-parts a)
    (multiple-value-bind (b0 b1 b2 b3)
	(qd-parts b)
      (let ((s0 (+ a0 b0))
	    (s1 (+ a1 b1))
	    (s2 (+ a2 b2))
	    (s3 (+ a3 b3)))
	(declare (double-float s0 s1 s2 s3))
	(let ((v0 (- s0 a0))
	      (v1 (- s1 a1))
	      (v2 (- s2 a2))
	      (v3 (- s3 a3)))
	  (let ((u0 (- s0 v0))
		(u1 (- s1 v1))
		(u2 (- s2 v2))
		(u3 (- s3 v3)))
	    (let ((w0 (- a0 u0))
		  (w1 (- a1 u1))
		  (w2 (- a2 u2))
		  (w3 (- a3 u3)))
	      (let ((u0 (- b0 v0))
		    (u1 (- b1 v1))
		    (u2 (- b2 v2))
		    (u3 (- b3 v3)))
		(let ((t0 (+ w0 u0))
		      (t1 (+ w1 u1))
		      (t2 (+ w2 u2))
		      (t3 (+ w3 u3)))
		  (multiple-value-bind (s1 t0)
		      (c::two-sum s1 t0)
		    (multiple-value-bind (s2 t0 t1)
			(three-sum s2 t0 t1)
		      (multiple-value-bind (s3 t0)
			  (three-sum2 s3 t0 t2)
			(declare (double-float t0))
			(setf t0 (+ t0 t1 t3))
			;; Renormalize
			(multiple-value-setq (s0 s1 s2 s3)
			  (renorm-5 s0 s1 s2 s3 t0))
			(%make-qd-d s0 s1 s2 s3)))))))))))))

(defun neg-qd (a)
  (declare (type %quad-double a))
  (multiple-value-bind (a0 a1 a2 a3)
      (qd-parts a)
    (%make-qd-d (- a0) (- a1) (- a2) (- a3))))

(defun sub-qd (a b)
  (declare (type %quad-double a b))
  (add-qd a (neg-qd b)))

(defun sub-qd-dd (a b)
  (declare (type %quad-double a)
	   (type double-double-float b))
  (add-qd-dd a (- b)))

(defun sub-qd-d (a b)
  (declare (type %quad-double a)
	   (type double-float b))
  (add-qd-d a (- b)))

(defun sub-d-qd (a b)
  (declare (type double-float a)
	   (type %quad-double b))
  ;; a - b = a + (-b)
  (add-d-qd a (neg-qd b)))
  

;; Works
;; (mul-qd-d (sqrt-qd (make-qd-dd 2w0 0w0)) 10d0) ->
;; 14.1421356237309504880168872420969807856967187537694807317667973799q0
;;
;; Clisp says
;; 14.142135623730950488016887242096980785696718753769480731766797379908L0
;;
(defun mul-qd-d (a b)
  "Multiply quad-double A with B"
  (declare (type %quad-double a)
	   (double-float b)
	   (optimize (speed 3)
		     (space 0)))
  (multiple-value-bind (p0 q0)
      (c::two-prod (qd-0 a) b)
    (multiple-value-bind (p1 q1)
	(c::two-prod (qd-1 a) b)
      (multiple-value-bind (p2 q2)
	  (c::two-prod (qd-2 a) b)
	(let ((p3 (* (qd-3 a) b))
	      (s0 p0))
	  #+nil
	  (format t "q0 p1 = ~A ~A~%" q0 p1)
	  (multiple-value-bind (s1 s2)
	      (c::two-sum q0 p1)
	    #+nil
	    (format t "s1 s2 = ~A ~A~%" s1 s2)
	    #+nil
	    (format t "s2 q1 p2 = ~A ~A ~A~%"
		    s2 q1 p2)
	    (multiple-value-bind (s2 q1 p2)
		(three-sum s2 q1 p2)
	      #+nil
	      (format t "s2,q1,p2 = ~A ~A ~A~%"
		      s2 q1 p2)
	      #+nil
	      (format t "q1 q2 p3 = ~A ~A ~A~%"
		      q1 q2 p3)
	      (multiple-value-bind (q1 q2)
		  (three-sum2 q1 q2 p3)
		#+nil
		(format t "q1 q2, p3 = ~A ~A ~A~%" q1 q2 p3)
		(let ((s3 q1)
		      (s4 (+ q2 p2)))
		  #+nil
		  (progn
		    (format t "~a~%" s0)
		    (format t "~a~%" s1)
		    (format t "~a~%" s2)
		    (format t "~a~%" s3)
		    (format t "~a~%" s4))
		  (multiple-value-bind (s0 s1 s2 s3)
		      (renorm-5 s0 s1 s2 s3 s4)
		    #+nil
		    (progn
		      (format t "~a~%" s0)
		      (format t "~a~%" s1)
		      (format t "~a~%" s2)
		      (format t "~a~%" s3)
		      (format t "~a~%" s4))
		    (%make-qd-d s0 s1 s2 s3)))))))))))

;; a0 * b0                        0
;;      a0 * b1                   1
;;      a1 * b0                   2
;;           a1 * b1              3
;;           a2 * b0              4
;;                a2 * b1         5
;;                a3 * b0         6
;;                     a3 * b1    7

;; Not working.
;; (mul-qd-dd (sqrt-qd (make-qd-dd 2w0 0w0)) 10w0) ->
;; 14.142135623730950488016887242097022172449805747901877456053837224q0
;;
;; But clisp says
;; 14.142135623730950488016887242096980785696718753769480731766797379908L0
;;                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;;
;; Running a test program using qd (2.1.210) shows that we get the
;; same wrong answer.
#+(or)
(defun mul-qd-dd (a b)
  (declare (type %quad-double a)
	   (double-double-float b)
	   (optimize (speed 3)))
  (multiple-value-bind (p0 q0)
      (c::two-prod (qd-0 a) (kernel:double-double-hi b))
    (multiple-value-bind (p1 q1)
	(c::two-prod (qd-0 a) (kernel:double-double-lo b))
      (multiple-value-bind (p2 q2)
	  (c::two-prod (qd-1 a) (kernel:double-double-hi b))
	(multiple-value-bind (p3 q3)
	    (c::two-prod (qd-1 a) (kernel:double-double-lo b))
	  (multiple-value-bind (p4 q4)
	      (c::two-prod (qd-2 a) (kernel:double-double-hi b))
	    (format t "p0, q0 = ~A ~A~%" p0 q0)
	    (format t "p1, q1 = ~A ~A~%" p1 q1)
	    (format t "p2, q2 = ~A ~A~%" p2 q2)
	    (format t "p3, q3 = ~A ~A~%" p3 q3)
	    (format t "p4, q4 = ~A ~A~%" p4 q4)
	    (multiple-value-setq (p1 p2 q0)
	      (three-sum p1 p2 q0))
	    (format t "p1 = ~A~%" p1)
	    (format t "p2 = ~A~%" p2)
	    (format t "q0 = ~A~%" q0)
	    ;; five-three-sum
	    (multiple-value-setq (p2 p3 p4)
	      (three-sum p2 p3 p4))
	    (format t "p2 = ~A~%" p2)
	    (format t "p3 = ~A~%" p3)
	    (format t "p4 = ~A~%" p4)
	    (multiple-value-setq (q1 q2)
	      (c::two-sum q1 q2))
	    (multiple-value-bind (s0 t0)
		(c::two-sum p2 q1)
	      (multiple-value-bind (s1 t1)
		  (c::two-sum p3 q2)
		(multiple-value-setq (s1 t0)
		  (c::two-sum s1 t0))
		(let ((s2 (+ t0 t1 p4))
		      (p2 s0)
		      (p3 (+ (* (qd-2 a)
				(kernel:double-double-hi b))
			     (* (qd-3 a)
				(kernel:double-double-lo b))
			     q3 q4)))
		  (multiple-value-setq (p3 q0 s1)
		    (three-sum2 p3 q0 s1))
		  (let ((p4 (+ q0 s2)))
		    (multiple-value-call #'%make-qd-d
		      (renorm-5 p0 p1 p2 p3 p4))))))))))))

;; a0 * b0                    0
;;      a0 * b1               1
;;      a1 * b0               2
;;           a0 * b2          3
;;           a1 * b1          4
;;           a2 * b0          5
;;                a0 * b3     6
;;                a1 * b2     7
;;                a2 * b1     8
;;                a3 * b0     9 

;; Works
;; (mul-qd (sqrt-qd (make-qd-dd 2w0 0w0)) (make-qd-dd 10w0 0w0)) ->
;; 14.1421356237309504880168872420969807856967187537694807317667973799q0
;;
;; Clisp says
;; 14.142135623730950488016887242096980785696718753769480731766797379908L0
(defun mul-qd (a b)
  (declare (type %quad-double a b)
	   (optimize (speed 3)
		     (space 0)))
  (multiple-value-bind (a0 a1 a2 a3)
      (qd-parts a)
    (multiple-value-bind (b0 b1 b2 b3)
	(qd-parts b)
      (multiple-value-bind (p0 q0)
	  (c::two-prod a0 b0)
	(multiple-value-bind (p1 q1)
	    (c::two-prod a0 b1)
	  (multiple-value-bind (p2 q2)
	      (c::two-prod a1 b0)
	    (multiple-value-bind (p3 q3)
		(c::two-prod a0 b2)
	      (multiple-value-bind (p4 q4)
		  (c::two-prod a1 b1)
		(multiple-value-bind (p5 q5)
		    (c::two-prod a2 b0)
		  ;; Start accumulation
		  (multiple-value-setq (p1 p2 q0)
		    (three-sum p1 p2 q0))

		  ;; six-three-sum of p2, q1, q2, p3, p4, p5
		  (multiple-value-setq (p2 q1 q2)
		    (three-sum p2 q1 q2))
		  (multiple-value-setq (p3 p4 p5)
		    (three-sum p3 p4 p5))
		  ;; Compute (s0,s1,s2) = (p2,q1,q2) + (p3,p4,p5)
		  (multiple-value-bind (s0 t0)
		      (c::two-sum p2 p3)
		    (multiple-value-bind (s1 t1)
			(c::two-sum q1 p4)
		      (let ((s2 (+ q2 p5)))
			(declare (double-float s2))
			(multiple-value-bind (s1 t0)
			    (c::two-sum s1 t0)
			  (declare (double-float s1))
			  (incf s2 (+ t0 t1))
			  ;; O(eps^3) order terms.  This is the sloppy
			  ;; multiplication version.  Should we use
			  ;; the precise version?  It's significantly
			  ;; more complex.
			  
			  (incf s1 (+ (* a0 b3)
				      (* a1 b2)
				      (* a2 b1)
				      (* a3 b0)
				      q0 q3 q4 q5))
			  #+nil
			  (format t "p0,p1,s0,s1,s2 = ~a ~a ~a ~a ~a~%"
				  p0 p1 s0 s1 s2)
			  (multiple-value-bind (p0 p1 s0 s1)
			      (renorm-5 p0 p1 s0 s1 s2)
			    (%make-qd-d p0 p1 s0 s1)))))))))))))))

;; This is the non-sloppy version.  I think this works just fine, but
;; since qd defaults to the sloppy multiplication version, we do the
;; same.
#+nil
(defun mul-qd (a b)
  (declare (type %quad-double a b)
	   (optimize (speed 3)))
  (multiple-value-bind (a0 a1 a2 a3)
      (qd-parts a)
    (multiple-value-bind (b0 b1 b2 b3)
	(qd-parts b)
      (multiple-value-bind (p0 q0)
	  (c::two-prod a0 b0)
	(multiple-value-bind (p1 q1)
	    (c::two-prod a0 b1)
	  (multiple-value-bind (p2 q2)
	      (c::two-prod a1 b0)
	    (multiple-value-bind (p3 q3)
		(c::two-prod a0 b2)
	      (multiple-value-bind (p4 q4)
		  (c::two-prod a1 b1)
		(declare (double-float q4))
		#+nil
		(progn
		  (format t"p0, q0 = ~a ~a~%" p0 q0)
		  (format t"p1, q1 = ~a ~a~%" p1 q1)
		  (format t"p2, q2 = ~a ~a~%" p2 q2)
		  (format t"p3, q3 = ~a ~a~%" p3 q3)
		  (format t"p4, q4 = ~a ~a~%" p4 q4))
		(multiple-value-bind (p5 q5)
		    (c::two-prod a2 b0)
		  #+nil
		  (format t"p5, q5 = ~a ~a~%" p5 q5)
		  ;; Start accumulation
		  (multiple-value-setq (p1 p2 q0)
		    (three-sum p1 p2 q0))
		  #+nil
		  (format t "p1 p2 q0 = ~a ~a ~a~%" p1 p2 q0)
		  ;; six-three-sum of p2, q1, q2, p3, p4, p5
		  (multiple-value-setq (p2 q1 q2)
		    (three-sum p2 q1 q2))
		  (multiple-value-setq (p3 p4 p5)
		    (three-sum p3 p4 p5))
		  ;; Compute (s0,s1,s2) = (p2,q1,q2) + (p3,p4,p5)
		  (multiple-value-bind (s0 t0)
		      (c::two-sum p2 p3)
		    (multiple-value-bind (s1 t1)
			(c::two-sum q1 p4)
		      (declare (double-float t1))
		      (let ((s2 (+ q2 p5)))
			(declare (double-float s2))
			(multiple-value-bind (s1 t0)
			    (c::two-sum s1 t0)
			  (declare (double-float s1))
			  (incf s2 (+ t0 t1))
			  (multiple-value-bind (p6 q6)
			      (c::two-prod a0 b3)
			    (multiple-value-bind (p7 q7)
				(c::two-prod a1 b2)
			      (multiple-value-bind (p8 q8)
				  (c::two-prod a2 b1)
				(multiple-value-bind (p9 q9)
				    (c::two-prod a3 b0)
				  (multiple-value-setq (q0 q3)
				    (c::two-sum q0 q3))
				  (multiple-value-setq (q4 q5)
				    (c::two-sum q4 q5))
				  (multiple-value-setq (p6 p7)
				    (c::two-sum p6 p7))
				  (multiple-value-setq (p8 p9)
				    (c::two-sum p8 p9))
				  ;; (t0,t1) = (q0,q3)+(q4,q5)
				  (multiple-value-setq (t0 t1)
				    (c::two-sum q0 q4))
				  (setf t1 (+ q3 q5))
				  ;; (r0,r1) = (p6,p7)+(p8,p9)
				  (multiple-value-bind (r0 r1)
				      (c::two-sum p6 p8)
				    (declare (double-float r1))
				    (incf r1 (+ p7 p9))
				    ;; (q3,q4) = (t0,t1)+(r0,r1)
				    (multiple-value-setq (q3 q4)
				      (c::two-sum t0 r0))
				    (incf q4 (+ t1 r1))
				    ;; (t0,t1)=(q3,q4)+s1
				    (multiple-value-setq (t0 t1)
				      (c::two-sum q3 s1))
				    (incf t1 q4)
				    ;; O(eps^4) terms
				    (incf t1
					  (+ (* a1 b3)
					     (* a2 b2)
					     (* a3 b1)
					     q6 q7 q8 q9
					     s2))
				    #+nil
				    (format t "p0,p1,s0,t0,t1 = ~a ~a ~a ~a ~a~%"
					    p0 p1 s0 t0 t1)
				    (multiple-value-call #'%make-qd-d
				      (renorm-5 p0 p1 s0 t0 t1))))))))))))))))))))

(defun sqr-qd (a)
  "Square A"
  (declare (type %quad-double a)
	   (optimize (speed 3)
		     (space 0)))
  (multiple-value-bind (p0 q0)
      (c::two-sqr (qd-0 a))
    (multiple-value-bind (p1 q1)
	(c::two-prod (* 2 (qd-0 a)) (qd-1 a))
      (multiple-value-bind (p2 q2)
	  (c::two-prod (* 2 (qd-0 a)) (qd-2 a))
	(multiple-value-bind (p3 q3)
	    (c::two-sqr (qd-1 a))
	  (multiple-value-setq (p1 q0)
	    (c::two-sum q0 p1))
	  (multiple-value-setq (q0 q1)
	    (c::two-sum q0 q1))
	  (multiple-value-setq (p2 p3)
	    (c::two-sum p2 p3))

	  (multiple-value-bind (s0 t0)
	      (c::two-sum q0 p2)
	    (declare (double-float t0))
	    (multiple-value-bind (s1 t1)
		(c::two-sum q1 p3)
	    (declare (double-float t1))
	      (multiple-value-setq (s1 t0)
		(c::two-sum s1 t0))
	      (incf t0 t1)

	      (multiple-value-setq (s1 t0)
		(c::quick-two-sum s1 t0))
	      (multiple-value-setq (p2 t1)
		(c::quick-two-sum s0 s1))
	      (multiple-value-setq (p3 q0)
		(c::quick-two-sum t1 t0))

	      (let ((p4 (* 2 (qd-0 a) (qd-3 a)))
		    (p5 (* 2 (qd-1 a) (qd-2 a))))
		(declare (double-float p4))
		(multiple-value-setq (p4 p5)
		  (c::two-sum p4 p5))
		(multiple-value-setq (q2 q3)
		  (c::two-sum q2 q3))

		(multiple-value-setq (t0 t1)
		  (c::two-sum p4 q2))

		(incf t1 (+ p5 q3))

		(multiple-value-setq (p3 p4)
		  (c::two-sum p3 t0))
		(incf p4 (+ q0 t1))

		(multiple-value-call #'%make-qd-d
		  (renorm-5 p0 p1 p2 p3 p4))))))))))
	      

(defun div-qd (a b)
  (declare (type %quad-double a b)
	   (optimize (speed 3)
		     (space 0)))
  (let ((a0 (qd-0 a))
	(b0 (qd-0 b)))
    (let* ((q0 (/ a0 b0))
	   (r (sub-qd a (mul-qd-d b q0)))
	   (q1 (/ (qd-0 r) b0)))
      (setf r (sub-qd r (mul-qd-d b q1)))
      (let ((q2 (/ (qd-0 r) b0)))
	(setf r (sub-qd r (mul-qd-d b q2)))
	(let ((q3 (/ (qd-0 r) b0)))
	  (multiple-value-bind (q0 q1 q2 q3)
	      (renorm-4 q0 q1 q2 q3)
	    (%make-qd-d q0 q1 q2 q3)))))))

;; Non-sloppy divide
#+(or)
(defun div-qd (a b)
  (declare (type %quad-double a b))
  (let ((a0 (qd-0 a))
	(b0 (qd-0 b)))
    (let* ((q0 (/ a0 b0))
	   (r (sub-qd a (mul-qd-d b q0)))
	   (q1 (/ (qd-0 r) b0)))
      (setf r (sub-qd r (mul-qd-d b q1)))
      (let ((q2 (/ (qd-0 r) b0)))
	(setf r (sub-qd r (mul-qd-d b q2)))
	(let ((q3 (/ (qd-0 r) b0)))
	  (setf r (sub-qd r (mul-qd-d b q3)))
	  (let ((q4 (/ (qd-0 r) b0)))
	  (multiple-value-bind (q0 q1 q2 q3)
	      (renorm-5 q0 q1 q2 q3 q4)
	    (%make-qd-d q0 q1 q2 q3))))))))

;; quad-double / double
(defun div-qd-d (a b)
  (declare (type %quad-double a)
	   (double-float b)
	   (optimize (speed 3)
		     (space 0)))
  ;; Compute approximate quotient using high order doubles, then
  ;; correct it 3 times using the remainder.  Analogous to long
  ;; division.
  (let ((q0 (/ (qd-0 a) b)))
    ;; Compute remainder a - q0 * b
    (multiple-value-bind (t0 t1)
	(c::two-prod q0 b)
      (let ((r (sub-qd-dd a (kernel:make-double-double-float t0 t1))))
	;; First correction
	(let ((q1 (/ (qd-0 r) b)))
	  (multiple-value-bind (t0 t1)
	      (c::two-prod q1 b)
	    (setf r (sub-qd-dd r (kernel:make-double-double-float t0 t1)))
	    ;; Second correction
	    (let ((q2 (/ (qd-0 r) b)))
	      (multiple-value-bind (t0 t1)
		  (c::two-prod q2 b)
		(setf r (sub-qd-dd r (kernel:make-double-double-float t0 t1)))
		;; Final correction
		(let ((q3 (/ (qd-0 r) b)))
		  (make-qd-d q0 q1 q2 q3))))))))))

;; Sloppy version
(defun div-qd-dd (a b)
  (declare (type %quad-double a)
	   (double-double-float b)
	   (optimize (speed 3)
		     (space 0)))
  (let* ((q0 (/ (qd-0 a) (kernel:double-double-hi b)))
	 (r (sub-qd-dd a (* b q0))))
    (let ((q1 (/ (qd-0 r) (kernel:double-double-hi b))))
      (setf r (sub-qd-dd r (* b q1)))
      (let ((q2 (/ (qd-0 r) (kernel:double-double-hi b))))
	(setf r (sub-qd-dd r (* b q2)))
	(let ((q3 (/ (qd-0 r) (kernel:double-double-hi b))))
	  (make-qd-d q0 q1 q2 q3))))))

(defun make-qd-dd (a0 a1)
  "Create a %quad-double from two double-double-floats"
  (declare (double-double-float a0 a1)
	   (optimize (speed 3) (space 0)))
  (make-qd-d (kernel:double-double-hi a0)
	     (kernel:double-double-lo a0)
	     (kernel:double-double-hi a1)
	     (kernel:double-double-lo a1)))


#-qd-inline
(declaim (ext:end-block))

(defun abs-qd (a)
  (declare (type %quad-double a))
  (if (minusp (qd-0 a))
      (neg-qd a)
      a))

;; a^n for an integer n
(defun npow (a n)
  (declare (type %quad-double a)
	   (fixnum n)
	   (optimize (speed 3)
		     (space 0)))
  (when (= n 0)
    (return-from npow (make-qd-d 1d0)))

  (let ((r a)
	(s (make-qd-d 1d0))
	(abs-n (abs n)))
    (declare (type (and fixnum unsigned-byte) abs-n)
	     (type %quad-double r s))
    (cond ((> abs-n 1)
	   ;; Binary exponentiation
	   (loop while (plusp abs-n)
	     do
	     (when (= 1 (logand abs-n 1))
	       (setf s (mul-qd s r)))
	     (setf abs-n (ash abs-n -1))
	     (when (plusp abs-n)
	       (setf r (sqr-qd r)))))
	  (t
	   (setf s r)))
    (if (minusp n)
	(div-qd (make-qd-d 1d0) s)
	s)))

;; The n'th root of a.  n is an positive integer and a should be
;; positive too.
(defun nroot-qd (a n)
  (declare (type %quad-double a)
	   (type (and fixnum unsigned-byte) n)
	   (optimize (speed 3)
		     (space 0)))
  ;; Use Newton's iteration to solve
  ;;
  ;; 1/(x^n) - a = 0
  ;;
  ;; The iteration becomes
  ;;
  ;; x' = x + x * (1 - a * x^n)/n
  ;;
  ;; Since Newton's iteration converges quadratically, we only need to
  ;; perform it twice.
  (let ((r (make-qd-d (expt (the (double-float (0d0)) (qd-0 a))
			    (- (/ (float n 1d0)))))))
    (declare (type %quad-double r))
    (flet ((term ()
	     (div-qd-d (mul-qd r
			       (add-qd-d (neg-qd (mul-qd a (npow r n)))
					 1d0))
		     (float n 1d0))))
    (setf r (add-qd r (term)))
    (setf r (add-qd r (term)))
    (setf r (add-qd r (term)))
    (div-qd (make-qd-d 1d0) r))))

(defun qd-> (a b)
  "A > B"
  (or (> (qd-0 a) (qd-0 b))
      (and (= (qd-0 a) (qd-0 b))
	   (or (> (qd-1 a) (qd-1 b))
	       (and (= (qd-1 a) (qd-1 b))
		    (or (> (qd-2 a) (qd-2 b))
			(and (= (qd-2 a) (qd-2 b))
			     (> (qd-3 a) (qd-3 b)))))))))

(defun zerop-qd (a)
  "Is A zero?"
  (declare (type %quad-double a))
  (zerop (qd-0 a)))

(defun onep-qd (a)
  "Is A equal to 1?"
  (declare (type %quad-double a))
  (and (= (qd-0 a) 1d0)
       (zerop (qd-1 a))
       (zerop (qd-2 a))
       (zerop (qd-3 a))))

(defun plusp-qd (a)
  "Is A positive?"
  (declare (type %quad-double a))
  (plusp (qd-0 a)))
	 
(defun minusp-qd (a)
  "Is A positive?"
  (declare (type %quad-double a))
  (minusp (qd-0 a)))

(defun qd-= (a b)
  (and (= (qd-0 a) (qd-0 b))
       (= (qd-1 a) (qd-1 b))
       (= (qd-2 a) (qd-2 b))
       (= (qd-3 a) (qd-3 b))))


(defun integer-decode-qd (q)
  (declare (type %quad-double q))
  (multiple-value-bind (hi-int hi-exp sign)
      (integer-decode-float (realpart q))
    (if (zerop (imagpart q))
	(values (ash hi-int 106) (- hi-exp 106) sign)
	(multiple-value-bind (lo-int lo-exp lo-sign)
	    (integer-decode-float (imagpart q))
	  (values (+ (* (* sign lo-sign) lo-int)
		     (ash hi-int (- hi-exp lo-exp)))
		  lo-exp
		  sign)))))

#+nil
(defun scale-float-qd (qd k)
  (make-qd-d (scale-float (qd-0 qd) k)
	     (scale-float (qd-1 qd) k)
	     (scale-float (qd-2 qd) k)
	     (scale-float (qd-3 qd) k)))

(defun scale-float-qd (qd k)
  (declare (type %quad-double qd)
	   (type (integer -1022 1022) k)
	   (optimize (speed 3) (space 0)))
  ;; (space 0) to get scale-double-float inlined
  (let ((scale (scale-float 1d0 k)))
    (%make-qd-d (* (qd-0 qd) scale)
		(* (qd-1 qd) scale)
		(* (qd-2 qd) scale)
		(* (qd-3 qd) scale))))

