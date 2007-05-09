;;; -*- Mode: lisp -*-

(in-package "QD")

(eval-when (:compile-toplevel)
  (setf *inline-expansion-limit* 1600))

;; We need some object that can hold 4 double-float numbers.  A
;; (complex double-double-float) is perfect for that because CMUCL can
;; handle them without consing.

(deftype %quad-double ()
  '(complex double-double-float))

;; QD-0, QD-1, QD-2, and QD-3 extract the various parts of a
;; quad-double.  QD-0 is the most significant part and QD-3 is the
;; least.
(declaim (inline qd-0 qd-1 qd-2 qd-3))
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

(declaim (inline %make-qd-d))
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


(declaim (inline qd-parts))
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
(declaim (inline three-sum))
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

(declaim (inline three-sum2))
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


#+nil
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
			  #+(or ppc sparc) sqrt-qd))

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

(declaim (inline renorm-4))
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

(declaim (inline renorm-5))
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

(declaim (inline make-qd-d))
(defun make-qd-d (a0 a1 a2 a3)
  "Create a %quad-double from four double-floats, appropriately
  normalizing the result from the four double-floats.
"
  (declare (double-float a0 a1 a2 a3))
  (multiple-value-bind (s0 s1 s2 s3)
      (renorm-4 a0 a1 a2 a3)
    (%make-qd-d s0 s1 s2 s3)))

;;;; Addition

;; Quad-double + double
(declaim (inline add-qd-d))
(defun add-qd-d (a b)
  "Add a quad-double A and a double-float B"
  (declare (type %quad-double a)
	   (double-float b)
	   (optimize (speed 3)))
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

(declaim (inline add-d-qd))
(defun add-d-qd (a b)
  (declare (double-float a)
	   (type %quad-double b)
	   (optimize (speed 3)))
  (add-qd-d b a))

(declaim (inline add-qd-dd))
(defun add-qd-dd (a b)
  "Add a quad-double A and a double-double B"
  (declare (type %quad-double a)
	   (double-double-float b)
	   (optimize (speed 3)))
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

(declaim (inline add-dd-qd))
(defun add-dd-qd (a b)
  (declare (double-double-float a)
	   (type %quad-double b)
	   (optimize (speed 3)))
  (add-qd-dd b a))

(declaim (inline add-qd))

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
	   (optimize (speed 3)))
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

(declaim (inline neg-qd))
(defun neg-qd (a)
  (declare (type %quad-double a))
  (multiple-value-bind (a0 a1 a2 a3)
      (qd-parts a)
    (%make-qd-d (- a0) (- a1) (- a2) (- a3))))

(declaim (inline sub-qd))
(defun sub-qd (a b)
  (declare (type %quad-double a b))
  (add-qd a (neg-qd b)))

(declaim (inline sub-qd-dd))
(defun sub-qd-dd (a b)
  (declare (type %quad-double a)
	   (type double-double-float b))
  (add-qd-dd a (- b)))

(declaim (inline sub-qd-d))
(defun sub-qd-d (a b)
  (declare (type %quad-double a)
	   (type double-float b))
  (add-qd-d a (- b)))

(declaim (inline sub-d-qd))
(defun sub-d-qd (a b)
  (declare (type double-float a)
	   (type %quad-double b))
  (sub-qd (make-qd-d a 0d0 0d0 0d0)
	  b))
  

;; Works
;; (mul-qd-d (sqrt-qd (make-qd-dd 2w0 0w0)) 10d0) ->
;; 14.1421356237309504880168872420969807856967187537694807317667973799q0
;;
;; Clisp says
;; 14.142135623730950488016887242096980785696718753769480731766797379908L0
;;
(declaim (inline mul-qd-d))
(defun mul-qd-d (a b)
  "Multiply quad-double A with B"
  (declare (type %quad-double a)
	   (double-float b)
	   (optimize (speed 3)))
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
(declaim (inline mul-qd-dd))

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
(declaim (inline mul-qd))

;; Works
;; (mul-qd (sqrt-qd (make-qd-dd 2w0 0w0)) (make-qd-dd 10w0 0w0)) ->
;; 14.1421356237309504880168872420969807856967187537694807317667973799q0
;;
;; Clisp says
;; 14.142135623730950488016887242096980785696718753769480731766797379908L0
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

(declaim (inline sqr-qd))
(defun sqr-qd (a)
  "Square A"
  (declare (type %quad-double a)
	   (optimize (speed 3)))
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
	      

(declaim (inline div-qd))
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
(declaim (inline div-qd-d))
(defun div-qd-d (a b)
  (declare (type %quad-double a)
	   (double-float b)
	   (optimize (speed 3)))
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
(declaim (inline div-qd-dd))
(defun div-qd-dd (a b)
  (declare (type %quad-double a)
	   (double-double-float b))
  (let* ((q0 (/ (qd-0 a) (kernel:double-double-hi b)))
	 (r (sub-qd-dd a (* b q0))))
    (let ((q1 (/ (qd-0 r) (kernel:double-double-hi b))))
      (setf r (sub-qd-dd r (* b q1)))
      (let ((q2 (/ (qd-0 r) (kernel:double-double-hi b))))
	(setf r (sub-qd-dd r (* b q2)))
	(let ((q3 (/ (qd-0 r) (kernel:double-double-hi b))))
	  (make-qd-d q0 q1 q2 q3))))))

(declaim (inline make-qd-dd))
(defun make-qd-dd (a0 a1)
  "Create a %quad-double from two double-double-floats"
  (declare (double-double-float a0 a1)
	   (optimize (speed 3) (space 0)))
  (make-qd-d (kernel:double-double-hi a0)
	     (kernel:double-double-lo a0)
	     (kernel:double-double-hi a1)
	     (kernel:double-double-lo a1)))


#+nil
(declaim (ext:end-block))

(declaim (inline abs-qd))
(defun abs-qd (a)
  (declare (type %quad-double a))
  (if (minusp (qd-0 a))
      (neg-qd a)
      a))

(defun sqrt-qd (a)
  (declare (type %quad-double a)
	   (optimize (speed 3) (space 0)))
  ;; Perform the following Newton iteration:
  ;;
  ;;  x' = x + (1 - a * x^2) * x / 2
  ;;
  ;; which converges to 1/sqrt(a).
  (when (= a 0)
    (return-from sqrt-qd #c(0w0 0w0)))

  (let* ((r (make-qd-d (/ (sqrt (the (double-float (0d0))
				  (qd-0 a)))) 0d0 0d0 0d0))
	 (half (make-qd-dd 0.5w0 0w0))
	 (h (mul-qd a half)))
    (declare (type %quad-double r))
    ;;(setf h (mul-qd-d a .5d0))
    (setf r (add-qd r (mul-qd r (sub-qd half (mul-qd h (sqr-qd r))))))
    (setf r (add-qd r (mul-qd r (sub-qd half (mul-qd h (sqr-qd r))))))
    (setf r (add-qd r (mul-qd r (sub-qd half (mul-qd h (sqr-qd r))))))
    (mul-qd r a)))

;; a^n for an integer n
(defun npow (a n)
  (declare (type %quad-double a)
	   (fixnum n))
  (when (= n 0)
    (return-from npow (make-qd-dd 1w0 0w0)))

  (let ((r a)
	(s (make-qd-dd 1w0 0w0))
	(abs-n (abs n)))
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
	(div-qd (make-qd-dd 1w0 0w0) s)
	s)))

;; The n'th root of a.  n is an positive integer and a should be
;; positive too.
(defun nroot-qd (a n)
  (declare (type %quad-double a)
	   (type (and fixnum unsigned-byte) n)
	   (optimize (speed 3)))
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
  (let ((r (make-qd-d (expt (qd-0 a) (- (/ (float n 1d0)))) 0d0 0d0 0d0)))
    (flet ((term ()
	     (div-qd (mul-qd r
			     (add-qd-d (neg-qd (mul-qd a (npow r n)))
				       1d0))
		     (make-qd-dd (float n 1w0) 0w0))))
    (setf r (add-qd r (term)))
    (setf r (add-qd r (term)))
    (setf r (add-qd r (term)))
    (div-qd (make-qd-dd 1w0 0w0) r))))

(declaim (inline qd->))
(defun qd-> (a b)
  "A > B"
  (or (> (qd-0 a) (qd-0 b))
      (and (= (qd-0 a) (qd-0 b))
	   (or (> (qd-1 a) (qd-1 b))
	       (and (= (qd-1 a) (qd-1 b))
		    (or (> (qd-2 a) (qd-2 b))
			(and (= (qd-2 a) (qd-2 b))
			     (> (qd-3 a) (qd-3 b)))))))))



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

(declaim (inline scale-float-qd))
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

(defun qd-to-digits (v &optional position relativep)
  ;; V is the number to be printed.  If RELATIVEP is NIL, POSITION is
  ;; the number of digits to the left of the decimal point where we
  ;; want to stop printing.  If RELATIVEP is non-NIL, POSITION is the
  ;; total number of digits we want printed.
  ;;
  ;; Two values are returned: k, and the digit string, without a
  ;; decimal point.  k is the index into the string, before which the
  ;; decimal point would go.
  (let ((print-base 10)			; B
	(float-radix 2)			; b
	(float-digits (* 4 53)) ; p
	(min-e lisp::double-float-min-e))
    (multiple-value-bind (f e)
	(integer-decode-qd v)
      (let ( ;; FIXME: these even tests assume normal IEEE rounding
	    ;; mode.  I wonder if we should cater for non-normal?
	    (high-ok (evenp f))
	    (low-ok (evenp f))
	    (result (make-array 50 :element-type 'base-char
				:fill-pointer 0 :adjustable t)))
	(labels ((scale (r s m+ m-)
		   ;; Keep increasing k until it's big enough
		   (do ((k 0 (1+ k))
			(s s (* s print-base)))
		       ((not (let ((test (+ r m+)))
			       (or (> test s)
				   (and high-ok (= test s)))))
			;; k is too big.  Decrease until
			(do ((k k (1- k))
			     (r r (* r print-base))
			     (m+ m+ (* m+ print-base))
			     (m- m- (* m- print-base)))
			    ((not (let ((test (* (+ r m+) print-base)))
				    (or (< test s)
					(and (not high-ok) (= test s)))))
			     ;; k is correct.  Generate the digits.
			     (values k (generate r s m+ m-)))))))
		 (generate (r s m+ m-)
		   (multiple-value-bind (d r)
		       (truncate (* r print-base) s)
		     (let ((m+ (* m+ print-base))
			   (m- (* m- print-base)))
		       (let ((tc1 (or (< r m-) (and low-ok (= r m-))))
			     (tc2 (let ((test (+ r m+)))
				    (or (> test s)
					(and high-ok (= test s))))))
			 (cond
			   ((and (not tc1) (not tc2))
			    (vector-push-extend (char lisp::*digits* d) result)
			    ;; FIXME sucky tail recursion.  This whole
			    ;; kaboodle should be DO*/LOOPified.
			    (generate r s m+ m-))
			   ;; pedantically keeping all the conditions
			   ;; in so that we can move them around.
			   ((and (not tc1) tc2)
			    (vector-push-extend (char lisp::*digits* (1+ d)) result)
			    result)
			   ((and tc1 (not tc2))
			    (vector-push-extend (char lisp::*digits* d) result)
			    result)
			   ((and tc1 tc2)
			    (vector-push-extend (char lisp::*digits*
						      (if (< (* r 2) s) d (1+ d)))
						result)
			    result)))))))
	  (let (r s m+ m-)
	    (if (>= e 0)
		(let* ((be (expt float-radix e))
		       (be1 (* be float-radix)))
		  (if (/= f (expt float-radix (1-
					       float-digits)))
		      (setf r (* f be 2)
			    s 2
			    m+ be
			    m- be)
		      (setf r (* f be1 2)
			    s (* float-radix 2)
			    m+ be1
			    m- be)))
		(if (or (= e min-e) 
			(/= f (expt float-radix (1-
						 float-digits))))
		    (setf r (* f 2)
			  s (* (expt float-radix (- e)) 2)
			  m+ 1
			  m- 1)
		    (setf r (* f float-radix 2)
			  s (* (expt float-radix (- 1 e)) 2)
			  m+ float-radix
			  m- 1)))
	    (when position
	      (when relativep
		;;(aver (> position 0))
		(do ((k 0 (1+ k))
		     ;; running out of letters here
		     (l 1 (* l print-base)))
		    ((>= (* s l) (+ r m+))
		     ;; k is now \hat{k}
		     (if (< (+ r (* s (/ (expt print-base (- k
							     position)) 2)))
			    (* s (expt print-base k)))
			 (setf position (- k position))
			 (setf position (- k position 1))))))
	      (let ((low (max m- (/ (* s (expt print-base
					       position)) 2)))
		    (high (max m+ (/ (* s (expt print-base
						position)) 2))))
		(when (<= m- low)
		  (setf m- low)
		  (setf low-ok t))
		(when (<= m+ high)
		  (setf m+ high)
		  (setf high-ok t))))
	    (scale r s m+ m-)))))))

(defun qd-print-exponent (x exp stream)
  (let ((*print-radix* nil))
    (format stream "q~D" exp)))

(defun qd-output-aux (x &optional (stream *standard-output*))
  (multiple-value-bind (e string)
      (qd-to-digits x)
    (when (minusp (float-sign (qd-0 x)))
      (write-char #\- stream))
    (cond ((< -3 e 8)
	   ;; Free format
	   (cond ((plusp e)
	      (write-string string stream :end (min (length string) e))
	      (dotimes (i (- e (length string)))
		(write-char #\0 stream))
	      (write-char #\. stream)
	      (write-string string stream :start (min (length string) e))
	      (when (<= (length string) e)
		(write-char #\0 stream))
	      (qd-print-exponent x 0 stream))
	     (t
	      (write-string "0." stream)
	      (dotimes (i (- e))
		(write-char #\0 stream))
	      (write-string string stream)
	      (qd-print-exponent x 0 stream))))
      (t
       ;; Exponential format
       (write-string string stream :end 1)
       (write-char #\. stream)
       (write-string string stream :start 1)
       ;; CLHS 22.1.3.1.3 says at least one digit must be printed
       ;; after the decimal point.
       (when (= (length string) 1)
	 (write-char #\0 stream))
       (qd-print-exponent x (1- e) stream)))))

;; Function that can be used with FORMAT ~/
(defun qd-format (stream arg colon-p at-sign-p &rest par)
  ;; We should do something with colon-p and at-sign-p
  (declare (ignore colon-p at-sign-p par))
  (qd-output-aux arg stream))

(defun read-qd (stream)
  (labels ((read-digits (s)
	     ;; Read a sequence of digits and return the decimal value
	     ;; and the character that terminated the sequence.
	     (let ((val 0)
		   (ch nil)
		   (count 0))
	       (loop
		   (setf ch (peek-char nil s nil :eof))
		   (cond ((eq ch :eof)
			  (return))
			 ((digit-char-p ch)
			  (read-char s)
			  (incf count)
			  (setf val (+ (digit-char-p ch)
				       (* 10 val))))
			 (t
			  (return))))
	       (values ch val count)))
	   (read-sign (s)
	     (let ((maybe-sign (peek-char t s nil :eof)))
	       (cond ((eql maybe-sign #\-)
		      (read-char s)
		      -1
		      )
		     ((eql maybe-sign #\+)
		      (read-char s)
		      +1)
		     ((and (characterp maybe-sign)
			   (digit-char-p maybe-sign))
		      +1)
		     ((eql maybe-sign #\.)
		      +1)
		     (t
		      maybe-sign))))
	   (read-exp (s)
	     (let ((exp-sign (read-sign s)))
	       (when (eq :eof exp-sign)
		 (return-from read-exp 0))
	       (when (eq :eof (peek-char t s nil :eof))
		 (error "Sign of exponent, but no value"))
	       (multiple-value-bind (char expo)
		   (read-digits s)
		 (declare (ignore char))
		 (* exp-sign expo))))
	   (make-float (sign int-part frac-part scale exp)
	     (declare (type (member -1 1) sign)
		      (type unsigned-byte int-part frac-part)
		      (fixnum scale exp))
	     #+(or)
	     (progn
	       (format t "sign      = ~A~%" sign)
	       (format t "int-part  = ~A~%" int-part)
	       (format t "frac-part = ~A~%" frac-part)
	       (format t "scale     = ~A~%" scale)
	       (format t "exp       = ~A~%" exp))
	     (let ((int (+ (* int-part (expt 10 scale))
			   frac-part))
		   (power (- exp scale)))
	       #+(or)
	       (format t "~A * ~A * 10^(~A)~%" sign int power)
	       (let* ((len (integer-length int)))
		 #+(or)
		 (format t "len = ~A~%" len)
		 (cond ((<= len 106)
			(let ((xx (make-qd-dd (* sign (float int 1w0)) 0w0))
			      (yy (npow (make-qd-dd 10w0 0w0)
					power)))
			  #+(or)
			  (progn
			    (format t "int = ~A~%" int)
			    (format t "fl  = ~A~%" (float int 1w0))
			    (format t "s   = ~A~%" sign)
			    (format t "sint = ~A~%" (* sign  (float int 1w0)))
			    (format t "~A~%" xx)
			    (format t "npow = ~A~%" yy))
			  (mul-qd xx yy)))
		       (t
			(let* ((hi (ldb (byte 106 (- len 106)) int))
			       (lo (ldb (byte 106 (- len 212)) int))
			       (xx (make-qd-dd (* sign (scale-float (float hi 1w0)
								    (- len 106)))
					       (* sign (scale-float (float lo 1w0)
								    (- len 106 106)))))
			       (yy (npow (make-qd-dd 10w0 0w0)
					 power)))
			  #+(or)
			  (progn
			    (format t "xx = ~A~%" xx)
			    (format t "   = ~/qd::qd-format/~%" xx)
			    (format t "yy = ~A~%" yy)
			    (format t "   = ~/qd::qd-format/~%" yy)
			    (format t "hi = ~X (~A)~%" hi
				    (scale-float (float hi 1w0)
						 (- len 106)))
			    (format t "lo = ~X (~A)~%" lo
				    (scale-float (float lo 1w0)
						 (- len 106 106)))
			    (format t "~/qd::qd-format/~%" (mul-qd xx yy)))
			  (mul-qd xx yy))))))))
    (let ((sign (read-sign stream))
	  (int-part 0)
	  (frac-part 0)
	  (exp 0)
	  (scale 0))
      (when (eq :eof (peek-char t stream nil :eof))
	(error "Sign, but no value"))
      (multiple-value-bind (char int)
	  (read-digits stream)
	(setf int-part int)
	(cond ((eql :eof char)
	       )
	      ((eql char #\.)
	       (read-char stream)
	       (multiple-value-bind (char val scale-val)
		   (read-digits stream)
		 (declare (ignore char))
		 (setf frac-part val)
		 (setf scale scale-val)
		 (let ((next (peek-char nil stream nil :eof)))
		   (when (equalp next #\q)
		     (read-char stream)
		     (let ((exp-sign (read-sign stream)))
		       (setf exp (read-exp stream))
		       (setf exp (* exp exp-sign)))))))
	      ((equalp char #\q)
	       (read-char stream)
	       (setf exp (read-exp stream))
	       ))
	(make-float sign int-part frac-part scale exp)))))

(defun qd-reader (stream subchar arg)
  (read-qd stream))

(set-dispatch-macro-character #\# #\Q #'qd-reader)
			      
