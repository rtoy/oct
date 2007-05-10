;;; -*- Mode: lisp -*-

(in-package "QD")

;; Compute to how many bits EST and TRUE are equal.  If they are
;; identical, return T.
(defun bit-accuracy (est true)
  (let ((diff (abs (qd-0 (sub-qd est true)))))
    (if (zerop diff)
	t
	(- (log diff 2d0)))))

;; Machin's formula for pi
#+nil
(defun atan-series (x)
  (let* ((d 1d0)
	 (eps (make-qd-d (scale-float 1d0 -212)
			0d0 0d0 0d0))
	 (tmp x)
	 (r (sqr-qd tmp))
	 (s1 (make-qd-dd 0w0 0w0))
	 (k 0)
	 (sign 1))
    (loop while (qd-> tmp eps) do
	  (incf k)
	  (setf s1
		(if (minusp sign)
		    (sub-qd s1 (div-qd tmp (make-qd-d d 0d0 0d0 0d0)))
		    (add-qd s1 (div-qd tmp (make-qd-d d 0d0 0d0 0d0)))))
	  (incf d 2d0)
	  (setf tmp (mul-qd tmp r))
	  (setf sign (- sign)))
    s1))

;; pi =
;; 3.1415926535897932384626433832795028841971693993751058209749445923078L0
(defun test2 ()
  ;; pi/4 = 4 * arctan(1/5) - arctan(1/239)
  ;;
  ;; Arctan is computed using the Taylor series
  ;;
  ;;   arctan(x) = x - x^3/3 + x^5/5 - x^7/7
  (flet ((atan-series (x)
	   (let* ((d 1d0)
		  (eps (make-qd-d (scale-float 1d0 -212)
				  0d0 0d0 0d0))
		  (tmp x)
		  (r (sqr-qd tmp))
		  (s1 (make-qd-dd 0w0 0w0))
		  (k 0)
		  (sign 1))
	     (loop while (qd-> tmp eps) do
		   (incf k)
		   (setf s1
			 (if (minusp sign)
			     (sub-qd s1 (div-qd tmp (make-qd-d d 0d0 0d0 0d0)))
			     (add-qd s1 (div-qd tmp (make-qd-d d 0d0 0d0 0d0)))))
		   (incf d 2d0)
		   (setf tmp (mul-qd tmp r))
		   (setf sign (- sign)))
	     s1)))
    (let* ((x1 (div-qd (make-qd-dd 1w0 0w0)
		       (make-qd-dd 5w0 0w0)))
	   (s1 (atan-series x1))
	   (x2 (div-qd (make-qd-dd 1w0 0w0)
		       (make-qd-dd 239w0 0w0)))
	   (s2 (atan-series x2))
	   (p (mul-qd-d (sub-qd (mul-qd-d s1 4d0)
				s2)
			4d0)))
      (format t "est: ~/qd::qd-format/~%" p)
      (format t "tru: ~/qd::qd-format/~%" +qd-pi+)
      (format t "err: ~/qd::qd-format/~%" (sub-qd p +qd-pi+))
      (format t "bits: ~A~%" (bit-accuracy p +qd-pi+))
      p)))

(defun test3 ()
  (declare (optimize (speed 3)))
  ;; Salamin-Brent Quadratic formula for pi
  (let* ((a (make-qd-dd 1w0 0w0))
	 (b (sqrt-qd (make-qd-dd 0.5w0 0w0)))
	 (s (make-qd-dd 0.5w0 0w0))
	 (m 1d0)
	 (p (div-qd (mul-qd-d (sqr-qd a) 2d0)
		    s)))
    (declare (double-float m))
    (dotimes (k 9)
      (setf m (* 2 m))
      (let* ((a-new (mul-qd-d (add-qd a b) .5d0))
	     (b-new (sqrt-qd (mul-qd a b)))
	     (s-new (sub-qd s
			    (mul-qd-d (sub-qd (sqr-qd a-new)
					      (sqr-qd b-new))
				      m))))
	(setf a a-new)
	(setf b b-new)
	(setf s s-new)
	(setf p (div-qd (mul-qd-d (sqr-qd a) 2d0)
			s))))
    (format t "est: ~/qd::qd-format/~%" p)
    (format t "tru: ~/qd::qd-format/~%" +qd-pi+)
    (format t "err: ~/qd::qd-format/~%" (sub-qd p +qd-pi+))
    (format t "bits: ~A~%" (bit-accuracy p +qd-pi+))
    p))

(defun test4 ()
  (declare (optimize (speed 3)))
  ;; Borwein Quartic formula for pi
  (let* ((a (sub-qd (make-qd-dd 6w0 0w0)
		    (mul-qd-d (sqrt-qd (make-qd-dd 2w0 0w0))
			      4d0)))
	 (y (sub-qd (sqrt-qd (make-qd-dd 2w0 0w0))
		    (make-qd-dd 1w0 0w0)))
	 (m (make-qd-dd 2w0 0w0))
	 (p (div-qd (make-qd-dd 1w0 0w0)
		    a)))
    (declare (double-float m))
    (dotimes (k 9)
      (setf m (* 4 m))
      (let ((r (nroot-qd (sub-qd (make-qd-dd 1w0 0w0)
				 (sqr-qd (sqr-qd y)))
			 4)))
	(setf y (div-qd (sub-qd (make-qd-dd 1w0 0w0)
				r)
			(add-qd (make-qd-dd 1w0 0w0)
				r)))
	(setf a (sub-qd (mul-qd a
				(sqr-qd (sqr-qd (add-qd-d y 1d0))))
			(mul-qd-d (mul-qd y
					  (add-qd-d (add-qd y (sqr-qd y))
						    1d0))
				  m)))
	(setf p (div-qd (make-qd-dd 1w0 0w0)
			a))))
    (format t "est: ~/qd::qd-format/~%" p)
    (format t "tru: ~/qd::qd-format/~%" +qd-pi+)
    (format t "err: ~/qd::qd-format/~%" (sub-qd p +qd-pi+))
    (format t "bits: ~A~%" (bit-accuracy p +qd-pi+))
    p))

;; e =
;; 2.718281828459045235360287471352662497757247093699959574966967627724L0
(defun test5 ()
  ;; Taylor series for e
  (let ((s (make-qd-dd 2w0 0w0))
	(tmp (make-qd-dd 1w0 0w0))
	(n 1d0)
	(delta 0d0)
	(i 0))
    (loop while (qd-> tmp (make-qd-dd 1w-100 0w0)) do
	  (incf i)
	  (incf n)
	  (setf tmp (div-qd tmp
			    (make-qd-dd (float n 1w0) 0w0)))
	  (setf s (add-qd s tmp)))
    (format t "est: ~/qd::qd-format/~%" s)
    (format t "tru: ~/qd::qd-format/~%" +qd-e+)
    (format t "err: ~/qd::qd-format/~%" (sub-qd s +qd-e+))
    (format t "bits: ~A~%" (bit-accuracy s +qd-e+))
    s))

;; log(2) =
;; 0.6931471805599453094172321214581765680755001343602552541206800094934L0
(defun test6 ()
  ;; Taylor series for log 2
  ;;
  ;; -log(1-x) = x + x^2/2 + x^3/3 + x^4/4 + ...
  ;;
  ;; with x = 1/2 to get log(1/2) = -log(2)
  (let ((s (make-qd-dd .5w0 0w0))
	(tt (make-qd-dd .5w0 0w0))
	(n 1d0)
	(i 0))
    (loop while (qd-> tt (make-qd-dd 1w-100 0w0)) do
	  (incf i)
	  (incf n)
	  (setf tt (mul-qd-d tt .5d0))
	  (setf s (add-qd s
			  (div-qd tt (make-qd-dd (float n 1w0) 0w0)))))
    (format t "est: ~/qd::qd-format/~%" s)
    (format t "tru: ~/qd::qd-format/~%" +qd-log2+)
    (format t "err: ~/qd::qd-format/~%" (sub-qd s +qd-log2+))
    (format t "bits: ~A~%" (bit-accuracy s +qd-log2+))
    s))

(defun test-atan (&optional (fun #'atan-qd))
  ;; Compute atan for known values

  ;; atan(1/sqrt(3)) = pi/6
  (let* ((arg (div-qd +qd-one+ (sqrt-qd (make-qd-d 3d0))))
	 (y (div-qd (funcall fun arg) +qd-pi+))
	 (true (div-qd +qd-one+ (make-qd-d 6d0))))
    (format t "atan(1/sqrt(3))/pi = ~/qd::qd-format/~%" y)
    (format t "1/6                = ~/qd::qd-format/~%" true)
    (format t "bits               = ~A~%"
	    (bit-accuracy y true)))
  ;; atan(sqrt(3)) = %pi/3
  (let* ((arg (sqrt-qd (make-qd-d 3d0)))
	 (y (div-qd (funcall fun arg) +qd-pi+))
	 (true (div-qd +qd-one+ (make-qd-d 3d0))))
    (format t "atan(sqrt(3))/pi   = ~/qd::qd-format/~%" y)
    (format t "1/6                = ~/qd::qd-format/~%" true)
    (format t "bits               = ~A~%"
	    (bit-accuracy y true)))
  ;; atan(1) = %pi/4
  (let* ((arg +qd-one+)
	 (y (div-qd (funcall fun arg) +qd-pi+))
	 (true (div-qd +qd-one+ (make-qd-d 4d0))))
    (format t "atan(1)/pi         = ~/qd::qd-format/~%" y)
    (format t "1/4                = ~/qd::qd-format/~%" true)
    (format t "bits               = ~A~%"
	    (bit-accuracy y true))))

(defun test-sin ()
  (let* ((arg (div-qd +qd-pi+ (make-qd-d 6d0)))
	 (y (sin-qd arg))
	 (true (make-qd-d 0.5d0)))
    (format t "sin(pi/6)      = ~/qd::qd-format/~%" y)
    (format t "1/2            = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  (let* ((arg (div-qd +qd-pi+ (make-qd-d 4d0)))
	 (y (sin-qd arg))
	 (true (sqrt-qd (make-qd-d 0.5d0))))
    (format t "sin(pi/4)      = ~/qd::qd-format/~%" y)
    (format t "1/sqrt(2)      = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  (let* ((arg (div-qd +qd-pi+ (make-qd-d 3d0)))
	 (y (sin-qd arg))
	 (true (div-qd (sqrt-qd (make-qd-d 3d0)) (make-qd-d 2d0))))
    (format t "sin(pi/3)      = ~/qd::qd-format/~%" y)
    (format t "sqrt(3)/2      = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  )

(defun test-tan (&optional (f #'tan-qd))
  (let* ((arg (div-qd +qd-pi+ (make-qd-d 6d0)))
	 (y (funcall f arg))
	 (true (div-qd +qd-one+ (sqrt-qd (make-qd-d 3d0)))))
    (format t "tan(pi/6)      = ~/qd::qd-format/~%" y)
    (format t "1/sqrt(3)      = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  (let* ((arg (div-qd +qd-pi+ (make-qd-d 4d0)))
	 (y (funcall f arg))
	 (true +qd-one+))
    (format t "tan(pi/4)      = ~/qd::qd-format/~%" y)
    (format t "1              = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  (let* ((arg (div-qd +qd-pi+ (make-qd-d 3d0)))
	 (y (funcall f arg))
	 (true (sqrt-qd (make-qd-d 3d0))))
    (format t "tan(pi/3)      = ~/qd::qd-format/~%" y)
    (format t "sqrt(3)        = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  )
    
(defun test-asin ()
  (let* ((arg (make-qd-d 0.5d0))
	 (y (asin-qd arg))
	 (true (div-qd +qd-pi+ (make-qd-d 6d0))))
    (format t "asin(1/2)      = ~/qd::qd-format/~%" y)
    (format t "pi/6           = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  (let* ((arg (sqrt-qd (make-qd-d 0.5d0)))
	 (y (asin-qd arg))
	 (true (div-qd +qd-pi+ (make-qd-d 4d0))))
    (format t "asin(1/sqrt(2))= ~/qd::qd-format/~%" y)
    (format t "pi/4           = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  (let* ((arg (div-qd (sqrt-qd (make-qd-d 3d0)) (make-qd-d 2d0)))
	 (y (asin-qd arg))
	 (true (div-qd +qd-pi+ (make-qd-d 3d0))))
    (format t "asin(sqrt(3)/2)= ~/qd::qd-format/~%" y)
    (format t "pi/3           = ~/qd::qd-format/~%" true)
    (format t "bits           = ~A~%"
	    (bit-accuracy y true)))
  )
    

(defun all-tests ()
  (test2)
  (test3)
  (test4)
  (test5)
  (test6)
  (test-sin)
  (test-asin)
  (dolist (f '(atan-qd cordic-atan-qd atan-double-qd))
    (test-atan f))
  (dolist (f '(tan-qd cordic-tan-qd))
    (test-tan f))
  )
