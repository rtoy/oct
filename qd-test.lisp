;;; -*- Mode: lisp -*-

(in-package "QD")

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
      (qd-output-aux p)
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
    (qd-output-aux p)
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
    (qd-output-aux p)
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
    (qd-output-aux s)
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
    (qd-output-aux s)
    s))
