;;; -*- Mode: lisp -*-

(in-package "QD")

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


(defun nint-qd (a)
  (let ((x0 (fround (qd-0 a)))
	(x1 0d0)
	(x2 0d0)
	(x3 0d0))
    (cond ((= x0 (qd-0 a))
	   ;; First double is already an integer
	   (setf x1 (fround (qd-1 a)))
	   (cond ((= x1 (qd-1 a))
		  ;; Second is an integer
		  (setf x2 (fround (qd-2 a)))
		  (cond ((= x2 (qd-2 a))
			 ;; Third is an integer
			 (setf x3 (fround (qd-3 a))))
			(t
			 (when (and (zerop (abs (- x2 (qd-2 a))))
				    (minusp (qd-3 a)))
			   (decf x2)))))
		 (t
		  (when (and (zerop (abs (- x1 (qd-1 a))))
			     (minusp (qd-2 a)))
		    (decf x1)))))
	  (t
	   (when (and (zerop (abs (- x0 (qd-0 a))))
		      (minusp (qd-1 a)))
	     (decf x0))))
    (multiple-value-bind (s0 s1 s2 s3)
	(renorm-4 x0 x1 x2 x3)
      (make-qd-d s0 s1 s2 s3))))

(defun ffloor-qd (a)
  (let ((x0 (ffloor (qd-0 a)))
	(x1 0d0)
	(x2 0d0)
	(x3 0d0))
    (cond ((= x0 (qd-0 a))
	   (setf x1 (ffloor (qd-1 a)))
	   (when (= x1 (qd-1 a))
	     (setf x2 (ffloor (qd-2 a)))
	     (when (= x2 (qd-2 a))
	       (setf x3 (ffloor (qd-3 a)))))
	   (make-qd-d x0 x1 x2 x3))
	  (t
	   (%make-qd-d x0 x1 x2 x3)))))
	
			 
(defun exp-qd (a)
  ;; Strategy:  Reduce the size of x by noting that
  ;;
  ;; exp(k*r+m) = exp(m) * exp(r)^k
  ;;
  ;; Thus, by choosing m to be a multiple of log(2) closest to x, we
  ;; can make |kr| < log(2)/2 = 0.3466.  Now we can set k = 256, so
  ;; that |r| <= 0.00136.
  ;;
  ;; Then
  ;;
  ;; exp(x) = exp(k*r+s*log(2)) = 2^s*(exp(r))^256
  ;;
  ;; We can use Taylor series to evaluate exp(r).

  (when (< (qd-0 a) -709)
    (return-from exp-qd (make-qd-dd 0w0 0w0)))

  (when (> (qd-0 a) 709)
    (error "exp-qd overflow"))

  (when (= a 0)
    (return-from exp-qd (make-qd-dd 1w0 0w0)))

  (let* ((k 256)
	 (z (truncate (qd-0 (nint-qd (div-qd a +qd-log2+)))))
	 (r1 (sub-qd a (mul-qd-d +qd-log2+ (float z 1d0))))
	 (r (div-qd (sub-qd a (mul-qd-d +qd-log2+ (float z 1d0)))
		    (make-qd-d (float k 1d0) 0d0 0d0 0d0)))
	 (p (div-qd (sqr-qd r) (make-qd-d 2d0 0d0 0d0 0d0)))
	 (s (add-qd-d (add-qd r p) 1d0))
	 (m 2d0))
    (loop
       (incf m)
       (setf p (mul-qd p r))
       (setf p (div-qd p (make-qd-d (float m 1d0) 0d0 0d0 0d0)))
       (setf s (add-qd s p))
       (unless (> (abs (qd-0 p)) (expt 2d0 -200))
	 (return)))

    (setf r (npow s k))
    (setf r (scale-float-qd r z))
    r))

(defun log-qd (a)
  ;; The Taylor series for log converges rather slowly.  Hence, this
  ;; routine tries to determine the root of the function
  ;;
  ;; f(x) = exp(x) - a
  ;;
  ;; using Newton iteration.  The iteration is
  ;;
  ;; x' = x - f(x) / f'(x)
  ;;    = x - (1 - a * exp(-x))
  ;;    = x + a * exp(-x) - 1
  ;;
  ;; Two iterations are needed.
  (when (= a 1)
    (return-from log-qd (make-qd-dd 0w0 0w0)))

  (when (minusp (qd-0 a))
    (error "log of negative"))

  (let ((x (make-qd-d (log (qd-0 a)) 0d0 0d0 0d0))
	(one (make-qd-dd 1w0 0w0)))
    (setf x (sub-qd-d (add-qd x (mul-qd a (exp-qd (neg-qd x))))
		    1d0))
    (setf x (sub-qd-d (add-qd x (mul-qd a (exp-qd (neg-qd x))))
		    1d0))
    (setf x (sub-qd-d (add-qd x (mul-qd a (exp-qd (neg-qd x))))
		    1d0))
    x))

;; sin(a) and cos(a) using Taylor series
;;
;; Assumes |a| <= pi/2048
(defun sincos-taylor (a)
  (declare (type %quad-double a))
  (let ((thresh (* +qd-eps+ (abs (qd-0 a)))))
    (when (= a 0)
      (return-from sincos-taylor
	(values (make-qd-dd 0w0 0w0)
		(make-qd-dd 1w0 0w0))))
    (let* ((x (neg-qd (sqr-qd a)))
	   (s a)
	   (p a)
	   (m 1d0))
      (loop
	 (setf p (mul-qd p x))
	 (incf m 2)
	 (setf p (div-qd-d p (* m (1- m))))
	 (setf s (add-qd s p))
	 ;;(format t "p = ~A~%" (qd-0 p))
	 (when (< (abs (qd-0 p)) thresh)
	   (return)))
      ;; cos(c) = sqrt(1-sin(c)^2).  This seems to work ok, even
      ;; though I would have expected some round-off errors in
      ;; computing this.  sqrt(1-x^2) is normally better computed as
      ;; sqrt(1-x)*sqrt(1+x) for small x.
      (values s (sqrt-qd (add-qd-d (neg-qd (sqr-qd s)) 1d0))))))

(defun drem-qd (a b)
  (declare (type %quad-double a b))
  (let ((n (nint-qd (div-qd a b))))
    (sub-qd a (mul-qd n b))))

(defun divrem-qd (a b)
  (declare (type %quad-double a b))
  (let ((n (nint-qd (div-qd a b))))
    (values n (sub-qd a (mul-qd n b)))))
  
(defun sin-qd (a)
  (declare (type %quad-double a))
  ;; To compute sin(x), choose integers a, b so that
  ;;
  ;; x = s + a * (pi/2) + b*(pi/1024)
  ;;
  ;; with |x| <= pi/2048.  Using a precomputed table of sin(k*pi/1024)
  ;; and cos(k*pi/1024), we can compute sin(x) from sin(s) and cos(s).
  ;;
  ;; sin(x) = sin(s+k*(pi/1024) + j*pi/2)
  ;;        = sin(s+k*(pi/1024))*cos(j*pi/2)
  ;;            + cos(s+k*(pi/1024))*sin(j*pi/2)
  ;;
  ;; sin(s+k*pi/1024) = sin(s)*cos(k*pi/1024)
  ;;                     + cos(s)*sin(k*pi/1024)
  ;;
  ;; cos(s+k*pi/1024) = cos(s)*cos(k*pi/1024)
  ;;                     - sin(s)*sin(k*pi/1024)
  (when (= a 0)
    (return-from sin-qd (make-qd-dd 0w0 0w0)))

  ;; Reduce modulo 2*pi
  (let ((r (drem-qd a +qd-2pi+)))
    ;; Now reduce by pi/2 and then by pi/1024 so that we obtain
    ;; numbers a, b, t
    (multiple-value-bind (j tmp)
	(divrem-qd r +qd-pi/2+)
      (let* ((j (truncate (qd-0 j)))
	     (abs-j (abs j)))
	(multiple-value-bind (k tmp)
	    (divrem-qd tmp +qd-pi/1024+)
	  (let* ((k (truncate (qd-0 k)))
		 (abs-k (abs k)))
	    (assert (<= abs-j 2))
	    (assert (<= abs-k 256))
	    ;; Compute sin(s) and cos(s)
	    (multiple-value-bind (sin-t cos-t)
		(sincos-taylor tmp)
	      (multiple-value-bind (s c)
		  (cond ((zerop abs-k)
			 (values sin-t cos-t))
			(t
			 ;; Compute sin(s+k*pi/1024), cos(s+k*pi/1024)
			 (let ((u (aref +qd-cos-table+ (1- abs-k)))
			       (v (aref +qd-sin-table+ (1- abs-k))))
			   (cond ((plusp k)
				  ;; sin(s) * cos(k*pi/1024)
				  ;; + cos(s) * sin(k*pi/1024)
				  ;;
				  ;; cos(s) * cos(k*pi/1024)
				  ;; - sin(s) * sin(k*pi/1024)
				  (values (add-qd (mul-qd u sin-t)
						  (mul-qd v cos-t))
					  (sub-qd (mul-qd u cos-t)
						  (mul-qd v sin-t))))
				 (t
				  ;; sin(s) * cos(k*pi/1024)
				  ;; - cos(s) * sin(|k|*pi/1024)
				  ;;
				  ;; cos(s) * cos(k*pi/1024)
				  ;; + sin(s) * sin(|k|*pi/1024)
				  (values (sub-qd (mul-qd u sin-t)
						  (mul-qd v cos-t))
					  (add-qd (mul-qd u cos-t)
						  (mul-qd v sin-t))))))))
		(format t "s = ~/qd::qd-format/~%" s)
		(format t "c = ~/qd::qd-format/~%" c)
		;; sin(x) =  sin(s+k*pi/1024) * cos(j*pi/2)
		;;         + cos(s+k*pi/1024) * sin(j*pi/2)
		(cond ((zerop abs-j)
		       ;; cos(j*pi/2) = 1, sin(j*pi/2) = 0
		       s)
		      ((= j 1)
		       ;; cos(j*pi/2) = 0, sin(j*pi/2) = 1
		       c)
		      ((= j -1)
		       ;; cos(j*pi/2) = 0, sin(j*pi/2) = -1
		       (neg-qd c))
		      (t
		       ;; cos(j*pi/2) = -1, sin(j*pi/2) = 0
		       (neg-qd s)))))))))))
		     
(defun cos-qd (a)
  ;; Just like sin-qd, but for cos.
  (declare (type %quad-double a))
  ;; To compute sin(x), choose integers a, b so that
  ;;
  ;; x = s + a * (pi/2) + b*(pi/1024)
  ;;
  ;; with |x| <= pi/2048.  Using a precomputed table of sin(k*pi/1024)
  ;; and cos(k*pi/1024), we can compute sin(x) from sin(s) and cos(s).
  ;;
  ;; sin(x) = sin(s+k*(pi/1024) + j*pi/2)
  ;;        = sin(s+k*(pi/1024))*cos(j*pi/2)
  ;;            + cos(s+k*(pi/1024))*sin(j*pi/2)
  ;;
  ;; sin(s+k*pi/1024) = sin(s)*cos(k*pi/1024)
  ;;                     + cos(s)*sin(k*pi/1024)
  ;;
  ;; cos(s+k*pi/1024) = cos(s)*cos(k*pi/1024)
  ;;                     - sin(s)*sin(k*pi/1024)
  (when (= a 0)
    (return-from cos-qd (make-qd-dd 1w0 0w0)))

  ;; Reduce modulo 2*pi
  (let ((r (drem-qd a +qd-2pi+)))
    ;; Now reduce by pi/2 and then by pi/1024 so that we obtain
    ;; numbers a, b, t
    (multiple-value-bind (j tmp)
	(divrem-qd r +qd-pi/2+)
      (let* ((j (truncate (qd-0 j)))
	     (abs-j (abs j)))
	(multiple-value-bind (k tmp)
	    (divrem-qd tmp +qd-pi/1024+)
	  (let* ((k (truncate (qd-0 k)))
		 (abs-k (abs k)))
	    (assert (<= abs-j 2))
	    (assert (<= abs-k 256))
	    ;; Compute sin(s) and cos(s)
	    (multiple-value-bind (sin-t cos-t)
		(sincos-taylor tmp)
	      (multiple-value-bind (s c)
		  (cond ((zerop abs-k)
			 (values sin-t cos-t))
			(t
			 ;; Compute sin(s+k*pi/1024), cos(s+k*pi/1024)
			 (let ((u (aref +qd-cos-table+ (1- abs-k)))
			       (v (aref +qd-sin-table+ (1- abs-k))))
			   (cond ((plusp k)
				  ;; sin(s) * cos(k*pi/1024)
				  ;; + cos(s) * sin(k*pi/1024)
				  ;;
				  ;; cos(s) * cos(k*pi/1024)
				  ;; - sin(s) * sin(k*pi/1024)
				  (values (add-qd (mul-qd u sin-t)
						  (mul-qd v cos-t))
					  (sub-qd (mul-qd u cos-t)
						  (mul-qd v sin-t))))
				 (t
				  ;; sin(s) * cos(k*pi/1024)
				  ;; - cos(s) * sin(|k|*pi/1024)
				  ;;
				  ;; cos(s) * cos(k*pi/1024)
				  ;; + sin(s) * sin(|k|*pi/1024)
				  (values (sub-qd (mul-qd u sin-t)
						  (mul-qd v cos-t))
					  (add-qd (mul-qd u cos-t)
						  (mul-qd v sin-t))))))))
		(format t "s = ~/qd::qd-format/~%" s)
		(format t "c = ~/qd::qd-format/~%" c)
		;; sin(x) =  sin(s+k*pi/1024) * cos(j*pi/2)
		;;         + cos(s+k*pi/1024) * sin(j*pi/2)
		(cond ((zerop abs-j)
		       ;; cos(j*pi/2) = 1, sin(j*pi/2) = 0
		       c)
		      ((= j 1)
		       ;; cos(j*pi/2) = 0, sin(j*pi/2) = 1
		       (neg-qd s))
		      ((= j -1)
		       ;; cos(j*pi/2) = 0, sin(j*pi/2) = -1
		       s)
		      (t
		       ;; cos(j*pi/2) = -1, sin(j*pi/2) = 0
		       (neg-qd c)))))))))))

;; Compute sin and cos of a
(defun sincos-qd (a)
  (declare (type %quad-double a))
  (when (= a 0)
    (return-from sincos-qd
      (values (make-qd-dd 0w0 0w0)
	      (make-qd-dd 1w0 0w0))))

  ;; Reduce modulo 2*pi
  (let ((r (drem-qd a +qd-2pi+)))
    ;; Now reduce by pi/2 and then by pi/1024 so that we obtain
    ;; numbers a, b, t
    (multiple-value-bind (j tmp)
	(divrem-qd r +qd-pi/2+)
      (let* ((j (truncate (qd-0 j)))
	     (abs-j (abs j)))
	(multiple-value-bind (k tmp)
	    (divrem-qd tmp +qd-pi/1024+)
	  (let* ((k (truncate (qd-0 k)))
		 (abs-k (abs k)))
	    (assert (<= abs-j 2))
	    (assert (<= abs-k 256))
	    ;; Compute sin(s) and cos(s)
	    (multiple-value-bind (sin-t cos-t)
		(sincos-taylor tmp)
	      (multiple-value-bind (s c)
		  (cond ((zerop abs-k)
			 (values sin-t cos-t))
			(t
			 ;; Compute sin(s+k*pi/1024), cos(s+k*pi/1024)
			 (let ((u (aref +qd-cos-table+ (1- abs-k)))
			       (v (aref +qd-sin-table+ (1- abs-k))))
			   (cond ((plusp k)
				  ;; sin(s) * cos(k*pi/1024)
				  ;; + cos(s) * sin(k*pi/1024)
				  ;;
				  ;; cos(s) * cos(k*pi/1024)
				  ;; - sin(s) * sin(k*pi/1024)
				  (values (add-qd (mul-qd u sin-t)
						  (mul-qd v cos-t))
					  (sub-qd (mul-qd u cos-t)
						  (mul-qd v sin-t))))
				 (t
				  ;; sin(s) * cos(k*pi/1024)
				  ;; - cos(s) * sin(|k|*pi/1024)
				  ;;
				  ;; cos(s) * cos(k*pi/1024)
				  ;; + sin(s) * sin(|k|*pi/1024)
				  (values (sub-qd (mul-qd u sin-t)
						  (mul-qd v cos-t))
					  (add-qd (mul-qd u cos-t)
						  (mul-qd v sin-t))))))))
		(format t "s = ~/qd::qd-format/~%" s)
		(format t "c = ~/qd::qd-format/~%" c)
		;; sin(x) =  sin(s+k*pi/1024) * cos(j*pi/2)
		;;         + cos(s+k*pi/1024) * sin(j*pi/2)
		(cond ((zerop abs-j)
		       ;; cos(j*pi/2) = 1, sin(j*pi/2) = 0
		       (values s c))
		      ((= j 1)
		       ;; cos(j*pi/2) = 0, sin(j*pi/2) = 1
		       (values c (neg-qd s)))
		      ((= j -1)
		       ;; cos(j*pi/2) = 0, sin(j*pi/2) = -1
		       (values (neg-qd c) s))
		      (t
		       ;; cos(j*pi/2) = -1, sin(j*pi/2) = 0
		       (values (neg-qd s)
			       (neg-qd c))))))))))))

  
(defun atan2-qd (y x)
  (declare (type %quad-double y x))
  ;; Instead of using Taylor series to compute atan, we instead use
  ;; Newton's iteration to solve the equation
  ;;
  ;;   sin(z) = y/r or cos(z) = x/r
  ;;
  ;; where r = sqrt(x^2+y^2)
  ;;
  ;; The iteration is
  ;;
  ;;   z' = z + (y - sin(z))/cos(z)       (for sin)
  ;;   z' = z + (x - cos(z))/sin(z)       (for cos)
  ;;
  ;; Here, x and y are normalized so that x^2 + y^2 = 1.
  ;;
  ;; If |x| > |y|, then the first iteration is used since the
  ;; denominator is larger.  Otherwise the second is used.
  (cond ((= x 0)
	 (cond ((= y 0)
		(error "atan2(0,0)"))
	       (t
		(return-from atan2-qd
		  (cond ((qd-> y (make-qd-dd 0w0 0w0))
			 +qd-pi/2+)
			(t
			 (neg-qd +qd-pi/2+)))))))
	((= y 0)
	 (return-from atan2-qd
	   (cond ((qd-> x (make-qd-dd 0w0 0w0))
		  (make-qd-dd 0w0 0w0))
		 (t
		  +qd-pi+)))))

  (when (= x y)
    (return-from atan2-qd
      (if (qd-> y (make-qd-dd 0w0 0w0))
	  +qd-pi/4+
	  +qd-3pi/4+)))

  (when (= x (neg-qd y))
    (return-from atan2-qd
      (if (qd-> y (make-qd-dd 0w0 0w0))
	  +qd-3pi/4+
	  (neg-qd +qd-pi/4+))))

  (let* ((r (sqrt-qd (add-qd (sqr-qd x)
			     (sqr-qd y))))
	 (xx (div-qd x r))
	 (yy (div-qd y r)))
    (format t "r  = ~/qd::qd-format/~%" r)
    (format t "xx = ~/qd::qd-format/~%" xx)
    (format t "yy = ~/qd::qd-format/~%" yy)
    
    ;; Compute double-precision approximation to atan
    (let ((z (%make-qd-d (atan (qd-0 y) (qd-0 x)) 0d0 0d0 0d0))
	  (sinz (make-qd-dd 0w0 0w0))
	  (cosz (make-qd-dd 0w0 0w0)))
      (cond ((qd-> xx yy)
	     ;; Newton iteration  z' = z + (y - sin(z))/cos(z)
	     (dotimes (k 3)
	       (multiple-value-setq (sinz cosz) (sincos-qd z))
	       (setf z (add-qd z (div-qd (sub-qd yy sinz)
					 cosz)))))
	    (t
	     ;; Newton iteration z' = z - (x - cos(z))/sin(z)
	     (dotimes (k 3)
	       (multiple-value-setq (sinz cosz) (sincos-qd z))
	       (format t "sinz ~/qd::qd-format/~%" sinz)
	       (format t "cosz ~/qd::qd-format/~%" cosz)
	       (setf z (sub-qd z (div-qd (sub-qd yy cosz)
					 sinz))))))
      z)))

(defun atan-d (y x)
  (let* ((r (abs (complex x y)))
	 (xx (/ x r))
	 (yy (/ y r)))
    (let ((z (atan (float y 1f0) (float x 1f0)))
	  (sinz 0d0)
	  (cosz 0d0))
      (format t "z = ~A~%" z)
      (cond ((> xx yy)
	     (dotimes (k 3)
	       (let ((sinz (sin z))
		     (cosz (cos z)))
		 (setf z (+ z (/ (- y sinz)
				 cosz))))))
	    (t
	     (dotimes (k 20)
	       (let ((sinz (sin z))
		     (cosz (cos z)))
		 (format t "sz, dz = ~A ~A~%" sinz cosz)
		 
		 (setf z (- z (/ (- x cosz)
				 sinz)))
		 (format t "z = ~A~%" z)))))
      z)))

(defvar *table*)
(defvar *ttable*)

(defun setup-cordic ()
  (let ((table (make-array 34))
	(ttable (make-array 34)))
    (setf (aref table 0) 1d0)
    (setf (aref table 1) 1d0)
    (setf (aref table 2) 1d0)
    (setf (aref ttable 0) (/ pi 4))
    (setf (aref ttable 1) (/ pi 4))
    (setf (aref ttable 2) (/ pi 4))
    (loop for k from 3 below 34 do
	 (setf (aref table k) (* 0.5d0 (aref table (1- k))))
	 (setf (aref ttable k) (atan (aref table k))))
    (setf *table* table)
    (setf *ttable* ttable)))


(defun cordic-rot (x y)
  (let ((z 0))
    (dotimes (k (length *table*))
      (cond ((plusp y)
	     (psetq x (+ x (* y (aref *table* k)))
		    y (- y (* x (aref *table* k))))
	     (incf z (aref *ttable* k)))
	    (t
	     (psetq x (- x (* y (aref *table* k)))
		    y (+ y (* x (aref *table* k))))
	     (decf z (aref *ttable* k)))
	    ))
    (values z x y)))

(defun atan2-d (y x)
  (multiple-value-bind (z dx dy)
      (cordic-rot x y)
    (let ((theta (/ dy dx)))
      (format t "theta = ~A~%" theta)
      (let ((corr (+ theta
		     (- (/ (expt theta 3)
			   3))
		     (/ (expt theta 5)
			5))))
	(format t "corr = ~A~%" corr)
	(+ z corr)))))

;; This is the basic CORDIC rotation.  Based on code from
;; http://www.voidware.com/cordic.htm.  The only difference between
;; this version and the typical CORDIC implementation is that the
;; first 3 rotations are all by pi/4.  This makes sense.  If the angle
;; is greater than pi/4, the rotations will reduce it to at most pi/4.
;; If the angle is less than pi/4, the 3 rotations by pi/4 will cause
;; us to end back at the same place.  (Should we try to be smarter?)
(defun cordic-rot-qd (x y)
  (declare (type %quad-double y x)
	   (optimize (speed 3)))
  (let* ((zero (%make-qd-d 0d0 0d0 0d0 0d0))
	 (z zero))
    (declare (type %quad-double zero z))
    (dotimes (k (length +atan-table+))
      (declare (fixnum k))
      (cond ((qd-> y zero)
	     (psetq x (add-qd x (mul-qd-d y (aref +atan-power-table+ k)))
		    y (sub-qd y (mul-qd-d x (aref +atan-power-table+ k))))
	     (setf z (add-qd z (aref +atan-table+ k))))
	    (t
	     (psetq x (sub-qd x (mul-qd-d y (aref +atan-power-table+ k)))
		    y (add-qd y (mul-qd-d x (aref +atan-power-table+ k))))
	     (setf z (sub-qd z (aref +atan-table+ k))))))
    (values z x y)))
  
(defun atan2-qd (y x)
  (declare (type %quad-double y x))
  ;; Use the CORDIC rotation to get us to a small angle.  Then use the
  ;; Taylor series for atan to finish the computation.
  (multiple-value-bind (z dx dy)
      (cordic-rot-qd x y)
    ;; Use Taylor series to finish off the computation
    (let* ((arg (div-qd dy dx))
	   (sq (neg-qd (sqr-qd arg)))
	   (sum (make-qd-d 1d0 0d0 0d0 0d0)))
      ;; atan(x) = x - x^3/3 + x^5/5 - ...
      ;;         = x*(1-x^2/3+x^4/5-x^6/7+...)
      (do ((k 3d0 (+ k 2d0))
	   (term sq))
	  ((< (abs (qd-0 term)) +qd-eps+))
	(setf sum (add-qd sum (div-qd-d term k)))
	(setf term (mul-qd term sq)))
      (setf sum (mul-qd arg sum))
      (values (add-qd z sum) z sum))))
