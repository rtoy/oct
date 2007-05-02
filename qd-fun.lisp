;;; -*- Mode: lisp -*-

(in-package "QD")

#-sparc
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

(defun expm1-qd (a)
  (declare (type %quad-double a))
  ;; D(x) = exp(x) - 1
  ;;
  ;; First, write x = s*log(2) + r*k where s is an integer and |r*k| <
  ;; log(2)/2.
  ;;
  ;; Then D(x) = D(s*log(2)+r*k) = 2^s*exp(r*k) - 1
  ;;           = 2^s*(exp(r*k)-1) - 1 + 2^s
  ;;           = 2^s*D(r*k)+2^s-1
  ;; But
  ;; exp(r*k) = exp(r)^k
  ;;          = (D(r) + 1)^k
  ;;
  ;; So
  ;; D(r*k) = (D(r) + 1)^k - 1
  ;;
  ;; For small r, D(r) can be computed using the Taylor series around
  ;; zero.  To compute D(r*k) = (D(r) + 1)^k - 1, we use the binomial
  ;; theorem to expand out the power and to exactly cancel out the -1
  ;; term, which is the source of inaccuracy.
  ;;
  ;; We want to have small r so the Taylor series converges quickly,
  ;; but that means k is large, which means the binomial expansion is
  ;; long.  We need to compromise.  Let use choose k = 8.  Then |r| <
  ;; log(2)/16 = 0.0433.  For this range, the Taylor series converges
  ;; to 212 bits of accuracy with about 28 terms.
  ;;
  ;;
  (flet ((taylor (x)
	   (declare (type %quad-double x))
	   ;; Taylor series for exp(x)-1
	   ;; = x+x^2/2!+x^3/3!+x^4/4!+...
	   ;; = x*(1+x/2!+x^2/3!+x^3/4!+...)
	   (let ((sum (make-qd-d 1d0 0d0 0d0 0d0))
		 (term (make-qd-d 1d0 0d0 0d0 0d0)))
	     (dotimes (k 28)
	       (setf term (div-qd-d (mul-qd term x) (float (+ k 2) 1d0)))
	       (setf sum (add-qd sum term)))
	     (mul-qd x sum)))
	 (binom (x)
	   (declare (type %quad-double x))
	   ;; (1+x)^8-1
	   ;; = x*(8 + 28*x + 56*x^2 + 70*x^3 + 56*x^4 + 28*x^5 + 8*x^6 + x^7)
	   ;; = x (x (x (x (x (x (x (x + 8) + 28) + 56) + 70) + 56) + 28) + 8)
	   (mul-qd
	    x
	    (add-qd-d
	     (mul-qd x
		     (add-qd-d
		      (mul-qd x
			      (add-qd-d
			       (mul-qd x
				       (add-qd-d
					(mul-qd x
						(add-qd-d
						 (mul-qd x
							 (add-qd-d
							  (mul-qd x
								  (add-qd-d x 8d0))
							  28d0))
						 56d0))
					70d0))
			       56d0))
		      28d0))
	     8d0)))
	 (arg-reduce (x)
	   (declare (type %quad-double x))
	   ;; Write x = s*log(2) + r*k where s is an integer and |r*k|
	   ;; < log(2)/2, and k = 8.
	   (let* ((s (truncate (qd-0 (nint-qd (div-qd a +qd-log2+)))))
		  (r*k (sub-qd x (mul-qd-d +qd-log2+ (float s 1d0))))
		  (r (div-qd-d r*k 8d0)))
	     (values s r))))
    (multiple-value-bind (s r)
	(arg-reduce a)
      (let* ((d (taylor r))
	     (dr (binom d)))
	(add-qd-d (scale-float-qd dr s)
		  (- (scale-float 1d0 s) 1))))))
    
	 
  
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

(defun log1p-qd (x)
  (declare (type %quad-double x))
  ;; Brent gives the following duplication formula for log1p(x) =
  ;; log(1+x):
  ;;
  ;; log1p(x) = 2*log1p(x/(1+sqrt(1+x)))
  ;;
  ;; So we apply the duplication formula until x is small enough, and
  ;; then use the series
  ;;
  ;; log(1+x) = 2*sum((x/(2+x))^(2*k+1)/(2*k+1),k,0,inf)
  ;;
  (cond ((> (abs (qd-0 x)) .005d0)
	 ;; log1p(x) = 2*log1p(x/(1+sqrt(1+x)))
	 (mul-qd-d (log1p-qd (div-qd x
				     (add-d-qd 1d0
					       (sqrt-qd (add-d-qd 1d0 x)))))
		   2d0))
	(t
	 ;; Use the series
	 (let* ((term (div-qd x (add-qd-d x 2d0)))
		(mult (sqr-qd term))
		(sum term))
	   (loop for k of-type double-float from 3d0 by 2d0
	      while (> (abs (qd-0 term)) +qd-eps+)
	      do
		(setf term (mul-qd term mult))
		(setf sum (add-qd sum (div-qd-d term k))))
	   (mul-qd-d sum 2d0)))))

(defun agm-qd (x y)
  (declare (type %quad-double x y))
  (let ((diff (qd-0 (abs-qd (sub-qd x y)))))
    (cond ((< diff +qd-eps+)
	   x)
	  (t
	   (let ((a-mean (div-qd-d (add-qd x y) 2d0))
		 (g-mean (sqrt-qd (mul-qd x y))))
	     (agm-qd a-mean g-mean))))))

(defun log-agm-qd (x)
  (declare (type %quad-double x))
  ;; log(x) ~ pi/2/agm(1,4/x)*(1+O(1/x^2))
  ;;
  ;; Need to make x >= 2^(d/2) to get d bits of precision.  We use
  ;;
  ;; log(2^k*x) = k*log(2)+log(x)
  ;;
  ;; to compute log(x).  log(2^k*x) is computed using AGM.
  ;;
  (multiple-value-bind (frac exp)
      (decode-float (qd-0 x))
    (declare (ignore frac))
    (cond ((>= exp 106)
	   ;; Big enough to use AGM
	   (div-qd +qd-pi/2+
		   (agm-qd (make-qd-d 1d0 0d0 0d0 0d0)
			   (div-qd (make-qd-d 4d0 0d0 0d0 0d0)
				   x))))
	  (t
	   ;; log(x) = log(2^k*x) - k * log(2)
	   (let* ((k (- 107 exp))
		  (big-x (scale-float-qd x k)))
	     (format t "exp = ~A~%" exp)
	     (format t "k = ~A~%" k)
	     (format t "big-x = ~A~%" big-x)
	     (sub-qd (log-agm-qd big-x)
		     (mul-qd-d +qd-log2+ (float k 1d0))))))))

(defun log-agm2-qd (x)
  (declare (type %quad-double x))
  ;; log(x) ~ pi/4/agm(theta2(q^4)^2,theta3(q^4)^2)
  ;;
  ;; where q = 1/x
  ;;
  ;; Need to make x >= 2^(d/36) to get d bits of precision.  We use
  ;;
  ;; log(2^k*x) = k*log(2)+log(x)
  ;;
  ;; to compute log(x).  log(2^k*x) is computed using AGM.
  ;;
  (multiple-value-bind (frac exp)
      (decode-float (qd-0 x))
    (declare (ignore frac))
    (cond ((>= exp 6)
	   ;; Big enough to use AGM
	   (let* ((q (div-qd (make-qd-d 1d0 0d0 0d0 0d0)
			     x))
		  (q^4 (npow q 4))
		  (q^8 (mul-qd q^4 q))
		  ;; theta2(q^4) = 2*q*(1+q^8+q^24)
		  ;;             = 2*q*(1+q^8+(q^8)^3)
		  (theta2 (mul-qd-d
			   (mul-qd
			    q
			    (add-qd-d
			     (add-qd q^8
				     (npow q^8 3))
			     1d0))
			   2d0))
		  ;; theta3(q^4) = 1+2*(q^4+q^16)
		  ;;             = 1+2*(q^4+(q^4)^4)
		  (theta3 (add-qd-d
			   (mul-qd-d
			    (add-qd q^4
				    (npow q^4 4))
			    2d0)
			   1d0)))
	   (div-qd +qd-pi/4+
		   (agm-qd (sqr-qd theta2)
			   (sqr-qd theta3)))))
	  (t
	   ;; log(x) = log(2^k*x) - k * log(2)
	   (let* ((k (- 7 exp))
		  (big-x (scale-float-qd x k)))
	     (format t "exp = ~A~%" exp)
	     (format t "k = ~A~%" k)
	     (format t "big-x = ~A~%" big-x)
	     (sub-qd (log-agm2-qd big-x)
		     (mul-qd-d +qd-log2+ (float k 1d0))))))))
	     
      

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
		;;(format t "s = ~/qd::qd-format/~%" s)
		;;(format t "c = ~/qd::qd-format/~%" c)
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
		#+nil
		(progn
		  (format t "s = ~/qd::qd-format/~%" s)
		  (format t "c = ~/qd::qd-format/~%" c))
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
		#+nil
		(progn
		  (format t "s = ~/qd::qd-format/~%" s)
		  (format t "c = ~/qd::qd-format/~%" c))
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
    #+nil
    (progn
      (format t "r  = ~/qd::qd-format/~%" r)
      (format t "xx = ~/qd::qd-format/~%" xx)
      (format t "yy = ~/qd::qd-format/~%" yy))
    
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
	     (format t "xx > yy~%")
	     (dotimes (k 5)
	       (let* ((sinz (sin z))
		      (cosz (cos z))
		      (delta (/ (- yy sinz)
				cosz)))
		 (format t "sz, dz = ~A ~A~%" sinz cosz)
		 (format t "delta  = ~A~%" delta)
		 (setf z (+ z delta))
		 (format t "z = ~A~%" z))))
	    (t
	     (dotimes (k 20)
	       (let ((sinz (sin z))
		     (cosz (cos z)))
		 (format t "sz, dz = ~A ~A~%" sinz cosz)
		 
		 (setf z (- z (/ (- xx cosz)
				 sinz)))
		 (format t "z = ~A~%" z)))))
      z)))

(defvar *table*)
(defvar *ttable*)
(defvar *cordic-scale*)

#+nil
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

(defun setup-cordic ()
  (let ((table (make-array 34))
	(ttable (make-array 34)))
    (setf (aref table 0) 4d0)
    (setf (aref table 1) 2d0)
    (setf (aref table 2) 1d0)
    (setf (aref ttable 0) (atan 4d0))
    (setf (aref ttable 1) (atan 2d0))
    (setf (aref ttable 2) (/ pi 4))
    (loop for k from 3 below 34 do
	 (setf (aref table k) (* 0.5d0 (aref table (1- k))))
	 (setf (aref ttable k) (atan (aref table k))))
    (setf *table* table)
    (setf *ttable* ttable)))

(defun setup-cordic ()
  (let ((table (make-array 34))
	(ttable (make-array 34))
	(scale 1d0))
    (loop for k from 0 below 34 do
	 (setf (aref table k) (scale-float 1d0 (- 2 k)))
	 (setf (aref ttable k) (atan (aref table k)))
	 (setf scale (* scale (cos (aref ttable k)))))
    (setf *table* table)
    (setf *ttable* ttable)
    (setf *cordic-scale* scale)))


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

(defun cordic-vec (z)
  (let ((x 1d0)
	(y 0d0)
	(scale 1d0))
    (dotimes (k 12 (length *table*))
      (setf scale (* scale (cos (aref *ttable* k))))
      (cond ((minusp z)
	     (psetq x (+ x (* y (aref *table* k)))
		    y (- y (* x (aref *table* k))))
	     (incf z (aref *ttable* k)))
	    (t
	     (psetq x (- x (* y (aref *table* k)))
		    y (+ y (* x (aref *table* k))))
	     (decf z (aref *ttable* k)))
	    ))
    (values x y z scale)))

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

(defun tan-d (r)
  (multiple-value-bind (x y z)
      (cordic-vec r)
    (setf x (* x *cordic-scale*))
    (setf y (* y *cordic-scale*))
    (format t "x = ~A~%" x)
    (format t "y = ~A~%" y)
    (format t "z = ~A~%" z)
    ;; Need to finish of the rotation
    (let ((st (sin z))
	  (ct (cos z)))
      (format t "st, ct = ~A ~A~%" st ct)
      (psetq x (- (* x ct) (* y st))
	     y (+ (* y ct) (* x st)))
      (format t "x = ~A~%" x)
      (format t "y = ~A~%" y)
      (/ y x)
      )))

(defun sin-d (r)
  (declare (type double-float r))
  (multiple-value-bind (x y z s)
      (cordic-vec r)
    
    ;; Need to finish the rotation
    (let ((st (sin z))
	  (ct (cos z)))
      (psetq x (- (* x ct) (* y st))
	     y (+ (* y ct) (* x st)))
      (* s y))))

;; This is the basic CORDIC rotation.  Based on code from
;; http://www.voidware.com/cordic.htm and
;; http://www.dspcsp.com/progs/cordic.c.txt.
;;
;; The only difference between this version and the typical CORDIC
;; implementation is that the first 3 rotations are all by pi/4.  This
;; makes sense.  If the angle is greater than pi/4, the rotations will
;; reduce it to at most pi/4.  If the angle is less than pi/4, the 3
;; rotations by pi/4 will cause us to end back at the same place.
;; (Should we try to be smarter?)
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

(defun cordic-atan2-qd (y x)
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
      (add-qd z sum))))

(defun atan-qd (y)
  (declare (type %quad-double y))
  (atan2-qd y (%make-qd-d 1d0 0d0 0d0 0d0)))

(defun asin-qd (a)
  (declare (type %quad-double a))
  (atan2-qd a (sqrt-qd (sub-qd (%make-qd-d 1d0 0d0 0d0 0d0)
			       (sqr-qd a)))))

(defun acos-qd (a)
  (declare (type %quad-double a))
  (atan2-qd (sqrt-qd (sub-qd (%make-qd-d 1d0 0d0 0d0 0d0)
			     (sqr-qd a)))
	    a))
  

(defun tan-qd (r)
  (declare (type %quad-double r))
  (multiple-value-bind (z x y)
      (cordic-vec-qd r)
    ;; Need to finish the rotation
    (multiple-value-bind (st ct)
	(sincos-taylor z)
      (psetq x (sub-qd (mul-qd x ct) (mul-qd y st))
	     y (add-qd (mul-qd y ct) (mul-qd x st)))
      (div-qd y x))))

(defun cordic-vec-qd (z)
  (declare (type %quad-double z)
	   (optimize (speed 3)))
  (let* ((x (%make-qd-d 1d0 0d0 0d0 0d0))
	 (y (%make-qd-d 0d0 0d0 0d0 0d0))
	 (zero (%make-qd-d 0d0 0d0 0d0 0d0))
	 )
    (declare (type %quad-double zero x y))
    (dotimes (k 30 (length +atan-table+))
      (declare (fixnum k)
	       (inline mul-qd-d sub-qd add-qd))
      (cond ((qd-> z zero)
	     (psetq x (sub-qd x (mul-qd-d y (aref +atan-power-table+ k)))
		    y (add-qd y (mul-qd-d x (aref +atan-power-table+ k))))
	     (setf z (sub-qd z (aref +atan-table+ k))))
	    (t
	     (psetq x (add-qd x (mul-qd-d y (aref +atan-power-table+ k)))
		    y (sub-qd y (mul-qd-d x (aref +atan-power-table+ k))))
	     (setf z (add-qd z (aref +atan-table+ k))))))
    (values z x y)))
  
(defun cordic-sin-qd (r)
  (declare (type %quad-double r))
  (multiple-value-bind (z x y)
      (cordic-vec-qd r)
    #+nil
    (progn
      (format t "~&x = ~/qd::qd-format/~%" x)
      (format t "~&y = ~/qd::qd-format/~%" y)
      (format t "~&z = ~/qd::qd-format/~%" z)
      (format t "~&s = ~/qd::qd-format/~%" s))
    ;; Need to finish the rotation
    (multiple-value-bind (st ct)
	(sincos-taylor z)
      #+nil
      (progn
	(format t "~&st = ~/qd::qd-format/~%" st)
	(format t "~&ct = ~/qd::qd-format/~%" ct)
	(format t "~&y  = ~/qd::qd-format/~%" (mul-qd +cordic-scale+ y)))

      (psetq x (sub-qd (mul-qd x ct) (mul-qd y st))
	     y (add-qd (mul-qd y ct) (mul-qd x st)))
      (mul-qd +cordic-scale+ y))))

(defun sinh-qd (a)
  (declare (type %quad-double a))
  ;; Hart et al. suggests sinh(x) = 1/2*(D(x) + D(x)/(D(x)+1))
  ;; where D(x) = exp(x) - 1.
  (let ((d (expm1-qd a)))
    (scale-float-qd (add-qd d
			    (div-qd d (add-qd-d d 1d0)))
		    -1)))

(defun cosh-qd (a)
  (declare (type %quad-double a))
  ;; cosh(x) = 1/2*(exp(x)+exp(-x))
  (let ((e (exp-qd a)))
    (scale-float-qd (add-qd e (div-qd (make-qd-d 1d0 0d0 0d0 0d0) e))
		    -1)))

(defun tanh-qd (a)
  (declare (type %quad-double a))
  ;; Hart et al. suggests tanh(x) = D(2*x)/(2+D(2*x))
  (let* ((a2 (mul-qd-d a 2d0))
	 (d (expm1-qd a2)))
    (div-qd d (add-qd-d d 2d0))))

(defun asinh-qd (a)
  (declare (type %quad-double a))
  ;; asinh(x) = log(x + sqrt(1+x^2))
  (log-qd (add-qd a
		  (sqrt-qd (add-qd-d (sqr-qd a)
				     1d0)))))

(defun acosh-qd (a)
  (declare (type %quad-double a))
  ;; acosh(x) = log(x + sqrt(x^2-1))
  (log-qd (add-qd a
		  (sqrt-qd (sub-qd-d (sqr-qd a)
				     1d0)))))

(defun atanh-qd (a)
  (declare (type %quad-double a))
  ;; atanh(x) = 1/2*log((1+x)/(1-x))
  ;;          = 1/2*log(1+(2*x)/(1-x))
  ;; This latter expression works better for small x
  (scale-float-qd (log-qd (div-qd (add-d-qd 1d0 a)
				  (sub-d-qd 1d0 a)))
		  -1))
  
  
