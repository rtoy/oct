(in-package "QD")

(defconstant +pi+
  (make-instance 'qd-real :value qdi:+qd-pi+))

(defun read-qd-real-or-complex (stream)
  (let ((c (peek-char t stream)))
    (cond ((char= c #\()
	   ;; Read a QD complex
	   (read-char stream)		; Skip the paren
	   (let ((real (read-qd stream))
		 (imag (read-qd stream)))
	     (unless (char= (peek-char t stream) #\))
	       (error "Illegal QD-COMPLEX number format"))
	     ;; Read closing paren
	     (read-char stream)
	     (make-instance 'qd-complex :real real :imag imag)))
	  (t
	   (make-instance 'qd-real :value (read-qd stream))))))
	
(defun qd-class-reader (stream subchar arg)
  (declare (ignore subchar arg))
  (read-qd-real-or-complex stream))

;; Yow!  We redefine the #q reader that is in qd-io.lisp to read in
;; and make a real qd-real float, instead of the hackish
;; %qd-real.
(set-dispatch-macro-character #\# #\Q #'qd-class-reader)

(defmethod add1 ((a number))
  (cl::1+ a))

(defmethod add1 ((a qd-real))
  (make-instance 'qd-real :value (add-qd-d (qd-value a) 1d0)))

(defmethod sub1 ((a number))
  (cl::1- a))

(defmethod sub1 ((a qd-real))
  (make-instance 'qd-real :value (sub-qd-d (qd-value a) 1d0)))

(declaim (inline 1+ 1-))

(defun 1+ (x)
  (add1 x))

(defun 1- (x)
  (sub1 x))

(defmethod two-arg-+ ((a qd-real) (b qd-real))
  (make-instance 'qd-real :value (add-qd (qd-value a) (qd-value b))))

(defmethod two-arg-+ ((a qd-real) (b real))
  (make-instance 'qd-real :value (add-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-+ ((a real) (b qd-real))
  (make-instance 'qd-real :value (add-d-qd (float a 1d0) (qd-value b))))

(defmethod two-arg-+ ((a number) (b number))
  (cl:+ a b))

(defun + (&rest args)
  (if (null args)
      0
      (do ((args (cdr args) (cdr args))
	   (res (car args)
		(two-arg-+ res (car args))))
	  ((null args) res))))

(defmethod two-arg-- ((a qd-real) (b qd-real))
  (make-instance 'qd-real :value (sub-qd (qd-value a) (qd-value b))))

(defmethod two-arg-- ((a qd-real) (b real))
  (make-instance 'qd-real :value (sub-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-- ((a real) (b qd-real))
  (make-instance 'qd-real :value (sub-d-qd (float a 1d0) (qd-value b))))

(defmethod two-arg-- ((a number) (b number))
  (cl:- a b))

(defmethod unary-minus ((a number))
  (cl:- a))

(defmethod unary-minus ((a qd-real))
  (make-instance 'qd-real :value (neg-qd (qd-value a))))

(defun - (number &rest more-numbers)
  (if more-numbers
      (do ((nlist more-numbers (cdr nlist))
	   (result number))
	  ((atom nlist) result)
         (declare (list nlist))
	 (setq result (two-arg-- result (car nlist))))
      (unary-minus number)))


(defmethod two-arg-* ((a qd-real) (b qd-real))
  (make-instance 'qd-real :value (mul-qd (qd-value a) (qd-value b))))

(defmethod two-arg-* ((a qd-real) (b real))
  (make-instance 'qd-real :value (mul-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-* ((a real) (b qd-real))
  (make-instance 'qd-real :value (mul-qd-d (qd-value b) (float a 1d0))))

(defmethod two-arg-* ((a number) (b number))
  (cl:* a b))

(defun * (&rest args)
  (if (null args)
      1
      (do ((args (cdr args) (cdr args))
	   (res (car args)
		(two-arg-* res (car args))))
	  ((null args) res))))

(defmethod two-arg-/ ((a qd-real) (b qd-real))
  (make-instance 'qd-real :value (div-qd (qd-value a) (qd-value b))))

(defmethod two-arg-/ ((a qd-real) (b real))
  (make-instance 'qd-real :value (div-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-/ ((a real) (b qd-real))
  (make-instance 'qd-real :value (div-qd (make-qd-d (float a 1d0))
					     (qd-value b))))

(defmethod two-arg-/ ((a number) (b number))
  (cl:/ a b))

(defmethod unary-divide ((a number))
  (cl:/ a))

(defmethod unary-divide ((a qd-real))
  (make-instance 'qd-real :value (div-qd +qd-one+ (qd-value a))))

(defun / (number &rest more-numbers)
  (if more-numbers
      (do ((nlist more-numbers (cdr nlist))
	   (result number))
	  ((atom nlist) result)
         (declare (list nlist))
	 (setq result (two-arg-/ result (car nlist))))
      (unary-divide number)))

(defmethod qlog ((a real) &optional b)
  (if b
      (cl:log a b)
      (cl:log a)))

(defmethod qlog ((a qd-real) &optional b)
  (make-instance 'qd-real
		 :value (if b
			    (/ (log-qd (qd-value a))
			       (if (realp b)
				   (cl:log b)
				   (log-qd (make-qd-d (float b 1d0)))))
			    (log-qd (qd-value a)))))

(declaim (inline log))
(defun log (a &optional b)
  (qlog a b))

(macrolet ((frob (name)
	     (let ((method-name (intern (concatenate 'string "Q" (symbol-name name))))
		   (cl-name (intern (symbol-name name) :cl))
		   (qd-name (intern (concatenate 'string (symbol-name name) "-QD"))))
	       `(progn
		 (defmethod ,method-name ((x number))
		   (,cl-name x))
		 (defmethod ,method-name ((x qd-real))
		   (make-instance 'qd-real :value (,qd-name (qd-value x))))
		 (declaim (inline ,name))
		 (defun ,name (x)
		   (,method-name x))))))
  (frob abs)
  (frob sqrt)
  (frob exp)
  (frob sin)
  (frob cos)
  (frob tan)
  (frob asin)
  (frob acos)
  (frob sinh)
  (frob cosh)
  (frob tanh)
  (frob asinh)
  (frob acosh)
  (frob atanh))

(macrolet ((frob (name &optional (type 'real))
	     (let ((method-name (intern (concatenate 'string "Q" (symbol-name name))))
		   (cl-name (intern (symbol-name name) :cl))
		   (qd-name (intern (concatenate 'string (symbol-name name) "-QD"))))
	       `(progn
		  (defmethod ,method-name ((x ,type))
		    (,cl-name x))
		  (defmethod ,method-name ((x qd-real))
		    (,qd-name (qd-value x)))
		  (declaim (inline ,name))
		  (defun ,name (x)
		    (,method-name x))))))
  (frob zerop)
  (frob plusp)
  (frob minusp))

(defmethod qatan ((y real) &optional x)
  (cl:atan y x))

(defmethod qatan ((y qd-real) &optional x)
  (make-instance 'qd-real
		 :value
		 (if x
		     (atan2-qd (qd-value y) (qd-value x))
		     (atan-qd (qd-value y)))))

(defun atan (y &optional x)
  (qatan y x))


(defmethod qexpt ((x number) (y number))
  (cl:expt x y))

(defmethod qexpt ((x qd-real) (y real))
  (exp (* y (log x))))

(defmethod qexpt ((x real) (y qd-real))
  (exp (* y (log x))))

(defmethod qexpt ((x qd-real) (y qd-real))
  ;; x^y = exp(y*log(x))
  (exp (* y (log x))))

(declaim (inline expt))
(defun expt (x y)
  (qexpt x y))



(macrolet
    ((frob (op)
       (let ((method (intern (concatenate 'string "TWO-ARG-" (symbol-name op))))
	     (cl-fun (find-symbol (symbol-name op) :cl))
	     (qd-fun (intern (concatenate 'string "QD-" (symbol-name op))
			     (find-package :qdi))))
	 `(progn
	    (defmethod ,method ((a real) (b real))
	      (,cl-fun a b))
	    (defmethod ,method ((a qd-real) (b real))
	      (,qd-fun (qd-value a) (make-qd-d (float b 1d0))))
	    (defmethod ,method ((a real) (b qd-real))
	      (,qd-fun (make-qd-d (float a 1d0)) (qd-value b)))))))
  (frob <)
  (frob >)
  (frob <=)
  (frob >=))

(defmethod two-arg-= ((a number) (b number))
  (cl:= a b))
(defmethod two-arg-= ((a qd-real) (b number))
  (if (realp b)
      (qdi::qd-= (qd-value a) (make-qd-d (float b 1d0)))
      nil))
(defmethod two-arg-= ((a number) (b qd-real))
  (if (realp a)
      (qdi::qd-= (make-qd-d (float a 1d0)) (qd-value b))
      nil))

(defun = (number &rest more-numbers)
  "Returns T if all of its arguments are numerically equal, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do ((nlist more-numbers (cdr nlist)))
      ((atom nlist) T)
    (declare (list nlist))
    (if (not (two-arg-= (car nlist) number))
	(return nil))))

(defun /= (number &rest more-numbers)
  "Returns T if no two of its arguments are numerically equal, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((head number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
    (declare (list nlist))
    (unless (do* ((nl nlist (cdr nl)))
		 ((atom nl) T)
	      (declare (list nl))
	      (if (two-arg-= head (car nl))
		  (return nil)))
      (return nil))))

(defun < (number &rest more-numbers)
  "Returns T if its arguments are in strictly increasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg-< n (car nlist))) (return nil))))

(defun > (number &rest more-numbers)
  "Returns T if its arguments are in strictly decreasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg-> n (car nlist))) (return nil))))

(defun <= (number &rest more-numbers)
  "Returns T if arguments are in strictly non-decreasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg-<= n (car nlist))) (return nil))))

(defun >= (number &rest more-numbers)
  "Returns T if arguments are in strictly non-increasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg->= n (car nlist))) (return nil))))

(defmethod qcomplex ((x real) &optional y)
  (cl:complex x (if y y 0)))

(defmethod qcomplex ((x qd-real) &optional y)
  (make-instance 'qd-complex
		 :real (qd-value x)
		 :imag (if y (qd-value y) +qd-zero+)))

(defun complex (x &optional (y 0))
  (qcomplex x y))

(defmethod qinteger-decode-float ((f float))
  (cl:integer-decode-float f))

(defmethod qinteger-decode-float ((f qd-real))
  (integer-decode-qd (qd-value f)))

(declaim (inline integer-decode-float))
(defun integer-decode-float (f)
  (qinteger-decode-float f))

(defmethod qdecode-float ((f float))
  (cl:decode-float f))

(defmethod qdecode-float ((f qd-real))
  (multiple-value-bind (frac exp s)
      (decode-float-qd (qd-value f))
    (values (make-instance 'qd-real :value frac)
	    exp
	    (make-instance 'qd-real :value  s))))

(declaim (inline decode-float))
(defun decode-float (f)
  (qdecode-float f))

(defmethod qscale-float ((f float) (n integer))
  (cl:scale-float f n))

(defmethod qscale-float ((f qd-real) (n integer))
  (make-instance 'qd-real :value (scale-float-qd (qd-value f) n)))

(declaim (inline scale-float))
(defun scale-float (f n)
  (qscale-float f n))

(defmethod qfloor ((x real) &optional y)
  (if y
      (cl:floor x y)
      (cl:floor x)))

(defmethod qfloor ((x qd-real) &optional y)
  (if (and y (/= y 1))
      (let ((f (qfloor (/ x y))))
	(values f
		(- x (* f y))))
      (let ((f (qdi::ffloor-qd (qd-value x))))
	(multiple-value-bind (int exp sign)
	    (integer-decode-qd f)
	  (values (ash (* sign int) exp)
		  (make-instance 'qd-real
				 :value (qd-value
					 (- x (make-instance 'qd-real
							     :value f)))))))))

(defun floor (x &optional y)
  (qfloor x y))

(defmethod qffloor ((x real) &optional y)
  (if y
      (cl:ffloor x y)
      (cl:ffloor x)))

(defmethod qffloor ((x qd-real) &optional y)
  (if (and y (/= y 1))
      (let ((f (qffloor (/ x y))))
	(values f
		(- x (* f y))))
      (let ((f (make-instance 'qd-real :value (qdi::ffloor-qd (qd-value x)))))
	(values f
		(- x f)))))

(defun ffloor (x &optional y)
  (qffloor x y))

(defun ceiling (x &optional y)
  (multiple-value-bind (f rem)
      (floor x y)
    (if (zerop rem)
	(values (+ f 1)
		rem)
	(values (+ f 1)
		(- rem 1)))))

(defun fceiling (x &optional y)
  (multiple-value-bind (f rem)
      (ffloor x y)
    (if (zerop rem)
	(values (+ f 1)
		rem)
	(values (+ f 1)
		(- rem 1)))))

(defun truncate (x &optional (y 1))
  (if (minusp x)
      (ceiling x y)
      (floor x y)))

(defun ftruncate (x &optional (y 1))
  (if (minusp x)
      (fceiling x y)
      (ffloor x y)))

(defmethod %unary-round ((x real))
  (cl::round x))

(defmethod %unary-round ((number qd-real))
  (multiple-value-bind (bits exp)
      (integer-decode-float number)
    (let* ((shifted (ash bits exp))
	   (rounded (if (and (minusp exp)
			     (oddp shifted)
			     (not (zerop (logand bits
						 (ash 1 (- -1 exp))))))
			(1+ shifted)
			shifted)))
      (if (minusp number)
	  (- rounded)
	  rounded))))

(defun round (number &optional (divisor 1))
  (if (eql divisor 1)
      (let ((r (%unary-round number)))
	(values r
		(- number r)))
      (multiple-value-bind (tru rem)
	  (truncate number divisor)
	(if (zerop rem)
	    (values tru rem)
	    (let ((thresh (/ (abs divisor) 2)))
	      (cond ((or (> rem thresh)
			 (and (= rem thresh) (oddp tru)))
		     (if (minusp divisor)
			 (values (- tru 1) (+ rem divisor))
			 (values (+ tru 1) (- rem divisor))))
		    ((let ((-thresh (- thresh)))
		       (or (< rem -thresh)
			   (and (= rem -thresh) (oddp tru))))
		     (if (minusp divisor)
			 (values (+ tru 1) (- rem divisor))
			 (values (- tru 1) (+ rem divisor))))
		    (t (values tru rem))))))))
