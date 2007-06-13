(in-package "QD")

(define-symbol-macro * cl:*)
(define-symbol-macro - cl:-)
(define-symbol-macro / cl:/)

(defclass qd-real ()
  ((qd :initform +qd-zero+
       :reader qd-value
       :initarg :value
       :type %quad-double)))

(defclass qd-complex ()
  ((real :initform +qd-zero+
	 :reader qd-real
	 :initarg :real
	 :type %quad-double)
   (imag :initform +qd-zero+
	 :reader qd-imag
	 :initarg :imag
	 :type %quad-double)))

#-cmu
(defmethod print-object ((qd qd-real) stream)
  (format stream "~/qdi::qd-format/" (qd-value qd)))

#+cmu
(defmethod print-object ((qd qd-real) stream)
  (let ((q (qd-value qd)))
    (if (or (ext:float-infinity-p (qd-0 q))
	    (ext:float-nan-p (qd-0 q)))
	(format stream "~/qdi::qd-format/" q)
	(format stream "#q~/qdi::qd-format/" q))))

(defmethod make-qd ((x real))
  (make-instance 'qd-real :value (make-qd-d (float x 1d0))))

(defmethod make-qd ((x qd-real))
  (make-instance 'qd-real :value (qd-value x)))
  
(defmethod print-object ((qd qd-complex) stream)
  (format stream "#q(~/qdi::qd-format/ ~/qdi::qd-format/)"
	  (qd-real qd)
	  (qd-imag qd)))

(defmethod qd-value ((x real))
  (make-qd-d (float x 1d0)))

(defmethod make-load-form ((qd qd-real) &optional environment)
  (declare (ignore environment))
  `(make-instance ',(class-of qd)
		  :value ',(qd-value qd)))

(defmethod make-load-form ((qd qd-complex) &optional environment)
  (declare (ignore environment))
  `(make-instance ',(class-of qd)
		  :real ',(qd-value (realpart qd))
		  :imag ',(qd-value (imagpart qd))))

(defmethod describe-object ((q qd-real) stream)
  (multiple-value-bind (q0 q1 q2 q3)
      (qd-parts (qd-value q))
    (format stream "~&~S is a quad-double with components ~
                    ~%  ~A, ~A, ~A, ~A~%"
	    q q0 q1 q2 q3)))
