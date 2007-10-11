;;;; -*- Mode: lisp -*-
;;;;
;;;; Copyright (c) 2007 Raymond Toy
;;;;
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

;;; Here are some simple timing tests, based on Yozo Hida's qd_timer
;;; test code. I've tried to make these versions time the same
;;; operations as Yozo's.

(in-package #:qd)

(defun time-add (&optional (n 100000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a (/ #q7))
	(b #q0))
    (declare (type qd-real a b))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (incf b a)))
    (format t "n = ~D~%" n)
    (format t "b = ~W~%" b)))

(defun time-mul (&optional (n 100000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a (+ 1 (/ (float n #q1))))
	(b #q1))
    (declare (type qd-real a b))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (setf b (* a b))))
    (format t "n = ~D~%" n)
    (format t "b = ~W~%" b)))
  
(defun time-mul (&optional (n 100000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a (+ 1 (/ (float n #q1))))
	(b #q1))
    (declare (type qd-real a b))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (setf b (* a b))))
    (format t "n = ~D~%" n)
    (format t "b = ~W~%" b)))

(defun time-div (&optional (n 100000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a (+ 1 (/ (float n #q1))))
	(b #q1))
    (declare (type qd-real a b))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (setf b (/ b a))))
    (format t "n = ~D~%" n)
    (format t "b = ~W~%" b)))

(defun time-sqrt (&optional (n 100000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a #q0)
	(b (+ 2 +pi+)))
    (declare (type qd-real a b))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (setf a (sqrt (+ a b)))))
    (format t "n = ~D~%" n)
    (format t "a = ~W~%" a)))

(defun time-sin (&optional (n 2000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a #q0)
	(b (/ +pi+ n))
	(c #q0))
    (declare (type qd-real a b c))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (incf a b)
	    (incf c (sin a))))
    (format t "n = ~D~%" n)
    (format t "a = ~W~%" a)))

(defun time-log (&optional (n 1000))
  (declare (fixnum n)
	   (optimize (speed 3)))
  (let ((a #q0)
	(c (exp #q-50.1))
	(d (exp (/ #q100.2 n))))
    (declare (type qd-real a c d))
    (time (dotimes (k n)
	    (declare (fixnum k))
	    (incf a (log c))
	    (setf c (* c d))))
    (format t "n = ~D~%" n)
    (format t "a = ~W~%" a)))

#||

Some test results.  These were all run on a Sun Blade 1500 using a 1.5
GHz Ultrasparc III.  I used the default configuration when compiling
qd, and used Sun's C++ compiler.  For the Lisp timing, I used CMUCL.

Executive summary:

Test	    Time
	qd	oct
----	-----------
add	0.036	0.09
mul	0.117	0.13
div	0.388	0.29
sqrt	0.142	0.11
sin	0.136	0.14
log	0.231	0.12

Times are in sec for the test.  The default number of iterations were
used.  Some of the results are a bit surprising.  I was expecting the
C++ code to be faster, and that's the case for add and mul.  But oct
with CMUCL was faster in div, sqrt, and log.  Note that oct uses a
different algorithm for log than qd, so that could explain the
difference.  Note also that CLOS dispatch is included in these timing
results.  CMUCL's CLOS dispatch is ok, but not great.


-------------------------------------------------------------------------------
The raw data:

The output from qd_timer -qd -v:

Timing qd_real
--------------

Timing addition...
n = 100000   t = 0.0362288
b = 1.428571e+04
100000 operations in 0.0362288 s.
  0.362288 us

Timing multiplication ...
n = 100000   t = 0.11686
b = 2.718268e+00
100000 operations in 0.11686 s.
  1.168602 us

Timing division ...
n = 100000   t = 0.388279
b = 0.367881
100000 operations in 0.388279 s.
  3.882788 us

Timing square root ...
n = 10000   t = 0.141866
a = 2.821980
10000 operations in 0.141866 s.
 14.186575 us

Timing sin ...
n = 2000   t = 0.136080
a = 3.141593
2000 operations in 0.136080 s.
 68.039791 us

Timing log ...
n = 1000   t = 0.230506
a = -50.100000
1000 operations in 0.230506 s.
230.506166 us

The output from CMUCL:

QD> (time-add)

; Evaluation took:
;   0.09 seconds of real time
;   0.1 seconds of user run time
;   0.0 seconds of system run time
;   147,285,856 CPU cycles
;   0 page faults and
;   7,200,016 bytes consed.
; 
n = 100000
b = #q14285.7142857142857142857142857142857142857142857142857142857142855q0
NIL
QD> (time-mul)

; Evaluation took:
;   0.13 seconds of real time
;   0.1 seconds of user run time
;   0.02 seconds of system run time
;   203,790,588 CPU cycles
;   0 page faults and
;   7,200,824 bytes consed.
; 
n = 100000
b = #q2.71826823717448966803506482442604644797444693267782286300915989397q0
NIL
QD> (time-div)

; Evaluation took:
;   0.29 seconds of real time
;   0.28 seconds of user run time
;   0.01 seconds of system run time
;   460,956,912 CPU cycles
;   0 page faults and
;   7,200,016 bytes consed.
; 
n = 100000
b = #q0.36788128056098406210328658773118942247132502490133718973918140856q0
NIL
QD> (time-sqrt 10000)

; Evaluation took:
;   0.11 seconds of real time
;   0.1 seconds of user run time
;   0.0 seconds of system run time
;   173,209,708 CPU cycles
;   0 page faults and
;   2,402,560 bytes consed.
; 
n = 10000
a = #q2.82198033014704783016853125515542796898998765943212617578596649019q0
NIL
QD> (time-sin)

; Evaluation took:
;   0.14 seconds of real time
;   0.14 seconds of user run time
;   0.0 seconds of system run time
;   213,378,476 CPU cycles
;   0 page faults and
;   3,105,800 bytes consed.
; 
n = 2000
a = #q3.14159265358979323846264338327950288419716939937510582097494459409q0
NIL
QD> (time-log)

; Evaluation took:
;   0.12 seconds of real time
;   0.12 seconds of user run time
;   0.01 seconds of system run time
;   192,187,304 CPU cycles
;   0 page faults and
;   1,621,792 bytes consed.
; 
n = 1000
a = #q-50.100000000000000000000000000000000000000000000000000000000208796q0
NIL
QD> 

---------------------------------------------
||#
