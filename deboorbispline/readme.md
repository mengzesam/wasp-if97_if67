
bisurfit.f is a total fortran file which is consist of other *.f files.

deboorbispline.h and deboorbispline.c are c interface functions.

bisplinetest.c is a test file

csurfit, cbspl and cbisp function are bivariate b-spline functions, they base on surfit.f, fpbspl.f and fpbisp.f
*which were programmed by Prof. P. Dierckx and are re-written with cpp.

*surfit_ is c interface function of surfit.f

*fpbspl.f and fpbisp.f author:

        Prof. P. Dierckx
	
		Department of Computer Science
		
		K.U.Leuven
		
		Celestijnenlaan 200 A
		
		B 3001, Heverlee, Belgium 
	
*please refer to http://www.netlib.org/dierckx/index.html for details
*Re-written by:

        mengzesen<mengzesam@126.com>

c and fortran sources compiling:

gcc cfile1.c cfile2 ffile1.f ffile2.f -o out.exe

There are two library files:deboorbispline.o and bisurfit.o,which are compiled with -O2

In cpp,you can combine them with your application by means of:
g++ bisurfit.o deboorbispline.o  yourcppfile.cpp -o yourapplication.out(or exe)
