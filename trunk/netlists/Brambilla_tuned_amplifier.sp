

.PARAM Af = 0.99
.PARAM Af_inv = 1.01010101
.PARAM Ar=0.02
.PARAM Ar_inv = 50
.PARAM Is=1.0e-15

*.option maxdeg = 9
*.option method = pade 1 1
*.option truntol = 1e-4
*.option outpath =  /home/egad/projects/HiSPICE/trunk/Hispice/results/
*.option Error_Control =  ERROR_BY_NSTEP
*.hb tones = 5XHz harms =10 intmodmax=0
.OP
*.hb 5Xhz harms =100
.tran 3e-10 4e-3


.print TRAN v(nVin)	v(nvo) 
*.print HB v(nVin)	v(nVo)


Vdd    ndd  0		20

******************
* Source for harmonic balance analysis
*Vin nVin  	       	0 	   	HB .9 0 1 1 0 0 
*******************


**************
* Source for transient analysis
*Vin     nVin 	        0		Pulse (0 1 0.1n 0.1n 0.1n 1n 5n)
Vin    nVin	       	0		sin(0 1.0 5e6 0 0 0) 
*************

C1     nVin  n1			33e-9
C2     n1    0		3.3e-9
R1     n1    ndd		8.2e3
R2     n1    0		4.7e3
R3     ndd   nvo		750.0	

xQ     n1    nvo		a1	npn
R4     a1     n3			750.0
Le     n3    0		220e-6

xCrystal    a1		0	xtal


.subckt     xtal 	nA	nB

Rm	    nA		nC	120
Cm	    C1		nB	0.165786399
Co	    nA		nB	5e-12
Lm	    nC		C1	6.11154981
.ends

.subckt npn nb nc ne
	GD1 n01 n02 cur = '(is*Ar_inv)*(exp((v(n01)-v(n02))*40) - 1)'  
	GD2 n01 n03 cur = '(Is*Af_inv)*(exp((v(n01)-v(n03))*40) - 1)'  

*	D1  n01 n02 5.0e-14 0.025
*	D2  n01 n03 1.010101010101010e-015 0.025

	Vdc n02 n04 0
	Vde n03 n05 0
	
	F1 Vde n04 n01  0.99 
	F2 Vdc n05 n01  0.02 
	
	*Ccx nb n04 0.1p
	*Ce n01 n05 1p
	*Cc n01 n04 0.8p
	Rb nb n01 100
  Rc n04 nc 10
  Re n05 ne 1

.ends

*.ic
*+ v(a1) =   6.209274727314458e+00
*+ v(n1) =   7.039477573282669e+00
*+ v(n3) =  0
*+ v(ndd) =   20.0000
*+ v(nvin) =   1.0000
*+ v(nvo) =   1.385281801992159e+01
*+ v(xq.n01) =   7.031198540317818e+00
*+ v(xq.n02) =   1.377085559352054e+01
*+ v(xq.n03) =     6.217553760284211e+00
*+ v(xq.n04) =  1.377085559352054e+01
*+ v(xq.n05) =    6.217553760284211e+00
*+ i(vdd) =    -9.776794155557879e-03
*+ i(Vin) =        0
*+ i(le) =    8.279032969752610e-03
*+ i(xq.vdc) =   -0.0000
*+ i(xq.vde) =    8.279032969752059e-03



*.END
