V1 n1 0 1
R1 n1 0 2
*Xa1 n2 n3 Mina2
*C1 n2 0 1
*Xb1 n3 0 Mina

R2 n3 0 1

.include ./netlists/include.sp

.subckt Mina n1 n2
R1 n1 n2 10
.ends

.subckt Mina2 n1 n2
L1 n1 n2 1
.ends

 
.tran 0.1 1
.print tran V(n3)