
.param Is = 1e-12

V1 n1 0 sin(0 1 1 0 0 0) 
*V1 0 n1 sin(1 1 90)
*R1 n1 0 1
*Xa1 n2 n3 Mina2
*C1 n2 0 1
*Xb1 n3 0 Mina

*F1 V1 n2 0 100.1
R2 n2 0 10
G1 n1 n2 CUR='Is*(exp((V(n1)-V(n2))*40)-1)'
*.include ./netlists/include.sp

*.subckt Mina n1 n2
*R1 n1 n2 10
*.ends

*.subckt Mina2 n1 n2
*L1 n1 n2 1
*.ends

 
.op
.tran 0.001 3
.print tran V(n1) V(n2)
*.option method = pade 1 1

.end
