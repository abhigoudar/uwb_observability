(* ::Package:: *)

p = {px,py,pz} ;
v = {vx,vy,vz} ;
q = {q0,q1,q2,q3} ;
ba = {bax,bay,baz}; 
bw = {bwx,bwy,bwz} ;
pt = {ptx,pty,ptz} ;
R = {{q0^2 + q1^2 - q2^2 - q3^2, 2*q1*q2 - 2*q0*q3, 2*q0*q2 + 2*q1*q3},
	 {2*q0*q3 + 2*q1*q2, q0^2-q1^2 + q2^2 - q3^2, 2*q2*q3 - 2*q0*q1},
	 {2*q1*q3 - 2*q0*q2, 2*q0*q1 + 2*q2*q3, q0^2 - q1^2 - q2^2 + q3^2}};
pwp = {pwpx,pwpy,pwpz};
pwq = {pwqx,pwqy,pwqz};
pwr = {pwrx,pwry,pwrz};

state = Join[p,v, q, ba, bw, pt];
zv = {0,0,0};
zm = {{0,0,0},{0,0,0},{0,0,0}};
zm4 = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
g = {0,0,9.8};
Xi = {{-q1, -q2, -q3},
      { q0, -q3,  q2},
      { q3,  q0, -q1},
      {-q2,  q1,  q0}};
      
f0 = Join[v, -R.ba + g, -0.5 * Xi.bw,zv,zv,zv];
f1 = Join[zm, R, zm4, zm, zm, zm];
f2 = Join[zm, zm, 0.5* Xi, zm, zm, zm];

deltapwp = pwp - R.pt-p;
deltapwq = pwq - R.pt-p;
deltapwr = pwr - R.pt-p;

L0y =  0.5 *{ deltapwp.deltapwp, deltapwq.deltapwq, deltapwr.deltapwr};
"Zeroth Order Lie Derivative"
StringForm["\t L0y Dimension=``",Dimensions[L0y]]
gradL0y = D[L0y, {state}];
StringForm["\t\t Dimension of gradL0y=``",Dimensions[gradL0y]]
(*StringForm["\t\t Rank of gradL0y=``",MatrixRank[gradL0y]]*)

"First Order Lie Derivative"

L1f0y = gradL0y.f0;
StringForm["\t L1f0y Dimension=``",Dimensions[L1f0y]]
gradL1f0y = D[L1f0y, {state}];
StringForm["\t\t Dimension of gradL1f0y=``",Dimensions[gradL1f0y]]
(*StringForm["\t\t Rank of gradL1f0y=``",MatrixRank[gradL1f0y]]*)

L1f2y = gradL0y.f2;
StringForm["\t L1f2y Dimension=``",Dimensions[L1f2y]]
gradL1f2y1 =  D[L1f2y[[All,1]], {state}];
gradL1f2y2 =  D[L1f2y[[All,2]], {state}];
gradL1f2y3 =  D[L1f2y[[All,3]], {state}];
gradL1f2y =Join[gradL1f2y1,gradL1f2y2,gradL1f2y3] ;
StringForm["\t\t Dimension of gradL1f2y=``",Dimensions[gradL1f2y]]

"Second Order Lie Derivative"

L2f0y = gradL1f0y.f0;
StringForm["\t L2f0y Dimension=``",Dimensions[L2f0y]]
gradL2f0y = D[L2f0y, {state}];
StringForm["\t\t Dimension of gradL2f0y=``",Dimensions[gradL2f0y]]

L2f2f0y = gradL1f0y.f2;
StringForm["\t L2f2f0y Dimension=``",Dimensions[L2f2f0y]]
gradL2f2f0y1 =  D[L2f2f0y[[All,1]], {state}];
gradL2f2f0y2 =  D[L2f2f0y[[All,2]], {state}];
gradL2f2f0y3 =  D[L2f2f0y[[All,3]], {state}];
gradL2f2f0y =Join[gradL2f2f0y1, gradL2f2f0y2, gradL2f2f0y3] ;
StringForm["\t\t Dimension of gradL2f2f0y=``",Dimensions[gradL2f2f0y]]




Obs = Join[gradL0y, gradL1f0y, gradL2f0y, gradL1f2y, gradL2f2f0y1];
StringForm["\t\t Original Observability matrix dimension=``",Dimensions[Obs]]
ObsMat = Drop[Obs, None, None];
(*ObsMat = ObsMat/.{pwpx->0,pwpy->0,pwpz->0,pwqx\[Rule]0,pwqy->0,pwqz->0,pwrx->5,pwry->5,pwrz->5}*)
MatrixPlot[ObsMat]
StringForm["\t\t Observability matrix dimension=``",Dimensions[ObsMat]]
"\t\t Calculating observability...``"
" Iteration 1"
StringForm["\t\t Rank=``",MatrixRank[ObsMat]]
"----------------------------------------------"
