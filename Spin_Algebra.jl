#%###################################
#%
#% Tensor algebra and Spin algebra  
#%
# A small library for quantum simlation of fermions written in Julia.
# Copyright (c) 2017-8 James Daniel Whitfield
# Dartmouth College, Department of Physics and Astronomy
#
# see parity-qsim.pdf for an introduction to the techniques used here.
#
######################
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
######################

#%###################################
#%We'll deal with systems of spin-half particles. These are governed by the
#%Pauli spin matrices.  Additionally, the raising and lowering spin operators are
#%useful.  We'll assign them all numbers and define algebra as a multiplication
#%table
#%
#% Pauli basis: I, X, Y, Z
#% Index basis: P_00, P_01, P_10, P_11
#%   with
#%  P_00 = (I+Z)/2
#%  P_11 = (I-Z)/2
#%  P_01 = S_+ = (X+iY )/2
#%  P_10 = S_- = (X-iY )/2
#%
#% Basic algebraic relations
#%   XY= iZ    YX= -iZ
#%   ZX= iY    XZ= -iY
#%   YZ= iX    ZY= -iX
#%  matrix    =  GetSpinAlgebra()
#%  MatrixID  =  MultiplyPauli(MatrixID, MatrixID)

########### Z BASIS ################
#Pauli Matrix basis
ii=[1. 1.0; 0.0  1.0];
x= [0  1.0; 1.0  0.0];
y= [0 -1im; 1im  0.0];
z= [1. 0.0; 0.0 -1.0];
#Index basis
p00=[1. 0.; 0. 0.];
p01=[0. 1.; 0. 0.];
p10=[0. 0.; 1. 0.];
p11=[0. 0.; 0. 1.];

# We use integer MatrixIDs for coding the algebra
I   =1;
X   =2;
Y   =3;
Z   =4;
P00 =5;
P01 =6;
P10 =7;
P11 =8;

function GetSpinAlgebra()

   SPIN_ALGEBRA=zeros(Complex64,8,8)
   #These are the MatrixIDs for coding algebra
   I   =1;
   X   =2;
   Y   =3;
   Z   =4;
   P00 =5;
   P01 =6;
   P10 =7;
   P11 =8;

   #ROW I.j
   #   [ 11  1X  1Y  1Z  1 P_00   1 P_01   1 P_10   1 P_11
   # = [ +1   X   Y   Z  +P_00    +P_01    +P_10    +P_11
   SPIN_ALGEBRA[I,I]=I;
   SPIN_ALGEBRA[I,X]=X;
   SPIN_ALGEBRA[I,Y]=Y;
   SPIN_ALGEBRA[I,Z]=Z;
   SPIN_ALGEBRA[I,P00]=P00;
   SPIN_ALGEBRA[I,P01]=P01;
   SPIN_ALGEBRA[I,P10]=P10;
   SPIN_ALGEBRA[I,P11]=P11;

   #ROW X.j
   #   [ X1  XX  XY  XZ  X P_00   X P_01   X P_10   X P_11]
   # = [ +X  +1 +iZ -iY  +P_10    +P_11    +P_00    +P_01 ]
   SPIN_ALGEBRA[X,I]	=    X;
   SPIN_ALGEBRA[X,X]	=    I;
   SPIN_ALGEBRA[X,Y]	=  1im*Z;
   SPIN_ALGEBRA[X,Z]	= -1im*Y;
   SPIN_ALGEBRA[X,P00]     =   P10;
   SPIN_ALGEBRA[X,P01]     =   P11;
   SPIN_ALGEBRA[X,P10]     =   P00;
   SPIN_ALGEBRA[X,P11]     =   P01;

   #ROW Y.j
   #   [  Y1  YX  YY  YZ   Y P_00   Y P_01   Y P_10   Y P_11 ]
   # = [  +Y -iZ  +1 +iX   +iP_10   +iP_11   -iP_00   -iP_01 ]
   SPIN_ALGEBRA[Y,I]	=    Y;
   SPIN_ALGEBRA[Y,X]	= -1im*Z;
   SPIN_ALGEBRA[Y,Y]	=    I;
   SPIN_ALGEBRA[Y,Z]	= +1im*X;
   SPIN_ALGEBRA[Y,P00]     = +1im*P10;
   SPIN_ALGEBRA[Y,P01]     = +1im*P11;
   SPIN_ALGEBRA[Y,P10]     = -1im*P00;
   SPIN_ALGEBRA[Y,P11]     = -1im*P01;

   #ROW Z.j
   #     Z1  ZX  ZY  ZZ   Z  P_00    Z  P_01    Z  P_10     Z  P_11 ]
   # = [ +Z +iY -iX  +1    +P_00      +P_01      -P_10       -P_11  ]
   SPIN_ALGEBRA[Z,I]	=    Z;
   SPIN_ALGEBRA[Z,X]       = +1im*Y;
   SPIN_ALGEBRA[Z,Y]	= -1im*X;
   SPIN_ALGEBRA[Z,Z]	=    I;
   SPIN_ALGEBRA[Z,P00] 	=   P00;
   SPIN_ALGEBRA[Z,P01]	=   P01;
   SPIN_ALGEBRA[Z,P10]	=  -P10;
   SPIN_ALGEBRA[Z,P11]	=  -P11;

   #ROW P00.j
   #   [ P_00 1   P_00 X   P_00 Y   P_00 Z    P_00 P_00   P_00 P_01   P_00 P_10   P_00 P_11
   # = [  +P_00    +P_01   -iP_01    +P_00      P_00         P_01         0          0
   SPIN_ALGEBRA[P00,I]	=    P00;
   SPIN_ALGEBRA[P00,X]	=    P01;
   SPIN_ALGEBRA[P00,Y]	= -1im*P01;
   SPIN_ALGEBRA[P00,Z]	=    P00;
   SPIN_ALGEBRA[P00,P00]   =    P00;
   SPIN_ALGEBRA[P00,P01]   =    P01;
   SPIN_ALGEBRA[P00,P10]   =     0;
   SPIN_ALGEBRA[P00,P11]   =     0;

   #ROW P01.j
   #   [  P_01 1   P_01 X   P_01 Y   P_01 Z   P_01 P_00   P_01 P_01   P_01 P_10   P_01 P_11
   # = [   +P_01    +P_00   +iP_00    -P_01       0            0         P_00       P_01
   SPIN_ALGEBRA[P01,I]	=    P01;
   SPIN_ALGEBRA[P01,X]	=    P00;
   SPIN_ALGEBRA[P01,Y]	= +1im*P00;
   SPIN_ALGEBRA[P01,Z]	=   -P01;
   SPIN_ALGEBRA[P01,P00]   =     0;
   SPIN_ALGEBRA[P01,P01]   =     0;
   SPIN_ALGEBRA[P01,P10]   =    P00;
   SPIN_ALGEBRA[P01,P11]   =    P01;

   #ROW P10.j
   #   [  P_10 1   P_10 X   P_10 Y   P_10 Z   P_10 P_00   P_10 P_01   P_10 P_10   P_10 P_11
   # = [   +P_10    +P_11   -iP_11    +P_10       P_10         P_11         0          0
   SPIN_ALGEBRA[P10,I]	=    P10;
   SPIN_ALGEBRA[P10,X]	=    P11;
   SPIN_ALGEBRA[P10,Y]	= -1im*P10;
   SPIN_ALGEBRA[P10,Z]	=   +P01;
   SPIN_ALGEBRA[P10,P00]   =    P10;
   SPIN_ALGEBRA[P10,P01]   =    P11;
   SPIN_ALGEBRA[P10,P10]   =     0;
   SPIN_ALGEBRA[P10,P11]   =     0;

   #ROW P11.j
   #   [  P_11 1   P_11 X   P_11 Y   P_11 Z   P_11 P_00   P_11 P_01   P_11 P_10   P_11 P_11
   # = [   +P_11    +P_10   +iP_10    -P_11       0            0         P_10       P_11
   SPIN_ALGEBRA[P11,I]	=    P11;
   SPIN_ALGEBRA[P11,X]	=    P10;
   SPIN_ALGEBRA[P11,Y]	= +1im*P10;
   SPIN_ALGEBRA[P11,Z]	=   -P11;
   SPIN_ALGEBRA[P11,P00]   =     0;
   SPIN_ALGEBRA[P11,P01]   =     0;
   SPIN_ALGEBRA[P11,P10]   =    P10;
   SPIN_ALGEBRA[P11,P11]   =    P11;

   return SPIN_ALGEBRA;
end

function MultiplyPauli(A,B)
  algebra=GetSpinAlgebra();
  return algebra[A,B];
end

#%###################################
#%          Tensor algebra          #
#%###################################
#% Spin chain multiplication
#%
#% We need both multiplication and addition of Pauli sums.   We will store the
#% spin chain operators as a list with a complex prefactor. Plain ASCII output
#% is required for the final output of the result.  The present
#% example is for small matrix size so we will explicitly store each tensor
#% contribution.
#%
#% PS    = PauliSum(PL1,PL2,...,PLn)
#% PList = PauliList(Prefactor,OperatorList)
#%
#% Pauli list functions
#%  PList        =  StandardForm(PList)
#%  PList        =  Minus(Plist)
#%  PList        =  MakeLocalOperator(MatrixID, tensor_index, num_spins)
#%  PList        =  MultiplyPauliList(PauliL::PauliList,PauliR::PauliList)
#%  PLists       =  MultiplyPauliLists(PauliSumL,PauliSumR)
#%  PList        =  Sum(PList)

function StandardForm(P::PauliList)

  #a little error catching to make sure MatrixIDs are okay
  if( 0!=sum(abs.(P.List).>8.0) || 0!=sum(abs.(P.List).<1.0))
    println("Warning invalid MatrixIDs in StandardForm function")
  end

  #pull StandardForm
  b=(P.List)./abs.(P.List);

  #if entry of list is zero, extracting the phase via A/|A| is invalid. Hence,
  ## !isnan ## command. All MatrixIDs should be 1,...,8 so this isn't necessary.
  b=prod( b[.!isnan.(b)] );

  P.Prefactor=P.Prefactor*b;

  #put everything as standard matrix IDs
  P.List=Complex.(round.(float(abs.(P.List))));

  return P;
end


function MakeLocalOperator(MatrixID, j, M)
  list=ones(M);
  list[j]=MatrixID;
  return PauliList(1,list)
end

function MultiplyPauliList(L::PauliList,R::PauliList)
  algebra=GetSpinAlgebra();

  out=zeros(Complex,max(length(L.List),length(R.List)));
  if(length(L.List)!=length(R.List))
    println("Warning: tensor products must be same length. Returning 0")
    return 0;
  end


  F=L.Prefactor*R.Prefactor;
  for k=1:length(L.List)
    out[k]=algebra[Int64(abs(L.List[k])),Int64(abs(R.List[k]))]
  end

  return StandardForm(PauliList(F,out));
end

function MultiplyPauliLists(aL,aR)
  #Here aL and aR are the arrays of PauliLists
  out=Array{PauliList}([]);
  for i=1:length(aL)
    for j=1:length(aR)
      if length(out)==0
        out=[MultiplyPauliList( aL[i], aR[j] )]
      else
        out=[out;MultiplyPauliList(aL[i],aR[j])]
      end
    end
  end
  return out;
end



#adds like list terms for simplification (S. Manski 2016)
#input: unsimplified PauliList
#output: PauliList with combined like terms
function Sum(arrayPL)

  arrayA=PauliList[];
  for i=1:length(arrayPL)
    for j=i+1:length(arrayPL)
      if( arrayPL[i].List == arrayPL[j].List && arrayPL[i].Prefactor!=0)
          arrayPL[i].Prefactor=arrayPL[i].Prefactor+arrayPL[j].Prefactor;
          arrayPL[j].Prefactor=0;
      end
    end
    if(arrayPL[i].Prefactor!=0)
      push!(arrayA,arrayPL[i]);
    end
  end
  return arrayA
end




# Wrapper functions
#  ans    =  Multiply(A,B)
#  output =  Print(A)


function Multiply(A,B)

  if typeof(A)<: Number
    #put the scalar second
    temp=A;
    A=B;
    B=temp;
  end


  if  typeof(A) <: Vector{PauliList} && typeof(B) <: Number
    for a in A
      a.Prefactor=a.Prefactor * B;
    end
    return A;
  end

  if typeof(A) <: PauliList && typeof(B) <: Number
    A.Prefactor=A.Prefactor*B;
    return A;
  end


  if(string(typeof(A))== "PauliList" && string(typeof(B))== "PauliList")
		return MultiplyPauliList(A,B)
	end

	if(   (     string(typeof(A))== "Array{PauliList,1}"  ||      string(typeof(A))== "Vector{PauliList}" )
           && (     string(typeof(B))== "Array{PauliList,1}"  ||      string(typeof(B))== "Vector{PauliList}" ))
		return MultiplyPauliLists(A,B)
	end

	println("Types: ",typeof(A)," , ",typeof(B))
	return A*B
end

#TODO OVERLOAD * TO USE Multiply
#import Base.*
#println(*(3,9))

