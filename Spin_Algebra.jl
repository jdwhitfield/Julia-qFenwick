#%###################################
#%
#% Tensor algebra and Spin algebra
#%
# A library for quantum simlation of fermions written in Julia.
# Copyright (c) 2017-9 James Daniel Whitfield
# Dartmouth College, Department of Physics and Astronomy
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
#% Tensor algebra and Spin algebra  #
#%###################################
#%   The Pauli list consists of a complex prefactor and an array listing the
#%   matrix ids in the Pauli string
#%
#%   Methods:
#%       PrintPauliList(PauliList)
#%       PrintPauliLists(PauliLists)

#Pauli list
temp = isdefined(Main,:PauliList)
if( temp == false )
  mutable struct PauliList
    Prefactor::Complex
    List::Array{Complex,1}
  end
end

"""PRINTPAULILIST(pauli_list) is a function to print Pauli list. Optional IOStream"""
function PrintPauliList(p::PauliList, io::IO = stdout)

  #correct the state
  P=StandardForm(p)

  #These are the MatrixIDs for coding the algebra
  I   =1;
  X   =2;
  Y   =3;
  Z   =4;
  P00 =5;
  P01 =6;
  P10 =7;
  P11 =8;

  """
  You can stick strings together (a process often called concatenation) using the multiply (*) operator:
    julia> "s" * "t"
    "st"
  If you've used other programming languages, you might expect to use the addition (+) operator:
  """

  #printing
  print(io,"(")
  if(real(P.Prefactor)>=0) #add plus sign to make it look better
    print(io,"+")
  end

  """
  # numbers can be converted to strings and formatted using printf
  @printf "e = %0.2f\n" e
  #> 2.718
  """
  print(io,P.Prefactor)

  print(io,") ")
  for i=1:length(P.List)
    if(P.List[i]==I)
      print(io,"  I  ")
    end
    if(P.List[i]==X)
      print(io,"  X  ")
    end
    if(P.List[i]==Y)
      print(io,"  Y  ")
    end
    if(P.List[i]==Z)
      print(io,"  Z  ")
    end
    if(P.List[i]==P00)
      print(io," P00 ")
    end
    if(P.List[i]==P01)
      print(io," P01 ")
    end
    if(P.List[i]==P10)
      print(io," P10 ")
    end
    if(P.List[i]==P11)
      print(io," P11 ")
    end
  end
  println(io,"")
  return
end

"""PRINTPAULILISTS(pauli_list) is a function to print Pauli list"""
function PrintPauliLists(arrayPL,io::IO=stdout)

  if(string(typeof(arrayPL))=="PauliList")
    arrayPL=[arrayPL];
  end

  if(length(arrayPL)==0)
    print(io,"0")
    return
  end

  if(string(typeof(arrayPL[1]))!="PauliList")
    println(io,"Warning: bad input; expected array of PauliLists in function Print_sum")
  end

  for k=1:length(arrayPL)
    if(k>1)
      print(io,"    +")
    else
      print(io,"op = ")
    end
    Print(arrayPL[k],io)
  end
end

"""Extended print function for Pauli lists and arrays of Pauli lists"""
function Print(P,io::IO = stdout)
  if(string(typeof(P))=="PauliList")
    PrintPauliList(P,io)
    return;
  end
  if( string(typeof(P))=="Array{PauliList,1}")
    PrintPauliLists(P,io)
    return;
  end
  if( string(typeof(P))=="Vector{PauliList}" )
    PrintPauliLists(P,io)
    return;
  end

  #println(typeof(P))
  println(io,P)
  return;
end


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
# i,x,y,z,p_00,p_01,p_10,p_11 = 1,2,3,4,5,6,7,8
I   =1;
X   =2;
Y   =3;
Z   =4;
P00 =5;
P01 =6;
P10 =7;
P11 =8;

"""This function returns the full ALGEBRA stored in terms of Matrix IDs
e.g.

SA=GetSpinAlgebra()

SA[P01,X] = 5     # MatrixID P00 = 5
SA[X,Y]   = 4im   # MatrixID Z   = 4
"""
function GetSpinAlgebra()

   SPIN_ALGEBRA=zeros(Complex,8,8)
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
#%  PList        =  MakeLocalOperator(MatrixID, tensor_index, num_spins)
#%  PList        =  MultiplyPauliList(PauliL::PauliList,PauliR::PauliList)
#%  PLists       =  MultiplyPauliLists(PauliSumL,PauliSumR)
#%  PList        =  SimplifySum(PList)

"""Pulls all phases to the prefactor so that only real Matrix IDs remain in list"""
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

"""MakeLocalOperator(MatrixID, j, M) makes a PauliList with MatrixID at index j amongst M qubits"""
function MakeLocalOperator(MatrixID, j, M)
  list=ones(M);
  list[j]=MatrixID;
  return PauliList(1,list)
end

"""Multiply function for two PauliLists"""
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

"""Multiply function for two arrays of PauliLists"""
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
  return SimplifySum(out);
end


#=
"""
SimplifySum adds like list terms for simplification (S. Manski 2016, Whitfield 2019)
input: unsimplified PauliList
output: PauliList with combined like terms
"""
function Sum!(A,B)

  # the deep copy is done to disconnect the input to the addition from the original object
  # otherwise the simplify sum will change the input values
  arrayPL=deepcopy(arrayPLgiven);
  #THIS IS A BOTTLE NECK IN MEMORY AND SPEED

  arrayA=PauliList[];

  for i=1:length(arrayPL)
    if arrayPL[i].Prefactor!=0
      #check if any other terms match
      for j=i+1:length(arrayPL)

        if( arrayPL[i].List == arrayPL[j].List )
            arrayPL[i].Prefactor=arrayPL[i].Prefactor+arrayPL[j].Prefactor;
            arrayPL[j].Prefactor=0;
        end

      end

    end
    if(arrayPL[i].Prefactor!=0)
      push!(arrayA,arrayPL[i]);
    end
  end
  return arrayA
end

#doesn't modify B but will modify A
function PlusEqual!(A,B)
  arrayA=PauliList[];

  for i=1:length(arrayPL)
    if arrayPL[i].Prefactor!=0
      #check if any other terms match
      for j=i+1:length(arrayPL)

        if( arrayPL[i].List == arrayPL[j].List )
            arrayPL[i].Prefactor=arrayPL[i].Prefactor+arrayPL[j].Prefactor;
            arrayPL[j].Prefactor=0;
        end

      end

    end
    if(arrayPL[i].Prefactor!=0)
      push!(arrayA,arrayPL[i]);
    end
  end
  return arrayA
end
=#
"""
Norm computes the L_1 coefficent norm of a PauliList or an array of PauliLists

The L_1 coefficent norm for \$h= \\sum a_i h_i\$  is given by \$|h|_c = \\sum |a_i|\$.
"""
function Norm(arrayPL::Array{PauliList,1})
  n = 0.0
  for i=1:length(arrayPL)
    n+=abs(arrayPL[i].Prefactor)
  end
  return n
end

function Norm(PL::PauliList)
  return abs(PL.Prefactor)
end

"""
SimplifySum adds like list terms for simplification (S. Manski 2016)
input: unsimplified PauliList
output: PauliList with combined like terms

The deepcopy of this function goes against the idea of Julia as a pass-by-
reference language.
"""
function SimplifySum(arrayPL::Array{PauliList,1})

  # the deep copy is done to disconnect the input to the addition from the original object
  # otherwise the simplify sum will change the input values

  #ref=deepcopy(arrayPL);

  arrayA=PauliList[];

  #to avoid copying and editing the incoming list we'll add like terms and then
  #when we combine them into the new array, we just the later summand altogether
  #effectively zeroing the value without modifying the array.
  skips=Int[]

  for i=1:length(arrayPL)
    #reset flag
    flag=0

    if i in skips || abs(arrayPL[i].Prefactor) < eps()
      continue
    end

    #check if any other terms match
    for j = i+1 : length(arrayPL)
      if arrayPL[i].List == arrayPL[j].List
        push!(skips,j)
        flag=1
        break
      end
    end

    if flag == 0

      push!(arrayA, arrayPL[i])

    else

      factor=arrayPL[i].Prefactor+arrayPL[skips[end]].Prefactor;
      if abs(factor)>eps()
        push!(arrayA, PauliList(factor,arrayPL[i].List))
      end

    end

    #=
    for k=1:length(ref)
      if (ref[k].List != arrayPL[k].List) || (ref[k].Prefactor!=arrayPL[k].Prefactor)
        println("warning")
        print("arraypl\n")
        println(arrayPL)
        print("ref\n")
        println(ref)
      end
    end
    =#

  end

  return arrayA
end

"""
SimplifySum adds like list terms for simplification (S. Manski 2016)
input: unsimplified PauliList
output: PauliList with combined like terms

The deepcopy of this function goes against the idea of Julia as a pass-by-
reference language.
"""
function SimplifySum_slow(arrayPLgiven)

  # the deep copy is done to disconnect the input to the addition from the original object
  # otherwise the simplify sum will change the input values
  arrayPL=deepcopy(arrayPLgiven);

  arrayA=PauliList[];

  for i=1:length(arrayPL)
    if arrayPL[i].Prefactor!=0

      #check if any other terms match
      for j=i+1:length(arrayPL)

        if( arrayPL[i].List == arrayPL[j].List )
            arrayPL[i].Prefactor=arrayPL[i].Prefactor+arrayPL[j].Prefactor;
            arrayPL[j].Prefactor=0;
        end
      end

    end
    if(arrayPL[i].Prefactor!=0)
      #push! is equivalent to Python's append. append! is equivalent to Python's extend.
      push!(arrayA,arrayPL[i])
    end
  end

  #=
  if arrayA == SimplifySum_slow(arrayPLgiven)
    println("Failure.\nSimplifySum\n")
    println(arrayA)
    println("SimplifySum2\n")
    println(SimplifySum_slow(arrayPLgiven))
  end
  =#

  return arrayA
end

"Returns maximum locality of an operator in PauliLists"
function Locality(arrayPL)
  if(string(typeof(arrayPL))=="PauliList")
    arrayPL=[arrayPL];
  end

  locality=-9;

  for k = 1 : length(arrayPL)
    pl=arrayPL[k];
    localityk=sum(pl.List .!=1)

    if(localityk > locality)
      locality=localityk
    end
  end

  return locality
end

# Wrapper functions
#  ans    =  Multiply(A,B)

"""
This function multiplies scalars, PauliLists, and vectors of PauliList without changing inputs
"""
function Multiply(Agiven,Bgiven)
  # the deep copy is done to disconnect the input to the addition from the original object
  # otherwise the simplify sum will change the input values
  A=deepcopy(Agiven)
  B=deepcopy(Bgiven)

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

  if( string(typeof(A))== "PauliList" && (     string(typeof(B))== "Array{PauliList,1}"  ||      string(typeof(B))== "Vector{PauliList}" ) )
		return MultiplyPauliList([A],B);
	end

  if( (     string(typeof(A))== "Array{PauliList,1}"  ||      string(typeof(A))== "Vector{PauliList}" ) && string(typeof(B))== "PauliList")
		return MultiplyPauliList(A,[B]);
	end

  if(string(typeof(A))== "PauliList" && string(typeof(B))== "PauliList")
		return MultiplyPauliList(A,B);
	end

	if(   (     string(typeof(A))== "Array{PauliList,1}"  ||      string(typeof(A))== "Vector{PauliList}" )
           && (     string(typeof(B))== "Array{PauliList,1}"  ||      string(typeof(B))== "Vector{PauliList}" ))
		return MultiplyPauliLists(A,B);
	end

	#println("Types: ",typeof(A)," , ",typeof(B))
	return A*B
end

#overloading operators and the print function
import Base.*
Base.:*(A::PauliList, B::PauliList) = Multiply(A,B);
Base.:*(A::Number, B::PauliList) = Multiply(A,B);
Base.:*(A::PauliList, B::Number) = Multiply(A,B);
Base.:*(A::Vector{PauliList}, B::Vector{PauliList}) = MultiplyPauliLists(A,B);
Base.:*(A::Array{PauliList,1}, B::Array{PauliList,1}) = MultiplyPauliLists(A,B);
Base.:*(A::Array{PauliList,1}, B::Vector{PauliList}) = MultiplyPauliLists(A,B);
Base.:*(A::Vector{PauliList}, B::Array{PauliList,1}) = MultiplyPauliLists(A,B);
Base.:*(A::Vector{PauliList}, B::Vector{PauliList}) = MultiplyPauliLists(A,B);

import Base.+
Base.:+(A::PauliList, B::PauliList) = SimplifySum([A;B])
Base.:+(A::Array{PauliList,1}, B::PauliList) = SimplifySum([A;B])
Base.:+(A::Array{PauliList,1}, B::Array{PauliList,1}) = SimplifySum([A;B])
Base.:+(A::PauliList, B::Array{PauliList,1}) = SimplifySum([A;B])
Base.:+(A::Vector{PauliList}, B::PauliList) = SimplifySum([A;B])
Base.:+(A::PauliList, B::Vector{PauliList}) = SimplifySum([A;B])

import Base.-
Base.:-(A::PauliList, B::PauliList) = SimplifySum(A+(-1)*B)
Base.:-(A::Vector{PauliList}, B::PauliList) = SimplifySum(A+(-1)*B)
Base.:-(A::PauliList, B::Vector{PauliList}) = SimplifySum(A+(-1)*B)
Base.:-(A::Array{PauliList,1}, B::Array{PauliList,1}) = SimplifySum(A+(-1)*B)
Base.:-(A::Array{PauliList,1}, B::PauliList) = SimplifySum(A+(-1)*B)
Base.:-(A::PauliList, B::Array{PauliList,1}) = SimplifySum(A+(-1)*B)

import Base.print
Base.:print(A::PauliList) = Print(A)
Base.:print(A::Vector{PauliList}) = Print(A)
Base.:print(A::Array{PauliList,1}) = Print(A)

import Base.println
Base.:println(A::PauliList) = Print(A) #+print("\n")
Base.:println(A::Vector{PauliList}) = Print(A) #+print("\n")
Base.:println(A::Array{PauliList,1}) = Print(A) #+print("\n")
