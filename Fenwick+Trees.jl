######################
#
#FENWICK+TREE.jl 
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
#%#           Parity trees          #
#%###################################
#%ParityTree(num_of_sites,Parents)
struct PTree
    M::Int64
    Parents::Array{Int64}
end

#%Parity tree functions
#%--
#%function Ancestors(      j,PT::PTree)
#%function Children(       j,PT::PTree)
#%function YoungerCousins( j,PT::PTree)
#%function DisjointRoots(  j,PT::PTree)
#%function Progeny(        j,PT::PTree)

function Ancestors(j,PT::PTree,verbose::Bool=false)
    A=[];
    cur=j;
    ctr=1;
    while ctr<1000
        if(verbose)
            println("ctr - ",ctr)
        end
        if(PT.Parents[cur]==-1)
            if(verbose)
                println("breaking val of cur=",cur)
            end
            break;
        else
            A=[A; PT.Parents[cur]];
            cur=PT.Parents[cur];
        end
        
        ctr=ctr+1
    
    end
    if ctr>PT.M
        if(verbose)
            println("ctr - ",ctr,"; PT.M - ",PT.M)
        end
        print("Warning, the loop didn't close automatically")
    end
    
    return A;
end

function Children(j,PT::PTree)
    C=[];
    for k=1:PT.M
        if(PT.Parents[k]==j)
            C=[C; k];
        end
    end
    return C;
end

function Progeny(j,tree::PTree)
  #children's children
  all_kids=[];

  ctr=0;
  g0kids=Children(j,tree);


  while(length(g0kids)>0)
    #Basis cases to illustrate the logic
    #
    # kids=Children(j)
    # grandkids=[];
    # for i=1:length(kids)
    #   child=kids[i];
    #   grandkids=[grandkids Children(child)];
    # end
    #
    # greatgrandkids=[];
    # for i=1:length(grandkids)
    #   grandchild=grandkids[i];
    #   greatgrandkids=[greatgrandkids Children(grandchild)];
    # end

    #g0 is the previous generation, g1 is the next generation
    g1kids=[];
    for i=1:length(g0kids)

      #pull child
      child0=g0kids[i];

      g1kids=append!(g1kids,Children(child0,tree));
    end

    #record this generation
    all_kids=append!(all_kids,g0kids);

    #check the next generation
    g0kids=g1kids;

    ctr=ctr+1;
    #prevent excessive nesting
    if(ctr>100)
      print("\n\tWarning: Excessive nesting, pruning\n")
      return;
    end
  end

  return all_kids;
end

function YoungerCousins(j,PT::PTree)
    #children with same common ancestor are cousins
    K=[];
    A=Ancestors(j,PT)
    for k=1:length(A)
        Ck=Children(A[k],PT);
        for kk=1:length(Ck)
            if(Ck[kk]<j)
                K=[K; Ck[kk]]
            end
        end
    end
    return K;
end

function DisjointRoots(j,PT::PTree)
    L=[];

    for k=j-1:-1:1
        if PT.Parents[k]==-1;
            L=[L; k];
        end
    end

    return L
end

#%###################################
#%       Parity tree structure      #
#%###################################
#%Generate JW, BK, and SBK. This is predicated on the manipulation of trees.
#%The trees are store such that each node has at most one parent. Their parent of
#%the kth node is stored in index [k] of a linear tree list.  If the node has no
#%parent, then a marker is stored instead of listing a parent.
#%
#% TREE MANIPULATION
#%--
#% tree   = ShiftTreeIndices(shift,tree)
#% tree   = CombineTrees(left_tree,right_tree)
#% matrix = BetaMatrix(tree)



function ShiftTreeIndices(shift,tree)

  shifted_tree=tree.+shift;

  #if it used to be empty it needs to reset
  shifted_tree[shifted_tree .== shift-1] .=-1;
  return shifted_tree
end

function CombineTrees(treeL,treeR)
  m=length(treeL)
  shifted_treeR=ShiftTreeIndices(m,treeR)
  return [treeL;shifted_treeR]
end

function BetaMatrix(PT::PTree)

  #make MxM zeros matrix
  B=zeros(Int64,PT.M,PT.M)

  #for each node in tree
  for j=1:PT.M

    #get progeny list
    pj=Progeny(j,PT)

    #mark appropriate matrix elements

    #use the state ordering from eq (24) of Seeley12

    #equation 24 of peter's paper has a strange ordering
    #based on b_i=\sum B_{ij} f_j instead of the usual change of basis formula
    # b_i=\sum B_{ji} f_j

    B[j,j]=1;
    for k=1:length(pj)
      B[PT.M+1-j,PT.M+1-pj[k]]=1;
    end

  end

  #return matrix
  return B;
end

#%######################################
#% Fenwick tree and related structures #
#%######################################
#%--
#% tree = MakeFenwickTree(M)
#% PT = JWTree(M)
#% PT = BKTree(M)
#% PT = SBKTree(M,stride_length)
#% PT = ParityTree(M)

function MakeFenwickTree(N)

    if(N==0)
      return [];
    end
    #global variables with respect to FENWICK
    tree=-ones(N,1)
    ctr=0;

    function FENWICK(L,R)
        #The recursive algorithm Vojta suggested in our paper


        if(ctr>100)
            #prevent excessive nesting
            print("\n\tWarning: Excessive nesting, pruning\n")
            return;
        end

        #algorithm is to connect from middle between L & R to right then recuse
        #on each half. Be careful about floor versus ceiling.

        if(L==R)            #if L=R then there is no middle, so we're done
            return;
        end
        n=Int64(floor((L+R)/2))

        #Connect middle to R
        #+1 because we count from 1 in this language (Julia)
        tree[n+1]=R+1;

        #branch
        ctr=ctr+1;
        FENWICK(L+0.0, n)
        FENWICK(1.0+n, R)

    end

    #run the recursive algorithm on (0,N-1) and it generates the BK tree over N
    #nodes.
    FENWICK(0.0,N-1.0)
    return tree
end

function JWTree(M)
    Parents=-ones(M,1)
    return PTree(M,Parents)
end

function BKTree(M)
  return PTree(M,MakeFenwickTree(M))
end

function SBKTree(M,L)
  #L=stride_length;

  #how many strides?
  n=M/L;
  n=n-n%1;
  #remainder
  r=M-n*L;

  sktree=[];
  #Now we just need to make n strides
  for i=1:n
    sktree=CombineTrees(MakeFenwickTree(L),sktree);
  end
  #plus the remaining nodes

  sktree=CombineTrees(MakeFenwickTree(r),sktree);

  return PTree(M,sktree);
end

function ParityTree(M)

  paritytree=-ones(M,1);

  for i=1:(M-1)
    paritytree[i]=i+1;
  end

  return PTree(M,paritytree);
end


#unit testing
pt=ParityTree(3)
ct=CombineTrees(pt.Parents,pt.Parents)
@assert Array([2;3;-1;5;6;-1])≈ Array(ct) "unit test for CombineTrees failed"


#%###################################
#%           Spin algebra           #
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
#%       PrintPauliList(PauliList)
#%       PrintPauliLists(PauliLists)
#%  PList        =  Sum(PList)

#Pauli list
temp = @isdefined PauliList
if( temp == false ) 
  struct PauliList
    Prefactor::Complex
    List::Array{Complex}
  end
end

function StandardForm(P::PauliList)

  #a little error catching to make sure MatrixIDs are okay
  if( 0!=sum(abs(P.List).>8.0) || 0!=sum(abs(P.List).<1.0))
    println("Warning invalid MatrixIDs in StandardForm function")
  end

  #pull StandardForm
  b=(P.List)./abs(P.List);

  #if entry of list is zero extracting the phase via A/|A| is invalid. Hence,
  ## !isnan ## command. All MatrixIDs should be 1,...,8 so this isn't necessary.
  b=prod( b[!isnan(b)] );

  P.Prefactor=P.Prefactor*b;

  #put everything as standard matrix IDs
  P.List=Complex.(round(abs(P.List)));

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

function PrintPauliList(P::PauliList)
  #correct the state
  P=StandardForm(P)

  #These are the MatrixIDs for coding the algebra
  I   =1;
  X   =2;
  Y   =3;
  Z   =4;
  P00 =5;
  P01 =6;
  P10 =7;
  P11 =8;

  #printing
  print("(")
  if(real(P.Prefactor)>=0) #add plus sign to make it look better
    print("+")
  end
  print(P.Prefactor)
  print(") ")
  for i=1:length(P.List)
    if(P.List[i]==I)
      print("  I  ")
    end
    if(P.List[i]==X)
      print("  X  ")
    end
    if(P.List[i]==Y)
      print("  Y  ")
    end
    if(P.List[i]==Z)
      print("  Z  ")
    end
    if(P.List[i]==P00)
      print(" P00 ")
    end
    if(P.List[i]==P01)
      print(" P01 ")
    end
    if(P.List[i]==P10)
      print(" P10 ")
    end
    if(P.List[i]==P11)
      print(" P11 ")
    end
  end
  println("")
  return
end

function PrintPauliLists(arrayPL)

  if(string(typeof(arrayPL))=="PauliList")
    arrayPL=[arrayPL];
  end

  if(string(typeof(arrayPL[1]))!="PauliList")
    println("Warning: bad input; expected array of PauliLists in function Print_sum")
  end

  for k=1:length(arrayPL)
    if(k>1)
      print("    +")
    else
      print("op = ")
    end
    Print(arrayPL[k])
  end
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

function Print(P)
  if(string(typeof(P))=="PauliList")
    PrintPauliList(P)
    return;
  end
  if( string(typeof(P))=="Array{PauliList,1}")
    PrintPauliLists(P)
    return;
  end
  if( string(typeof(P))=="Vector{PauliList}" )
    PrintPauliLists(P)
    return;
  end
  println(typeof(P))
  println(P)
  return;
end




#%###################################
#%       Fermionic operators        #
#%###################################
#%Given a Fenwick tree structure, generate operators. To simplify, we use the
#%Majorana basis which aligns better with the Pauli basis of spin-\frac12 states.
#%
#%FUNCTIONS:
#% PL = c_op(     idx       ,parity_tree)
#% PL = d_op(     idx       ,parity_tree)
#% PL = a_op(     idx       ,parity_tree)
#% PL = ad_op(    idx       ,parity_tree)
#% PS =   nj(     idx       ,parity_tree)
#% PS =  tij( idx_i,idx_j   ,parity_tree)

#basic definitions
function c_op(j,PT::PTree)
  # c_j = Z_{DisjointRoots} Z_{Children} Z_{YCousins} X_j X_{Ancestors}

  #These are the MatrixIDs for coding the algebra
  I   =1;
  X   =2;
  Y   =3;
  Z   =4;
  P00 =5;
  P01 =6;
  P10 =7;
  P11 =8;
  ops=ones(Complex,PT.M);

  #Z_{DisjointRoots}
  roots=DisjointRoots(j,PT)
  ops[roots]=Z;

  #Z_children
  children=Children(j,PT)
  ops[children]=Z;

  #Z_cousins
  cousins=YoungerCousins(j,PT)
  ops[cousins]=Z;

  #X_j
  ops[j]=X;

  #X_ancestors
  elders=Ancestors(j,PT)
  ops[elders]=X;

  return PauliList(1,ops)
end

function d_op(j,PT::PTree)
  # d_j = Z_{DisjointRoots}  Z_{YCousins} Y_j X_{Ancestors}

  #These are the MatrixIDs for coding the algebra
  #Would be nice to have these as globals e.g. in a header file
  I   =1;
  X   =2;
  Y   =3;
  Z   =4;
  P00 =5;
  P01 =6;
  P10 =7;
  P11 =8;

  ops=ones(Complex,PT.M);

  #Z_{DisjointRoots}
  roots=DisjointRoots(j,PT)
  ops[roots]=Z;

  #Z_cousins
  cousins=YoungerCousins(j,PT)
  ops[cousins]=Z;

  #Y_j
  ops[j]=Y;

  #X_ancestors
  elders=Ancestors(j,PT)
  ops[elders]=X;

  return PauliList(1,ops)
end

function a_op(j,f::PTree)
  return Multiply(.5,[c_op(j,f),Multiply(1im,d_op(j,f))])
end

function ad_op(j,f::PTree)
  return Multiply(.5,[c_op(j,f),Multiply(-1im,d_op(j,f))])
end

#One body terms

#n_j = a_j^\dag a_j =
function nj(j,f::PTree)
  #a_j^\dag a_j&              =\frac12 (\mathbf 1 + i c_jd_j)

  I=PauliList((1/2.),ones(f.M))
  cd=MultiplyPauliList(c_op(j,f),d_op(j,f))
  icd=PauliList((1im/2.)*cd.Prefactor,cd.List)
  return [I,icd]
end

tree=BKTree(4);

# t_ij = a_i^\dag a_j&+a_j^\dag a_i =\frac i2 (c_id_j+c_jd_i)
function Tij(i,j,f::PTree)

  factor=1im/2;

  term1_PL= MultiplyPauliList(c_op(i,f),d_op(j,f))
  term1_PL.Prefactor=term1_PL.Prefactor*factor;
  # Print(c_op(i,f))
  # Print(d_op(j,f))
  # Print(term1_PL)

  term2_PL= MultiplyPauliList(c_op(j,f),d_op(i,f))
  term2_PL.Prefactor=term2_PL.Prefactor*factor;

  return [term1_PL,term2_PL]
end

#Two body terms

# W_ij = a_n^\dag  a_m^\dag a_m a_n
#      = \frac14 (\mathbf 1 + i c_md_m)(\mathbf 1 + i c_nd_n)
function Wij(j,k,f::PTree)
  I=PauliList((1/4.),ones(f.M))

  icd_j= Multiply(.25im,Multiply(c_op(j,f),d_op(j,f)))
  icd_k= Multiply(.25im,Multiply(c_op(k,f),d_op(k,f)))


  return [I;icd_j;icd_k;Multiply(4,Multiply(icd_j,icd_k))]
end

#M_ijkl  = a_i^\dag a_j^\dag a_k a_l + a_l^\dag a_k^\dag a_j a_i
function Wijkl(i,j,k,l,tree)
   if(i==l && j==k)
     return Wij(i,j,tree);
   end
   return Sum([Multiply(ad_op(i,tree),Multiply(ad_op(j,tree),Multiply(a_op(k,tree),a_op(l,tree))));Multiply(ad_op(l,tree),Multiply(ad_op(k,tree),Multiply(a_op(j,tree),a_op(i,tree))))])

end



#HERE'S TEST CASES FROM arXiv:1208.5986
tree=BKTree(4);

println("eq. 68")
Print(Wij(1,2,tree))

println("eq. 69")
Print(Wij(3,4,tree))

println("eq. 70")
Print(Wij(1,4,tree))

println("eq. 74")
Print(Tij(1,3,tree))

println("eq. 76")
Print(Tij(2,4,tree))

println("eq. 77")
Print(Wijkl(1,4,2,3,tree))

println("eq. 78")
Print(Wijkl(1,2,4,3,tree))

#HYDROGEN INTEGRALS AT SEPARATION @ 1.401 Bohr
h00=-1.252477;
h11=-1.252477;

h22=-0.475934;
h33=-0.475934;

h0110=0.674493;
h1001=0.674493;

h2332=0.697397;
h3223=0.697397;

h0220=0.663472;
h0330=0.663472;
h1221=0.663472;
h1331=0.663472;
h2002=0.663472;
h3003=0.663472;
h2112=0.663472;
h3113=0.663472;

h0202=0.181287;
h1313=0.181287;
h2130=0.181287;
h2310=0.181287;
h0312=0.181287;
h0132=0.181287;


#one-body Hamiltonian
H1=Sum([Multiply(h00,nj(1,tree));
        Multiply(h11,nj(2,tree));
        Multiply(h22,nj(3,tree));
        Multiply(h33,nj(4,tree))]);
#two-body Hamiltonian
 H2=Sum([
   Multiply(h0110,Wijkl(1,2,2,1,tree));
   Multiply(h2332,Wijkl(3,4,4,3,tree));
   Multiply(h0330,Wijkl(1,4,4,1,tree));
   Multiply(h1221,Wijkl(2,3,3,2,tree));
   Multiply((h0220-h0202),Wijkl(1,3,3,1,tree));
   Multiply((h1331-h1313),Wijkl(2,4,4,2,tree));
   Multiply((h0132),Wijkl(1,2,4,3,tree));
   Multiply(h0312,Wijkl(1,4,2,3,tree))
   ])
println("eq. 76")
Print(Sum([H1;H2]))
