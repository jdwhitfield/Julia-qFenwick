######################
#
#FENWICK+TREES.jl 
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
temp = Base.isdefined(Base,:PTree)
if( temp == false ) 
  mutable struct PTree
      M::Int64
      Parents::Array{Int64}
  end
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
#% Tensor algebra and Spin algebra  #
#%###################################


#%       PrintPauliList(PauliList)
#%       PrintPauliLists(PauliLists)

#Pauli list
temp = Base.isdefined(Base,:PauliList)
if( temp == false ) 
  mutable struct PauliList
    Prefactor::Complex
    List::Array{Complex}
  end
end

function PrintPauliList(p::PauliList)
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

include("Spin_Algebra.jl")

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
  ops[roots].=Z;

  #Z_children
  children=Children(j,PT)
  ops[children].=Z;

  #Z_cousins
  cousins=YoungerCousins(j,PT)
  ops[cousins].=Z;

  #X_j
  ops[j]=X;

  #X_ancestors
  elders=Ancestors(j,PT)
  ops[elders].=X;

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
  ops[roots].=Z;

  #Z_cousins
  cousins=YoungerCousins(j,PT)
  ops[cousins].=Z;

  #Y_j
  ops[j]=Y;

  #X_ancestors
  elders=Ancestors(j,PT)
  ops[elders].=X;

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



#HERE'S TEST CASES FROM https://arxiv.org/pdf/1208.5986
tree=BKTree(4);

println("eq. 68")
Print(Wij(1,2,tree))

println("eq. 69")
Print(Wij(3,4,tree))

println("eq. 70")
Print(Wij(1,4,tree))

println("eq. 74")
Print(Multiply(ad_op(1,tree),a_op(3,tree)))

println("eq. 76")
Print(
  Multiply(ad_op(1,tree),
   Multiply(a_op(3,tree),
    Multiply(ad_op(2,tree),
            a_op(4,tree)
            ))))

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
println("eq. 79")
Print(Sum([H1;H2]))
