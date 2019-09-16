######################
#
#Fermions.jl
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
#%  Get Spin algebra and Tree lib   #
#%###################################

include("Spin_Algebra.jl")
include("Trees.jl")

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
#% PS =  Tij( idx_i,idx_j   ,parity_tree)


#basic definitions
"""
c_j = a_j + a_j^+
c_j = Z_{DisjointRoots} Z_{Children} Z_{YCousins} X_j X_{Ancestors}
"""
function c_op(j,PT::PTree)
  #

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

  #ops[roots]=Z;
  # FAILS WITH ERROR
  #
  """MethodError: no method matching setindex_shape_check(::Int64, ::Int64)
  Closest candidates are:
     setindex_shape_check(!Matched::AbstractArray{#s57,1} where #s57, ::Integer) at indices.jl:193
     setindex_shape_check(!Matched::AbstractArray{#s57,1} where #s57, ::Integer, !Matched::Integer) at indices.jl:196
     setindex_shape_check(!Matched::AbstractArray{#s57,2} where #s57, ::Integer, !Matched::Integer) at indices.jl:200
     ...

  Stacktrace:
    [1] macro expansion at ./multidimensional.jl:641 [inlined]
    [2] _unsafe_setindex!(::IndexLinear, ::Array{Complex,1}, ::Int64, ::Array{Any,1}) at ./multidimensional.jl:636
    [3] c_op(::Int64, ::PTree) at ./multidimensional.jl:631
    [4] top-level scope at In[59]:1"""
  #do it by hand
  #for k=1:length(roots)
  #   ops[roots[k]]=Z;
  #end

  ops[roots].=Z;

  """#Z_{DisjointRoots}
  roots=DisjointRoots(j,PT)
  ops[roots].=Z;"""

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
  ops[elders].=X

  return PauliList(1,ops)
end

"d_j = -i (a_j - a_j^+)"
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

"a_j following fermion algebra"
function a_op(j,f::PTree)
  return Multiply(.5,[c_op(j,f),Multiply(1im,d_op(j,f))])
end

"a_j^+ following fermion algebra"
function ad_op(j,f::PTree)
  return Multiply(.5,[c_op(j,f),Multiply(-1im,d_op(j,f))])
end

#One body terms
"""
n_j = a_j^+ a_j
n_j = a_j^+ a_j = (1 + i c_jd_j) /2
"""
function nj(j,f::PTree)
  #

  I=PauliList((1/2.),ones(f.M))
  cd=MultiplyPauliList(c_op(j,f),d_op(j,f))
  icd=PauliList((1im/2.)*cd.Prefactor,cd.List)
  return [I,icd]
end


"t_ij = a_i^+ a_j + a_j^+ a_i = i (c_i d_j+c_j d_i)/2"
function Tij(i,j,f::PTree)

  factor=1im/2;

  term1_PL= MultiplyPauliList(c_op(i,f),d_op(j,f))
  term1_PL.Prefactor=term1_PL.Prefactor*factor;
  # Print(c_op(i,f))
  # Print(d_op(j,f))
  # Print(term1_PL)

  term2_PL= MultiplyPauliList(c_op(j,f),d_op(i,f))
  term2_PL.Prefactor=term2_PL.Prefactor*factor;

  return SimplifySum([term1_PL,term2_PL])
end

#Two body terms

""" W_ij = a_n^+  a_m^+ a_m a_n
      =  ( 1 + i c_m d_m)( 1 + i c_n d_n) / 4"""
function Wij(j,k,f::PTree)
  I=PauliList((1/4.),ones(f.M))

  icd_j= Multiply(.25im,Multiply(c_op(j,f),d_op(j,f)))
  icd_k= Multiply(.25im,Multiply(c_op(k,f),d_op(k,f)))


  return SimplifySum([I;icd_j;icd_k;Multiply(4,Multiply(icd_j,icd_k))])
end

"W_ijkl  = a_i^+ a_j^+ a_k a_l + a_l^+ a_k^+ a_j a_i"
function Wijkl(i,j,k,l,tree::PTree)
   if(i==l && j==k)
     return Wij(i,j,tree);
   end
   return SimplifySum([Multiply(ad_op(i,tree),Multiply(ad_op(j,tree),Multiply(a_op(k,tree),a_op(l,tree))));Multiply(ad_op(l,tree),Multiply(ad_op(k,tree),Multiply(a_op(j,tree),a_op(i,tree))))])

end

function Hubbard(L::Int, tree ;  snake=false, dim=2,  periodic=false,printlayout=false)
  H= PauliList[]
  idx(i::Int,j::Int)=-9


  if !snake
    zigzagidx(i::Int,j::Int) = i+(j-1)*L;
    idx=zigzagidx
  else
    snakeidx(i::Int,j::Int) = j%2 * ( i+ (j-1)*L ) + (j+1)%2 * ( j*L -i +1 );
    idx=snakeidx
  end

  for i=1:L
    for j=1:L
      if printlayout
        print(idx(i,j)," ")
      end

      if i<L
        #print("connecting ",i," ",j," to ",i+1," ",j,"\n")
        H+=Tij(idx(i,j),idx(i+1,j),tree);
      elseif periodic
        #print("connecting ",i," ",j," to ", 1 ," ",j,"\n")
        H+=Tij(idx(i,j),idx(1,j),tree);
      end

      if j<L
        #print("connecting ",i," ",j," to ",i," ",j+1,"\n")
        H+=Tij(idx(i,j),idx(i,j+1),tree);
      elseif periodic
        #print("connecting ",i," ",j," to ", i ," ",1,"\n")
        H+=Tij(idx(i,j),idx(i,1),tree);
      end

      #print("onsite interaction at ",i," ",j,"\n")
      H+=( 2*(rand(0:1)-.5) )*  nj(idx(i,j),tree);

      #print(j%2 * ( i+ (j-1)*L ) + (j+1)%2 * ( j*L -i +1 ),"\t")
      #print(i+(j-1)*L,"\t")
    end
    if printlayout
      print("\n")
    end
  end
  return H
end