#%###################################
#%
#%  Tree structures
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

#   set expandtab
#   set tabstop=2
#   set shiftwidth=2

using Random
using GraphRecipes
using Plots

#%###################################
#%#           Parity trees          #
#%###################################

#% PTree(num_of_sites,Parents)
temp = isdefined(Main,:PTree)
if( temp == false )
  """M the number of nodes; Parents[j] is the parent of j or -1 if orphaned.
  Tree structure is enforced by constructor."""
  mutable struct PTree
      M::Int64
      Parents::Array{Int64,1}
      #How to check the parents as a constructor, the following didn't work
      #PTree(M,Parents) = new(M,CheckParents(Parents))
  end
end

"converts an array to a tree. Directly modifies given array"
function CheckParents!(parents::Array{Int,1})
  B=parents

  for k = 1:length(B)


    if ( B[k] > k || B[k] == -1)
      continue #everything fine
    end

    #can't be your own parent
    if B[k] == k
      B[k] = -1
      continue
    end

    #k   :  1   2   3   4   5
    #B(k): -1  -1   1   2   5

    #k   :  1   2   3   4   5   6   7   8   9  10
    #B(k): -1   3   9   6  -1   5   9  10  -1   1

    if B[k]<k
      proposed_child = k
      proposed_parent= B[k]

      if B[proposed_parent] == -1 #then we can swap them

        B[min(proposed_parent,proposed_child)]=max(proposed_parent,proposed_child)
        B[max(proposed_parent,proposed_child)]=-1
        continue

      else #can't swap orientation
           #instead swap direction
        delta = k - B[k];
        p = min(k+delta, length(B));
        B[k] = p;
      end
      continue
    end

  end
end

"Converts a beta matrix to a parent array"
function BetaToParents(B::Array{})
  # Seeley12 defined everything backward such that B[1,1] corresponds to
  # occupancy n_M being included in node x_M. This makes the indexing
  # complicated.  We have indexed the columns and rows of matrix B using
  # separate variables. We count from bottom which is why row is computed
  # using -j and the column is computed using -nodek.

  M=size(B,1);

  #list of oprhans
  parents=-1*ones(Int,M);

  for nodek = 1:M
    col=M+1-nodek;
    for j=1:M-nodek
      row=M+1-nodek-j;
      if( B[row,col]==1 )
        parents[nodek]=j+nodek;
        break;
      end
    end
  end

  return parents


end

"Make sure the parent list is valid, sort node so that roots have highest index"
function CheckTree(pt::PTree)
  for k=1:pt.M

    #can't be your own parent
    if(pt.Parents[k]==k)
      pt.Parents[k]=-1;
    end

    if(pt.Parents[k]>k)
      k-pt.Parents[k]
    end

  end
end

#%Parity tree functions
#%--
#%function Ancestors(      j,PT::PTree)
#%function Children(       j,PT::PTree)
#%function YoungerCousins( j,PT::PTree)
#%function DisjointRoots(  j,PT::PTree)
#%function Progeny(        j,PT::PTree)

"The ancestors of a tree node are the collection of all parents of parents of the node."
function Ancestors(j::Int,PT::PTree,verbose::Bool=false)
    A=Int[];
    cur::Int=j;
    ctr=1;
    while ctr<1000
        if(verbose)
            println("ctr - ",ctr)
        end
        if(PT.Parents[cur] .== -1)
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

"Returns immediate child of a node"
function Children(j::Int,PT::PTree)
    C=Int[];
    for k=1:PT.M
        if(PT.Parents[k]==j)
            C=[C; k];
        end
    end
    return C;
end

"Return set of children's children"
function Progeny(j,tree::PTree)
  all_kids=Int[];

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
      print(all_kids)
      return;
    end
  end

  return all_kids;
end

"""We consider cousins as the children (no progeny) of all ancestors of a node.
Note that this differs from the familial definition of the term; for instance
this definition includes siblings and discludes second cousins. The set of
younger cousin is the collection of cousins with lower node index."""
function YoungerCousins(j,PT::PTree)
    #children of ancestors younger than j
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

"Returns nodes that have no parents and lower node index"
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
#% pt     = ShiftTree(shift, ptree)
#% pt     = CombineTrees(pt_L, pt_R)
#% parent_array   = ShiftTreeIndices(shift,parent_array)
#% parent_array   = CombineTreeIndices(left_parent_array,right_parent_array)
#% matrix = BetaMatrix(tree)

"Shift the tree indicies of an array, while preserving the orphan flag"
function ShiftTreeIndices(shift,parent_array)

  shifted_tree = parent_array.+shift;

  #if it used to be empty it needs to reset
  shifted_tree[shifted_tree .== shift-1] .=-1;
  return shifted_tree
end

"Shifts all indicies of a parity tree by a constant offset"
function ShiftTree(shift,pt::PTree)
  treeL = PTree(shift,-ones(Int,shift))
  return CombineTrees(treeL,pt)
end

"Combines two parity tree"
function CombineTrees(treeL::PTree,treeR::PTree)

  M=treeL.M+treeR.M;

  return PTree(M, CombineTreeIndices(treeL.Parents,treeR.Parents))
end

"Combine two parity tree arrays"
function CombineTreeIndices(treeL,treeR)
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
#% PT = RandomTree(M)

function MakeFenwickTree(N::Int)
    #clean input
    if(N==0)
      return [];
    end

    N=convert(Int,N);



    #global variables with respect to FENWICK
    tree=-ones(Int,N)
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

function JWTree(M::Int)
    Parents=-ones(Int, M)
    return PTree(M,Parents)
end

function BKTree(M::Int)
  return PTree(M,MakeFenwickTree(M))
end

function SBKTree(M::Int,L::Int)
  #L=stride_length;

  #how many strides?
  n=floor(Int,M/L);

  #remainding nodes
  r=convert(Int, M-n*L);

  sktree=[];
  #Now we just need to make n strides
  for i=1:n
    sktree=CombineTreeIndices(MakeFenwickTree(L),sktree);
  end
  #plus the remaining nodes

  sktree=CombineTreeIndices(MakeFenwickTree(r),sktree);

  return PTree(M,sktree);
end

function ParityTree(M)

  paritytree=-ones(Int,M);

  for i=1:(M-1)
    paritytree[i]=i+1;
  end
  
  return PTree(M,paritytree);
end

"""
RandomTree(M::Int,[p_orphan,seed_val])

Generate a random binary tree structure with M nodes
There is a optional second parameter controlling the probability of a node
being an oprhan.  The third optional parameter is the seed value for the
random number generator.
"""
function RandomTree(M; p_orphan=-9,seed_val=-9)

  if seed_val != -9
    Random.seed!(seed_val)
  end

  if p_orphan<0
    p_orphan=rand()
  end

  tree=-ones(Int,M);

  for i=1:(M-1)
    if(rand() < p_orphan)
      continue
    end

    #pick any elder
    tree[i]=rand(i+1:M)

  end

  return PTree(M,tree);
end

#####################################

function PlotTree(pt::PTree)
  src=Int[]
  tgt=Int[]

  for m=1:pt.M
    if(pt.Parents[m]>0)
      src=[src;m]
      tgt=[tgt;convert(Int,pt.Parents[m])]
    else
      src=[src;m]
      tgt=[tgt;m]
    end
  end
  graphplot(src,tgt, names=1:pt.M, method=:arcdiagram, nodesize=min(3,pt.M))
end
