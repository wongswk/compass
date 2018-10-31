#library(stringr)
#atom_abbr<-read.csv("table_of_atom_abbrev.csv", header=TRUE)
#pdb<-read.pdb("TR928truth.pdb")


calscco <-function(residn, vdangle, pdb){
  ####----Load in pdb file, extract necessary material, convert to data frame----
  pdb_atom<-pdb$atom
  pdb_atom_orig_df<-as.data.frame(pdb_atom)
  pdb_atom_df<-as.data.frame(pdb_atom)
  
  #####----Find out what atom type the resid is ----
  atyperows<-which(pdb_atom$resno == residn, arr.ind = TRUE)
  atyperow<-atyperows[1]
  atype<-pdb_atom[atyperow, 5]
  
  #####load in atomdep list which contains information on the names of the atoms, coordinates, 
  #bond angles, bond length, and which atoms are included in calculating the side chain 
  #source("atomdep_list.R")
  
  
  #####find one letter abbreviation for resid----
  place_a<-which(atom_abbr==atype, arr.ind=TRUE)
  r_place_a<-as.numeric(place_a[1:1])
  rabbr1<-as.character(atom_abbr[r_place_a, 3]) #resid abbreviation 
  
  #####go to atomdep and find which atoms are needed: names----
  xnew<-atomdeps[[rabbr1]]$names; xnew
  
  ####Run calco on each of the atoms----
  for(j in 1:length(xnew)){
    
    #find names of atoms 
    elements<-atomdeps[[rabbr1]]$matx[j,]
    
    #get coordinates from pdb file 
    #search in pdb_atom$eleno by the user imput number
    #locate resid 90 and N and then grab xyz for it 
    xpdbca1<-as.double(pdb_atom_df[pdb_atom_df$resno==residn & pdb_atom_df$elety==elements[1], c('x', 'y', 'z')])
    xpdbca2<-as.double(pdb_atom_df[pdb_atom_df$resno==residn & pdb_atom_df$elety==elements[2], c('x', 'y', 'z')])
    xpdbca3<-as.double(pdb_atom_df[pdb_atom_df$resno==residn & pdb_atom_df$elety==elements[3], c('x', 'y', 'z')])
    
    #turn grabbed coordinates into a matrix 
    xa_matrix<-(matrix(c(xpdbca1,xpdbca2, xpdbca3), nrow=3, byrow=TRUE))
    
    #grad bangle and length from atomdep_list 
    xbangle<-atomdeps[[rabbr1]]$bangle        
    xblength<-atomdeps[[rabbr1]]$blength      
    
    #run calco
    #sourceCpp("calco_final.cpp")
    result<-calCo(xa_matrix, xblength[j], xbangle[j], vdangle[j])
    
    #for next iteration, use the new coordinates from the prev
    #modify pdb file 
    pdb_atom_df[pdb_atom_df$resno==residn & pdb_atom_df$elety==xnew[j], c('x', 'y', 'z')]<-result
    
    
    print(result)
  }
}


#####Testing-----
# one<-c(115.674,  65.563,   9.602)
# two<-c(114.548,  66.114,  10.305 )
# three<-c(114.137,  65.222,  11.470 )
# four<-c(112.645,  65.804,  12.297)
# Torsion(one, two, three, four)
# test_fxn1<-side_chain(residn=90, vdangle=-175.2707, pdb)
# 
# #ser 185
# c1<-c(116.988,  69.666,  21.693)
# c2<-c(118.370,  69.294,  21.334)
# c3<-c(119.045,  68.675,  22.569)
# c4<-c( 118.076,  68.170,  23.480)
# Torsion(c1, c2, c3, c4)
# f3<-side_chain(resid=185, vdangle=22.51443, pdb)
# 
# 
# #multiple atoms
# #asp 154
# #atom 1: CG: N, CA, CB
# ma1N<-c(119.977,  55.051,  -0.562 )
# ma2CA<-c(119.120,  56.075,  -1.143)
# ma3CB<-c(118.397,  55.443,  -2.335)
# ma4CGorig<-c(117.567,  56.427,  -3.105) #result line 1
# Torsion(ma1N, ma2CA, ma3CB, ma4CGorig)
# #atom2: od1: CA, CB, CG
# od1<-c(117.680,  57.645,  -2.842) #result line 2 
# Torsion(ma2CA, ma3CB, ma4CGorig, od1)
# #atom3: od2: od1, cb, cg 
# od2<-c(116.790,  55.974,  -3.974) #result line 3 
# Torsion(od1, ma3CB, ma4CGorig, od2)
# f4<-side_chain(resid=154, vdangle=c(-174.9018, 9.66114, -179.5561), pdb)


