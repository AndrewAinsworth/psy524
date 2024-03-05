* Encoding: UTF-8.
TITLE Setting your working directory.
cd 'c:\temp'.

TITLE Open data file.
get file 'samplecancorr.sav'.

TITLE Run the correlations between the variables and save the corr matrix.

  CORRELATIONS
  /VARIABLES=ts tc bs bc
  /matrix out ('corr.sav').

MATRIX.
print /TITLE="Reading in the Matrix File and Defining the Pieces of the XX and YY matrices".

MGET  /file='corr.sav'/type=corr.
print CR. 
compute xx=CR(1:2,1:2).
print xx.
compute yy=CR(3:4,3:4).
print yy.
compute xy=CR(1:2,3:4).
print xy.
compute yx=CR(3:4,1:2).
print yx.

print /TITLE="Creating the R matrix from the XX and YY matrices".
compute r=(inv(yy))*yx*(inv(xx))*xy.
print r.

print /TITLE="Computing the Eigenvalues/Vectors from R".
call svd(r,eigvec,eigval,w).
print eigval.
print eigvec.
print w.

print /TITLE="Computing the Canonical Correlation".
compute cc=sqrt(eigval).
print cc. 

print /TITLE="Normalized matrix of YY eigenvectors".
call eigen(yy,eigvecy,eigvaly).
print GSCH(eigvec).

print /TITLE="The sqrt below is a special matrix sqrt that spss doesn't have a function for so here it is from R".
compute msqrtyy={0.99967373, 0.02554261; 0.02554261, 0.99967373}.
compute timsqrty=t(inv(msqrtyy)).
print timsqrty.

print /TITLE="B_y is found by multiplying timsqrty by a normalized (GSCH function) matrix of the eigenvalues".
compute Bhaty = GSCH(eigvec).
print Bhaty.
compute Bsuby = timsqrty * Bhaty.
print Bsuby.

get y/ variables = Zbs Zbc.
print y. 

compute yscores = y * Bsuby.
print yscores.

END MATRIX.