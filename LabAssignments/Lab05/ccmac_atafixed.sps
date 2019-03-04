preserve.                                                                       
set printback=off.                                                              
define cancorr (set1  =!charend('/')                                            
               /set2  =!charend('/')                                            
               /debug =!charend('/') !DEFAULT ('N')                             
               /KEEPSC=!charend('/') !DEFAULT ('Y') 
               /PRCOR =!charend('/') !DEFAULT (25 ) ).                                                                                                        
preserve.
!IF ( !DEBUG !EQ 'N') !THEN
set printback=off mprint off.                                                   
!ELSE
set printback on mprint on.
!IFEND .

*---------------------------------------------------------------------------.
* Save the original file for later retrieval.
*---------------------------------------------------------------------------.
!IF (!DEBUG !EQ 'N') !THEN 
SET RESULTS ON.
DO IF $CASENUM=1.
PRINT / "NOTE: ALL OUTPUT INCLUDING ERROR MESSAGES HAVE BEEN TEMPORARILY" 
      / "SUPPRESSED. IF YOU EXPERIENCE UNUSUAL BEHAVIOR THEN RERUN THIS" 
      / "MACRO WITH AN ADDITIONAL ARGUMENT /DEBUG='Y'."  
      / "BEFORE DOING THIS YOU SHOULD RESTORE YOUR DATA FILE."
      / "THIS WILL FACILITATE FURTHER DIAGNOSTICS OF ANY PROBLEMS".
END IF.
!IFEND .
save outfile='C:\Temp\cc__temp1.sav'.
DATASET NAME temp_copy.

*---------------------------------------------------------------------------.
* Compute the correlation matrix and pass information to MATRIX.
*---------------------------------------------------------------------------.

* DEFAULT:  SET RESULTS AND ERRORS OFF TO SUPPRESS CORRELATION PIVOT TABLE *.
!IF (!DEBUG='N') !THEN
set results off errors off.
!IFEND 
corr variables=!set1 !set2 /missing=listwise/matrix out(*).                                                                              
set errors on results listing.                                                                      

*---------------------------------------------------------------------------.
* Get correlations and compute basic quantities needed for analysis.
*---------------------------------------------------------------------------.
* SET command placed to prevent exceeding MXLOOPS 40 default * .                

SET MXLOOPS=199 MITERATE 199.
matrix.
get r /variables=!set1/file=*.
compute p1=ncol(r).
get r /file=* /names=varname/variables=!set1 !set2.
compute p2=ncol(r)-p1.
compute nx1=varname(1:p1).
compute nv=p1+p2.
compute nx2=varname((p1+1):nv).
compute rr=r(4:(nv+3),1:nv).
compute ns=r(3,1).
compute r11=rr(1:p1,1:p1).
compute r22=rr((p1+1):nv,(p1+1):nv).
compute r12=rr(1:p1,(p1+1):nv).
compute d1=r(2,1:p1).
compute d2=r(2,(p1+1):nv).
compute d1=mdiag(d1).
compute d2=mdiag(d2).
compute s1=d1*r11*d1.
compute s12=d1*r12*d2.
compute s2=d2*r22*d2.
compute d1=inv(d1).
compute d2=inv(d2).

*---------------------------------------------------------------------------.
* Print correlations.
*---------------------------------------------------------------------------.

do if (p1 <= !PRCOR and p2 <= !PRCOR ).
print r11 /format "f7.4"/title 'Correlations for Set-1'                         
          /space 2 /rnames=nx1 /cnames=nx1.
print r22 /format "f7.4"/title 'Correlations for Set-2'
          /space 2 /rnames=nx2 /cnames=nx2.
print r12 /format "f7.4"/title 'Correlations Between Set-1 and Set-2 '
         /space 2 /rnames=nx1 /cnames=nx2.
end if.  
                                                                             
*---------------------------------------------------------------------------.
* R1 and r2 are cholesky decomp of r11 and r22.
*---------------------------------------------------------------------------.

compute r1=chol(r11).
compute r2=chol(r22).

*---------------------------------------------------------------------------.
* R1_inv and r2_inv are inverse of r1 and r2.
*---------------------------------------------------------------------------.

compute r1_inv=inv(r1).
compute r2_inv=inv(r2).

*---------------------------------------------------------------------------.
* compute omega matrix.
*---------------------------------------------------------------------------.

do if (p1 le p2).
compute omega=t(r1_inv)*r12*r2_inv.
else.
compute omega=t(r2_inv)*t(r12)*r1_inv.
end if.

*---------------------------------------------------------------------------.
* SVD computes the singular value decomposition of omega.
*---------------------------------------------------------------------------.

call svd(omega,u,lambda,v).

*---------------------------------------------------------------------------.
* Create a list of names for use later in labels .
*---------------------------------------------------------------------------.

!LET !@=!NULL !LET !@1=!NULL !LET !@2=!NULL                                     
  !DO !N= 1 !TO 199                                                             
      !LET !@=!CONCAT(!@,!QUOTE(!N),",")                                        
      !LET !@1=!CONCAT(!@1,!QUOTE(!CONCAT('CV1-',!N)),",")                      
      !LET !@2=!CONCAT(!@2,!QUOTE(!CONCAT('CV2-',!N)),",")                      
  !DOEND                                                                        
  !LET !@=!CONCAT(!@,!QUOTE(@@))                                                
  !LET !@1=!CONCAT(!@1,!QUOTE(@@))                                              
  !LET !@2=!CONCAT(!@2,!QUOTE(@@)).                                                                               
Compute num={!@}.

*---------------------------------------------------------------------------.
* Lambda stores the canonical correlations. Print them now.
*---------------------------------------------------------------------------.

print diag(lambda)/format "f8.3"/title 'Canonical Correlations'                 
          /space 2/rnames=num.

compute dlam=diag(lambda).

*---------------------------------------------------------------------------.
* Compute the eigenvalues and test of remaining canonical correlations.
*---------------------------------------------------------------------------.

compute eign=(1 &/ (1-dlam &**2)) - 1.
compute wlam=1 &/ (1+eign).
compute n=nrow(wlam).
compute wilk=wlam.
compute df=wlam.
compute sig=wlam.
compute bart2=wlam.
compute tem=1.
loop  #l=1 to n.
+  compute tem=tem*wlam(n-#l+1).
+  compute df(n-#l+1)=(p1-n+#l)*(p2-n+#l).
+  compute dof=df(n-#l+1).
+  compute bart2(n-#l+1)=-(ns-0.5*(p1+p2+3))*ln(tem).
+  compute chi=bart2(n-#l+1).
+  compute sig(n-#l+1)=1-chicdf(chi,dof).
+  compute wilk(n-#l+1)=tem.
end loop.
compute test={wilk,bart2,df,sig}.
print test /format "f8.3"/title 'Test that remaining correlations are zero:'
  /space 2/rnames=num
  /cnames={"Wilk's ","Chi-SQ","  DF  ","  Sig."}.

*---------------------------------------------------------------------------.
* Compute and print the standardized canonical coefficients for set-1.
*---------------------------------------------------------------------------.

do if (p1 le p2).
compute a=r1_inv*u.
else.
compute a=r1_inv*v.
end if.
do if (p2 lt p1).
compute a=a(:,1:p2).
end if.
print a/format "f8.3"/title 'Standardized Canonical Coefficients for Set-1'     
      /space 2/rnames=nx1/cnames=num.

*---------------------------------------------------------------------------.
* Compute and print the unstandardized canonical coefficients for set-1.
*---------------------------------------------------------------------------.

compute a1=d1*a.
print a1/format "f8.3"/title 'Raw Canonical Coefficients for Set-1'             
      /space 2/rnames=nx1/cnames=num.

*---------------------------------------------------------------------------.
* Compute and print the standardized canonical coefficients for set-2.
*---------------------------------------------------------------------------.

do if (p1 le p2).
compute b=r2_inv*v.
else.
compute b=r2_inv*u.
end if.
do if (p1 le p2).
compute b=b(:,1:p1).
end if.
print b /format "f8.3"/title 'Standardized Canonical Coefficients for Set-2'    
      /space 2/rnames=nx2/cnames=num.

*---------------------------------------------------------------------------.
* Compute and print the unstandardized canonical coefficients for set-2.
*---------------------------------------------------------------------------.

compute b1=d2*b.
print b1/format "f8.3"/title 'Raw Canonical Coefficients for Set-2'             
      /space 2/rnames=nx2/cnames=num.

*---------------------------------------------------------------------------.
* Compute and print the canonical loadings for set-1.
*---------------------------------------------------------------------------.

compute tem=d1*s1*a1.
print tem /format "f8.3"/title 'Canonical Loadings for Set-1'
   /space 2/rnames=nx1/cnames=num.

*---------------------------------------------------------------------------.
* Compute the redundancy index as the proportion of variance in set-1
* explained by its own canonical variates.
*---------------------------------------------------------------------------.

compute f1=cssq(tem)/p1.
compute f1=t(f1).

*---------------------------------------------------------------------------.
* Compute and print cross loadings for set-1.
*---------------------------------------------------------------------------.

compute tem=d1*s12*b1.
print tem /format "f8.3"/title 'Cross Loadings for Set-1'                       
   /space 2/rnames=nx1/cnames=num.
                                                                               
*---------------------------------------------------------------------------.
* Compute the redundancy index as the proportion of variance in set-1          
* explained by the set-2 canonical variates.
*---------------------------------------------------------------------------.                                                                               

compute cs3=cssq(tem)/p1.
compute cs3=t(cs3).

*---------------------------------------------------------------------------.
* Compute and print cross loadings for set-2.
*---------------------------------------------------------------------------.

compute tem=d2*s2*b1.
print tem /format "f8.3"/title 'Canonical Loadings for Set-2'                   
   /space 2/rnames=nx2/cnames=num.

*---------------------------------------------------------------------------.
* Compute the redundancy index as the proportion of variance in set-2          
* explained by its own canonical variates.
*---------------------------------------------------------------------------.

compute f2=cssq(tem)/p2.
compute f2=t(f2).

*---------------------------------------------------------------------------.
* Compute and print cross loadings for set-2.
*---------------------------------------------------------------------------.

compute tem=d2*t(s12)*a1.
print tem /format "f8.3"/title 'Cross Loadings for Set-2'                       
   /space 2/rnames=nx2/cnames=num.

*---------------------------------------------------------------------------.
* Compute the redundancy index as the proportion of variance in set-2
* explained by the set-1 canonical variates.
*---------------------------------------------------------------------------.

compute cs4=cssq(tem)/p2.
compute cs4=t(cs4).

*---------------------------------------------------------------------------.
* Print redundancy analysis results.
*---------------------------------------------------------------------------.

compute c1={!@1}.
compute c2={!@2}.

print /title '            Redundancy Analysis:' /space 2.
print f1/format "f15.3"
  /title 'Proportion of Variance of Set-1 Explained by Its Own Can. Var.'
    /space 2/rnames=c1/cnames= {"Prop Var"}.
print cs3/format "f15.3"
  /title 'Proportion of Variance of Set-1 Explained by Opposite Can.Var.'
    /space 2/rnames=c2/cnames= {"Prop Var"}.
print f2/format "f15.3"
  /title 'Proportion of Variance of Set-2 Explained by Its Own Can. Var.'
    /space 2/rnames=c2/cnames= {"Prop Var"}.
print cs4/format "f15.3"
  /title 'Proportion of Variance of Set-2 Explained by Opposite Can. Var.'
    /space 2/rnames=c1/cnames= {"Prop Var"}.

*---------------------------------------------------------------------------.
* Create files for use in calculation of canonical scores.
*---------------------------------------------------------------------------.

SAVE {P1,P2} / OUTFILE 'C:\Temp\CC__SIZE.SAV'.
SAVE {T(A1),T(B1)} / OUTFILE 'C:\Temp\CC__AB.SAV' .
END MATRIX.
*---------------------------------------------------------------------------.
* Create a file with variable names and set number variable
*---------------------------------------------------------------------------.

SET MESSAGES OFF RESULTS OFF.
SELECT IF $CASENUM=1.
DO REPEAT V=!SET1.
COMPUTE V=1.
END REPEAT.
DO REPEAT V=!SET2.
COMPUTE V=2.
END REPEAT.
STRING VARNAME (A8).
COMPUTE VARNAME='SET_NUM'.
FLIP VARIABLES !SET1 !SET2 / NEWNAMES=VARNAME .
COMPUTE VARSEQ=1.
SPLIT FILE BY SET_NUM.
CREATE VARSEQ=CSUM(VARSEQ).
SAVE OUTFILE 'C:\Temp\CC_NAMES.SAV'.
GET FILE 'C:\Temp\CC__SIZE.SAV' .

*---------------------------------------------------------------------------.
* Set up required information to create compute statements for scoring.
*---------------------------------------------------------------------------.

WRITE OUTFILE 'C:\Temp\CC__AB.INC'
   /'STRING @NMA001 TO @NMA',COL1 (N3), ' (A8)'
   /'       @NMB001 TO @NMB',COL2 (N3), ' (A8)'
   /'VECTOR @NMA= @NMA001 TO @NMA',COL1 (N3)
   /'       @NMB= @NMB001 TO @NMB',COL2 (N3)
   /'COMPUTE N_A='COL1 (N3)
   /'COMPUTE N_B='COL2 (N3)
   /'IF (SET_NUM=1) @NMA(VARSEQ)=CASE_LBL'
   /'IF (SET_NUM=2) @NMB(VARSEQ)=CASE_LBL'
   /'COMPUTE @=1'
   /'AGGREGATE OUTFILE "C:\Temp\CC__SPRD.SAV" / BREAK @'
   / ' / N_A=MAX(N_A) / N_B=MAX(N_B)'
   / ' / @NMA001 TO @NMA',COL1 (N3) '=MAX (@NMA001 TO @NMA',COL1 (N3),')'
   / ' / @NMB001 TO @NMB',COL2 (N3) '=MAX (@NMB001 TO @NMB',COL2 (N3),')'
   / 'GET FILE "C:\Temp\CC__AB.SAV"'
   / 'COMPUTE @=1'
   / 'MATCH FILES FILE * / TABLE "C:\Temp\CC__SPRD.SAV"/BY @'
   / 'VECTOR @NMA= @NMA001 TO @NMA',COL1 (N3)
   / '       @NMB= @NMB001 TO @NMB',COL2 (N3)
   / '       COEF= COL1 TO @'.
EXECUTE.

GET FILE 'C:\Temp\CC_NAMES.SAV'.
INCLUDE FILE 'C:\Temp\CC__AB.INC'.
SET PRINTBACK OFF.

*---------------------------------------------------------------------------.
* Write out the compute statements for scoring.
*---------------------------------------------------------------------------.

STRING @SCNM@ (A8).
STRING @OLDNM@ (A8).
COMPUTE @SCNM@=CONCAT('S1_CV',STRING($CASENUM,N3)).
WRITE OUTFILE 'C:\Temp\CC__.INC' /'COMPUTE ',@SCNM@ ,'= 0'.
LOOP CC@@@ = 1 TO N_A.
COMPUTE @OLDNM@=@NMA(CC@@@).
COMPUTE @COEF@ =COEF(CC@@@).
WRITE OUTFILE 'C:\Temp\CC__.INC' / ' +',@COEF@ (comma18.16),' * ',@OLDNM@ .
END LOOP.
COMPUTE @SCNM@=CONCAT('S2_CV',STRING($CASENUM,N3)).
WRITE OUTFILE 'C:\Temp\CC__.INC' /'COMPUTE ',@SCNM@ ,'= 0'.
LOOP CC@@@=1 TO N_B.
COMPUTE @OLDNM@=@NMB(CC@@@).
COMPUTE @COEF@ =COEF(CC@@@+N_A).
WRITE OUTFILE 'C:\Temp\CC__.INC' / ' +',@COEF@ (comma18.16),' * ',@OLDNM@ .
END LOOP.
EXECUTE.

*---------------------------------------------------------------------------.
* Get the original data and run the scoring program.
*---------------------------------------------------------------------------.

DATASET ACTIVATE temp_copy.
INCLUDE FILE 'C:\Temp\CC__.INC' .
ERASE FILE 'C:\Temp\CC__SIZE.SAV' .
ERASE FILE 'C:\Temp\CC__AB.INC'.
ERASE FILE 'C:\Temp\CC_NAMES.SAV'.
ERASE FILE 'C:\Temp\CC__AB.SAV'.
ERASE FILE "C:\Temp\CC__SPRD.SAV" 
!IF (!KEEPSC ='N') !THEN
ERASE FILE 'C:\Temp\CC__.INC'.
!ELSE

DO IF ($CASENUM=1).
SET RESULTS ON.
PRINT /'The canonical scores have been written to the active file.'
      /'Also, a file containing an SPSS Scoring program has been written'
      /'To use this file GET a system file with the SAME variables'
      /'Which were used in the present analysis.  Then use an INCLUDE command'
      /'to run the scoring program.'
      /'For example :' /
      /'GET FILE anotherfilename'
      /'INCLUDE FILE "C:\Temp\CC__.INC".'
      /'EXECUTE.'.
END IF.
EXECUTE.

!IFEND.
RESTORE.
!enddefine.
RESTORE.
