COMMENT in the INCLUDE line you need to put the specific path to the the Canonical Correlation macro.
COMMENT if you save it to your temp file on the C drive it will work as written.  
COMMENT you must have the ats_data set open in order for this to work.  
COMMENT set1 includes the Xs and set2 includes the Ys.

INCLUDE 'C:\Temp\ccmac_atafixed.sps'.
CANCORR set1 =ts, tc/
	     set2 =bs, bc.