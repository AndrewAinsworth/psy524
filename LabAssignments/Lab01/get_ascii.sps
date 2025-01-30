* Encoding: UTF-8.
COMMENT for the CD command below, replace with the file directory for your anova4.dat data file. 
cd 'C:\Users\ata20315\Box\Classes\Psy524\Spring24\psy524-master\LabAssignments\Lab01'. 
data list fixed file="anova4.dat"
 / gender 1 beforet1 2-4 beforet2 5-6 beforet3 7-9 beforet4 10-11 aftert1 12-14 aftert2 15-17 aftert3 18-20 after4 21-23.
LIST.

SAVE OUTFILE='anova4.sav'.


