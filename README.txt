READ_ME - snp_annotation.py

author: Jurrien de Jong
date: 16-10-2021
files: snp_annotation.py, data.msf, seq.txt

Quick Note:
	If you are looking for some test data, you can find it at the usage section!

Description of the code:
	The program snp_annotation.py will handle a few input arguments:
    -i : Input file in .msf ( so the MSA alignment )
    -seq : The sequence ( DNA or Protein ) to be aligned with the MSA ( for this instance a text file )
    -pos : The position in the above sequence where a SNP needs to be implemented
    -SNP : The SNP as a string; like "A"
  When provided with arguments, the program will read the .msf file and safe the alignment in a variable. 
  After obtaining a proper protein sequence, the SNP will be implemented at the given position ( if valid ) by another function.
  The alignment itself will be done by the function 'get_align_score'. It will compare each column of the MSA with the given sequence,
  and give a score between 0 and the length of the sequence. So 10 matches out of 50 will give a score of 10 / 50. Finally,
  the scores will be processed by the function 'write_results'. It will calculate a percentage score between 0 and 100,
  and print a corresponding message. Afterward it will end the program.

  How to interpret the score:
    * 0 - 10 % similarity:
        This is a pretty bad SNP, it has a lot of consequences! The AA might be pretty preserved!

    * 10 - 40 % similarity:
        This SNP has some bad affects, but is not terrible.

    * 40 - 80 % similarity:
        The SNP has very few bad effect.

    * 80 - 100 % similarity:
        This SNP has a neutral effect!

Usage:
	Please open your terminal and head to the folder containing the files, if you cannot get to it
	please visit the links in the support section.

	snp_annotation.py can be run with the following syntax:
	
	python3 SNP_Annotation.py -i [input file in .msf format ] -seq [.txt file with a sequence] -pos [position of SNP] -snp [SNP as a string]
	
  Here are a few examples:

python3 snp_annotation.py -i data.msf -seq seq.txt -pos 0 -snp "A"
  
python3 snp_annotation.py -i data.msf -seq seq.txt -pos 0 -snp "T"

Support:
	Please view these pages if you are having trouble with cd, running the command or classes
	
	
	https://docs.python.org/3/tutorial/classes.html
	
	https://linuxize.com/post/linux-cd-command/
	
	https://realpython.com/run-python-scripts/
  
  For further questions, please send me an email:
  ju.de.jong@st.hanze.nl
