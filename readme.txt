-----------------------------------------------------------
Dependencies
-----------------------------------------------------------
1. Chromosome3D is implemented to run in Linux environment.
   It was tested in "x86_64 GNU/Linux" OS.
2. Perl v5.10.1 was used during development and testing;
   But it should run with other versions of perl as well.
3. CNS suite

-----------------------------------------------------------
Installation
-----------------------------------------------------------
1. Install CNS suite
   1.1. To download CNS suite, provide your academic profile related 
        information at http://cns-online.org/cns_request/. An email
        with (a) link to download, (b) login, and (c) password
        will be sent to you. Follow the link, possibly
        http://cns-online.org/download/, and download 
        CNS suite "cns_solve_1.3_all_intel-mac_linux.tar.gz".
   1.2. Unzip
        $ tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz
   1.3. Change directory to cns_solve
        $ cd cns_solve_1.3
   1.4. Unhide the file '.cns_solve_env_sh'
        $ mv .cns_solve_env_sh cns_solve_env.sh
   1.5. Edit 'cns_solve_env.sh' and 'cns_solve_env' to replace
        '_CNSsolve_location_' with CNS installation directory.
        For instance, if your CNS installation path is
        '/home/user/programs/cns_solve_1.3' replace
        '_CNSsolve_location_' with this path
   1.6. Install
        $ make install
   1.7. Increase the value for ‘nrestraints’ (maximum number of restraints it can take)
        Change the code at line 60 of the module ‘cns_solve_1.3/modules/nmr/readdata’
		$ vim cns_solve_1.3/modules/nmr/readdata
		- change 20000 to 200000 (by adding a zero)
   1.8. Test CNS installation
        $ source cns_solve_env.sh
        $ cd test 
        $ ../bin/run_tests -tidy *.inp
 
2. Download chromosome3D_v1.0.tar.gz
   3.1 Download chromosome3D_v1.0.tar.gz if you don't have it.
   3.2 Untar
       $ tar zxvf chromosome3D_v1.0.tar.gz
   3.3 Change directory to chromosome3D_v1.0
       $ cd chromosome3D_v1.0

3. Change variable values in the chromosome3D.pl file
   3.1 Change the path of the variable $cns_suite 
       to CNS installation directory
   3.2 Make it executable
       $chmod +x chromosome3D.pl
  
4. Test Chromosome3D
   4.1 Execute "perl ./chromosome3D.pl" or "./chromosome3D.pl"
       It should print the usage information.
   4.2 Test using an example
       $ ./chromosome3D.pl -if "./input/chr22_1mb_matrix.txt" -o "./output/chr22_1mb"
   4.3 Test by building structures for all the chromosomes at 1MB and 500KB
       $ ./test.sh
 
-----------------------------------------------------------
Please cite:
"Chromosome3D: a Distance Geometry Method for Reconstructing Three-Dimensional 
Chromosomal Structures from Hi-C Interaction Frequency Data" (submitted),
Badri Adhikari^, Tuan Trieu^, Jianlin Cheng
^These authors contributed equally to this work 

-----------------------------------------------------------
chengji@missouri.edu (PI)
