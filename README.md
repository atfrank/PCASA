# PCASA: Protein Cα-based Solvent Accessibilities
  
- Predicts residue solvent accessibilities based on Cα coordinates.
- <img src="http://www.sciweavers.org/tex2img.php?eq=SASA_%7Bi%7D%20%3D%20%5Calpha_%7Bi%7D%20-%20%5Csum_%7B%28j%20%20%5Cin%20%20r_%7Bij%7D%20%20%5Cleq%20%20r_%7Bcut%7D%29%7D%20%5Cbeta_%7Bij%7D%20r_%7Bij%7D%5E%7B%5Cgamma%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="SASA_{i} = \alpha_{i} - \sum_{(j  \in  r_{ij}  \leq  r_{cut})} \beta_{ij} r_{ij}^{\gamma}" width="218" height="43" />

## Install
```shell
$ cd /path/to/CASA/
$ make clean
$ make 
```

## Usage manual
```shell
$ bin/casa -h
Usage: casa [-options] <PDBfile>
Options:
         -sasafile : input file with comparison SASA (string)
              format >> white space-separated >> 3-letter AA code, residue number, SASA 
         -header : print header information (flag)
         -predict : run in prediction mode (flag)
         -cut : cutoff distances; only contacts with the cutoff are considered (int)
         -trj : trajectory file (string)
         -skip : skip rate for reading trajectory (integer)
         -start : frame at which to start reading trajectory (integer)
         -stop : frame at which to stop reading trajectoryframe (integer)
         -identification : ID tag used in output (string)
```

## Examples
```shell
$ # predict solvent accessibilities from a coordinate file 
$ bin/casa -predict tests/file.pdb (PDB format)
$
$ # predict chemical shifts from a trajectory file (DCD format) 
$ bin/casa -predict -trj tests/file.dcd tests/file.pdb
```
bin/casa -predict tests/file.pdb

## Output
### format
_trajectory-frame, id-tag, atomname, residue-number, residue-name, reference-sasa, target-sasa, predicted-sasa_
### example
```shell
$ bin/casa -predict -sasafile tests/file.dat tests/file.pdb
frame ID atomname resid resname refSASA targetSASA predSASA
1 test CA 1 LYS 206.142 198.653 146.896
1 test CA 2 SER 150.044 90.561 90.1534
1 test CA 3 ALA 137.505 86.5428 73.1248
1 test CA 4 LYS 206.142 77.8203 72.1502
1 test CA 5 ASP 186.735 65.3621 76.9402
1 test CA 6 ALA 137.505 57.1637 50.7573
1 test CA 7 LEU 189.262 8.87559 30.5331
1 test CA 8 LEU 189.262 33.5205 27.855
1 test CA 9 LEU 189.262 96.2658 45.6996
1 test CA 10 TRP 270.484 23.5421 107.336
1 testing CA 11 CYS 153.447 0.000153999 2.1017
1 testing CA 12 GLN 220.199 55.2656 38.3275
1 testing CA 13 MET 208.696 113.255 85.332
1 testing CA 14 LYS 206.142 65.1948 57.4488
1 testing CA 15 THR 173.03 4.51759 7.755
1 testing CA 16 ALA 137.505 73.8344 50.0637
1 testing CA 17 GLY 101.908 84.0864 51.4855
1 testing CA 18 TYR 265.457 28.7181 84.9654
1 testing CA 19 PRO 178.326 95.2985 79.2481
1 testing CA 20 ASN 184.135 76.9341 78.2766
1 testing CA 21 VAL 185.402 13.5256 25.9715
1 testing CA 22 ASN 184.135 110.19 55.9439
1 testing CA 23 ILE 211.18 4.55099 33.9554
1 testing CA 24 HIS 183.775 125.466 77.3049
1 testing CA 25 ASN 184.135 48.5255 60.2922
1 testing CA 26 PHE 261.103 1.22718 22.3274
1 testing CA 27 THR 173.03 29.5924 5.7757
1 testing CA 28 THR 173.03 58.0605 27.1849
1 testing CA 29 SER 150.044 13.5907 24.9031
1 testing CA 30 TRP 270.484 0.000361735 19.7456
1 testing CA 31 ARG 236.902 66.9004 0
1 testing CA 32 ASP 186.735 13.1151 0
1 testing CA 33 GLY 101.908 0.000104317 0
1 testing CA 34 MET 208.696 1.27006 0
1 testing CA 35 ALA 137.505 0.000155516 0
1 testing CA 36 PHE 261.103 0.000267815 0
1 testing CA 37 ASN 184.135 0.00051426 0
1 testing CA 38 ALA 137.505 0.000155628 0
1 testing CA 39 LEU 189.262 0.000490252 0
1 testing CA 40 ILE 211.18 6.09256e-05 0
1 testing CA 41 HIS 183.775 36.2151 5.70499
1 testing CA 42 LYS 206.142 103.993 39.2768
1 testing CA 43 HIS 183.775 57.291 70.5003
1 testing CA 44 ARG 236.902 74.5092 71.1028
1 testing CA 45 PRO 178.326 84.8675 50.1049
1 testing CA 46 ASP 186.735 110.4 95.4533
1 testing CA 47 LEU 189.262 25.5779 49.4521
1 testing CA 48 ILE 211.18 9.93742 38.7266
1 testing CA 49 ASP 186.735 75.6502 85.4309
1 testing CA 50 PHE 261.103 8.50226 62.2941
1 testing CA 51 ASP 186.735 121.303 102.209
1 testing CA 52 LYS 206.142 155.843 127.431
1 testing CA 53 LEU 189.262 11.7579 35.4757
1 testing CA 54 LYS 206.142 134.867 67.5259
1 testing CA 55 LYS 206.142 102.651 66.8599
1 testing CA 56 SER 150.044 104.178 72.7155
1 testing CA 57 ASN 184.135 71.5457 66.3684
1 testing CA 58 ALA 137.505 10.0498 0
1 testing CA 59 HIS 183.775 87.4487 33.9589
1 testing CA 60 TYR 265.457 139.655 72.5834
1 testing CA 61 ASN 184.135 0.000173838 0
1 testing CA 62 LEU 189.262 0 0
1 testing CA 63 GLN 220.199 79.562 35.755
1 testing CA 64 ASN 184.135 28.5225 22.8425
1 testing CA 65 ALA 137.505 0 0
1 testing CA 66 PHE 261.103 0.0266165 0
1 testing CA 67 ASN 184.135 41.9797 37.4739
1 testing CA 68 LEU 189.262 33.6049 29.7702
1 testing CA 69 ALA 137.505 0.00024102 4.40694
1 testing CA 70 GLU 216.958 56.6934 59.9704
1 testing CA 71 GLN 220.199 149.57 97.7012
1 testing CA 72 HIS 183.775 101.584 87.3304
1 testing CA 73 LEU 189.262 9.37555 61.6599
1 testing CA 74 GLY 101.908 62.3747 41.6254
1 testing CA 75 LEU 189.262 18.3699 41.7409
1 testing CA 76 THR 173.03 88.9944 63.2823
1 testing CA 77 LYS 206.142 89.6859 53.1617
1 testing CA 78 LEU 189.262 95.4654 65.469
1 testing CA 79 LEU 189.262 6.64622 50.8683
1 testing CA 80 ASP 186.735 68.1461 70.6404
1 testing CA 81 PRO 178.326 14.9921 8.49702
1 testing CA 82 GLU 216.958 129.157 82.0767
1 testing CA 83 ASP 186.735 94.3001 58.7114
1 testing CA 84 ILE 211.18 0.832089 0
1 testing CA 85 SER 150.044 32.1735 6.13299
1 testing CA 86 VAL 185.402 76.7603 59.4181
1 testing CA 87 ASP 186.735 113.969 90.1513
1 testing CA 88 HIS 183.775 154.211 74.6631
1 testing CA 89 PRO 178.326 5.75296 22.3515
1 testing CA 90 ASP 186.735 84.9023 52.1474
1 testing CA 91 GLU 216.958 64.9065 60.6884
1 testing CA 92 LYS 206.142 169.512 89.418
1 testing CA 93 SER 150.044 47.4317 32.5154
1 testing CA 94 ILE 211.18 0.000266664 20.4314
1 testing CA 95 ILE 211.18 21.9915 26.6108
1 testing CA 96 THR 173.03 64.0118 62.0393
1 testing CA 97 TYR 265.457 0.000379593 40.3207
1 testing CA 98 VAL 185.402 0.000770054 12.6376
1 testing CA 99 VAL 185.402 45.6866 50.1855
1 testing CA 100 THR 173.03 30.2165 45.7425
1 testing CA 101 TYR 265.457 0.000366015 27.2282
1 testing CA 102 TYR 265.457 84.1802 73.7535
1 testing CA 103 HIS 183.775 79.53 77.6874
1 testing CA 104 TYR 265.457 54.4673 91.1056
1 testing CA 105 PHE 261.103 19.1687 68.0807
1 testing CA 106 SER 150.044 54.3165 69.3892
1 testing CA 107 LYS 206.142 142.239 124.314
1 testing CA 108 MET 208.696 197.075 155.057


