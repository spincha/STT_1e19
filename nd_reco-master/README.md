# nd_reco

* Input : Edep-Sim root file

* Output: ROOT tree file

# how to run

1. cmake .
2. make
3.
```
./smearBMgo -g 870 -v both -i inputEdepSimFile.root -o out.root
```
or check the detailed usage by
```
./smearBMgo -h
```
# how to read output

> root readEdep_smeared.C
