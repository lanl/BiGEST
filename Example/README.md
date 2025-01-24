# BiGEST Example


For this example, we will use the fasta file included, [Chlorella sorokiniana UTEX 1230](https://phycocosm.jgi.doe.gov/Chloso1230_1_1/Chloso1230_1_1.home.html). 



The command to run BiGEST for this example file is:

```
cd BiGEST/
python3 BiGEST.py --db db/BiGEST_DB -f Example/Chloso1230_1_1_AssemblyScaffolds_2023-12-07.fasta.gz -o Example/output -a True
```

This will run rpstblastn and antiSMASH on our fasta file and create the BiGEST output files, with a final output of a GenBank file with the sections of this genome: Chlorella sorokiniana UTEX 1230. 

If you have already run blast and antismash, then you can run:

```
python3 BiGEST.py -i Example/blast_results/Chloso1230_BiGEST_CDD_rpstblastn.txt -f Example/Chloso1230_1_1_AssemblyScaffolds_2023-12-07.fasta.gz -o Example/output -a Example/antismash_results/Chloso1230_1_1_AssemblyScaffolds_2023-12-07.fasta.gz/1230_1_1_AssemblyScaffolds_2023-12-07.gbk
```


As the BiGEST script runs, it will give status updates on which contig it is on. If an antiSMASH flag is given, it will compile those results after processing each contig, so a delay at the end of the script is normal. 
For Chloso1230_1_1_AssemblyScaffolds_2023-12-07.fasta, the output should look like this: 


![Example output for BiGESt](https://gitlab.lanl.gov/ladriani/bigest/-/blob/main/images/example_output.png)

