The parameter ```-m``` is optional. Adds a marge of 100 for the variant's breakpoints (when you compare variants with the Truth file) 


work1: number of barcodes in each variant.
```python work1.py -vcf donnees1/SVs/HSapiensChr1_Simulated/candidateSV_inversion.vcf -bam donnees1/BAM/HSapiensChr1_Simulated/possorted_bam.bam -t donnees1/SVs/HSapiensChr1_Simulated/Truth -m```

work2: number of isolated barcodes in each variant.
```python work2.py -vcf donnees1/SVs/HSapiensChr1_Simulated/candidateSV_inversion.vcf -bam donnees1/BAM/HSapiensChr1_Simulated/possorted_bam.bam -bci NewBCIs/HSapiensChr1_Simulated/possorted_bam.bci -t donnees1/SVs/HSapiensChr1_Simulated/Truth -m```


work3: number of common barcodes between the left-region of breakpoint and the right one.
```python work3.py -vcf donnees1/SVs/Ecoli_Simulated/candidateSV_inversion.vcf -bam donnees1/BAM/Ecoli_Simulated/possorted_bam.bam -t donnees1/SVs/Ecoli_Simulated/Truth -m```
