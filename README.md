# Stupid Simple Structural Variant View

A two-step process that can help visualize the coverage near a variant across multiple BAMs.

Step 1: Create a coverage file for each variant in a VCF.  (NOTE: use [GQT](https://github.com/ryanlayer/gqt), [bcftools](https://github.com/samtools/bcftools), [bedtools](https://github.com/arq5x/bedtools2) etc to filter your full VCF down to a more interesting/managable subset.)

    $ bcftools view test/test.vcf.gz | ./sv_depth.py test/test.ped 1000
    
`sv_depth.py` creates a covrage file  "var_ID.txt" for each SV in the VCF, where ID is the variant ID (3rd column of the variant record).  In this case, `test/test.vcf.gz` has two variants, `DEL_pindel_18715` and  `UW_VH_5456`, and so it creates `var_DEL_pindel_18715.txt` and `var_UW_VH_5456.txt`.  The two parameters to `sv_depth.py` are the ped file (where the 8th column is a path to the BAM files) and a the distance up and down stream of the the SV that is inclued in the coverage file.

Step 2: Visualize these coverages with `spark.py`, which takes an exome file, a map from transcrip to gene name, and an ouputfile.
  
    $ cat var_DEL_pindel_18715.txt \
        | ./spark.py \
            -e data/h37_ensemble_exons.bed.gz \
            -n data/h37_ensemble_exons.togenename.txt \
            -o var_DEL_pindel_18715.png

![var_DEL_pindel_18715.png](test/var_DEL_pindel_18715.png)

    $ cat var_UW_VH_5456.txt \
        | ./spark.py \
            -e data/h37_ensemble_exons.bed.gz \
            -n data/h37_ensemble_exons.togenename.txt \
            -o var_UW_VH_5456.png

![var_UW_VH_5456.png](test/var_UW_VH_5456.png)

```
gqt query -i 15-0017322.15-0017321.15-0017320.15-0017319.15-0017318.15-0017317.15-0017316.15-0017315.x_btu356_LCR-hs37d5.MT.hs37d5.phix.bed.non_ref.IHH.vcf.gz -d ihh.ped.db -p "phenotype==2" -g "HET HOM_ALT" -p "phenotype==1" -g "HOM_REF" -v | bedtools intersect -header -a genes.bed -b stdin | sort -u
1      	109296900      	109296900      	STXBP3
1      	175091590      	175091645      	TNN
13     	113499583      	113500053      	ATP11A
15     	58850831       	58850831       	LIPC
15     	78418614       	78423877       	CIB2
2      	182443195      	182445136      	CERKL
2      	189119618      	189120981      	LINC01090
2      	196999182      	196999182      	STK17B
3      	183876861      	183877361      	DVL3
4      	124630701      	124630701      	LINC01091
4      	187614520      	187614843      	FAT1
7      	105160926      	105160926      	PUS7
7      	105161909      	105161909      	PUS7
7      	136834298      	136834750      	LOC349160
8      	645443 	645708 	ERICH1
```
