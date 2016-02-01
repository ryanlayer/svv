# SVView

A simple tool that can help visualize the coverage near a variant across multiple BAMs.


    $ gqt query -v -i target.vcf.gz -d target.ped.db \
      -p "phenotype == 2" -g "HET HOM_ALT" \
      -p "phenotype == 1" -g "HOM_REF" \
    | bedtools intersect -u -header -wa -a stdin -b data/h37_ensemble_exons.bed.gz \
    | bcftools view -i "SVTYPE=='DEL' && SVLEN>-500000" \
    | ~/scratch/svv/sv_depth.py target.ped 10000
    
  This will create a covrage file named "var_ID.txt".
  
    $ cat var_1.txt | spark.py \
      -e data/h37_ensemble_exons.bed.gz \
      -n data/h37_ensemble_exons.togenename.txt \
      -o var_1.png
      
