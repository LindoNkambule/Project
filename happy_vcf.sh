export HGREF=/path/to/reference.fa
happy="/path/to/hap.py-build/bin/hap.py"
golden_vcf="/path/to/golden.vcf"
for filename in *.vcf; do
  vcf=${filename%.*}
  mkdir ${filename%.*}_HAPPY
  $happy $golden_vcf $filename -o $vcf
  rm *.gz *.tbi *.json
  mv $vcf* ${filename%.*}_HAPPY/
done
