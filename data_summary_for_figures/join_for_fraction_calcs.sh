
file="${1:-all_methyl_call_counts.csv}"
header='calls,analysis,reference,aligner,method,filtering,mass,miscalls,joinkey,clinvar_lines'

echo $header
csvjoin -H -c 8 \
  <(grep ',calls,' $file | grep ',nofilt,' | tr -d '\r' | awk -F, '{print $0","$3"_"$4"_"$5"_"$7}') \
  <(grep ',miscalls,' $file | tr -d '\r' | awk -F, '{print $0","$3"_"$4"_"$5"_"$7}') \
  <(grep ',clinvar_miscalls,' $file | tr -d '\r' | awk -F, '{print $0","$3"_"$4"_"$5"_"$7}') 2>/dev/null | #hides warnings about dup headers
  tail -n +2 | #to remove generic header 
  tr -d '\r' | awk -F, -v OFS=',' '{print $1,$2,$3,$4,$5,$6,$7,$9,$8,$16}'
