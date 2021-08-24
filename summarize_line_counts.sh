dir=methyl_calls
echo lines,analysis,reference,aligner,method,filtering,mass
for file in `find "$dir" -name '*.tsv' ! -name '*called_genes*' ! -name '*clinvar_genes*'`; do
  # for quick testing:
  # cnt=`head "$file" | wc -l | cut -f 1 -d ' '`
  cnt=`wc -l "$file"| cut -f 1 -d ' '`
  basename=$(basename $file)
  if [ $(dirname $file) == 'methyl_calls/clinvar_calls' ]; then # to prevent name conflict
    basename=$(echo $basename | sed -r 's/calls/clinvar_calls/')
  fi
  #echo $file
  printf "$cnt,$basename\n" | \
  sed -r 's/(\.*)_(t2t|grch38)/\1,\2/' |\
  sed -r 's/(t2t|grch38)_(bismark|bwameth)_(emseq|wgbs)_((filt|nofilt)_)?([0-9]*ng).*/\1,\2,\3,\5,\6/' |\
  sed 's/,,/,combined,/'
done

for file in `find "$dir" -name '*called_genes*.tsv' -o -name '*clinvar_genes*.tsv'`; do
  cnt=`cut -f 10 $file | sed "s/.*gene_id=//" | sed "s/;.*//" | sort | uniq | wc -l | cut -f 1 -d ' '`
  basename=$(basename $file)
  printf "$cnt,$basename\n" | \
  sed -r 's/(\.*)_(t2t|grch38)/\1,\2/' |\
  sed -r 's/(t2t|grch38)_(bismark|bwameth)_(emseq|wgbs)_((filt|nofilt)_)?([0-9]*ng).*/\1,\2,\3,\5,\6/' |\
  sed 's/,,/,combined,/'
done