for tt in `sed 's/\t/:/' scratch/test_data/mOTUs.txt`;
do
echo $tt | cut -f1 -d":";
fs=`echo $tt | cut -f2 -d":"| sed 's#;#.faa scratch/test_data/proteoms/#g'`;
name=`echo $tt | cut -f1 -d":"`
mOTUlizer/bin/__main__.py --genome2cog_only -n $name --faas scratch/test_data/proteoms/${fs}.faa >> test_data/cog_sets/${name}.json;
done

for tt in `sed 's/\t/:/' scratch/test_data/mOTUs.txt`;
do
  name=`echo $tt | cut -f1 -d":"`
  mOTUlizer/bin/__main__.py --ls -n ${name}.ls --cogs test_data/cog_sets/${name}.json > /dev/null
  for s in `seq 30 3 95`
  do
#    mOTUlizer/bin/__main__.py -s $s -n ${name}.$s --cogs test_data/cog_sets/${name}.json > /dev/null
  mOTUlizer/bin/__main__.py --rs -n ${name}.${s}.rd --cogs test_data/cog_sets/${name}.json > /dev/null
  done
done  2>> optimising.txt
