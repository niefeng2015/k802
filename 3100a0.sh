rm -rf hashcat.potfile
rm -rf hashcat.potfile00000

./hashcat -m 3100 -a 0 3100.hash chenxi.dict
