wget -O compound.xml "http://biocyc.org/xmlquery?[x:x%3C-ecoli^^compounds]&detail=full"
printf "compounds table done\n"
wget -O reaction.xml "http://biocyc.org/xmlquery?[x:x%3C-ecoli^^reactions]&detail=full"
printf "reactions table done\n"
wget -O pathway.xml "http://biocyc.org/xmlquery?[x:x%3C-ecoli^^pathways]&detail=full"
printf "pathways table done\n"
wget -O protein.xml "http://biocyc.org/xmlquery?[x:x%3C-ecoli^^proteins]&detail=full"
printf "proteins table done\n"