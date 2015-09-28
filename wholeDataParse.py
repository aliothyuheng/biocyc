import csv
from bs4 import BeautifulSoup
import time

start_time = time.time()
#compound table
soup = BeautifulSoup(open("fly_compound_full.xml"), "xml")
compound = soup.find_all("Compound")
temp = [None]*10
compound_result = []
temp_left = []
temp_right = []
reaction_left = []
reaction_right = []
temp_syn = []
syn = []

for row in compound:
    if row.has_attr('ID'):
        temp[0] = row['ID'].strip()
        temp[0] = temp[0].decode("utf-8")
        temp[0] = temp[0][4:]
        temp_left.append(temp[0])
        temp_right.append(temp[0])
        temp_syn.append(temp[0])
        if row.find("molecule"):
            molecule = row.find("molecule")
            if molecule.find("formula"):
                temp[2] = molecule.find("formula")['concise'].strip()
                temp[2] = temp[2].encode("utf-8")
            if molecule.find("string"):
                temp[4] = molecule.find("string").string.strip()
                temp[4] = temp[4].encode("utf-8")
        if row.find("appears-in-left-side-of"):
            left = row.find("appears-in-left-side-of")
            for item in left.find_all("Reaction"):
                temp_left.append(item['frameid'].strip().encode("utf-8"))
        if row.find("appears-in-right-side-of"):
            right = row.find("appears-in-right-side-of")
            for item in right.find_all("Reaction"):
                temp_right.append(item['frameid'].strip().encode("utf-8"))
        if row.find("dblink"):
            db = row.find_all("dblink")
            for item in db:
                db_type = item.find("dblink-db")
                db_id = item.find("dblink-oid")
                if db_type.string.strip() == "CAS":
                    temp[6] = db_id.string.strip()
                    temp[6] = temp[6].encode("utf-8")
                elif db_type.string.strip() == "LIGAND-CPD":
                    temp[3] = db_id.string.strip()
                    temp[3] = temp[3].encode("utf-8")
        if row.find("inchi"):
            temp[9] = row.find("inchi").string.strip()
            temp[9] = temp[9].encode("utf-8")
        if row.find("synonym"):
            synonym = row.find_all("synonym")
            for item in synonym:
                temp_syn.append(item.string.strip().encode("utf-8"))
        if row.find("molecular-weight"):
            mw = float(row.find("molecular-weight").string.strip())
            temp[8] = mw
            temp[1] = int(mw)
        if row.find("common-name"):
            temp[5] = row.find("common-name").string.strip()
            temp[5] = temp[5].encode("utf-8")
        temp[7] = 1
        syn.append(temp_syn)
        reaction_left.append(temp_left)
        reaction_right.append(temp_right)
        compound_result.append(temp)
    temp = [None]*10
    temp_syn = []
    temp_left = []
    temp_right = []

outFile = open("compound.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(compound_result)
outFile.close()

outFile = open("compound_left.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_left)
outFile.close()

outFile = open("compound_right.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_right)
outFile.close()

outFile = open("synonym.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(syn)
outFile.close()

print("compound table done")
print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()

#parse reaction
soup = BeautifulSoup(open("fly_reaction_full.xml"), "xml")
reaction = soup.find_all("Reaction")
temp = []
reaction_result = []
temp_left = []
temp_right = []
reaction_left = []
reaction_right = []
temp_syn = []
reaction_synonym = []
temp_path = []
reaction_pathway = []
temp_enzyme = []
reaction_enzyme = []

for row in reaction:
    if row.has_attr('ID'):
        temp.append(row['ID'].strip())
        temp[0] = temp[0].encode("utf-8")
        temp[0] = temp[0][4:]
        if row.find("ec-number"):
            temp_enzyme.append(temp[0])
            ec_number = row.find("ec-number").text.strip()
            temp_enzyme.append(ec_number.encode("utf-8"))
        if row.find("enzymatic-reaction"):
            enzymatic_reaction = row.find_all("Enzymatic-Reaction")
            if len(temp_enzyme) == 0:
                temp_enzyme.append(temp[0])
                temp_enzyme.append("")
            for item in enzymatic_reaction:
                if item.find("Protein"):
                    protein = item.find("Protein")
                    temp_enzyme.append(protein['frameid'].strip().encode("utf-8"))
        if row.find("reaction-direction"):
            reaction_direction = row.find("reaction-direction").string.strip()
            temp.append(reaction_direction.encode("utf-8"))
        if row.find("left"):
            left = row.find_all("left")
            for item in left:
                if item.find("Compound"):
                    temp_left.append(temp[0])
                    cpd = item.find("Compound")['frameid']
                    cpd = cpd.strip()
                    temp_left.append(cpd.encode("utf-8"))
                    if item.find("coefficient"):
                        temp_left.append(item.find("coefficient").string.strip())
                    else:
                        temp_left.append(1)
                if temp_left:
                    reaction_left.append(temp_left)
                    temp_left = []
        if row.find("right"):
            right = row.find_all("right")
            for item in right:
                if item.find("Compound"):
                    temp_right.append(temp[0])
                    cpd = item.find("Compound")['frameid']
                    cpd = cpd.strip()
                    temp_right.append(cpd.encode("utf-8"))
                    if item.find("coefficient"):
                        temp_right.append(item.find("coefficient").string.strip())
                    else:
                        temp_right.append(1)
                if temp_right:
                    reaction_right.append(temp_right)
                    temp_right = []
        if row.find("common-name"):
            temp_syn.append(temp[0])
            name = row.find("common-name").string.strip()
            temp_syn.append(name.encode("utf-8"))
        if row.find("synonym"):
            synonym = row.find_all("synonym")
            if len(temp_syn) == 0:
                temp_syn.append(temp[0])
                temp_syn.append("")
            for item in synonym:
                syn_name = item.string.strip()
                temp_syn.append(syn_name.encode("utf-8"))
        if row.find("in-pathway"):
            temp_path.append(temp[0])
            pathway = row.find_all("in-pathway")
            for item in pathway:
                if item.find("Pathway"):
                    pathway_id = item.find("Pathway")['frameid'].strip()
                    temp_path.append(pathway_id.encode("utf-8"))
        if len(temp) == 1:
            temp.append("")
            temp.append(1)
        else:
            temp.append(1)
        reaction_result.append(temp)
        temp = []
        if temp_enzyme:
            reaction_enzyme.append(temp_enzyme)
            temp_enzyme = []
        if temp_syn:
            reaction_synonym.append(temp_syn)
            temp_syn = []
        if temp_path:
            reaction_pathway.append(temp_path)
            temp_path = []


outFile = open("reaction.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_result)
outFile.close()

outFile = open("reaction_left.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_left)
outFile.close()

outFile = open("reaction_right.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_right)
outFile.close()

outFile = open("reaction_synonym.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_synonym)
outFile.close()

outFile = open("reaction_pathway.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_pathway)
outFile.close()

outFile = open("reaction_enzyme.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(reaction_enzyme)
outFile.close()

print("reaction table done")
print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()

#parse pathway
soup = BeautifulSoup(open("fly_pathway_full.xml"), "xml")
enzyme = soup.find_all("Pathway")
temp = [None]*2
pathway_result = []

for row in enzyme:
    if row.has_attr('ID'):
        temp[0] = row['ID'].strip()
        temp[0] = temp[0].decode("utf-8")
        temp[0] = temp[0][4:]
        if row.find("common-name"):
            temp[1] = row.find("common-name").string.strip().encode("utf-8")
        else:
            temp[1] = ''
        pathway_result.append(temp)
        temp = [None]*2

outFile = open("pathway.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(pathway_result)
outFile.close()

print("pathway table done")
print("--- %s seconds ---" % (time.time() - start_time))

#parse enzyme
soup = BeautifulSoup(open("fly_enzyme_full.xml"), "xml")
enzyme = soup.find_all("Protein")
temp = []
temp_syn = []
enzyme_result = []
syn_result = []
temp_go = []
go_result = []
temp_gene = []
gene_result = []
temp_location = []
location_result = []

for row in enzyme:
    if row.has_attr('ID'):
        temp.append(row['ID'].strip())
        temp[0] = temp[0].decode("utf-8")
        temp[0] = temp[0][4:]
        enzyme_result.append(temp)
        if row.find("common-name"):
            temp_syn.append(temp[0])
            temp_syn.append(row.find("common-name").string.strip().encode("utf-8"))
        if row.find("synonym"):
            if len(temp_syn) == 0:
                temp_syn.append(temp[0])
                temp_syn.append("")
            synonym = row.find_all("synonym")
            for item in synonym:
                temp_syn.append(item.text.strip().encode("utf-8"))
        if row.find("has-go-term"):
            temp_go.append(temp[0])
            go = row.find_all("has-go-term")
            for item in go:
                temp_go.append(item.find("GO-Term")['frameid'].strip().encode("utf-8"))
            go_result.append(temp_go)
            temp_go = []
        if row.find("Gene"):
            temp_gene.append(temp[0])
            temp_gene.append(row.find("Gene")['frameid'].strip().encode("utf-8"))
            gene_result.append(temp_gene)
            temp_gene = []
        if row.find("location"):
            location = row.find("cco")
            temp_location.append(temp[0])
            temp_location.append(location['frameid'].strip().encode("utf-8"))
            location_result.append(temp_location)
            temp_location = []
        if temp_syn:
            syn_result.append(temp_syn)
            temp_syn = []
        temp = []

outFile = open("enzyme.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(enzyme_result)
outFile.close()

outFile = open("enzyme_synonym.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(syn_result)
outFile.close()

outFile = open("enzyme_go.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(go_result)
outFile.close()

outFile = open("enzyme_gene.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(gene_result)
outFile.close()

outFile = open("enzyme_location.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(location_result)
outFile.close()

print("enzyme table done")
print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()

#parse gene
soup = BeautifulSoup(open("fly_gene_full.xml"), "xml")
enzyme = soup.find_all("Gene")
temp = [None]*2
gene_result = []

for row in enzyme:
    if row.has_attr('ID'):
        temp[0] = row['ID'].strip()
        temp[0] = temp[0].decode("utf-8")
        temp[0] = temp[0][4:]
        if row.find("common-name"):
            temp[1] = row.find("common-name").string.strip().encode("utf-8")
        else:
            temp[1] = ''
        if row.find("synonym"):
            synonym = row.find_all("synonym")
            for item in synonym:
                temp.append(item.string.strip().encode("utf-8"))
        gene_result.append(temp)
        temp = [None]*2

outFile = open("gene.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(gene_result)
outFile.close()

#parse id
soup = BeautifulSoup(open("fly_compound_full.xml"), "xml")
compound = soup.find_all("Compound")
temp = [None]*3
ID = []

for row in compound:
    if row.has_attr('ID'):
        if row.find("dblink"):
            db = row.find_all("dblink")
            for item in db:
                temp[0] = row['ID'].strip()
                temp[0] = temp[0].decode("utf-8")
                temp[0] = temp[0][4:]
                db_type = item.find("dblink-db")
                db_id = item.find("dblink-oid")
                temp[1] = db_type.string.strip().encode("utf-8")
                temp[2] = db_id.string.strip().encode("utf-8")
                ID.append(temp)
                temp = [None]*3

outFile = open("compound_ExRef.csv", "wb")
writer = csv.writer(outFile)
writer.writerows(ID)
outFile.close()