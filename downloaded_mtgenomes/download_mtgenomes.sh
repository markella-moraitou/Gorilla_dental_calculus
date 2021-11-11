#Look for gorilla mitochondrial genomes to include in the haplotype network
#First - download mt genomes from: Van der Valk (2019) 'Significant loss of mitochondrial diversity within the last century due to extinction of peripheral populations in eastern gorillas': mitochondrial genome sequences are available at Genbank under accession numbers MH177628 - MH177754
#Second - search NCBI using the query "gorilla mitochondrion
 

#Get accession numbers and some metadata for the genomes of the Van der Valk paper
#echo -e "accession\tindividual\torganism\tgenome\tcompleteness" > gorilla_mt_genomes_to_download.txt

#for i in {177628..177754}
#do
#    id="MH"$i
    #Get some info
#    echo $id | esummary -db nuccore > summary.temp
#    accession=$(cat summary.temp | grep "AccessionVersion" | cut -d ">" -f 2 | cut -d "<" -f 1 )
#    individual=$(cat summary.temp | grep "Title" | cut -d " " -f 8)
#    organism=$(cat summary.temp | grep "Organism" | cut -d ">" -f 2 | cut -d "<" -f 1 )
#    genome=$(cat summary.temp | grep "Genome" | cut -d ">" -f 2 | cut -d "<" -f 1 )
#    completeness=$(cat summary.temp | grep "Completeness" | cut -d ">" -f 2 | cut -d "<" -f 1 )
#    echo -e "${accession}\t${individual}\t${organism}\t${genome}\t${completeness}" >> gorilla_mt_genomes_to_download.txt
#done

#rm summary.temp


#Search nuccore for the term "gorilla mitochondrion" - save accession numbers and some metadata

#Save accession numbers that appear when using looking for gorilla mitochondria in the NCBI database
esearch -db nuccore -query "gorilla mitochondrion" | esummary | grep "AccessionVersion" | cut -d ">" -f 2 | cut -d "<" -f 1 > gorilla_mt_accession_numbers.txt

#For each accession number get some key info
cat gorilla_mt_accession_numbers.txt | while read i
do
    esummary -id $i -db nuccore > summary.temp
    accession=$(cat summary.temp | grep "AccessionVersion" | cut -d ">" -f 2 | cut -d "<" -f 1 )
    #If this accession was already added in the previous step, skip
    grep -qP $accession gorilla_mt_genomes_to_download.txt && continue
    #Else, continue adding metadata
    individual=$(cat summary.temp | grep "Title" | cut -d " " -f 8)
    organism=$(cat summary.temp | grep "Organism" | cut -d ">" -f 2 | cut -d "<" -f 1 )
    genome=$(cat summary.temp | grep "Genome" | cut -d ">" -f 2 | cut -d "<" -f 1 )
    completeness=$(cat summary.temp | grep "Completeness" | cut -d ">" -f 2 | cut -d "<" -f 1 )
    echo -e "${accession}\t${individual}\t${organism}\t${genome}\t${completeness}" >> gorilla_mt_genomes_to_download.txt
done

rm summary.temp
rm gorilla_mt_accession_numbers.txt 

#Filter table to keep only mitochondrial genomes from gorillas
head -1 gorilla_mt_genomes_to_download.txt > gorilla_mt_genomes_to_download_filtered.txt
awk -v FS="\t" '$3 ~ "gorilla" && $4=="mitochondrion" {print $0}' gorilla_mt_genomes_to_download.txt >> gorilla_mt_genomes_to_download_filtered.txt

#Replace original list with the filtered one
mv gorilla_mt_genomes_to_download_filtered.txt gorilla_mt_genomes_to_download.txt

#Rename e.g. eastern lowland to Grauer's gorilla
sed -i "s/eastern mountain gorilla/Gbb/g" gorilla_mt_genomes_to_download.txt
sed -i "s/eastern lowland gorilla/Gbg/g" gorilla_mt_genomes_to_download.txt
sed -i "s/western lowland gorilla/Ggg/g" gorilla_mt_genomes_to_download.txt
sed -i "s/western gorilla/Gg/g" gorilla_mt_genomes_to_download.txt
sed -i "s/eastern gorilla/Gb/g" gorilla_mt_genomes_to_download.txt

#Many genomes are assigned to a species (e.g. Gb), but the individual ID reveals the subspecies (e.g. GBG_89)
#Use this info to update the organism column
awk '$2~"MG" && $3=="Gb" {print $0}' gorilla_mt_genomes_to_download.txt | sed "s/Gb/Gbb/g" >>  gorilla_mt_genomes_to_download_temp.txt
awk '$2~"GBG" && $3=="Gb" {print $0}' gorilla_mt_genomes_to_download.txt | sed "s/Gb/Gbg/g" >>  gorilla_mt_genomes_to_download_temp.txt

#The add the rest of entries that contain subspecies info
awk '$3=="Gbb" {print $0}' gorilla_mt_genomes_to_download.txt >>  gorilla_mt_genomes_to_download_temp.txt
awk '$3=="Gbg" {print $0}' gorilla_mt_genomes_to_download.txt >>  gorilla_mt_genomes_to_download_temp.txt
awk '$3=="Ggg" {print $0}' gorilla_mt_genomes_to_download.txt >>  gorilla_mt_genomes_to_download_temp.txt

mv gorilla_mt_genomes_to_download_temp.txt gorilla_mt_genomes_to_download.txt

#Download fasta sequences based on accession list
tail -n +2 gorilla_mt_genomes_to_download.txt | while read i
do
    #if available, use individual name to name the file
    name=$(echo $i | awk '{print $1$3}')
    echo $name
    epost -db nuccore -id $(echo $i | awk '{print $1}') | efetch -format fasta  > ${name}.fasta
    gzip ${name}.fasta
done
