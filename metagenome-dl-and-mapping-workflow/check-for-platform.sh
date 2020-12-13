rm -rf platform-info.tsv

printf "\n"

for ERR in $(cat $1)
do

    sample=$(grep "${ERR}" $2 | cut -f 1)

    printf "Checking on: ${sample}\r"

    esearch -db sra -query "${ERR}" | efetch > info.tmp

    if grep -q "Illumina" info.tmp; then
        platform="Illumina"
    elif grep -q "454" info.tmp; then
        platform="454"
    else
        platform="not sure"
    fi

    printf "${sample}\t${ERR}\t${platform}\n" >> platform-info.tsv

    rm info.tmp

done

printf "\r\n\n"
